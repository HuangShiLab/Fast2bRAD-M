use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::process::Command;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
use rayon::prelude::*;
use fxhash::FxHasher;
use std::hash::Hasher;
use tracing;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;
use crate::types::{DigestStats, InputType, QualityControl};

const BATCH_SIZE: usize = 10000;
const CHANNEL_BUFFER: usize = 16;

pub type Hash = u64;

/// Compute hash of the canonical (lexicographically smaller of forward/RC) sequence.
/// Uses a fixed stack buffer — zero heap allocation. tag_length must be ≤ 40.
#[inline]
fn canonical_hash(seq: &[u8]) -> Hash {
    let mut rc_buf = [0u8; 40];
    let len = seq.len();
    for i in 0..len {
        rc_buf[i] = match seq[len - 1 - i] {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            x => x,
        };
    }
    let rc = &rc_buf[..len];
    let canonical = if seq <= rc { seq } else { rc };
    let mut hasher = FxHasher::default();
    hasher.write(canonical);
    hasher.finish()
}

// [Optimization] RawRecord struct adjustments
// 1. id changed to Vec<u8> to avoid UTF-8 validation and String conversion overhead during parsing
// 2. All fields retain capacity for memory reuse
#[derive(Debug, Clone)]
struct RawRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>, // Use empty Vec to represent None, avoiding Option unwrap overhead
}

impl RawRecord {
    fn new() -> Self {
        Self {
            id: Vec::with_capacity(64),
            seq: Vec::with_capacity(150),
            qual: Vec::with_capacity(150),
        }
    }

    // [Core optimization] Memory reuse logic
    // No new memory allocation; data is copied directly into existing buffers
    fn populate_from(&mut self, rec: &SequenceRecord) {
        self.id.clear();
        self.id.extend_from_slice(rec.id());

        self.seq.clear();
        // [Fix] rec.seq() returns Cow<[u8]>; extend_from_slice requires &[u8]
        // Add & to borrow the Cow and use Deref to auto-convert to &[u8]
        self.seq.extend_from_slice(&rec.seq());

        self.qual.clear();
        if let Some(q) = rec.qual() {
            self.qual.extend_from_slice(q);
        }
    }
}

struct WriteTask {
    hash: Hash,
    id_str: String,
}

/// Reference-genome record carrying the contig id and the tag's position, so the
/// writer can emit both the position-less database file and the optional
/// `contig|pos` file from the same stream.
struct GenomeTask {
    hash: Hash,
    contig: String,
    pos: usize,
}

#[derive(Args, Debug, Clone)]
pub struct ExtractArgs {
    // ── Input (choose one of --list or -i) ──
    /// Batch mode: genome/sample list file (TSV: name<TAB>path1[<TAB>path2]). Mutually exclusive with -i
    #[arg(short = 'l', long = "list", conflicts_with = "input", help_heading = "Input")]
    pub genome_list: Option<PathBuf>,
    /// Single mode: one or two FASTQ/FASTA files (PE: R1 R2). Mutually exclusive with --list
    #[arg(short = 'i', long = "input", num_args = 1..=2, conflicts_with = "list", help_heading = "Input")]
    pub input: Vec<PathBuf>,
    /// Input type: 1=reference genome, 2=shotgun reads (SE/PE), 3=single 2bRAD tags
    #[arg(short = 't', long = "type", help_heading = "Input")]
    pub input_type: u8,
    /// Enzyme name (e.g. BcgI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", help_heading = "Input")]
    pub enzyme_site: String,

    // ── Output ──
    /// Output directory
    #[arg(long = "od", help_heading = "Output")]
    pub output_dir: PathBuf,
    /// Output file prefix (only used in single mode with -i)
    #[arg(long = "op", num_args = 1, help_heading = "Output")]
    pub output_prefix: Vec<String>,
    /// Reference genome only (-t 1): also write `{prefix}.{enzyme}.pos.iibdb`
    /// recording each tag's position as `contig|offset` (offset = distance of the
    /// tag's first base from its contig's first base). The default position-less
    /// database is still written either way.
    #[arg(long = "record-pos", help_heading = "Output")]
    pub record_pos: bool,

    // ── Quality Control (for sample reads, -t 2/3) ──
    /// Enable quality control filtering (yes/no)
    #[arg(long = "qc", default_value = "yes", help_heading = "Quality Control")]
    pub quality_control: String,
    /// Maximum allowed N-base ratio per read (reads exceeding this are discarded)
    #[arg(short = 'n', long, default_value = "0.08", help_heading = "Quality Control")]
    pub max_n: f64,
    /// Minimum base quality score (Phred)
    #[arg(short = 'q', long, default_value = "30", help_heading = "Quality Control")]
    pub min_quality: u8,
    /// Minimum percentage of bases that must pass the quality threshold
    #[arg(short = 'p', long, default_value = "80", help_heading = "Quality Control")]
    pub min_quality_percent: u8,
    /// Quality score encoding base (33=Phred+33/Sanger, 64=Phred+64)
    #[arg(short = 'b', long, default_value = "33", help_heading = "Quality Control")]
    pub quality_base: u8,

    // ── PEAR Merging (only for PE reads, -t 2) ──
    /// Enable PEAR merging for paired-end reads (yes/no). Significantly slower when enabled
    #[arg(long = "use-pear", default_value = "no", help_heading = "PEAR Merging (PE only)")]
    pub use_pear: String,
    /// Path to PEAR executable
    #[arg(long = "pe", default_value = "pear", help_heading = "PEAR Merging (PE only)")]
    pub pear_bin: String,
    /// Threads per PEAR process
    #[arg(long = "pc", default_value = "1", help_heading = "PEAR Merging (PE only)")]
    pub pear_threads: usize,

    // ── Performance ──
    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4", help_heading = "Performance")]
    pub threads: usize,
}

pub fn run(args: ExtractArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global();

    if let Some(genome_list) = args.genome_list.clone() {
        return run_batch_mode(args, &genome_list);
    }
    run_single_sample(args)
}

fn run_single_sample(args: ExtractArgs) -> Result<()> {
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let input_type = InputType::from_u8(args.input_type)
        .ok_or_else(|| anyhow!("Invalid input type: {}", args.input_type))?;

    let qc = QualityControl {
        enabled: args.quality_control.eq_ignore_ascii_case("yes"),
        max_n: args.max_n,
        min_quality: args.min_quality,
        min_quality_percent: args.min_quality_percent,
        quality_base: args.quality_base,
    };

    io_utils::ensure_directory(&args.output_dir)?;

    match input_type {
        InputType::ReferenceGenome => extract_reference_genome(&args, enzyme)?,
        InputType::ShotgunMetagenome => extract_shotgun(&args, enzyme, &qc)?,
        InputType::Single2bRAD => extract_single_tag(&args, enzyme, &qc)?,
    }
    Ok(())
}

fn run_batch_mode(base_args: ExtractArgs, genome_list: &std::path::Path) -> Result<()> {
    use std::io::BufRead;
    tracing::info!("### Batch processing mode: {}", genome_list.display());
    let file = File::open(genome_list)?;
    let reader = std::io::BufReader::new(file);
    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 { continue; }
        let sample_name = fields[0].to_string();
        let input1 = PathBuf::from(fields[1]);
        let input2 = if fields.len() > 2 && !fields[2].is_empty() { Some(PathBuf::from(fields[2])) } else { None };
        samples.push((sample_name, input1, input2));
    }

    let pb = ProgressBar::new(samples.len() as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner} {pos}/{len}").unwrap());

    samples.into_par_iter().for_each(|(sample_name, input1, input2)| {
        let mut sample_args = base_args.clone();
        sample_args.genome_list = None;
        sample_args.input = vec![input1];
        if let Some(in2) = input2 { sample_args.input.push(in2); }
        sample_args.output_prefix = vec![sample_name.clone()];

        match run_single_sample(sample_args) {
            Ok(_) => {},
            Err(e) => tracing::error!("Sample {} processing failed: {}", sample_name, e),
        }
        pb.inc(1);
    });
    pb.finish();
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
     if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme); }
     if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme); } }
     bail!("Unknown enzyme: {}", site)
}

// ==========================================
// [Core optimization] General-purpose pipeline reader
// ==========================================
// Spawns a background thread to read the file and uses recycle_rx to receive
// used Batches for memory reuse
type BatchData = (Vec<RawRecord>, usize); // (Buffer container, number of valid records)

fn spawn_reader_thread(input_path: PathBuf) -> (
    mpsc::Receiver<Result<BatchData>>,
    mpsc::Sender<BatchData>,
    thread::JoinHandle<()>
) {
    // work_tx: sends filled data to the consumer
    let (work_tx, work_rx) = mpsc::sync_channel::<Result<BatchData>>(CHANNEL_BUFFER);
    // recycle_tx: consumer returns used containers to the producer
    let (recycle_tx, recycle_rx) = mpsc::channel::<BatchData>();

    let handle = thread::spawn(move || {
        let mut reader = match parse_fastx_file(&input_path) {
            Ok(r) => r,
            Err(e) => {
                let _ = work_tx.send(Err(anyhow!(e).context("Failed to open input file")));
                return;
            }
        };

        loop {
            // 1. Obtain a Batch container (prefer recycled; allocate new if none available)
            let (mut batch, _) = recycle_rx.try_recv().unwrap_or_else(|_| {
                let mut v = Vec::with_capacity(BATCH_SIZE);
                for _ in 0..BATCH_SIZE { v.push(RawRecord::new()); }
                (v, 0)
            });

            // 2. Fill data
            let mut count = 0;
            let mut exhausted = false;

            for i in 0..BATCH_SIZE {
                match reader.next() {
                    Some(Ok(rec)) => {
                        // [Memory reuse] No new allocation here; reuses the Vec inside batch[i]
                        batch[i].populate_from(&rec);
                        count += 1;
                    },
                    Some(Err(e)) => {
                        let _ = work_tx.send(Err(anyhow!(e).context("Fastx parse error")));
                        return;
                    },
                    None => {
                        exhausted = true;
                        break;
                    }
                }
            }

            // 3. Send data
            if count > 0 {
                if work_tx.send(Ok((batch, count))).is_err() {
                    break; // consumer disconnected
                }
            } else {
                // If no data was read this round and input is exhausted, stop sending
                break;
            }

            if exhausted {
                break;
            }
        }
    });

    (work_rx, recycle_tx, handle)
}

// ========== Type 1: Reference genome ==========

fn extract_reference_genome(args: &ExtractArgs, enzyme: &'static Enzyme) -> Result<()> {
    let input_path = args.input[0].clone();
    let prefix = args.output_prefix[0].clone();
    // Default (position-less) database — always written.
    let output_path = args.output_dir.join(format!("{}.{}.iibdb", prefix, enzyme.name));
    // Optional companion file recording each tag's `contig|pos`.
    let pos_output_path = if args.record_pos {
        Some(args.output_dir.join(format!("{}.{}.pos.iibdb", prefix, enzyme.name)))
    } else {
        None
    };

    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<GenomeTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path).context("Failed to create output file")?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        // Second writer (only when --record-pos) records id = "contig|pos".
        let mut pos_writer = match pos_output_path {
            Some(ref p) => {
                let f = File::create(p).context("Failed to create position output file")?;
                Some(BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, f))
            }
            None => None,
        };
        while let Ok(batch) = write_rx.recv() {
            for task in batch {
                io_utils::write_binary_record(&mut writer, task.hash, &task.contig)?;
                if let Some(pw) = pos_writer.as_mut() {
                    let id = format!("{}|{}", task.contig, task.pos);
                    io_utils::write_binary_record(pw, task.hash, &id)?;
                }
            }
        }
        Ok(())
    });

    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);

    let input_sequences = Arc::new(AtomicUsize::new(0));
    let total_tags = Arc::new(AtomicUsize::new(0));
    let auto_numbered = Arc::new(AtomicUsize::new(0));
    // Contigs are read in order; track the running base so each record gets a
    // stable 1-based ordinal for auto-numbering.
    let mut contig_base = 0usize;

    // Main thread consumer loop
    while let Ok(result) = work_rx.recv() {
        let (batch, count) = result?;

        // &batch[..count] is a slice; only valid records are processed
        process_genome_batch(
            &batch[..count],
            enzyme,
            &write_tx,
            &input_sequences,
            &total_tags,
            contig_base,
            &auto_numbered,
        )?;
        contig_base += count;

        // [Recycle] Return the whole container to the reader after processing
        let _ = recycle_tx.send((batch, 0));
    }

    drop(write_tx);
    let _ = reader_handle.join();
    writer_handle.join().unwrap()?;

    let auto_n = auto_numbered.load(Ordering::Relaxed);
    if auto_n > 0 {
        tracing::warn!(
            "Genome {}: {} contig(s) had no sequence ID in the FASTA header; they were auto-numbered as contig<N> (N = 1-based contig order).",
            prefix,
            auto_n
        );
    }

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", prefix, enzyme.name));
    let stats = DigestStats {
        sample_id: prefix.clone(),
        enzyme: enzyme.name.to_string(),
        input_sequences: input_sequences.load(Ordering::Relaxed),
        tag_count: total_tags.load(Ordering::Relaxed),
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    Ok(())
}

fn process_genome_batch(
    batch: &[RawRecord], // changed to slice
    enzyme: &Enzyme,
    tx: &mpsc::SyncSender<Vec<GenomeTask>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
    contig_base: usize,
    auto_numbered: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    // enumerate() over the indexed parallel iterator yields each record's position
    // within the batch, giving a stable global contig ordinal (contig_base + i).
    let results: Vec<GenomeTask> = batch.par_iter().enumerate().flat_map(|(local_idx, record)| {
        // record.seq is already a Vec<u8>; convert to uppercase.
        let mut sequence = record.seq.clone();
        sequence.make_ascii_uppercase();

        let positions_iter = enzyme.find_all_tags(&sequence);

        // Contig id = first whitespace token of the header. If the header has no
        // id (empty/whitespace only), auto-number it and count it for the warning.
        let id_utf8 = String::from_utf8_lossy(&record.id);
        let contig = match id_utf8.split_whitespace().next() {
            Some(tok) => tok.to_string(),
            None => {
                auto_numbered.fetch_add(1, Ordering::Relaxed);
                format!("contig{}", contig_base + local_idx + 1)
            }
        };

        let mut tasks = Vec::new();
        for (pos, len) in positions_iter {
            if pos + len > sequence.len() { continue; }
            let tag_seq = &sequence[pos..pos + len];
            let hash = canonical_hash(tag_seq);
            tasks.push(GenomeTask { hash, contig: contig.clone(), pos });
        }
        tasks
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() {
        tx.send(results).context("Failed to send results to writer")?;
    }
    Ok(())
}

// ========== Type 2: Shotgun ==========

fn extract_shotgun(args: &ExtractArgs, enzyme: &'static Enzyme, qc: &QualityControl) -> Result<()> {
    // [Change] Core logic: decide whether to merge based on use_pear, or build the list of files to process
    let inputs_to_process = if args.input.len() == 2 && args.use_pear.eq_ignore_ascii_case("yes") {
        tracing::info!("Run PEAR merging (use-pear=yes) ...");
        // Merging produces one new file; only process that file
        vec![run_pear_and_combine(args, enzyme)?]
    } else {
        if args.input.len() == 2 {
            tracing::info!("Skip PEAR merging (use-pear=no), process paired-end files sequentially ...");
        }
        // Process the original input files directly (may be 1 or 2)
        args.input.clone()
    };

    let output_path = args.output_dir.join(format!("{}.{}.iibsp", args.output_prefix[0], enzyme.name));

    // Create writer thread: tags extracted from all input files are written to this single file
    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path)?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(batch) = write_rx.recv() {
            for task in batch { io_utils::write_binary_record(&mut writer, task.hash, &task.id_str)?; }
        }
        Ok(())
    });

    // Shared statistics: accumulated across files
    let input_sequences = Arc::new(AtomicUsize::new(0));
    let tag_count = Arc::new(AtomicUsize::new(0));

    // Iterate over all files to process (1 file in PEAR mode, up to 2 otherwise)
    for input_path in inputs_to_process {
        tracing::info!("Extracting file: {}", input_path.display());
        let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);

        // Pipeline loop for a single file
        while let Ok(result) = work_rx.recv() {
            let (batch, count) = result?;
            process_shotgun_batch(&batch[..count], enzyme, qc, &write_tx, &input_sequences, &tag_count)?;
            let _ = recycle_tx.send((batch, 0));
        }

        // Wait for the current file's reader thread to finish before processing the next one
        let _ = reader_handle.join();
    }

    // All files processed; close the write channel and wait for the writer thread
    drop(write_tx);
    writer_handle.join().unwrap()?;

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", args.output_prefix[0], enzyme.name));
    let stats = DigestStats {
        sample_id: args.output_prefix[0].clone(),
        enzyme: enzyme.name.to_string(),
        input_sequences: input_sequences.load(Ordering::Relaxed),
        tag_count: tag_count.load(Ordering::Relaxed),
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    Ok(())
}

fn process_shotgun_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    let results: Vec<WriteTask> = batch.par_iter().flat_map(|record| {
        // QC check
        if !qc.check_n(&record.seq) { return Vec::new(); }
        if !record.qual.is_empty() {
            if !qc.check_quality(&record.qual) { return Vec::new(); }
        }

        let mut sequence = record.seq.clone(); // Clone for modification (uppercase)
        sequence.make_ascii_uppercase();

        let positions = enzyme.find_all_tags(&sequence);
        if positions.is_empty() { return Vec::new(); }

        let id_utf8 = String::from_utf8_lossy(&record.id);
        // Fastq IDs: the ID is everything before the first space
        let seq_id = id_utf8.split_whitespace().next().unwrap_or(&id_utf8);

        let mut tasks = Vec::with_capacity(positions.len());
        for (i, (pos, len)) in positions.iter().enumerate() {
            let tag_seq = &sequence[*pos..*pos + len];
            let tag_hash = canonical_hash(tag_seq);
            let id_str = format!("{}_tag{}", seq_id, i + 1);
            tasks.push(WriteTask { hash: tag_hash, id_str });
        }
        tasks
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

fn run_pear_and_combine(args: &ExtractArgs, enzyme: &Enzyme) -> Result<PathBuf> {
    let r1 = &args.input[0];
    let r2 = &args.input[1];
    let prefix = &args.output_prefix[0];
    let base = args.output_dir.join(format!("{}.{}", prefix, enzyme.name));
    let status = Command::new(&args.pear_bin)
        .args(["-f", r1.to_str().unwrap(), "-r", r2.to_str().unwrap(), "-e", "-o", base.to_str().unwrap(), "-j", &args.pear_threads.to_string()])
        .status()?;
    if !status.success() { bail!("PEAR failed"); }
    let pear_fastq = args.output_dir.join(format!("{}.{}.pear.fastq", prefix, enzyme.name));
    {
        let mut out = File::create(&pear_fastq)?;
        for suffix in [".assembled.fastq", ".unassembled.forward.fastq", ".unassembled.reverse.fastq"] {
             let p = args.output_dir.join(format!("{}.{}{}", prefix, enzyme.name, suffix));
             if p.exists() { std::io::copy(&mut File::open(&p)?, &mut out)?; std::fs::remove_file(p)?; }
        }
    }
    let discarded = args.output_dir.join(format!("{}.{}.discarded.fastq", prefix, enzyme.name));
    if discarded.exists() { std::fs::remove_file(discarded)?; }
    Ok(pear_fastq)
}

// ========== Type 3: Single tag ==========
fn extract_single_tag(args: &ExtractArgs, enzyme: &'static Enzyme, qc: &QualityControl) -> Result<()> {
    let input_path = args.input[0].clone();
    let output_path = args.output_dir.join(format!("{}.{}.iibsp", args.output_prefix[0], enzyme.name));

    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path)?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(batch) = write_rx.recv() { for task in batch { io_utils::write_binary_record(&mut writer, task.hash, &task.id_str)?; } }
        Ok(())
    });

    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);

    let input_sequences = Arc::new(AtomicUsize::new(0));
    let enzyme_reads = Arc::new(AtomicUsize::new(0));
    let qc_passed = Arc::new(AtomicUsize::new(0));

    while let Ok(result) = work_rx.recv() {
        let (batch, count) = result?;
        process_single_tag_batch(&batch[..count], enzyme, qc, &write_tx, &input_sequences, &enzyme_reads, &qc_passed)?;
        let _ = recycle_tx.send((batch, 0));
    }

    drop(write_tx);
    let _ = reader_handle.join();
    writer_handle.join().unwrap()?;

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", args.output_prefix[0], enzyme.name));
    let seqs = input_sequences.load(Ordering::Relaxed);
    let passed = qc_passed.load(Ordering::Relaxed);
    let mut stat_file = File::create(&stat_path)?;
    writeln!(stat_file, "sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent")?;
    let percent = if seqs > 0 { (passed as f64 / seqs as f64) * 100.0 } else { 0.0 };
    writeln!(stat_file, "{}\t{}\t{}\t{}\t{}\t{:.2}%",
        args.output_prefix[0], enzyme.name, seqs, enzyme_reads.load(Ordering::Relaxed), passed, percent)?;

    Ok(())
}

fn process_single_tag_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    count_seq: &AtomicUsize,
    count_enz: &AtomicUsize,
    count_qc: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);
    let results: Vec<WriteTask> = batch.par_iter().filter_map(|record| {
        // Clone seq for safety and logical consistency.
        // For performance, a reference could be used if logic permits,
        // but the original logic includes truncate(50) which mutates the data.
        let mut sequence = record.seq.clone();
        if sequence.len() > 50 { sequence.truncate(50); }
        sequence.make_ascii_uppercase();

        for pattern in enzyme.patterns {
            if sequence.len() < enzyme.tag_length { break; }
            for offset in 0..sequence.len() {
                if offset + enzyme.tag_length > sequence.len() { break; }
                let window = &sequence[offset..offset + enzyme.tag_length];
                if pattern.matches(window) {
                    let mut pass = true;
                    if !qc.check_n(window) { pass = false; }
                    if pass {
                        if !record.qual.is_empty() {
                            if offset + enzyme.tag_length <= record.qual.len() {
                                let qs = &record.qual[offset..offset + enzyme.tag_length];
                                if !qs.is_empty() && !qc.check_quality(qs) { pass = false; }
                            }
                        }
                    }
                    if pass {
                        let hash = canonical_hash(window);
                        // Only convert ID string when confirmed to pass
                        let id_str = String::from_utf8_lossy(&record.id).to_string();
                        return Some((true, WriteTask { hash, id_str }));
                    } else {
                        // Enzyme site matched but QC failed
                        return Some((false, WriteTask { hash: 0, id_str: String::new() }));
                    }
                }
            }
        }
        None
    }).map(|(passed, task)| {
        if passed { (1, 1, Some(task)) } else { (1, 0, None) }
    }).collect::<Vec<_>>().into_iter().fold(Vec::new(), |mut acc, (enz, qc_pass, task_opt)| {
        count_enz.fetch_add(enz, Ordering::Relaxed);
        count_qc.fetch_add(qc_pass, Ordering::Relaxed);
        if let Some(t) = task_opt { acc.push(t); }
        acc
    });
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

