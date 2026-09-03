use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex};
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
use crate::types::{InputType, QualityControl};

const BATCH_SIZE: usize = 10000;
const CHANNEL_BUFFER: usize = 16;
/// Soft cap on the sequence bytes carried by a single batch. Batching purely by
/// record count is fine for reads but not for reference genomes, where one
/// record is a whole contig: a 10 000-record batch could hold an entire genome
/// and, with `CHANNEL_BUFFER` batches in flight, peak at many times the genome
/// size. Filling stops as soon as this many bases are buffered, bounding
/// in-flight memory to roughly `CHANNEL_BUFFER * BATCH_MAX_BYTES`.
const BATCH_MAX_BYTES: usize = 8 * 1024 * 1024;

/// Largest tag length `canonical_hash` can reverse-complement in its stack
/// buffer. Enforced by a debug assertion there and by a test over every enzyme.
pub const MAX_TAG_LENGTH: usize = 64;

/// Prefix used to name contigs whose FASTA header carries no id. Deliberately
/// unlikely to collide with a real header (a plain `contig1` collides with the
/// very common real contig name `contig1`, which made two different contigs
/// indistinguishable in `--record-pos` output).
const UNNAMED_CONTIG_PREFIX: &str = "unnamed_2bRAD_contig";

pub type Hash = u64;

/// Compute hash of the canonical (lexicographically smaller of forward/RC) sequence.
/// Uses a fixed stack buffer — zero heap allocation.
#[inline]
fn canonical_hash(seq: &[u8]) -> Hash {
    debug_assert!(
        seq.len() <= MAX_TAG_LENGTH,
        "tag of {} bases exceeds canonical_hash's {}-byte buffer",
        seq.len(),
        MAX_TAG_LENGTH
    );
    let mut rc_buf = [0u8; MAX_TAG_LENGTH];
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

/// Reference-genome record carrying the contig id and the tag's position, so the
/// writer can emit both the position-less database file and the optional
/// `contig|pos` file from the same stream.
struct GenomeTask {
    hash: Hash,
    contig: String,
    pos: usize,
}

/// Which bases the quality filters look at for read input (`-t 2` / `-t 3`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QcScope {
    /// Filter on the whole read: a read with too many N or too many low-quality
    /// bases anywhere is dropped, tags and all. This is what `-t 2` has always
    /// done and mirrors a pre-filter applied to the FASTQ as a whole.
    Read,
    /// Filter on the tag window only: bases far away from the enzyme site do not
    /// disqualify the tag. This is what `-t 3` has always done, and it is the
    /// more sensible choice for long shotgun reads.
    Tag,
}

impl QcScope {
    fn parse(value: &str, input_type: InputType) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "read" => Ok(QcScope::Read),
            "tag" => Ok(QcScope::Tag),
            // Default per input type: keep each path's historical behaviour.
            "auto" => Ok(match input_type {
                InputType::Single2bRAD => QcScope::Tag,
                _ => QcScope::Read,
            }),
            other => bail!("Invalid --qc-scope: {} (expected read, tag or auto)", other),
        }
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExtractArgs {
    // ── Input (choose one of --list or -i) ──
    /// Batch mode: genome/sample list file (TSV: name<TAB>path1[<TAB>path2]). Mutually exclusive with -i
    #[arg(
        short = 'l',
        long = "list",
        // Accepted spellings used by the docs and helper scripts.
        alias = "genome-list",
        alias = "batch",
        conflicts_with = "input",
        required_unless_present = "input",
        help_heading = "Input"
    )]
    pub genome_list: Option<PathBuf>,
    /// Single mode: one or two FASTQ/FASTA files (PE: R1 R2). Mutually exclusive with --list
    #[arg(
        short = 'i',
        long = "input",
        num_args = 1..=2,
        conflicts_with = "genome_list",
        required_unless_present = "genome_list",
        help_heading = "Input"
    )]
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
    /// Output file prefix (required in single mode with -i; in batch mode the list's first column is used)
    #[arg(
        long = "op",
        num_args = 1,
        required_unless_present = "genome_list",
        help_heading = "Output"
    )]
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
    /// Which bases the QC filters apply to: read=whole read, tag=tag window only,
    /// auto=read for -t 2 and tag for -t 3 (historical per-type behaviour)
    #[arg(long = "qc-scope", default_value = "auto", help_heading = "Quality Control")]
    pub qc_scope: String,
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
    /// Number of parallel threads. In batch mode this is the number of samples
    /// processed concurrently; each one additionally runs a reader and a writer
    /// I/O thread (plus --pc threads if PEAR is enabled), so the process holds
    /// roughly THREADS*(1+2) OS threads
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

/// Validate the single-sample arguments up front. These used to be indexed
/// blindly (`args.input[0]`, `args.output_prefix[0]`), which panicked with a
/// bare "index out of bounds" whenever `--op` was forgotten or the args were
/// built programmatically (e.g. from `pipeline`).
fn single_sample_inputs(args: &ExtractArgs) -> Result<(&[PathBuf], &str)> {
    if args.input.is_empty() {
        bail!("No input file given: pass -i/--input (one or two files), or -l/--list for batch mode");
    }
    if args.input.len() > 2 {
        bail!("-i/--input takes at most two files (SE: one; PE: R1 R2), got {}", args.input.len());
    }
    let prefix = args
        .output_prefix
        .first()
        .map(|s| s.as_str())
        .filter(|s| !s.is_empty())
        .ok_or_else(|| anyhow!("--op (output prefix) is required in single-sample mode (-i)"))?;
    Ok((&args.input, prefix))
}

/// Every file this run may create, so a failed run can clean up after itself
/// instead of leaving a truncated or zero-byte database behind — downstream
/// `build-*-db`/`quantify` cannot tell such a file from a genuinely tag-free
/// sample.
fn expected_outputs(args: &ExtractArgs, prefix: &str, enzyme: &Enzyme, input_type: InputType) -> Vec<PathBuf> {
    let base = format!("{}.{}", prefix, enzyme.name);
    let mut paths = vec![args.output_dir.join(format!("{}.stat.tsv", base))];
    match input_type {
        InputType::ReferenceGenome => {
            paths.push(args.output_dir.join(format!("{}.iibdb", base)));
            paths.push(args.output_dir.join(format!("{}.pos.iibdb", base)));
        }
        InputType::ShotgunMetagenome | InputType::Single2bRAD => {
            paths.push(args.output_dir.join(format!("{}.iibsp", base)));
            // PEAR intermediates, in case the run died between PEAR and cleanup.
            for suffix in [
                ".pear.fastq",
                ".assembled.fastq",
                ".unassembled.forward.fastq",
                ".unassembled.reverse.fastq",
                ".discarded.fastq",
            ] {
                paths.push(args.output_dir.join(format!("{}{}", base, suffix)));
            }
        }
    }
    paths
}

fn remove_outputs(paths: &[PathBuf]) {
    for path in paths {
        if path.exists() {
            if let Err(e) = std::fs::remove_file(path) {
                tracing::warn!("Failed to remove incomplete output {}: {}", path.display(), e);
            }
        }
    }
}

fn run_single_sample(args: ExtractArgs) -> Result<()> {
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let input_type = InputType::from_u8(args.input_type)
        .ok_or_else(|| anyhow!("Invalid input type: {} (expected 1, 2 or 3)", args.input_type))?;
    let (inputs, prefix) = single_sample_inputs(&args)?;
    let qc_scope = QcScope::parse(&args.qc_scope, input_type)?;

    if args.record_pos && input_type != InputType::ReferenceGenome {
        tracing::warn!(
            "--record-pos only applies to reference genomes (-t 1); ignoring it for -t {}",
            args.input_type
        );
    }

    let qc = QualityControl {
        enabled: args.quality_control.eq_ignore_ascii_case("yes"),
        max_n: args.max_n,
        min_quality: args.min_quality,
        min_quality_percent: args.min_quality_percent,
        quality_base: args.quality_base,
    };

    io_utils::ensure_directory(&args.output_dir)?;

    let result = match input_type {
        InputType::ReferenceGenome => extract_reference_genome(&args, enzyme, inputs, prefix),
        InputType::ShotgunMetagenome => extract_shotgun(&args, enzyme, &qc, qc_scope, inputs, prefix),
        InputType::Single2bRAD => extract_single_tag(&args, enzyme, &qc, qc_scope, inputs, prefix),
    };

    if result.is_err() {
        remove_outputs(&expected_outputs(&args, prefix, enzyme, input_type));
    }
    result
}

fn run_batch_mode(base_args: ExtractArgs, genome_list: &std::path::Path) -> Result<()> {
    use std::io::BufRead;
    tracing::info!("### Batch processing mode: {}", genome_list.display());
    let file = File::open(genome_list)
        .with_context(|| format!("Failed to open list file: {}", genome_list.display()))?;
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
    if samples.is_empty() {
        bail!("No usable entries in list file {}", genome_list.display());
    }

    let pb = ProgressBar::new(samples.len() as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner} {pos}/{len}").unwrap());

    // A sample that fails must not be silently reduced to a log line: its
    // (already removed) output would otherwise be read downstream as a sample
    // with zero tags. Collect the failures and fail the whole run.
    let failures: Mutex<Vec<(String, String)>> = Mutex::new(Vec::new());

    samples.into_par_iter().for_each(|(sample_name, input1, input2)| {
        let mut sample_args = base_args.clone();
        sample_args.genome_list = None;
        sample_args.input = vec![input1];
        if let Some(in2) = input2 { sample_args.input.push(in2); }
        sample_args.output_prefix = vec![sample_name.clone()];

        match run_single_sample(sample_args) {
            Ok(_) => {},
            Err(e) => {
                tracing::error!("Sample {} processing failed: {:#}", sample_name, e);
                failures.lock().unwrap().push((sample_name, format!("{:#}", e)));
            }
        }
        pb.inc(1);
    });
    pb.finish();

    let mut failures = failures.into_inner().unwrap();
    if !failures.is_empty() {
        failures.sort();
        let detail = failures
            .iter()
            .map(|(name, err)| format!("  {}: {}", name, err))
            .collect::<Vec<_>>()
            .join("\n");
        bail!(
            "{} sample(s) failed during extraction (their partial output files were removed):\n{}",
            failures.len(),
            detail
        );
    }
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

fn spawn_reader_thread(input_path: PathBuf, max_batch_bytes: usize) -> (
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
                let _ = work_tx.send(Err(anyhow!(e)
                    .context(format!("Failed to open input file: {}", input_path.display()))));
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
            let mut bytes = 0usize;
            let mut exhausted = false;

            for i in 0..BATCH_SIZE {
                match reader.next() {
                    Some(Ok(rec)) => {
                        // [Memory reuse] No new allocation here; reuses the Vec inside batch[i]
                        batch[i].populate_from(&rec);
                        bytes += batch[i].seq.len();
                        count += 1;
                        // Stop early on very long records (contigs) to bound
                        // the memory held by the in-flight batches.
                        if bytes >= max_batch_bytes { break; }
                    },
                    Some(Err(e)) => {
                        let _ = work_tx.send(Err(anyhow!(e)
                            .context(format!("Fastx parse error in {}", input_path.display()))));
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

/// Consume every batch from one reader, recycling containers even when the
/// callback fails. Takes `work_rx` by value so it is dropped on return, which
/// unblocks a reader thread that is waiting on a full channel.
fn drain_batches<F>(
    work_rx: mpsc::Receiver<Result<BatchData>>,
    recycle_tx: &mpsc::Sender<BatchData>,
    mut f: F,
) -> Result<()>
where
    F: FnMut(&[RawRecord]) -> Result<()>,
{
    while let Ok(item) = work_rx.recv() {
        let (batch, count) = item?;
        let outcome = f(&batch[..count]);
        let _ = recycle_tx.send((batch, 0));
        outcome?;
    }
    Ok(())
}

/// Stream one file through `f`, always joining the reader thread afterwards.
fn stream_file<F>(path: &Path, max_batch_bytes: usize, f: F) -> Result<()>
where
    F: FnMut(&[RawRecord]) -> Result<()>,
{
    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(path.to_path_buf(), max_batch_bytes);
    let result = drain_batches(work_rx, &recycle_tx, f);
    drop(recycle_tx);
    let _ = reader_handle.join();
    result
}

/// Stream two mate files in lockstep, handing `f` the matching slices of R1 and
/// R2. Both readers use the same record-count batching (no byte cap), so batch
/// boundaries line up as long as the files hold the same number of reads — and
/// if they do not, that is reported rather than silently mis-pairing.
fn stream_paired<F>(path1: &Path, path2: &Path, mut f: F) -> Result<()>
where
    F: FnMut(&[RawRecord], &[RawRecord]) -> Result<()>,
{
    let (rx1, recycle1, handle1) = spawn_reader_thread(path1.to_path_buf(), usize::MAX);
    let (rx2, recycle2, handle2) = spawn_reader_thread(path2.to_path_buf(), usize::MAX);

    let mismatch = || {
        anyhow!(
            "Paired input files hold different numbers of reads: {} vs {}",
            path1.display(),
            path2.display()
        )
    };

    let result = (|| -> Result<()> {
        loop {
            let item1 = rx1.recv();
            let item2 = rx2.recv();
            match (item1, item2) {
                (Ok(a), Ok(b)) => {
                    let (batch1, count1) = a?;
                    let (batch2, count2) = b?;
                    let outcome = if count1 != count2 {
                        Err(mismatch())
                    } else {
                        f(&batch1[..count1], &batch2[..count2])
                    };
                    let _ = recycle1.send((batch1, 0));
                    let _ = recycle2.send((batch2, 0));
                    outcome?;
                }
                (Err(_), Err(_)) => break,
                // One file ran out first. Surface that file's own read error if
                // it had one, rather than blaming the length mismatch.
                (Ok(a), Err(_)) => {
                    a?;
                    return Err(mismatch());
                }
                (Err(_), Ok(b)) => {
                    b?;
                    return Err(mismatch());
                }
            }
        }
        Ok(())
    })();

    drop(rx1);
    drop(rx2);
    drop(recycle1);
    drop(recycle2);
    let _ = handle1.join();
    let _ = handle2.join();
    result
}

/// Join a writer thread, turning a panic into an error instead of a second panic.
fn join_writer(handle: thread::JoinHandle<Result<()>>) -> Result<()> {
    match handle.join() {
        Ok(result) => result,
        Err(_) => Err(anyhow!("Writer thread panicked")),
    }
}

// ========== Type 1: Reference genome ==========

fn extract_reference_genome(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    inputs: &[PathBuf],
    prefix: &str,
) -> Result<()> {
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
        writer.flush()?;
        if let Some(pw) = pos_writer.as_mut() { pw.flush()?; }
        Ok(())
    });

    let input_sequences = Arc::new(AtomicUsize::new(0));
    let total_bases = Arc::new(AtomicUsize::new(0));
    let total_tags = Arc::new(AtomicUsize::new(0));
    let auto_numbered = Arc::new(AtomicUsize::new(0));

    let consumer_result = (|| -> Result<()> {
        // Contigs are read in order; track the running base so each record gets
        // a stable 1-based ordinal for auto-numbering, across all input files.
        let mut contig_base = 0usize;
        for input_path in inputs {
            tracing::info!("Digesting: {}", input_path.display());
            stream_file(input_path, BATCH_MAX_BYTES, |batch| {
                process_genome_batch(
                    batch,
                    enzyme,
                    &write_tx,
                    &input_sequences,
                    &total_bases,
                    &total_tags,
                    contig_base,
                    &auto_numbered,
                )?;
                contig_base += batch.len();
                Ok(())
            })?;
        }
        Ok(())
    })();

    drop(write_tx);
    let writer_result = join_writer(writer_handle);
    consumer_result?;
    writer_result?;

    let auto_n = auto_numbered.load(Ordering::Relaxed);
    if auto_n > 0 {
        tracing::warn!(
            "Genome {}: {} contig(s) had no sequence ID in the FASTA header; they were auto-numbered as {}<N> (N = 1-based contig order).",
            prefix,
            auto_n,
            UNNAMED_CONTIG_PREFIX
        );
    }

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", prefix, enzyme.name));
    io_utils::write_genome_stats(
        &stat_path,
        prefix,
        enzyme.name,
        input_sequences.load(Ordering::Relaxed),
        total_bases.load(Ordering::Relaxed),
        total_tags.load(Ordering::Relaxed),
    )?;

    Ok(())
}

fn process_genome_batch(
    batch: &[RawRecord], // changed to slice
    enzyme: &Enzyme,
    tx: &mpsc::SyncSender<Vec<GenomeTask>>,
    count_seq: &AtomicUsize,
    count_bases: &AtomicUsize,
    count_tag: &AtomicUsize,
    contig_base: usize,
    auto_numbered: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);
    count_bases.fetch_add(batch.iter().map(|r| r.seq.len()).sum::<usize>(), Ordering::Relaxed);

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
                format!("{}{}", UNNAMED_CONTIG_PREFIX, contig_base + local_idx + 1)
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

// ========== Sample tag output (-t 2 / -t 3) ==========

/// Writer thread for the id-less sample tag stream. Read ids are not stored:
/// nothing downstream consumes them, and writing them made R1/R2 collide on
/// identical Illumina names (see io_utils' format notes).
fn spawn_sample_writer(output_path: PathBuf) -> (mpsc::SyncSender<Vec<Hash>>, thread::JoinHandle<Result<()>>) {
    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<Hash>>(CHANNEL_BUFFER);
    let handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        io_utils::write_sample_tag_header(&mut writer)?;
        while let Ok(batch) = write_rx.recv() {
            for hash in batch {
                io_utils::write_sample_tag(&mut writer, hash)?;
            }
        }
        writer.flush()?;
        Ok(())
    });
    (write_tx, handle)
}

/// Collect every tag hash carried by one read, applying QC at `scope`.
fn read_tag_hashes(
    record: &RawRecord,
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    out: &mut Vec<Hash>,
) {
    // Uppercase *before* QC: `check_n` counts N bases, and a lower-case `n`
    // used to slip through here while the -t 3 path (which uppercases first)
    // rejected the very same read.
    let mut sequence = record.seq.clone();
    sequence.make_ascii_uppercase();

    if scope == QcScope::Read {
        if !qc.check_n(&sequence) { return; }
        if !qc.check_quality(&record.qual) { return; }
    }

    for (pos, len) in enzyme.find_all_tags(&sequence) {
        let window = &sequence[pos..pos + len];
        if scope == QcScope::Tag {
            if !qc.check_n(window) { continue; }
            // Reads whose quality string is shorter than the tag (or absent,
            // i.e. FASTA input) simply skip the quality check.
            let qual_window = record.qual.get(pos..pos + len).unwrap_or(&[]);
            if !qc.check_quality(qual_window) { continue; }
        }
        out.push(canonical_hash(window));
    }
}

/// Merge the tags of one mate pair, dropping tags seen on both mates. Overlapping
/// mates cover the same physical fragment, so the same tag appearing in R1 and R2
/// is one observation, not two — counting it twice biased quantification.
fn merge_pair_tags(mut from_r1: Vec<Hash>, from_r2: Vec<Hash>) -> Vec<Hash> {
    let r1_len = from_r1.len();
    for hash in from_r2 {
        if !from_r1[..r1_len].contains(&hash) {
            from_r1.push(hash);
        }
    }
    from_r1
}

// ========== Type 2: Shotgun ==========

fn extract_shotgun(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    inputs: &[PathBuf],
    prefix: &str,
) -> Result<()> {
    let use_pear = inputs.len() == 2 && args.use_pear.eq_ignore_ascii_case("yes");
    // With PEAR, merging produces a single file that is processed on its own.
    let merged = if use_pear {
        tracing::info!("Run PEAR merging (use-pear=yes) ...");
        Some(run_pear_and_combine(args, enzyme, inputs, prefix)?)
    } else {
        if inputs.len() == 2 {
            tracing::info!("Skip PEAR merging (use-pear=no), process R1/R2 as pairs ...");
        }
        None
    };

    let output_path = args.output_dir.join(format!("{}.{}.iibsp", prefix, enzyme.name));
    let (write_tx, writer_handle) = spawn_sample_writer(output_path);

    // Shared statistics: accumulated across files
    let input_sequences = Arc::new(AtomicUsize::new(0));
    let tag_count = Arc::new(AtomicUsize::new(0));

    let consumer_result = (|| -> Result<()> {
        if let Some(merged_path) = merged.as_ref() {
            tracing::info!("Extracting file: {}", merged_path.display());
            return stream_file(merged_path, BATCH_MAX_BYTES, |batch| {
                process_shotgun_batch(batch, enzyme, qc, scope, &write_tx, &input_sequences, &tag_count)
            });
        }
        if inputs.len() == 2 {
            // Paired: read R1/R2 in lockstep so a fragment seen by both mates
            // contributes its tags once.
            return stream_paired(&inputs[0], &inputs[1], |b1, b2| {
                process_shotgun_pair_batch(b1, b2, enzyme, qc, scope, &write_tx, &input_sequences, &tag_count)
            });
        }
        for input_path in inputs {
            tracing::info!("Extracting file: {}", input_path.display());
            stream_file(input_path, BATCH_MAX_BYTES, |batch| {
                process_shotgun_batch(batch, enzyme, qc, scope, &write_tx, &input_sequences, &tag_count)
            })?;
        }
        Ok(())
    })();

    drop(write_tx);
    let writer_result = join_writer(writer_handle);

    // The merged FASTQ is a scratch file; remove it whether or not we succeeded.
    if let Some(merged_path) = merged.as_ref() {
        if let Err(e) = std::fs::remove_file(merged_path) {
            tracing::warn!("Failed to remove PEAR intermediate {}: {}", merged_path.display(), e);
        }
    }

    consumer_result?;
    writer_result?;

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", prefix, enzyme.name));
    io_utils::write_read_stats(
        &stat_path,
        prefix,
        enzyme.name,
        input_sequences.load(Ordering::Relaxed),
        tag_count.load(Ordering::Relaxed),
    )?;

    Ok(())
}

fn process_shotgun_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    tx: &mpsc::SyncSender<Vec<Hash>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    let results: Vec<Hash> = batch.par_iter().flat_map(|record| {
        let mut tags = Vec::new();
        read_tag_hashes(record, enzyme, qc, scope, &mut tags);
        tags
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

fn process_shotgun_pair_batch(
    batch1: &[RawRecord],
    batch2: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    tx: &mpsc::SyncSender<Vec<Hash>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch1.len() + batch2.len(), Ordering::Relaxed);

    let results: Vec<Hash> = batch1.par_iter().zip(batch2.par_iter()).flat_map(|(r1, r2)| {
        let mut tags1 = Vec::new();
        read_tag_hashes(r1, enzyme, qc, scope, &mut tags1);
        let mut tags2 = Vec::new();
        read_tag_hashes(r2, enzyme, qc, scope, &mut tags2);
        merge_pair_tags(tags1, tags2)
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

fn run_pear_and_combine(
    args: &ExtractArgs,
    enzyme: &Enzyme,
    inputs: &[PathBuf],
    prefix: &str,
) -> Result<PathBuf> {
    let r1 = &inputs[0];
    let r2 = &inputs[1];
    let base = args.output_dir.join(format!("{}.{}", prefix, enzyme.name));
    // Paths are passed as OsStr, so non-UTF-8 paths work instead of panicking.
    let output = Command::new(&args.pear_bin)
        .arg("-f").arg(r1)
        .arg("-r").arg(r2)
        .arg("-e")
        .arg("-o").arg(&base)
        .arg("-j").arg(args.pear_threads.to_string())
        .output()
        .with_context(|| format!("Failed to run PEAR ({})", args.pear_bin))?;
    if !output.status.success() {
        bail!(
            "PEAR failed ({}) for {} / {}:\n{}",
            output.status,
            r1.display(),
            r2.display(),
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }
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

/// What one single-tag read yielded, so the per-read statistics stay exact.
enum SingleTagOutcome {
    /// No enzyme site anywhere in the read.
    NoSite,
    /// Site found but the tag failed QC.
    QcFailed,
    /// Site found and the tag passed QC.
    Passed(Hash),
}

fn extract_single_tag(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    inputs: &[PathBuf],
    prefix: &str,
) -> Result<()> {
    let output_path = args.output_dir.join(format!("{}.{}.iibsp", prefix, enzyme.name));
    let (write_tx, writer_handle) = spawn_sample_writer(output_path);

    let input_sequences = Arc::new(AtomicUsize::new(0));
    let enzyme_reads = Arc::new(AtomicUsize::new(0));
    let qc_passed = Arc::new(AtomicUsize::new(0));
    let tag_count = Arc::new(AtomicUsize::new(0));

    let consumer_result = (|| -> Result<()> {
        if inputs.len() == 2 {
            // Both mates carry the same 2bRAD tag; pair them so it is written once.
            tracing::info!(
                "Extracting paired files: {} + {}",
                inputs[0].display(),
                inputs[1].display()
            );
            return stream_paired(&inputs[0], &inputs[1], |b1, b2| {
                process_single_tag_pair_batch(
                    b1, b2, enzyme, qc, scope, &write_tx,
                    &input_sequences, &enzyme_reads, &qc_passed, &tag_count,
                )
            });
        }
        for input_path in inputs {
            tracing::info!("Extracting file: {}", input_path.display());
            stream_file(input_path, BATCH_MAX_BYTES, |batch| {
                process_single_tag_batch(
                    batch, enzyme, qc, scope, &write_tx,
                    &input_sequences, &enzyme_reads, &qc_passed, &tag_count,
                )
            })?;
        }
        Ok(())
    })();

    drop(write_tx);
    let writer_result = join_writer(writer_handle);
    consumer_result?;
    writer_result?;

    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", prefix, enzyme.name));
    let seqs = input_sequences.load(Ordering::Relaxed);
    let passed = qc_passed.load(Ordering::Relaxed);
    let mut stat_file = File::create(&stat_path)
        .with_context(|| format!("Failed to write sample statistics: {}", stat_path.display()))?;
    writeln!(stat_file, "sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\ttag_count\tpercent")?;
    let percent = if seqs > 0 { (passed as f64 / seqs as f64) * 100.0 } else { 0.0 };
    writeln!(stat_file, "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}%",
        prefix, enzyme.name, seqs, enzyme_reads.load(Ordering::Relaxed), passed,
        tag_count.load(Ordering::Relaxed), percent)?;

    Ok(())
}

/// Find the first (left-most) tag in a single-tag read and QC it.
fn single_tag_outcome(
    record: &RawRecord,
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
) -> SingleTagOutcome {
    let mut sequence = record.seq.clone();
    sequence.make_ascii_uppercase();

    if scope == QcScope::Read {
        if !qc.check_n(&sequence) { return SingleTagOutcome::QcFailed; }
        if !qc.check_quality(&record.qual) { return SingleTagOutcome::QcFailed; }
    }

    // Delegate to Enzyme::find_first_tag, which tries each pattern (fwd
    // then rev) in order and returns the left-most match for the first
    // pattern that has any hit at all — this mirrors `Single_Lable`'s
    // per-site, leftmost, match-then-`last` behaviour in
    // 2bRADExtraction.pl, and works uniformly for both fixed-byte
    // (`Exact`) and IUPAC-degenerate (`Degenerate`) enzymes.
    //
    // The whole read is searched. This used to be truncated to the first 50
    // bases, which silently dropped every tag starting past base 50-tag_length
    // (offset 18 for a 32 bp BcgI tag) — including all tags in untrimmed
    // 100/150 bp reads.
    let (offset, len) = match enzyme.find_first_tag(&sequence) {
        Some(hit) => hit,
        None => return SingleTagOutcome::NoSite,
    };
    let window = &sequence[offset..offset + len];

    if scope == QcScope::Tag {
        if !qc.check_n(window) { return SingleTagOutcome::QcFailed; }
        let qual_window = record.qual.get(offset..offset + len).unwrap_or(&[]);
        if !qc.check_quality(qual_window) { return SingleTagOutcome::QcFailed; }
    }

    SingleTagOutcome::Passed(canonical_hash(window))
}

/// Fold per-read outcomes into the statistics counters, returning the passing hashes.
fn tally_single_tag(
    outcomes: impl IntoIterator<Item = SingleTagOutcome>,
    count_enz: &AtomicUsize,
    count_qc: &AtomicUsize,
    out: &mut Vec<Hash>,
) {
    let mut enz = 0usize;
    let mut passed = 0usize;
    for outcome in outcomes {
        match outcome {
            SingleTagOutcome::NoSite => {}
            SingleTagOutcome::QcFailed => enz += 1,
            SingleTagOutcome::Passed(hash) => {
                enz += 1;
                passed += 1;
                out.push(hash);
            }
        }
    }
    count_enz.fetch_add(enz, Ordering::Relaxed);
    count_qc.fetch_add(passed, Ordering::Relaxed);
}

#[allow(clippy::too_many_arguments)]
fn process_single_tag_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    tx: &mpsc::SyncSender<Vec<Hash>>,
    count_seq: &AtomicUsize,
    count_enz: &AtomicUsize,
    count_qc: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    let outcomes: Vec<SingleTagOutcome> = batch
        .par_iter()
        .map(|record| single_tag_outcome(record, enzyme, qc, scope))
        .collect();

    let mut results = Vec::with_capacity(outcomes.len());
    tally_single_tag(outcomes, count_enz, count_qc, &mut results);

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_single_tag_pair_batch(
    batch1: &[RawRecord],
    batch2: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    scope: QcScope,
    tx: &mpsc::SyncSender<Vec<Hash>>,
    count_seq: &AtomicUsize,
    count_enz: &AtomicUsize,
    count_qc: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch1.len() + batch2.len(), Ordering::Relaxed);

    // Per-read statistics still count both mates; only the written tags are
    // de-duplicated, because both mates describe one fragment.
    let per_pair: Vec<(SingleTagOutcome, SingleTagOutcome, Vec<Hash>)> = batch1
        .par_iter()
        .zip(batch2.par_iter())
        .map(|(r1, r2)| {
            let o1 = single_tag_outcome(r1, enzyme, qc, scope);
            let o2 = single_tag_outcome(r2, enzyme, qc, scope);
            let mut tags1 = Vec::new();
            if let SingleTagOutcome::Passed(h) = o1 { tags1.push(h); }
            let mut tags2 = Vec::new();
            if let SingleTagOutcome::Passed(h) = o2 { tags2.push(h); }
            (o1, o2, merge_pair_tags(tags1, tags2))
        })
        .collect();

    let mut results = Vec::with_capacity(per_pair.len());
    let mut enz = 0usize;
    let mut passed = 0usize;
    for (o1, o2, tags) in per_pair {
        for outcome in [o1, o2] {
            match outcome {
                SingleTagOutcome::NoSite => {}
                SingleTagOutcome::QcFailed => enz += 1,
                SingleTagOutcome::Passed(_) => { enz += 1; passed += 1; }
            }
        }
        results.extend(tags);
    }
    count_enz.fetch_add(enz, Ordering::Relaxed);
    count_qc.fetch_add(passed, Ordering::Relaxed);

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `canonical_hash` reverse-complements into a fixed stack buffer, so no
    /// enzyme may define a tag longer than it.
    #[test]
    fn every_enzyme_tag_fits_the_canonical_hash_buffer() {
        for enzyme in crate::enzymes::ENZYMES {
            assert!(
                enzyme.tag_length <= MAX_TAG_LENGTH,
                "{} has tag_length {} > MAX_TAG_LENGTH {}",
                enzyme.name,
                enzyme.tag_length,
                MAX_TAG_LENGTH
            );
        }
    }

    #[test]
    fn canonical_hash_is_strand_agnostic() {
        let fwd = b"ACGTTGCAAACCGGTTACGTACGTACGTACGT";
        let rc: Vec<u8> = fwd
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                x => *x,
            })
            .collect();
        assert_eq!(canonical_hash(fwd), canonical_hash(&rc));
    }

    #[test]
    fn mate_pair_tags_are_deduplicated() {
        assert_eq!(merge_pair_tags(vec![1, 2], vec![2, 3]), vec![1, 2, 3]);
        assert_eq!(merge_pair_tags(vec![], vec![7]), vec![7]);
        assert_eq!(merge_pair_tags(vec![7], vec![7]), vec![7]);
    }
}
