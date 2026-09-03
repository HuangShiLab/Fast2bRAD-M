use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow, bail};
use clap::{Parser, Subcommand};
use needletail::parse_fastx_file;
use tracing;

use f2brad_core::enzymes::{enzyme_by_id, enzyme_by_name};
use f2brad_core::extract::Hash;
use f2brad_core::io_utils::open_compact_reader;
use f2brad_host::genotype::{canonicalize, load_host_db, HostDb};

#[derive(Parser, Debug)]
#[command(name = "f2brad-holo", version, about = "One-pass holo driver for fast2bRAD")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Classify each fragment as host, microbial, both, or neither and report
    /// the host fraction.
    Classify(ClassifyArgs),
}

#[derive(Parser, Debug)]
struct ClassifyArgs {
    /// Host tag database TSV from `f2brad-host build-db`
    #[arg(short = 'd', long = "host-db", required = true)]
    host_db: PathBuf,

    /// Microbial compact database (.iibdb) from `f2brad-m build-*-db`
    #[arg(short = 'm', long = "microbe-db", required = true)]
    microbe_db: PathBuf,

    /// Read 1 FASTQ/FASTA (may be gzip-compressed)
    #[arg(short = '1', long = "r1", required = true)]
    r1: PathBuf,

    /// Optional read 2 FASTQ/FASTA (paired-end)
    #[arg(short = '2', long = "r2")]
    r2: Option<PathBuf>,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    output_dir: PathBuf,

    /// Maximum Hamming distance for host tag matching
    #[arg(long = "host-max-mismatch", default_value = "2")]
    host_max_mismatch: usize,

    /// Number of parallel threads (currently unused; reserved for future batching)
    #[arg(short = 'j', long = "threads", default_value = "4")]
    threads: usize,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();
    match cli.command {
        Commands::Classify(args) => classify::run(args),
    }
}

mod classify {
    use super::*;

    pub fn run(args: ClassifyArgs) -> Result<()> {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();

        let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
            enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
        } else {
            enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
        };

        std::fs::create_dir_all(&args.output_dir)
            .with_context(|| format!("Failed to create output directory: {}", args.output_dir.display()))?;

        tracing::info!("Loading host DB from {}", args.host_db.display());
        let host_db = load_host_db(&args.host_db)?;
        let host_index = f2brad_host::genotype::HostTagIndex::new(&host_db.loci, args.host_max_mismatch);
        tracing::info!("Loaded {} host tags", host_db.loci.len());

        tracing::info!("Loading microbial DB from {}", args.microbe_db.display());
        let microbe_hashes = load_microbe_hashes(&args.microbe_db)?;
        tracing::info!("Loaded {} unique microbial tag hashes", microbe_hashes.len());

        let mut stats = ClassifyStats::default();

        if let Some(r2) = &args.r2 {
            classify_paired(&args.r1, r2, enzyme, &host_db, &host_index, &microbe_hashes, &mut stats)?;
        } else {
            classify_single(&args.r1, enzyme, &host_db, &host_index, &microbe_hashes, &mut stats)?;
        }

        // Report and write stats.
        let total_classified = stats.host_only + stats.microbe_only + stats.both;
        let host_fraction = if total_classified > 0 {
            (stats.host_only + stats.both) as f64 / total_classified as f64
        } else {
            0.0
        };

        tracing::info!("Classified fragments: host_only={}, microbe_only={}, both={}, neither={}",
            stats.host_only, stats.microbe_only, stats.both, stats.neither);
        tracing::info!("Host fraction (of classified fragments): {:.4}", host_fraction);

        let out_path = args.output_dir.join("holo_classify.tsv");
        let file = File::create(&out_path)
            .with_context(|| format!("Failed to create output: {}", out_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "metric\tvalue")?;
        writeln!(writer, "input_fragments\t{}", stats.fragments)?;
        writeln!(writer, "host_only\t{}", stats.host_only)?;
        writeln!(writer, "microbe_only\t{}", stats.microbe_only)?;
        writeln!(writer, "both\t{}", stats.both)?;
        writeln!(writer, "neither\t{}", stats.neither)?;
        writeln!(writer, "host_fraction\t{:.6}", host_fraction)?;
        writer.flush()?;
        tracing::info!("Wrote {}", out_path.display());

        Ok(())
    }

    #[derive(Default)]
    struct ClassifyStats {
        fragments: usize,
        host_only: usize,
        microbe_only: usize,
        both: usize,
        neither: usize,
    }

    fn load_microbe_hashes(path: &PathBuf) -> Result<HashSet<Hash>> {
        let mut reader = open_compact_reader(path)?;
        let mut hashes = HashSet::new();
        while let Some((hash, _gcf_index)) = reader.next_record()? {
            hashes.insert(hash);
        }
        Ok(hashes)
    }

    fn classify_single(
        path: &PathBuf,
        enzyme: &f2brad_core::enzymes::Enzyme,
        host_db: &HostDb,
        host_index: &f2brad_host::genotype::HostTagIndex,
        microbe_hashes: &HashSet<Hash>,
        stats: &mut ClassifyStats,
    ) -> Result<()> {
        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open reads: {}", path.display()))?;

        while let Some(record) = reader.next() {
            let record = record.with_context(|| format!("Failed to read record from {}", path.display()))?;
            let seq = record.seq();
            let qual = record.qual();
            let (is_host, is_microbe) = classify_fragment(seq.as_ref(), qual, enzyme, host_db, host_index, microbe_hashes);
            update_stats(is_host, is_microbe, stats);
            stats.fragments += 1;
        }

        Ok(())
    }

    fn classify_paired(
        path1: &PathBuf,
        path2: &PathBuf,
        enzyme: &f2brad_core::enzymes::Enzyme,
        host_db: &HostDb,
        host_index: &f2brad_host::genotype::HostTagIndex,
        microbe_hashes: &HashSet<Hash>,
        stats: &mut ClassifyStats,
    ) -> Result<()> {
        let mut reader1 = parse_fastx_file(path1)
            .with_context(|| format!("Failed to open reads: {}", path1.display()))?;
        let mut reader2 = parse_fastx_file(path2)
            .with_context(|| format!("Failed to open reads: {}", path2.display()))?;

        loop {
            let rec1 = reader1.next();
            let rec2 = reader2.next();
            match (rec1, rec2) {
                (None, None) => break,
                (Some(r1), Some(r2)) => {
                    let r1 = r1.with_context(|| format!("Failed to read record from {}", path1.display()))?;
                    let r2 = r2.with_context(|| format!("Failed to read record from {}", path2.display()))?;
                    let (h1, m1) = classify_fragment(r1.seq().as_ref(), r1.qual(), enzyme, host_db, host_index, microbe_hashes);
                    let (h2, m2) = classify_fragment(r2.seq().as_ref(), r2.qual(), enzyme, host_db, host_index, microbe_hashes);
                    update_stats(h1 || h2, m1 || m2, stats);
                    stats.fragments += 1;
                }
                _ => bail!(
                    "Paired input files have different numbers of reads: {} vs {}",
                    path1.display(),
                    path2.display()
                ),
            }
        }

        Ok(())
    }

    fn classify_fragment(
        seq: &[u8],
        qual: Option<&[u8]>,
        enzyme: &f2brad_core::enzymes::Enzyme,
        _host_db: &HostDb,
        host_index: &f2brad_host::genotype::HostTagIndex,
        microbe_hashes: &HashSet<Hash>,
    ) -> (bool, bool) {
        let mut is_host = false;
        let mut is_microbe = false;

        // Check all tag windows in the read; if any matches host/microbe, mark it.
        let tags = enzyme.find_all_tags(seq);
        for (offset, len) in tags {
            let tag_seq = &seq[offset..offset + len];
            let canonical = canonicalize(tag_seq);

            if host_index.find(&canonical, 2).is_some() {
                is_host = true;
            }

            // Microbial matching is exact hash lookup.
            use f2brad_core::extract::canonical_hash;
            let hash = canonical_hash(&canonical);
            if microbe_hashes.contains(&hash) {
                is_microbe = true;
            }
        }

        let _ = qual; // quality not used for classification
        (is_host, is_microbe)
    }

    fn update_stats(is_host: bool, is_microbe: bool, stats: &mut ClassifyStats) {
        match (is_host, is_microbe) {
            (true, true) => stats.both += 1,
            (true, false) => stats.host_only += 1,
            (false, true) => stats.microbe_only += 1,
            (false, false) => stats.neither += 1,
        }
    }
}
