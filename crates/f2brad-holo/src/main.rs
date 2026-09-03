use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};
use clap::{Parser, Subcommand};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use tracing;

use f2brad_core::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use f2brad_core::extract::Hash;
use f2brad_core::io_utils::{open_compact_reader, write_sample_tag_header};
use f2brad_host::genotype::{
    add_match_to_pileup, canonicalize, extract_matches, load_host_db, write_bimbam, write_vcf,
    HostDb, Pileup, ReadMatch,
};

#[derive(Parser, Debug)]
#[command(name = "f2brad-holo", version, about = "One-pass holo driver for fast2bRAD")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// One-pass classification + host genotyping + microbial profiling.
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

    /// Optional microbial database directory. If supplied, microbial tag hashes
    /// are mapped to GCF/taxon abundances and written to `microbe_counts.tsv`.
    /// Expected contents: `{enzyme}.{level}.iibdb` and
    /// `abfh_classify_with_speciename.txt.gz`.
    #[arg(long = "microbe-db-dir")]
    microbe_db_dir: Option<PathBuf>,

    /// Optional microbial mask file: one canonical tag sequence per line.
    /// Tags in this list are removed from the microbial hash set before matching.
    #[arg(long = "microbe-mask")]
    microbe_mask: Option<PathBuf>,

    /// Read 1 FASTQ/FASTA (may be gzip-compressed). Required unless --sample-list is given.
    #[arg(short = '1', long = "r1")]
    r1: Option<PathBuf>,

    /// Optional read 2 FASTQ/FASTA (paired-end)
    #[arg(short = '2', long = "r2")]
    r2: Option<PathBuf>,

    /// Sample list TSV: sample_name<TAB>r1_path[<TAB>r2_path].
    /// When provided, --r1/--r2/--sample-name are ignored and all samples are
    /// processed in parallel subdirectories of --output.
    #[arg(short = 'l', long = "sample-list")]
    sample_list: Option<PathBuf>,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    output_dir: PathBuf,

    /// Sample name used in output files. Required unless --sample-list is given.
    #[arg(long = "sample-name", default_value = "sample")]
    sample_name: String,

    /// Maximum Hamming distance for host tag matching
    #[arg(long = "host-max-mismatch", default_value = "2")]
    host_max_mismatch: usize,

    /// Minimum base quality (Phred) for genotype pileup
    #[arg(short = 'q', long = "min-qual", default_value = "20")]
    min_qual: u8,

    /// Minimum per-locus depth to emit a genotype record
    #[arg(long = "min-depth", default_value = "4")]
    min_depth: usize,

    /// Taxonomy level for microbial counts (kingdom..strain)
    #[arg(short = 't', long = "taxonomy", default_value = "species")]
    taxonomy_level: String,

    /// Write a sample tag stream (`sample.iibsp.gz`) for downstream `f2brad-m quantify`
    #[arg(long = "output-iibsp", default_value = "true")]
    output_iibsp: bool,

    /// Number of parallel threads
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

        let max_mismatch = args.host_max_mismatch.max(1).min(3);

        tracing::info!("Loading host DB from {}", args.host_db.display());
        let host_db = Arc::new(load_host_db(&args.host_db, max_mismatch)?);
        tracing::info!("Loaded {} host tags", host_db.loci.len());

        tracing::info!("Loading microbial DB from {}", args.microbe_db.display());
        let (microbe_hashes, hash_to_gcfs) = load_microbe_db(&args.microbe_db)?;
        tracing::info!("Loaded {} microbial tag hashes", microbe_hashes.len());

        let microbe_hashes = if let Some(mask_path) = &args.microbe_mask {
            tracing::info!("Loading microbial mask from {}", mask_path.display());
            let mask = load_microbe_mask(mask_path)?;
            tracing::info!("Applying microbial mask: {} tag(s)", mask.len());
            microbe_hashes.difference(&mask).copied().collect()
        } else {
            microbe_hashes
        };
        tracing::info!("{} microbial hashes remain after masking", microbe_hashes.len());
        let microbe_hashes = Arc::new(microbe_hashes);
        let hash_to_gcfs = Arc::new(hash_to_gcfs);

        let tax_level_idx = taxonomy_level_index(&args.taxonomy_level)?;
        let taxonomy: Option<HashMap<String, Vec<String>>> = if let Some(db_dir) = &args.microbe_db_dir {
            let tax_path = db_dir.join("abfh_classify_with_speciename.txt.gz");
            if tax_path.exists() {
                Some(load_taxonomy(&tax_path)?)
            } else {
                tracing::warn!("Taxonomy file not found at {}; taxon counts disabled", tax_path.display());
                None
            }
        } else {
            None
        };
        let taxonomy = Arc::new(taxonomy);

        let samples = if let Some(list_path) = &args.sample_list {
            parse_sample_list(list_path)?
        } else {
            let r1 = args.r1.as_ref()
                .ok_or_else(|| anyhow!("--r1 is required when --sample-list is not provided"))?;
            vec![SampleEntry {
                name: args.sample_name.clone(),
                r1: r1.clone(),
                r2: args.r2.clone(),
            }]
        };

        if samples.is_empty() {
            bail!("No samples to process");
        }

        tracing::info!("Processing {} sample(s) in parallel", samples.len());

        let min_qual = args.min_qual;
        let min_depth = args.min_depth;
        let taxonomy_level = args.taxonomy_level.clone();
        let output_iibsp = args.output_iibsp;
        let output_dir = args.output_dir.clone();

        let results: Vec<Result<()>> = samples
            .par_iter()
            .map(|sample| {
                let sample_output = output_dir.join(&sample.name);
                std::fs::create_dir_all(&sample_output)
                    .with_context(|| format!("Failed to create sample directory: {}", sample_output.display()))?;
                process_one_sample(
                    &sample.name,
                    &sample.r1,
                    sample.r2.as_ref(),
                    enzyme,
                    &host_db,
                    &microbe_hashes,
                    &hash_to_gcfs,
                    taxonomy.as_ref().as_ref(),
                    tax_level_idx,
                    max_mismatch,
                    min_qual,
                    min_depth,
                    &taxonomy_level,
                    output_iibsp,
                    &sample_output,
                )
            })
            .collect();

        for (sample, result) in samples.iter().zip(results.iter()) {
            match result {
                Ok(_) => tracing::info!("Sample {} completed", sample.name),
                Err(e) => tracing::error!("Sample {} failed: {}", sample.name, e),
            }
        }

        // Return first error, if any.
        for result in results {
            result?;
        }

        Ok(())
    }

    struct SampleEntry {
        name: String,
        r1: PathBuf,
        r2: Option<PathBuf>,
    }

    fn parse_sample_list(path: &PathBuf) -> Result<Vec<SampleEntry>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open sample list: {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut samples = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 2 {
                bail!("Invalid sample list line {}: expected at least 2 columns", i + 1);
            }
            samples.push(SampleEntry {
                name: parts[0].to_string(),
                r1: PathBuf::from(parts[1]),
                r2: parts.get(2).map(|s| PathBuf::from(s)),
            });
        }
        Ok(samples)
    }

    fn process_one_sample(
        sample_name: &str,
        r1: &PathBuf,
        r2: Option<&PathBuf>,
        enzyme: &Enzyme,
        host_db: &HostDb,
        microbe_hashes: &HashSet<Hash>,
        hash_to_gcfs: &HashMap<Hash, Vec<String>>,
        taxonomy: Option<&HashMap<String, Vec<String>>>,
        tax_level_idx: usize,
        max_mismatch: usize,
        min_qual: u8,
        min_depth: usize,
        taxonomy_level: &str,
        output_iibsp: bool,
        output_dir: &PathBuf,
    ) -> Result<()> {
        let mut stats = ClassifyStats::default();
        let mut pileup = Pileup::new(host_db.loci.len(), host_db.loci[0].canonical.len());
        let mut gcf_counts: HashMap<String, usize> = HashMap::new();
        let mut taxon_counts: HashMap<String, usize> = HashMap::new();

        let iibsp_path = output_dir.join(format!("{}.iibsp.gz", sample_name));
        let mut iibsp_sink: std::io::Sink = std::io::sink();
        let mut iibsp_encoder: Option<GzEncoder<BufWriter<File>>> = if output_iibsp {
            let iibsp_file = File::create(&iibsp_path)
                .with_context(|| format!("Failed to create iibsp output: {}", iibsp_path.display()))?;
            let iibsp_buf = BufWriter::with_capacity(128 * 1024, iibsp_file);
            let mut enc = GzEncoder::new(iibsp_buf, Compression::default());
            write_sample_tag_header(&mut enc)?;
            Some(enc)
        } else {
            None
        };

        let iibsp_writer: &mut dyn std::io::Write = if let Some(ref mut enc) = iibsp_encoder {
            enc
        } else {
            &mut iibsp_sink
        };

        if let Some(r2) = r2 {
            classify_paired(
                r1,
                r2,
                enzyme,
                host_db,
                microbe_hashes,
                hash_to_gcfs,
                taxonomy,
                tax_level_idx,
                max_mismatch,
                min_qual,
                &mut stats,
                &mut pileup,
                &mut gcf_counts,
                &mut taxon_counts,
                iibsp_writer,
            )?;
        } else {
            classify_single(
                r1,
                enzyme,
                host_db,
                microbe_hashes,
                hash_to_gcfs,
                taxonomy,
                tax_level_idx,
                max_mismatch,
                min_qual,
                &mut stats,
                &mut pileup,
                &mut gcf_counts,
                &mut taxon_counts,
                iibsp_writer,
            )?;
        }

        if let Some(enc) = iibsp_encoder {
            enc.finish()?;
            if output_iibsp {
                tracing::info!("[{}] Wrote {}", sample_name, iibsp_path.display());
            }
        }

        // Host genotype outputs.
        let vcf_path = output_dir.join("genotypes.vcf");
        write_vcf(&vcf_path, &host_db.loci, &pileup, min_depth)?;
        tracing::info!("[{}] Wrote {}", sample_name, vcf_path.display());

        let bimbam_path = output_dir.join("dosages.bimbam");
        write_bimbam(&bimbam_path, &host_db.loci, &pileup, min_depth)?;
        tracing::info!("[{}] Wrote {}", sample_name, bimbam_path.display());

        // Microbial count outputs.
        if !gcf_counts.is_empty() {
            let gcf_path = output_dir.join("gcf_counts.tsv");
            write_count_table(&gcf_path, "gcf", &gcf_counts)?;
        }
        if !taxon_counts.is_empty() {
            let tax_path = output_dir.join(format!("{}_counts.tsv", taxonomy_level));
            write_count_table(&tax_path, taxonomy_level, &taxon_counts)?;
        }

        // Classification summary.
        let total_classified = stats.host_only + stats.microbe_only + stats.both;
        let host_fraction = if total_classified > 0 {
            (stats.host_only + stats.both) as f64 / total_classified as f64
        } else {
            0.0
        };

        let out_path = output_dir.join("holo_classify.tsv");
        let file = File::create(&out_path)
            .with_context(|| format!("Failed to create output: {}", out_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "metric\tvalue")?;
        writeln!(writer, "sample\t{}", sample_name)?;
        writeln!(writer, "input_fragments\t{}", stats.fragments)?;
        writeln!(writer, "host_only\t{}", stats.host_only)?;
        writeln!(writer, "microbe_only\t{}", stats.microbe_only)?;
        writeln!(writer, "both\t{}", stats.both)?;
        writeln!(writer, "neither\t{}", stats.neither)?;
        writeln!(writer, "host_fraction\t{:.6}", host_fraction)?;
        writeln!(writer, "microbial_hashes_observed\t{}", stats.microbe_hashes_observed)?;
        writer.flush()?;
        tracing::info!("[{}] Wrote {}", sample_name, out_path.display());

        Ok(())
    }

    #[derive(Default)]
    struct ClassifyStats {
        fragments: usize,
        host_only: usize,
        microbe_only: usize,
        both: usize,
        neither: usize,
        microbe_hashes_observed: usize,
    }

    /// Load a microbial .iibdb compact database. Returns the set of all tag hashes
    /// and a map from hash to the GCF id(s) carrying that tag.
    fn load_microbe_db(path: &PathBuf) -> Result<(HashSet<Hash>, HashMap<Hash, Vec<String>>)> {
        let mut reader = open_compact_reader(path)?;
        let gcf_table: Vec<String> = reader.gcf_table().to_vec();
        let mut hashes = HashSet::new();
        let mut hash_to_gcfs: HashMap<Hash, Vec<String>> = HashMap::new();
        while let Some((hash, gcf_index)) = reader.next_record()? {
            hashes.insert(hash);
            let gcf_id = gcf_table
                .get(gcf_index as usize)
                .cloned()
                .unwrap_or_else(|| format!("gcf_index_{}", gcf_index));
            hash_to_gcfs.entry(hash).or_default().push(gcf_id);
        }
        Ok((hashes, hash_to_gcfs))
    }

    /// Load a microbial mask file: one canonical sequence per line. Convert each
    /// sequence to the canonical hash used by the microbial database.
    fn load_microbe_mask(path: &PathBuf) -> Result<HashSet<Hash>> {
        use f2brad_core::extract::canonical_hash;
        let file = File::open(path)
            .with_context(|| format!("Failed to open microbe mask: {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut mask = HashSet::new();
        for line in reader.lines() {
            let line = line?;
            let seq = line.trim().as_bytes();
            if seq.is_empty() {
                continue;
            }
            let canonical = canonicalize(seq);
            let hash = canonical_hash(&canonical);
            mask.insert(hash);
        }
        Ok(mask)
    }

    fn taxonomy_level_index(level: &str) -> Result<usize> {
        let valid = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"];
        valid.iter().position(|&x| x.eq_ignore_ascii_case(level))
            .ok_or_else(|| anyhow!("Invalid taxonomy level: {}", level))
    }

    /// Load GCF -> taxonomy mapping from the compressed classify file.
    fn load_taxonomy(path: &PathBuf) -> Result<HashMap<String, Vec<String>>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open taxonomy file: {}", path.display()))?;
        let reader = BufReader::new(GzDecoder::new(file));
        let mut taxonomy = HashMap::new();
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() < 2 {
                continue;
            }
            let gcf = parts[0].to_string();
            let ranks: Vec<String> = parts[1..].iter().map(|s| s.to_string()).collect();
            taxonomy.insert(gcf, ranks);
        }
        Ok(taxonomy)
    }

    fn classify_single(
        path: &PathBuf,
        enzyme: &f2brad_core::enzymes::Enzyme,
        host_db: &HostDb,
        microbe_hashes: &HashSet<Hash>,
        hash_to_gcfs: &HashMap<Hash, Vec<String>>,
        taxonomy: Option<&HashMap<String, Vec<String>>>,
        tax_level_idx: usize,
        max_mismatch: usize,
        min_qual: u8,
        stats: &mut ClassifyStats,
        pileup: &mut Pileup,
        gcf_counts: &mut HashMap<String, usize>,
        taxon_counts: &mut HashMap<String, usize>,
        iibsp: &mut dyn Write,
    ) -> Result<()> {
        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open reads: {}", path.display()))?;

        let mut fragment_hashes: HashSet<Hash> = HashSet::new();
        while let Some(record) = reader.next() {
            let record = record.with_context(|| format!("Failed to read record from {}", path.display()))?;
            let seq = record.seq();
            let qual = record.qual();
            fragment_hashes.clear();
            let fragment_tags = classify_fragment(
                seq.as_ref(), qual, enzyme, host_db, microbe_hashes, hash_to_gcfs,
                taxonomy, tax_level_idx, max_mismatch, min_qual, pileup, gcf_counts,
                taxon_counts, &mut fragment_hashes,
            )?;
            for hash in &fragment_hashes {
                iibsp.write_all(&hash.to_le_bytes())?;
            }
            stats.microbe_hashes_observed += fragment_hashes.len();
            update_stats(fragment_tags.is_host, fragment_tags.is_microbe, stats);
            stats.fragments += 1;
        }

        Ok(())
    }

    fn classify_paired(
        path1: &PathBuf,
        path2: &PathBuf,
        enzyme: &f2brad_core::enzymes::Enzyme,
        host_db: &HostDb,
        microbe_hashes: &HashSet<Hash>,
        hash_to_gcfs: &HashMap<Hash, Vec<String>>,
        taxonomy: Option<&HashMap<String, Vec<String>>>,
        tax_level_idx: usize,
        max_mismatch: usize,
        min_qual: u8,
        stats: &mut ClassifyStats,
        pileup: &mut Pileup,
        gcf_counts: &mut HashMap<String, usize>,
        taxon_counts: &mut HashMap<String, usize>,
        iibsp: &mut dyn Write,
    ) -> Result<()> {
        let mut reader1 = parse_fastx_file(path1)
            .with_context(|| format!("Failed to open reads: {}", path1.display()))?;
        let mut reader2 = parse_fastx_file(path2)
            .with_context(|| format!("Failed to open reads: {}", path2.display()))?;

        let mut n_pairs = 0usize;
        let mut fragment_hashes: HashSet<Hash> = HashSet::new();
        loop {
            let rec1 = reader1.next();
            let rec2 = reader2.next();
            match (rec1, rec2) {
                (None, None) => break,
                (Some(r1), Some(r2)) => {
                    let r1 = r1.with_context(|| format!("Failed to read record from {}", path1.display()))?;
                    let r2 = r2.with_context(|| format!("Failed to read record from {}", path2.display()))?;

                    fragment_hashes.clear();
                    let t1 = classify_fragment(
                        r1.seq().as_ref(), r1.qual(), enzyme, host_db, microbe_hashes,
                        hash_to_gcfs, taxonomy, tax_level_idx, max_mismatch, min_qual, pileup,
                        gcf_counts, taxon_counts, &mut fragment_hashes,
                    )?;
                    let t2 = classify_fragment(
                        r2.seq().as_ref(), r2.qual(), enzyme, host_db, microbe_hashes,
                        hash_to_gcfs, taxonomy, tax_level_idx, max_mismatch, min_qual, pileup,
                        gcf_counts, taxon_counts, &mut fragment_hashes,
                    )?;
                    for hash in &fragment_hashes {
                        iibsp.write_all(&hash.to_le_bytes())?;
                    }
                    stats.microbe_hashes_observed += fragment_hashes.len();

                    update_stats(t1.is_host || t2.is_host, t1.is_microbe || t2.is_microbe, stats);
                    stats.fragments += 1;
                    n_pairs += 1;
                }
                _ => bail!(
                    "Paired input files have different numbers of reads: {} vs {}",
                    path1.display(),
                    path2.display()
                ),
            }
        }

        tracing::info!("Processed {} paired-end fragments", n_pairs);
        Ok(())
    }

    struct FragmentTags {
        is_host: bool,
        is_microbe: bool,
    }

    fn classify_fragment(
        seq: &[u8],
        qual: Option<&[u8]>,
        enzyme: &f2brad_core::enzymes::Enzyme,
        host_db: &HostDb,
        microbe_hashes: &HashSet<Hash>,
        hash_to_gcfs: &HashMap<Hash, Vec<String>>,
        taxonomy: Option<&HashMap<String, Vec<String>>>,
        tax_level_idx: usize,
        max_mismatch: usize,
        min_qual: u8,
        pileup: &mut Pileup,
        gcf_counts: &mut HashMap<String, usize>,
        taxon_counts: &mut HashMap<String, usize>,
        fragment_hashes: &mut HashSet<Hash>,
    ) -> Result<FragmentTags> {
        use f2brad_core::extract::canonical_hash;

        let mut is_host = false;
        let mut is_microbe = false;

        // Host matches: use genotype.rs extraction, then add best match per locus to pileup.
        let host_matches = extract_matches(seq, qual, enzyme, host_db, max_mismatch);
        if !host_matches.is_empty() {
            is_host = true;
        }
        let mut best_host: HashMap<usize, ReadMatch> = HashMap::new();
        for m in host_matches {
            best_host
                .entry(m.locus_idx)
                .and_modify(|existing| {
                    if m.dist < existing.dist {
                        *existing = m.clone();
                    }
                })
                .or_insert(m);
        }
        for m in best_host.values() {
            add_match_to_pileup(m, min_qual, host_db, pileup);
        }

        // Microbial matches: collect observed hashes, deduplicate per fragment, count once.
        let tags = enzyme.find_all_tags(seq);
        for (offset, len) in tags {
            let tag_seq = &seq[offset..offset + len];
            let canonical = canonicalize(tag_seq);
            let hash = canonical_hash(&canonical);
            if microbe_hashes.contains(&hash) && fragment_hashes.insert(hash) {
                is_microbe = true;
                if let Some(gcfs) = hash_to_gcfs.get(&hash) {
                    for gcf in gcfs {
                        *gcf_counts.entry(gcf.clone()).or_insert(0) += 1;
                        if let Some(tax) = taxonomy {
                            if let Some(ranks) = tax.get(gcf) {
                                if let Some(rank_name) = ranks.get(tax_level_idx) {
                                    *taxon_counts.entry(rank_name.clone()).or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(FragmentTags { is_host, is_microbe })
    }

    fn update_stats(is_host: bool, is_microbe: bool, stats: &mut ClassifyStats) {
        match (is_host, is_microbe) {
            (true, true) => stats.both += 1,
            (true, false) => stats.host_only += 1,
            (false, true) => stats.microbe_only += 1,
            (false, false) => stats.neither += 1,
        }
    }

    fn write_count_table(path: &PathBuf, label: &str, counts: &HashMap<String, usize>) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create count table: {}", path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{}\tcount", label)?;
        let mut items: Vec<_> = counts.iter().collect();
        items.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        for (name, count) in items {
            writeln!(writer, "{}\t{}", name, count)?;
        }
        writer.flush()?;
        tracing::info!("Wrote {}", path.display());
        Ok(())
    }
}
