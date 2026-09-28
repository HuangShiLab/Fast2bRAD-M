use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Args;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use tracing;

use f2brad_core::enzymes::{enzyme_by_id, enzyme_by_name};
use f2brad_core::extract::{canonical_hash, Hash};

#[derive(Args, Debug)]
pub struct DigestArgs {
    /// Reference genome FASTA/FASTQ (may be gzip-compressed)
    #[arg(short = 'i', long = "input", required = true)]
    pub input: PathBuf,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    pub enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    pub output_dir: PathBuf,

    /// Maximum Hamming distance to consider a tag non-unique (default: 2)
    #[arg(long = "max-mismatch", default_value = "2")]
    pub max_mismatch: usize,

    /// Optional BED file of CpG islands. If provided, each tag is annotated
    /// with whether it overlaps a CpG island.
    #[arg(long = "cpg-islands")]
    pub cpg_islands: Option<PathBuf>,

    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

#[derive(Debug, Clone)]
struct Locus {
    contig: String,
    pos: usize,
    strand: char,
    seq: Vec<u8>,
    canonical: Vec<u8>,
    hash: Hash,
    gc_frac: f64,
    cpg_count: usize,
    in_cpg_island: bool,
}

pub fn run(args: DigestArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();

    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("Failed to create output directory: {}", args.output_dir.display()))?;

    tracing::info!("Digesting {} with enzyme {} (tag length {})", args.input.display(), enzyme.name, enzyme.tag_length);

    let mut reader = parse_fastx_file(&args.input)
        .with_context(|| format!("Failed to open input: {}", args.input.display()))?;

    let mut loci: Vec<Locus> = Vec::new();
    let mut contig_count = 0usize;
    let mut total_bases = 0usize;

    while let Some(record) = reader.next() {
        let record = record.with_context(|| "Failed to read FASTA record")?;
        let full_header = String::from_utf8_lossy(record.id());
        let id = full_header.split_whitespace().next().unwrap_or("unnamed").to_string();
        let seq = record.seq();
        let seq_bytes = seq.as_ref();
        contig_count += 1;
        total_bases += seq_bytes.len();

        let tags = enzyme.find_all_tags(seq_bytes);
        for (offset, len) in tags {
            let tag_seq = &seq_bytes[offset..offset + len];
            let fwd = tag_seq.to_vec();
            let rc = reverse_complement(&fwd);
            let (canonical, strand) = if fwd <= rc {
                (fwd, '+')
            } else {
                (rc, '-')
            };
            let hash = canonical_hash(&canonical);
            let (gc_frac, cpg_count) = sequence_composition(&canonical);
            loci.push(Locus {
                contig: id.clone(),
                pos: offset,
                strand,
                seq: tag_seq.to_vec(),
                canonical,
                hash,
                gc_frac,
                cpg_count,
                in_cpg_island: false,
            });
        }
    }

    tracing::info!("Found {} sites in {} contigs ({} bases)", loci.len(), contig_count, total_bases);

    // Optional CpG-island annotation.
    if let Some(bed_path) = &args.cpg_islands {
        tracing::info!("Loading CpG-island BED: {}", bed_path.display());
        let islands = load_cpg_islands(bed_path)?;
        annotate_cpg_islands(&mut loci, &islands, args.threads);
        let island_count = loci.iter().filter(|l| l.in_cpg_island).count();
        tracing::info!("Tags overlapping CpG islands: {} / {} ({:.2}%)", island_count, loci.len(),
            if loci.is_empty() { 0.0 } else { island_count as f64 / loci.len() as f64 * 100.0 });
    }

    // Uniqueness analysis: a tag is usable if no other genomic tag lies within
    // `max_mismatch` Hamming distance of it. Uses the pigeonhole principle for
    // speed: split the tag into (d+1) parts and only compare tags that share an
    // exact part.
    let max_mismatch = args.max_mismatch.max(1).min(3);
    tracing::info!("Computing uniqueness at Hamming distance <= {}", max_mismatch);
    let unique = compute_unique_mask(&loci, max_mismatch, args.threads);
    let unique_count = unique.iter().filter(|&&b| b).count();
    let unique_frac = if loci.is_empty() { 0.0 } else { unique_count as f64 / loci.len() as f64 };
    tracing::info!("Unique tags: {} / {} ({:.2}%)", unique_count, loci.len(), unique_frac * 100.0);

    // Write per-locus TSV.
    let sites_path = args.output_dir.join(format!("{}.{}.sites.tsv", enzyme.name, args.max_mismatch));
    {
        let file = File::create(&sites_path)
            .with_context(|| format!("Failed to create sites file: {}", sites_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "contig\tpos\tstrand\tseq\tcanonical\thash\tgc_frac\tcpg_count\tcpg_island\tunique")?;
        for (locus, is_unique) in loci.iter().zip(unique.iter()) {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:016x}\t{:.4}\t{}\t{}\t{}",
                locus.contig,
                locus.pos,
                locus.strand,
                String::from_utf8_lossy(&locus.seq),
                String::from_utf8_lossy(&locus.canonical),
                locus.hash,
                locus.gc_frac,
                locus.cpg_count,
                if locus.in_cpg_island { 1 } else { 0 },
                if *is_unique { 1 } else { 0 }
            )?;
        }
    }

    // Write summary statistics.
    let stats_path = args.output_dir.join(format!("{}.stat.tsv", enzyme.name));
    {
        let avg_gc = if loci.is_empty() {
            0.0
        } else {
            loci.iter().map(|l| l.gc_frac).sum::<f64>() / loci.len() as f64
        };
        let avg_cpg = if loci.is_empty() {
            0.0
        } else {
            loci.iter().map(|l| l.cpg_count).sum::<usize>() as f64 / loci.len() as f64
        };
        let island_count = loci.iter().filter(|l| l.in_cpg_island).count();
        let file = File::create(&stats_path)
            .with_context(|| format!("Failed to create stats file: {}", stats_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "enzyme\ttag_length\tcontigs\ttotal_bases\tsites\tunique_sites\tunique_fraction\tmax_mismatch\tavg_gc_frac\tavg_cpg_count\tcpg_island_tags")?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.4}\t{:.4}\t{}",
            enzyme.name,
            enzyme.tag_length,
            contig_count,
            total_bases,
            loci.len(),
            unique_count,
            unique_frac,
            max_mismatch,
            avg_gc,
            avg_cpg,
            island_count
        )?;
    }

    tracing::info!("Wrote {} and {}", sites_path.display(), stats_path.display());
    Ok(())
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            x => x,
        })
        .collect()
}

fn sequence_composition(seq: &[u8]) -> (f64, usize) {
    let len = seq.len();
    if len == 0 {
        return (0.0, 0);
    }
    let mut gc = 0usize;
    let mut cpg = 0usize;
    for &b in seq {
        if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
            gc += 1;
        }
    }
    for window in seq.windows(2) {
        if (window[0] == b'C' || window[0] == b'c') && (window[1] == b'G' || window[1] == b'g') {
            cpg += 1;
        }
    }
    (gc as f64 / len as f64, cpg)
}

#[derive(Debug, Clone)]
struct BedInterval {
    contig: String,
    start: usize,
    end: usize,
}

fn load_cpg_islands(path: &PathBuf) -> Result<Vec<BedInterval>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open CpG-island BED: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let contig = parts[0].to_string();
        let start: usize = parts[1].parse().with_context(|| format!("Invalid BED start: {}", parts[1]))?;
        let end: usize = parts[2].parse().with_context(|| format!("Invalid BED end: {}", parts[2]))?;
        intervals.push(BedInterval { contig, start, end });
    }
    intervals.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));
    Ok(intervals)
}

fn annotate_cpg_islands(loci: &mut [Locus], islands: &[BedInterval], threads: usize) {
    // Build a map from contig to sorted intervals for fast lookup.
    let mut by_contig: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for island in islands {
        by_contig.entry(island.contig.clone()).or_default().push((island.start, island.end));
    }
    for intervals in by_contig.values_mut() {
        intervals.sort_by_key(|&(s, _)| s);
    }

    let _ = threads; // rayon global pool is configured by the caller
    loci.par_iter_mut().for_each(|locus| {
        let tag_start = locus.pos;
        let tag_end = locus.pos + locus.canonical.len();
        if let Some(intervals) = by_contig.get(&locus.contig) {
            locus.in_cpg_island = intervals
                .binary_search_by(|&(s, e)| {
                    if e <= tag_start {
                        std::cmp::Ordering::Less
                    } else if s >= tag_end {
                        std::cmp::Ordering::Greater
                    } else {
                        std::cmp::Ordering::Equal
                    }
                })
                .is_ok();
        }
    });
}

/// Compute a boolean mask: true if the tag at this locus has no other tag in
/// the genome within `max_mismatch` Hamming distance.
fn compute_unique_mask(loci: &[Locus], max_mismatch: usize, threads: usize) -> Vec<bool> {
    if loci.is_empty() {
        return Vec::new();
    }
    let len = loci[0].canonical.len();
    let parts = max_mismatch + 1;
    let part_size = (len + parts - 1) / parts;

    // Build pigeonhole buckets: for each part index and part sequence, store
    // the indices of loci whose canonical tag has that exact part.
    let mut buckets: Vec<HashMap<Vec<u8>, Vec<usize>>> = vec![HashMap::new(); parts];
    for (i, locus) in loci.iter().enumerate() {
        for p in 0..parts {
            let start = p * part_size;
            let end = ((p + 1) * part_size).min(len);
            if start >= len {
                continue;
            }
            let part = locus.canonical[start..end].to_vec();
            buckets[p].entry(part).or_default().push(i);
        }
    }

    // For each locus, collect candidate indices from its buckets, then check
    // full Hamming distance. Parallelise over loci.
    let unique: Vec<bool> = loci
        .par_iter()
        .enumerate()
        .map(|(i, locus)| {
            let seq = &locus.canonical;
            let mut candidates: Vec<usize> = Vec::new();
            for p in 0..parts {
                let start = p * part_size;
                let end = ((p + 1) * part_size).min(len);
                if start >= len {
                    continue;
                }
                let part = &seq[start..end];
                if let Some(idxs) = buckets[p].get(part) {
                    candidates.extend(idxs.iter().copied());
                }
            }
            candidates.sort_unstable();
            candidates.dedup();

            for &j in &candidates {
                if i == j {
                    continue;
                }
                if hamming_distance(seq, &loci[j].canonical) <= max_mismatch {
                    return false;
                }
            }
            true
        })
        .collect();

    let _ = threads; // rayon uses global thread pool configured earlier
    unique
}

fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}
