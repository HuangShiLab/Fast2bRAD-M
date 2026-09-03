use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Args;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use tracing;

use f2brad_core::enzymes::{enzyme_by_id, enzyme_by_name};
use f2brad_core::extract::Hash;

#[derive(Args, Debug)]
pub struct CrossArgs {
    /// Human tag TSV produced by `f2brad-host digest` (sites file)
    #[arg(short = 't', long = "human-tags", required = true)]
    pub human_tags: PathBuf,

    /// Microbial genome list TSV: species_name<TAB>genome_path
    #[arg(short = 'l', long = "genome-list", required = true)]
    pub genome_list: PathBuf,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    pub enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    pub output_dir: PathBuf,

    /// Maximum Hamming distance to call a collision (default: 2)
    #[arg(long = "max-mismatch", default_value = "2")]
    pub max_mismatch: usize,

    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

#[derive(Debug, Clone)]
struct HumanTag {
    canonical: Vec<u8>,
    hash: Hash,
}

#[derive(Debug, Default, Clone)]
struct CollisionRecord {
    species: String,
    distance: usize,
    human_hash: Hash,
    microbe_canonical: Vec<u8>,
}

pub fn run(args: CrossArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();

    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("Failed to create output directory: {}", args.output_dir.display()))?;

    tracing::info!("Loading human tags from {}", args.human_tags.display());
    let human_tags = load_human_tags(&args.human_tags)?;
    tracing::info!("Loaded {} human tags", human_tags.len());

    let max_mismatch = args.max_mismatch.max(1).min(3);
    tracing::info!("Building human tag index for Hamming distance <= {}", max_mismatch);
    let index = HumanTagIndex::new(&human_tags, max_mismatch);

    tracing::info!("Reading microbial genome list from {}", args.genome_list.display());
    let genomes = read_genome_list(&args.genome_list)?;
    tracing::info!("Loaded {} genome entries", genomes.len());

    // Process genomes in parallel. Each thread collects collisions for its subset.
    let results: Vec<Vec<CollisionRecord>> = genomes
        .par_iter()
        .map(|(species, path)| {
            process_microbial_genome(species, path, &index, max_mismatch, enzyme)
        })
        .collect();

    // Aggregate collisions.
    let mut species_counts: HashMap<String, [usize; 4]> = HashMap::new();
    let mut human_mask: HashSet<Hash> = HashSet::new();
    let collisions_path = args.output_dir.join(format!("collisions.{}.{}.tsv", enzyme.name, max_mismatch));
    {
        let file = File::create(&collisions_path)
            .with_context(|| format!("Failed to create collisions file: {}", collisions_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "human_hash\tspecies\tdistance\tmicrobe_canonical")?;
        for records in &results {
            for rec in records {
                human_mask.insert(rec.human_hash);
                let counts = species_counts.entry(rec.species.clone()).or_default();
                counts[rec.distance] += 1;
                writeln!(
                    writer,
                    "{:016x}\t{}\t{}\t{}",
                    rec.human_hash,
                    rec.species,
                    rec.distance,
                    String::from_utf8_lossy(&rec.microbe_canonical)
                )?;
            }
        }
    }

    // Write species summary.
    let summary_path = args.output_dir.join(format!("species_summary.{}.{}.tsv", enzyme.name, max_mismatch));
    {
        let file = File::create(&summary_path)
            .with_context(|| format!("Failed to create summary file: {}", summary_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "species\tcollisions_0\tcollisions_1\tcollisions_2\tcollisions_total")?;
        let mut species_vec: Vec<_> = species_counts.iter().collect();
        species_vec.sort_by(|a, b| b.1[0..=max_mismatch].iter().sum::<usize>().cmp(&a.1[0..=max_mismatch].iter().sum::<usize>()));
        for (species, counts) in species_vec {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                species,
                counts[0],
                counts[1],
                counts[2],
                counts[0..=max_mismatch].iter().sum::<usize>()
            )?;
        }
    }

    // Write human masking list.
    let mask_path = args.output_dir.join(format!("human_mask.{}.{}.txt", enzyme.name, max_mismatch));
    {
        let file = File::create(&mask_path)
            .with_context(|| format!("Failed to create mask file: {}", mask_path.display()))?;
        let mut writer = BufWriter::new(file);
        for hash in &human_mask {
            writeln!(writer, "{:016x}", hash)?;
        }
    }

    tracing::info!(
        "Found {} colliding human tags across {} species",
        human_mask.len(),
        species_counts.len()
    );
    tracing::info!(
        "Wrote {}, {}, {}",
        collisions_path.display(),
        summary_path.display(),
        mask_path.display()
    );

    Ok(())
}

fn load_human_tags(path: &PathBuf) -> Result<Vec<HumanTag>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open human tags: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut tags = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 7 {
            continue;
        }
        let canonical = parts[4].as_bytes().to_vec();
        let hash = u64::from_str_radix(parts[5], 16)
            .with_context(|| format!("Invalid hash on line {}: {}", i + 1, parts[5]))?;
        tags.push(HumanTag { canonical, hash });
    }
    Ok(tags)
}

fn read_genome_list(path: &PathBuf) -> Result<Vec<(String, PathBuf)>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open genome list: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut genomes = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let species = parts[0].to_string();
        let path = PathBuf::from(parts[1]);
        genomes.push((species, path));
    }
    Ok(genomes)
}

fn process_microbial_genome(
    species: &str,
    path: &PathBuf,
    index: &HumanTagIndex,
    max_mismatch: usize,
    enzyme: &f2brad_core::enzymes::Enzyme,
) -> Vec<CollisionRecord> {
    let mut collisions = Vec::new();
    let mut reader = match parse_fastx_file(path) {
        Ok(r) => r,
        Err(e) => {
            tracing::warn!("Failed to open {}: {}", path.display(), e);
            return collisions;
        }
    };

    while let Some(record) = reader.next() {
        let record = match record {
            Ok(r) => r,
            Err(e) => {
                tracing::warn!("Failed to read record from {}: {}", path.display(), e);
                continue;
            }
        };
        let seq = record.seq();
        let seq_bytes = seq.as_ref();
        let tags = enzyme.find_all_tags(seq_bytes);
        for (offset, len) in tags {
            let tag_seq = &seq_bytes[offset..offset + len];
            let canonical = canonicalize(tag_seq);
            if let Some((distance, human_hash)) = index.find_collision(&canonical, max_mismatch) {
                collisions.push(CollisionRecord {
                    species: species.to_string(),
                    distance,
                    human_hash,
                    microbe_canonical: canonical,
                });
            }
        }
    }
    collisions
}

fn canonicalize(seq: &[u8]) -> Vec<u8> {
    let fwd = seq.to_vec();
    let rc = reverse_complement(&fwd);
    if fwd <= rc { fwd } else { rc }
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

/// Index for finding human tags within a given Hamming distance of a query.
struct HumanTagIndex {
    /// tag length
    tag_len: usize,
    /// part size for pigeonhole
    part_size: usize,
    /// buckets[p][part_sequence] = list of human tag indices
    buckets: Vec<HashMap<Vec<u8>, Vec<usize>>>,
    /// canonical sequences
    canonical: Vec<Vec<u8>>,
    /// hashes
    hashes: Vec<Hash>,
}

impl HumanTagIndex {
    fn new(tags: &[HumanTag], max_mismatch: usize) -> Self {
        let tag_len = tags[0].canonical.len();
        let parts = max_mismatch + 1;
        let part_size = (tag_len + parts - 1) / parts;
        let mut buckets: Vec<HashMap<Vec<u8>, Vec<usize>>> = vec![HashMap::new(); parts];
        for (i, tag) in tags.iter().enumerate() {
            for p in 0..parts {
                let start = p * part_size;
                let end = ((p + 1) * part_size).min(tag_len);
                if start >= tag_len {
                    continue;
                }
                let part = tag.canonical[start..end].to_vec();
                buckets[p].entry(part).or_default().push(i);
            }
        }
        Self {
            tag_len,
            part_size,
            buckets,
            canonical: tags.iter().map(|t| t.canonical.clone()).collect(),
            hashes: tags.iter().map(|t| t.hash).collect(),
        }
    }

    /// Find the closest human tag within `max_mismatch` of `query`. Returns
    /// (distance, human_hash) for the first match found at the smallest
    /// distance. Exact matches are returned first.
    fn find_collision(&self, query: &[u8], max_mismatch: usize) -> Option<(usize, Hash)> {
        // Try exact match first.
        let parts = max_mismatch + 1;
        let mut candidates: Vec<usize> = Vec::new();
        for p in 0..parts {
            let start = p * self.part_size;
            let end = ((p + 1) * self.part_size).min(self.tag_len);
            if start >= self.tag_len {
                continue;
            }
            let part = &query[start..end];
            if let Some(idxs) = self.buckets[p].get(part) {
                candidates.extend(idxs.iter().copied());
            }
        }
        candidates.sort_unstable();
        candidates.dedup();

        // Check distances in ascending order to return the closest match.
        let mut best: Option<(usize, Hash)> = None;
        for &idx in &candidates {
            let dist = hamming_distance(query, &self.canonical[idx]);
            if dist <= max_mismatch {
                if best.is_none() || dist < best.unwrap().0 {
                    best = Some((dist, self.hashes[idx]));
                    if dist == 0 {
                        break;
                    }
                }
            }
        }
        best
    }
}

fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}
