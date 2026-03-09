use std::fs::File;
use std::path::{Path, PathBuf};
use std::hash::Hasher;
use std::sync::mpsc;
use std::thread;
use std::io::BufWriter;

use anyhow::{Result, bail};
use clap::Args;
// Remove flate2 imports
// use flate2::write::GzEncoder;
// use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use indicatif::ProgressBar;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use tracing;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;

pub type Hash = u64;

struct WriteTask {
    hash: Hash,
    id: String,
}

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

#[derive(Args, Debug)]
pub struct BuildQualDbArgs {
    // ── Required ──
    /// Genome list with taxonomy (TSV: GCF_id<TAB>taxonomy...<TAB>fasta_path)
    #[arg(short = 'l', long = "list", help_heading = "Required")]
    pub genome_list: PathBuf,
    /// Enzyme name (e.g. BcgI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", help_heading = "Required")]
    pub enzyme_site: String,
    /// Taxonomy level(s): kingdom|phylum|class|order|family|genus|species|strain, comma-separated or "all"
    #[arg(short = 't', long = "type", help_heading = "Required")]
    pub taxonomy_levels: String,
    /// Output directory
    #[arg(short = 'o', long = "output", help_heading = "Required")]
    pub output_dir: PathBuf,

    // ── Tag Source (mutually exclusive, both optional) ──
    /// Pre-built enzyme intermediate file. If provided, skip genome digestion. Mutually exclusive with --pre-digested-dir
    #[arg(short = 'e', long = "enzyme-file", conflicts_with = "pre_digested_dir", help_heading = "Tag Source (choose at most one)")]
    pub enzyme_file: Option<PathBuf>,
    /// Directory containing per-genome pre-digested .iibdb files. Mutually exclusive with -e
    #[arg(long = "pre-digested-dir", conflicts_with = "enzyme_file", help_heading = "Tag Source (choose at most one)")]
    pub pre_digested_dir: Option<PathBuf>,

    // ── Options ──
    /// Remove redundant (shared) tags across taxa (yes/no)
    #[arg(short = 'r', long = "remove-redundant", default_value = "yes", help_heading = "Options")]
    pub remove_redundant: String,

    // ── Performance ──
    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4", help_heading = "Performance")]
    pub threads: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TaxonomyLevel {
    Kingdom = 1, Phylum = 2, Class = 3, Order = 4, Family = 5, Genus = 6, Species = 7, Strain = 8,
}

impl TaxonomyLevel {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "kingdom" => Some(TaxonomyLevel::Kingdom),
            "phylum" => Some(TaxonomyLevel::Phylum),
            "class" => Some(TaxonomyLevel::Class),
            "order" => Some(TaxonomyLevel::Order),
            "family" => Some(TaxonomyLevel::Family),
            "genus" => Some(TaxonomyLevel::Genus),
            "species" => Some(TaxonomyLevel::Species),
            "strain" => Some(TaxonomyLevel::Strain),
            _ => None,
        }
    }
    pub fn name(&self) -> &'static str {
        match self {
            TaxonomyLevel::Kingdom => "kingdom",
            TaxonomyLevel::Phylum => "phylum",
            TaxonomyLevel::Class => "class",
            TaxonomyLevel::Order => "order",
            TaxonomyLevel::Family => "family",
            TaxonomyLevel::Genus => "genus",
            TaxonomyLevel::Species => "species",
            TaxonomyLevel::Strain => "strain",
        }
    }
    pub fn all_levels() -> Vec<Self> {
        vec![TaxonomyLevel::Kingdom, TaxonomyLevel::Phylum, TaxonomyLevel::Class, TaxonomyLevel::Order, TaxonomyLevel::Family, TaxonomyLevel::Genus, TaxonomyLevel::Species, TaxonomyLevel::Strain]
    }
}

#[derive(Debug, Clone)]
pub struct GenomeRecord {
    pub gcf_id: String,
    pub taxonomy: Vec<String>,
    pub genome_path: Option<PathBuf>,
}

pub fn run(args: BuildQualDbArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;
    let remove_redundant = args.remove_redundant.eq_ignore_ascii_case("yes");

    io_utils::ensure_directory(&args.output_dir)?;
    tracing::info!("Reading genome taxonomy list ...");
    let genomes = read_genome_list(&args.genome_list)?;
    tracing::info!("Total {} genomes", genomes.len());

    let enzyme_file = if let Some(ref file) = args.enzyme_file {
        tracing::info!("Using pre-digested file: {}", file.display());
        file.clone()
    } else if let Some(ref dir) = args.pre_digested_dir {
        tracing::info!("Merging pre-digested files from directory in parallel: {}", dir.display());
        let output_file = args.output_dir.join(format!("{}.enzyme.iibdb", enzyme.name));
        merge_pre_digested_files(&genomes, enzyme, dir, &output_file)?;
        output_file
    } else {
        tracing::info!("Digesting genomes and generating hashes (Parallel Binary) ...");
        let output_file = args.output_dir.join(format!("{}.enzyme.iibdb", enzyme.name));
        digest_genomes(&genomes, enzyme, &output_file)?;
        output_file
    };

    // [Optimization] Load enzyme file into memory once, reuse for all levels
    tracing::info!("Loading enzyme records into memory ...");
    let enzyme_records = load_enzyme_records(&enzyme_file)?;
    tracing::info!("Loaded {} enzyme records", enzyme_records.len());

    for level in &levels {
        tracing::info!("\n========== Building {}-level database (Hash mode) ==========", level.name());
        build_database_for_level(&enzyme_records, enzyme, &args.output_dir, *level, &genomes, remove_redundant)?;
    }
    tracing::info!("\nAll done!");
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
    if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme); }
    if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme); } }
    bail!("Unknown enzyme: {}", site)
}

fn parse_taxonomy_levels(levels_str: &str) -> Result<Vec<TaxonomyLevel>> {
    if levels_str.eq_ignore_ascii_case("all") { return Ok(TaxonomyLevel::all_levels()); }
    let mut levels = Vec::new();
    for part in levels_str.split(',') {
        let level = TaxonomyLevel::from_str(part.trim()).ok_or_else(|| anyhow::anyhow!("Invalid taxonomy level: {}", part))?;
        levels.push(level);
    }
    if levels.is_empty() { bail!("At least one taxonomy level must be specified"); }
    Ok(levels)
}

fn read_genome_list(path: &Path) -> Result<Vec<GenomeRecord>> {
    let file = File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut genomes = Vec::new();
    let mut is_gtdb_format = false;
    let mut first_data_line = true;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if first_data_line {
            first_data_line = false;
            if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                is_gtdb_format = true;
                if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; }
            }
        }
        if is_gtdb_format {
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord { gcf_id: extract_gcf_id(parts[0]), taxonomy: parse_gtdb_taxonomy(parts[1])?, genome_path: None });
        } else {
            if parts.len() < 9 { continue; }
            let genome_path = if parts.len() > 9 { Some(PathBuf::from(parts[9])) } else { None };
            genomes.push(GenomeRecord { gcf_id: parts[0].to_string(), taxonomy: parts[1..9].iter().map(|s| s.to_string()).collect(), genome_path });
        }
    }
    Ok(genomes)
}

fn extract_gcf_id(filename: &str) -> String {
    let name = filename.split('/').last().unwrap_or(filename);
    let name_clean = if let Some(pos) = name.find("_genomic") { &name[..pos] } else { name };
    if name_clean.starts_with("GCF_") || name_clean.starts_with("GCA_") {
        let parts: Vec<&str> = name_clean.split('_').collect();
        if parts.len() >= 2 { return format!("{}_{}", parts[0], parts[1]); }
    }
    name_clean.to_string()
}

fn parse_gtdb_taxonomy(gtdb_str: &str) -> Result<Vec<String>> {
    let parts: Vec<&str> = gtdb_str.split(';').collect();
    let mut taxonomy = Vec::new();
    for part in parts.iter() {
        if let Some(pos) = part.find("__") { taxonomy.push(part[pos+2..].to_string()); } else { taxonomy.push(part.to_string()); }
    }
    while taxonomy.len() < 8 {
        if let Some(last) = taxonomy.last() { taxonomy.push(format!("{}_strain", last)); } else { taxonomy.push("unknown".to_string()); }
    }
    Ok(taxonomy)
}

fn merge_pre_digested_files(genomes: &[GenomeRecord], enzyme: &'static Enzyme, pre_digested_dir: &Path, output_file: &Path) -> Result<()> {
    let (sender, receiver) = mpsc::channel::<Vec<WriteTask>>();
    let output_file_buf = output_file.to_path_buf();

    let writer_thread = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_file_buf)?;
        let buf_writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        let mut writer = zstd::Encoder::new(buf_writer, 3)?;
        while let Ok(tasks) = receiver.recv() {
            for task in tasks {
                io_utils::write_binary_record(&mut writer, task.hash, &task.id)?;
            }
        }
        writer.finish()?;
        Ok(())
    });

    let pb = ProgressBar::new(genomes.len() as u64);
    genomes.par_iter().for_each_with(sender, |s, genome| {
        let genome_id = genome.gcf_id.split('.').take(2).collect::<Vec<_>>().join(".");
        // Retain compatibility: check for .iibdb first
        let patterns = [format!("{}.{}.iibdb", genome_id, enzyme.name), format!("{}.{}.iibdb", genome.gcf_id, enzyme.name), format!("{}.{}.iibsp", genome_id, enzyme.name)];
        for pattern in &patterns {
            let test_path = pre_digested_dir.join(pattern);
            if test_path.exists() {
                if let Ok(tasks) = read_binary_content(&test_path, &genome.gcf_id) { s.send(tasks).expect("Failed to send tasks"); }
                break;
            }
        }
        pb.inc(1);
    });
    pb.finish();
    writer_thread.join().unwrap()?;
    Ok(())
}

fn read_binary_content(file_path: &Path, gcf_id: &str) -> Result<Vec<WriteTask>> {
    let mut tasks = Vec::new();
    let mut reader = io_utils::open_binary_reader(file_path)?;
    let mut id_buffer = String::with_capacity(128);
    while let Some(hash) = reader.next_record_reuse(&mut id_buffer)? {
        let (scaffold, pos) = if let Some((s, p)) = id_buffer.rsplit_once(':') { (s, p) } else { (id_buffer.as_str(), "0") };
        let new_id = format!("{}|0|{}|{}|0|-", gcf_id, scaffold, pos);
        tasks.push(WriteTask { hash, id: new_id });
    }
    Ok(tasks)
}

fn digest_genomes(genomes: &[GenomeRecord], enzyme: &'static Enzyme, output_file: &Path) -> Result<()> {
    let (sender, receiver) = mpsc::channel::<Vec<WriteTask>>();
    let output_file_buf = output_file.to_path_buf();

    let writer_thread = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_file_buf)?;
        let buf_writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        let mut writer = zstd::Encoder::new(buf_writer, 3)?;
        while let Ok(tasks) = receiver.recv() {
            for task in tasks {
                io_utils::write_binary_record(&mut writer, task.hash, &task.id)?;
            }
        }
        writer.finish()?;
        Ok(())
    });

    let pb = ProgressBar::new(genomes.len() as u64);
    genomes.par_iter().for_each_with(sender, |s, genome| {
        if let Some(genome_path) = &genome.genome_path {
            if genome_path.exists() {
                let mut tasks = Vec::new();
                let mut tag_index = 0usize;
                if let Ok(mut reader) = parse_fastx_file(genome_path) {
                    while let Some(record) = reader.next() {
                        if let Ok(record) = record {
                             let seq_id = std::str::from_utf8(record.id()).unwrap_or("seq").split_whitespace().next().unwrap_or("seq").to_string();
                            let mut sequence = record.seq().to_vec();
                            sequence.make_ascii_uppercase();
                            let positions = enzyme.find_all_tags(&sequence);
                            for (pos, len) in positions {
                                tag_index += 1;
                                let tag_seq = &sequence[pos..pos + len];
                                let hash_val = canonical_hash(tag_seq);
                                let id_str = format!("{}|{}|{}|{}|0|-", genome.gcf_id, tag_index, seq_id, pos + 1);
                                tasks.push(WriteTask { hash: hash_val, id: id_str });
                            }
                        }
                    }
                }
                if !tasks.is_empty() { s.send(tasks).expect("Failed to send digest tasks"); }
            }
        }
        pb.inc(1);
    });
    pb.finish();
    writer_thread.join().unwrap()?;
    Ok(())
}

// ================== Single-pass multi-level optimization ==================
// Previously the enzyme file was re-read from disk 2× per taxonomy level
// (once for tag taxonomy collection, once for output). With 8 levels, that's
// 16 disk reads of a potentially multi-GB file.
// Now we load the enzyme records into memory once and reuse them for all levels.

/// In-memory enzyme record: only stores the hash and GCF ID we actually need.
struct EnzymeRecord {
    hash: Hash,
    gcf_id: String,
}

fn load_enzyme_records(enzyme_file: &Path) -> Result<Vec<EnzymeRecord>> {
    let mut reader = io_utils::open_binary_reader(enzyme_file)?;
    let mut records = Vec::new();
    let mut id_buffer = String::with_capacity(256);
    while let Some(hash) = reader.next_record_reuse(&mut id_buffer)? {
        let gcf_id = id_buffer.split('|').next().unwrap_or("");
        if !gcf_id.is_empty() {
            records.push(EnzymeRecord { hash, gcf_id: gcf_id.to_string() });
        }
    }
    Ok(records)
}

fn build_database_for_level(
    records: &[EnzymeRecord],
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<()> {
    // Build taxonomy lookup for this level
    let mut gcf_to_taxonomy: FxHashMap<&str, String> = FxHashMap::default();
    for genome in genomes {
        if level as usize > genome.taxonomy.len() { bail!("Taxonomy level index out of range"); }
        let taxonomy_str = genome.taxonomy[0..level as usize].join("\t");
        gcf_to_taxonomy.insert(&genome.gcf_id, taxonomy_str);
    }

    // Pass 1 (in-memory): collect tag taxonomy and genome tag counts
    tracing::info!("Step 1: Collecting tag taxonomy information ...");
    let mut tag_taxonomy: FxHashMap<Hash, FxHashSet<String>> = FxHashMap::default();
    let mut genome_tags: FxHashMap<&str, FxHashMap<Hash, usize>> = FxHashMap::default();

    for rec in records {
        if let Some(taxonomy) = gcf_to_taxonomy.get(rec.gcf_id.as_str()) {
            tag_taxonomy.entry(rec.hash).or_default().insert(taxonomy.clone());
            if remove_redundant {
                *genome_tags.entry(rec.gcf_id.as_str()).or_default().entry(rec.hash).or_insert(0) += 1;
            }
        }
    }

    // Pass 2 (in-memory): output unique tags with compact + zstd format
    tracing::info!("Step 2: Identifying unique tags and writing database ...");
    let output_path = output_dir.join(format!("{}.{}.iibdb", enzyme.name, level.name()));
    let file = File::create(&output_path)?;
    let buf_writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);

    let gcf_ids: Vec<&str> = genomes.iter().map(|g| g.gcf_id.as_str()).collect();
    let gcf_to_index: FxHashMap<&str, u32> = gcf_ids.iter().enumerate()
        .map(|(i, id)| (*id, i as u32)).collect();
    let mut compact_writer = io_utils::CompactDatabaseWriter::new(buf_writer, &gcf_ids)?;
    let mut unique_counts: FxHashMap<&str, usize> = FxHashMap::default();

    for rec in records {
        if !gcf_to_taxonomy.contains_key(rec.gcf_id.as_str()) { continue; }

        let is_unique = tag_taxonomy.get(&rec.hash).map_or(false, |t| t.len() == 1);

        let is_unique = if is_unique && remove_redundant {
            genome_tags.get(rec.gcf_id.as_str())
                .and_then(|counts| counts.get(&rec.hash))
                .map_or(true, |&count| count <= 1)
        } else { is_unique };

        if is_unique {
            if let Some(&idx) = gcf_to_index.get(rec.gcf_id.as_str()) {
                *unique_counts.entry(rec.gcf_id.as_str()).or_insert(0) += 1;
                compact_writer.write_record(rec.hash, idx)?;
            }
        }
    }
    compact_writer.finish()?;

    tracing::info!("  Output database: {}", output_path.display());
    tracing::info!("  Contains unique tags for {} genomes", unique_counts.len());
    Ok(())
}
