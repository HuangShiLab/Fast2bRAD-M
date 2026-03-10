use anyhow::{bail, Context, Result, anyhow};
use clap::Parser;
// use flate2::write::GzEncoder; // Removed
// use flate2::Compression;      // Removed
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use indicatif::ProgressBar;
use needletail::parse_fastx_file;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::thread;
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

#[derive(Parser, Debug)]
pub struct BuildQuanDbArgs {
    // ── Required ──
    /// Genome list for this sample (TSV, typically sdb.list from find-genome)
    #[arg(short = 'l', long = "list", help_heading = "Required")]
    pub genome_list: PathBuf,
    /// Enzyme name (e.g. BcgI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", help_heading = "Required")]
    pub enzyme_site: String,
    /// Taxonomy level(s): kingdom|phylum|class|order|family|genus|species|strain, comma-separated or "all"
    #[arg(short = 't', long = "taxonomy", help_heading = "Required")]
    pub taxonomy_levels: String,
    /// Output directory
    #[arg(short = 'o', long = "output", help_heading = "Required")]
    pub output_dir: PathBuf,

    // ── Tag Source (provide exactly one) ──
    /// Pre-built enzyme intermediate file (typically from build-qual-db). Mutually exclusive with --pre-digested-dir
    #[arg(short = 'e', long = "enzyme-file", conflicts_with = "pre_digested_dir", help_heading = "Tag Source (provide one)")]
    pub enzyme_file: Option<PathBuf>,
    /// Directory containing per-genome pre-digested .iibdb files. Mutually exclusive with -e
    #[arg(long = "pre-digested-dir", conflicts_with = "enzyme_file", help_heading = "Tag Source (provide one)")]
    pub pre_digested_dir: Option<PathBuf>,

    // ── Options ──
    /// Remove redundant (shared) tags across taxa (yes/no)
    #[arg(short = 'r', long = "remove-redundant", default_value = "no", help_heading = "Options")]
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
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "kingdom" => Ok(TaxonomyLevel::Kingdom), "phylum" => Ok(TaxonomyLevel::Phylum), "class" => Ok(TaxonomyLevel::Class),
            "order" => Ok(TaxonomyLevel::Order), "family" => Ok(TaxonomyLevel::Family), "genus" => Ok(TaxonomyLevel::Genus),
            "species" => Ok(TaxonomyLevel::Species), "strain" => Ok(TaxonomyLevel::Strain), _ => bail!("Invalid taxonomy level"),
        }
    }
    fn as_str(&self) -> &'static str {
        match self {
            TaxonomyLevel::Kingdom => "kingdom", TaxonomyLevel::Phylum => "phylum", TaxonomyLevel::Class => "class",
            TaxonomyLevel::Order => "order", TaxonomyLevel::Family => "family", TaxonomyLevel::Genus => "genus",
            TaxonomyLevel::Species => "species", TaxonomyLevel::Strain => "strain",
        }
    }
}

struct GenomeRecord {
    gcf_id: String,
    taxonomy: Vec<String>,
}

pub fn run(args: BuildQuanDbArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };

    tracing::info!("Reading genome taxonomy list ...");
    let (genome_records, _) = read_genome_list(&args.genome_list, &levels)?;
    tracing::info!("Total {} genomes", genome_records.len());

    std::fs::create_dir_all(&args.output_dir)?;
    let remove_redundant = args.remove_redundant.to_lowercase() == "yes";

    let intermediate_enzyme_file = digest_genomes_to_intermediate_file(
        &genome_records, enzyme, &args.output_dir, args.enzyme_file.as_ref(), args.pre_digested_dir.as_ref()
    )?;

    for level in &levels {
        tracing::info!("\n========== Building {}-level database (Hash mode) ==========", level.as_str());
        build_database_for_level(&intermediate_enzyme_file, enzyme, &args.output_dir, *level, &genome_records, remove_redundant)?;
    }
    tracing::info!("\nAll done!");
    Ok(())
}

fn parse_taxonomy_levels(levels_str: &str) -> Result<Vec<TaxonomyLevel>> {
    if levels_str == "all" {
        Ok(vec![TaxonomyLevel::Kingdom, TaxonomyLevel::Phylum, TaxonomyLevel::Class, TaxonomyLevel::Order, TaxonomyLevel::Family, TaxonomyLevel::Genus, TaxonomyLevel::Species, TaxonomyLevel::Strain])
    } else {
        levels_str.split(',').map(|s| TaxonomyLevel::from_str(s.trim())).collect()
    }
}

fn read_genome_list(list_path: &Path, levels: &[TaxonomyLevel]) -> Result<(Vec<GenomeRecord>, FxHashMap<String, usize>)> {
    let file = File::open(list_path)?;
    let reader = BufReader::new(file);
    let mut genomes = Vec::new();
    let mut taxonomy_levels_map = FxHashMap::default();
    let mut is_gtdb_format = false;
    let mut first_data_line = true;
    for level in levels { taxonomy_levels_map.insert(level.as_str().to_string(), *level as usize); }
    for line in reader.lines() {
        let line = line?;
        let trimmed_line = line.trim();
        if trimmed_line.is_empty() || trimmed_line.starts_with('#') { continue; }
        let parts: Vec<&str> = trimmed_line.split('\t').collect();
        if first_data_line {
            first_data_line = false;
            if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") { is_gtdb_format = true; if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; } }
        }
        if is_gtdb_format {
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord { gcf_id: extract_gcf_id(parts[0].trim()), taxonomy: parse_gtdb_taxonomy(parts[1])? });
        } else {
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord { gcf_id: parts[0].trim().to_string(), taxonomy: parts[1..].iter().map(|s| s.trim().to_string()).collect() });
        }
    }
    Ok((genomes, taxonomy_levels_map))
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

/// Returns (path, was_generated). `was_generated` is false when using an existing file via `-e`.
fn digest_genomes_to_intermediate_file(genomes: &[GenomeRecord], enzyme: &'static Enzyme, output_dir: &Path, enzyme_file: Option<&PathBuf>, pre_digested_dir: Option<&PathBuf>) -> Result<PathBuf> {
    if let Some(existing_file) = enzyme_file {
        tracing::info!("Using pre-digested file: {}", existing_file.display());
        return Ok(existing_file.clone());
    }
    if let Some(dir) = pre_digested_dir {
        tracing::info!("Merging pre-digested files from directory in parallel: {}", dir.display());
        let output_file = output_dir.join(format!("{}.enzyme.iibdb", enzyme.name));
        merge_pre_digested_files(genomes, enzyme, dir, &output_file)?;
        return Ok(output_file);
    }
    bail!("Please provide -e or --pre-digested-dir");
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
        let patterns = [format!("{}.{}.iibdb", genome_id, enzyme.name), format!("{}.{}.iibsp", genome.gcf_id, enzyme.name), format!("{}.{}.iibdb", genome.gcf_id, enzyme.name)];
        for pattern in &patterns {
            let test_path = pre_digested_dir.join(pattern);
            if test_path.exists() {
                if let Ok(tasks) = process_file_content(&test_path, &genome.gcf_id) { s.send(tasks).expect("Failed"); }
                break;
            }
        }
        pb.inc(1);
    });
    pb.finish();
    writer_thread.join().unwrap()?;
    Ok(())
}

fn process_file_content(path: &Path, gcf_id: &str) -> Result<Vec<WriteTask>> {
    let mut tasks = Vec::new();
    let path_str = path.to_string_lossy();
    if path_str.ends_with(".fa.gz") || path_str.ends_with(".fq.gz") {
         let mut reader = parse_fastx_file(path)?;
         while let Some(record) = reader.next() {
            let record = record.context("Parse fail")?;
            let header_str = std::str::from_utf8(record.id()).unwrap_or("");
            let (scaffold, pos) = if let Some(idx) = header_str.rfind(':') { (&header_str[..idx], &header_str[idx+1..]) } else { (header_str, "0") };
            let seq_bytes = record.seq();
            let hash_val = if let Ok(val) = std::str::from_utf8(&seq_bytes).unwrap_or("").trim().parse::<u64>() { val } else {
                let mut seq_vec = seq_bytes.to_vec();
                seq_vec.make_ascii_uppercase();
                canonical_hash(&seq_vec)
            };
            let new_id = format!("{}|0|{}|{}|0|0", gcf_id, scaffold, pos);
            tasks.push(WriteTask { hash: hash_val, id: new_id });
         }
    } else {
        let mut reader = io_utils::open_binary_reader(path)?;
        let mut id_buffer = String::with_capacity(128); // Optimization: reuse buffer
        while let Some(hash) = reader.next_record_reuse(&mut id_buffer)? {
             let new_id = if id_buffer.starts_with("GCF_") || id_buffer.starts_with("GCA_") {
                 id_buffer.clone()
             } else {
                 let (scaffold, pos) = if let Some((s, p)) = id_buffer.rsplit_once(':') { (s, p) } else { (id_buffer.as_str(), "0") };
                 format!("{}|0|{}|{}|0|0", gcf_id, scaffold, pos)
             };
             tasks.push(WriteTask { hash, id: new_id });
        }
    }
    Ok(tasks)
}

// ================== Low-memory streaming approach ==================
// Instead of loading all enzyme records into memory (~2+ GB), we stream
// from disk 2–3 times per taxonomy level. This trades I/O for memory,
// reducing peak usage by ~78%.
//
// Pass 1: Determine per-tag uniqueness (taxonomy interned to u32 indices)
// Pass 2: (only if remove_redundant) Find within-genome duplicate tags
// Pass 3: Write unique tags to compact database

fn build_database_for_level(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<()> {
    // Build lookup tables: GCF ID → genome index, genome index → taxonomy index
    let gcf_to_idx: FxHashMap<&str, u32> = genomes.iter().enumerate()
        .map(|(i, g)| (g.gcf_id.as_str(), i as u32)).collect();

    let mut taxonomy_to_idx: FxHashMap<String, u32> = FxHashMap::default();
    let mut gcf_tax: Vec<u32> = Vec::with_capacity(genomes.len());
    for genome in genomes {
        let end_index = std::cmp::min(level as usize, genome.taxonomy.len());
        let taxonomy_str = genome.taxonomy[0..end_index].join("\t");
        let next_id = taxonomy_to_idx.len() as u32;
        let tax_idx = *taxonomy_to_idx.entry(taxonomy_str).or_insert(next_id);
        gcf_tax.push(tax_idx);
    }

    // ── Pass 1: determine per-tag uniqueness (stream from disk) ──
    // tag_status[hash] = taxonomy_idx   → tag belongs to exactly one taxon
    // tag_status[hash] = u32::MAX       → tag is shared across 2+ taxa
    tracing::info!("Step 1: Collecting tag taxonomy information ...");
    let mut tag_status: FxHashMap<Hash, u32> = FxHashMap::default();
    let mut processed_gcfs: FxHashSet<u32> = FxHashSet::default();
    {
        let mut reader = io_utils::open_binary_reader(enzyme_file)?;
        let mut buf = String::with_capacity(256);
        while let Some(hash) = reader.next_record_reuse(&mut buf)? {
            let gcf_id = buf.split('|').next().unwrap_or("");
            if let Some(&gi) = gcf_to_idx.get(gcf_id) {
                processed_gcfs.insert(gi);
                let tax = gcf_tax[gi as usize];
                tag_status.entry(hash)
                    .and_modify(|v| { if *v != tax && *v != u32::MAX { *v = u32::MAX; } })
                    .or_insert(tax);
            }
        }
    }
    if !genomes.is_empty() {
        let percent = (processed_gcfs.len() * 100) / genomes.len();
        tracing::info!("  Genomes covered: {}/{} ({}%)", processed_gcfs.len(), genomes.len(), percent);
    }
    drop(processed_gcfs);

    // ── Pass 2 (only when remove_redundant): find within-genome duplicate tags ──
    // Only tracks tags already known to be unique to one taxon, keeping the sets small.
    let dup_set: FxHashSet<(u32, Hash)> = if remove_redundant {
        tracing::info!("Step 2: Checking within-genome tag redundancy ...");
        let mut seen: FxHashSet<(u32, Hash)> = FxHashSet::default();
        let mut dup: FxHashSet<(u32, Hash)> = FxHashSet::default();
        {
            let mut reader = io_utils::open_binary_reader(enzyme_file)?;
            let mut buf = String::with_capacity(256);
            while let Some(hash) = reader.next_record_reuse(&mut buf)? {
                // Skip tags already known to be shared across taxa
                if tag_status.get(&hash).map_or(true, |&v| v == u32::MAX) { continue; }
                let gcf_id = buf.split('|').next().unwrap_or("");
                if let Some(&gi) = gcf_to_idx.get(gcf_id) {
                    let key = (gi, hash);
                    if !seen.insert(key) {
                        dup.insert(key);
                    }
                }
            }
        }
        // Release the larger `seen` set; keep only the much smaller `dup` set
        drop(seen);
        dup
    } else {
        FxHashSet::default()
    };

    // ── Pass 3: write unique tags to compact database ──
    let step = if remove_redundant { 3 } else { 2 };
    tracing::info!("Step {}: Writing unique tags to database ...", step);
    let output_path = output_dir.join(format!("{}.{}.iibdb", enzyme.name, level.as_str()));
    let file = File::create(&output_path)?;
    let buf_writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);

    let gcf_ids: Vec<&str> = genomes.iter().map(|g| g.gcf_id.as_str()).collect();
    let mut compact_writer = io_utils::CompactDatabaseWriter::new(buf_writer, &gcf_ids)?;
    let mut unique_genome_count: FxHashSet<u32> = FxHashSet::default();

    {
        let mut reader = io_utils::open_binary_reader(enzyme_file)?;
        let mut buf = String::with_capacity(256);
        while let Some(hash) = reader.next_record_reuse(&mut buf)? {
            // Skip tags shared across taxa
            if tag_status.get(&hash).map_or(true, |&v| v == u32::MAX) { continue; }
            let gcf_id = buf.split('|').next().unwrap_or("");
            if let Some(&gi) = gcf_to_idx.get(gcf_id) {
                let dominated = remove_redundant && dup_set.contains(&(gi, hash));
                if !dominated {
                    unique_genome_count.insert(gi);
                    compact_writer.write_record(hash, gi)?;
                }
            }
        }
    }
    compact_writer.finish()?;

    tracing::info!("  Output database: {}", output_path.display());
    tracing::info!("  Contains unique tags for {} genomes", unique_genome_count.len());
    Ok(())
}
