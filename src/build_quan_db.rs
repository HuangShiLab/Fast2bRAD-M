use anyhow::{bail, Context, Result, anyhow};
use clap::Parser;
// use flate2::write::GzEncoder; // Removed
// use flate2::Compression;      // Removed
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader, BufWriter, Write}; // Import BufWriter
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

pub fn hash_bytes(bytes: &[u8]) -> Hash {
    let mut hasher = FxHasher::default();
    hasher.write(bytes);
    hasher.finish()
}

fn get_canonical_sequence(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let complement = match b {
            b'A' | b'a' => b'T', b'T' | b't' => b'A', b'C' | b'c' => b'G', b'G' | b'g' => b'C', b'N' | b'n' => b'N', x => x,
        };
        rc.push(complement);
    }
    if seq <= rc.as_slice() { seq.to_vec() } else { rc }
}

#[derive(Parser, Debug)]
pub struct BuildQuanDbArgs {
    #[arg(short = 'l', long = "list")]
    pub genome_list: PathBuf,
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,
    #[arg(short = 't', long = "taxonomy")]
    pub taxonomy_levels: String,
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,
    #[arg(short = 'e', long = "enzyme-file")]
    pub enzyme_file: Option<PathBuf>,
    #[arg(long = "pre-digested-dir")]
    pub pre_digested_dir: Option<PathBuf>,
    #[arg(short = 'r', long = "remove-redundant", default_value = "no")]
    pub remove_redundant: String,
    #[arg(short = 'j', long = "threads", default_value = "4")]
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
        // Fix: BufWriter
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(tasks) = receiver.recv() {
            for task in tasks {
                io_utils::write_binary_record(&mut writer, task.hash, &task.id)?;
            }
        }
        writer.flush()?;
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
                let canonical = get_canonical_sequence(&seq_vec);
                hash_bytes(&canonical)
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

fn build_database_for_level(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<()> {
    tracing::info!("Step 1: Collecting tag taxonomy information ...");
    let mut gcf_to_taxonomy = FxHashMap::default();
    for genome in genomes {
        let end_index = std::cmp::min(level as usize, genome.taxonomy.len());
        let taxonomy_str = genome.taxonomy[0..end_index].join("\t");
        gcf_to_taxonomy.insert(genome.gcf_id.clone(), taxonomy_str);
    }
    let (tag_taxonomy, genome_tags) = collect_tag_taxonomies(enzyme_file, &gcf_to_taxonomy, genomes, remove_redundant)?;
    tracing::info!("Step 2: Identifying unique tags and writing database ...");
    identify_and_output_unique_tags(enzyme_file, enzyme, output_dir, level, &tag_taxonomy, &genome_tags, remove_redundant)?;
    Ok(())
}

type TagTaxonomyMap = FxHashMap<Hash, FxHashSet<String>>;
type GenomeTagCountMap = FxHashMap<String, FxHashMap<Hash, usize>>;

fn collect_tag_taxonomies(
    enzyme_file: &Path,
    gcf_to_taxonomy: &FxHashMap<String, String>,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<(TagTaxonomyMap, GenomeTagCountMap)> {
    let mut tag_taxonomy: TagTaxonomyMap = FxHashMap::default();
    let mut genome_tags: GenomeTagCountMap = FxHashMap::default();
    let mut reader = io_utils::open_binary_reader(enzyme_file)?;
    let mut processed_gcfs = FxHashSet::default();
    let mut id_buffer = String::with_capacity(256); // Reuse buffer

    while let Some(hash_val) = reader.next_record_reuse(&mut id_buffer)? {
        let mut parts = id_buffer.split('|');
        let gcf_id = parts.next().unwrap_or("");
        if gcf_id.is_empty() { continue; }
        if !gcf_to_taxonomy.contains_key(gcf_id) { continue; }
        let taxonomy = gcf_to_taxonomy.get(gcf_id).unwrap();
        processed_gcfs.insert(gcf_id.to_string());
        tag_taxonomy.entry(hash_val).or_insert_with(FxHashSet::default).insert(taxonomy.clone());
        if remove_redundant {
            *genome_tags.entry(gcf_id.to_string()).or_insert_with(FxHashMap::default).entry(hash_val).or_insert(0) += 1;
        }
    }
    let percent = (processed_gcfs.len() * 100) / genomes.len();
    tracing::info!("  Genomes covered: {}/{} ({}%)", processed_gcfs.len(), genomes.len(), percent);
    Ok((tag_taxonomy, genome_tags))
}

fn identify_and_output_unique_tags(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    tag_taxonomy: &TagTaxonomyMap,
    genome_tags: &GenomeTagCountMap,
    remove_redundant: bool,
) -> Result<()> {
    let output_path = output_dir.join(format!("{}.{}.iibdb", enzyme.name, level.as_str()));
    let file = File::create(&output_path)?;
    // Fix: BufWriter
    let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
    let mut reader = io_utils::open_binary_reader(enzyme_file)?;
    let mut unique_counts: FxHashMap<String, usize> = FxHashMap::default();
    let mut id_buffer = String::with_capacity(256); // Reuse buffer

    while let Some(hash_val) = reader.next_record_reuse(&mut id_buffer)? {
        let mut parts = id_buffer.split('|');
        let gcf_id = parts.next().unwrap_or("");
        let tag_index = parts.next().unwrap_or("0");
        let scaffold_id = parts.next().unwrap_or("scaffold");
        let pos = parts.next().unwrap_or("0");
        if gcf_id.is_empty() { continue; }
        let mut is_unique = false;
        if let Some(taxonomies) = tag_taxonomy.get(&hash_val) {
            if taxonomies.len() == 1 { is_unique = true; }
        }
        if is_unique && remove_redundant {
            if let Some(genome_tag_counts) = genome_tags.get(gcf_id) {
                if let Some(&count) = genome_tag_counts.get(&hash_val) { if count > 1 { is_unique = false; } }
            }
        }
        if is_unique {
            *unique_counts.entry(gcf_id.to_string()).or_insert(0) += 1;
            let new_id = format!("{}|{}|{}|{}|0|1", gcf_id, tag_index, scaffold_id, pos);
            io_utils::write_binary_record(&mut writer, hash_val, &new_id)?;
        }
    }
    writer.flush()?;

    tracing::info!("  Output database: {}", output_path.display());
    tracing::info!("  Contains unique tags for {} genomes", unique_counts.len());
    Ok(())
}
