use anyhow::{bail, Result, anyhow};
use clap::Parser;
use fxhash::FxHashMap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tracing;

use crate::enzymes::{enzyme_by_id, enzyme_by_name};
use crate::io_utils;

// Interned string type: Arc<str> allows cheap cloning (atomic refcount increment)
// instead of full heap allocation for each String::clone()
type Istr = Arc<str>;

#[derive(Parser, Debug)]
pub struct QuantifyArgs {
    // ── Input ──
    /// Sample list file (TSV: sample_name<TAB>path_to.iibsp)
    #[arg(short = 'l', long = "list", help_heading = "Input")]
    pub sample_list: PathBuf,
    /// Database directory (containing {enzyme}.{level}.iibdb and classify file)
    #[arg(short = 'd', long = "database", help_heading = "Input")]
    pub database_dir: PathBuf,
    /// Taxonomy level: kingdom|phylum|class|order|family|genus|species|strain
    #[arg(short = 't', long = "taxonomy", help_heading = "Input")]
    pub taxonomy_level: String,
    /// Enzyme name (e.g. BcgI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", help_heading = "Input")]
    pub enzyme_site: String,

    // ── Output ──
    /// Output directory
    #[arg(short = 'o', long = "output", help_heading = "Output")]
    pub output_dir: PathBuf,

    // ── Filtering ──
    /// G-score threshold: taxa with G-score below this are excluded (0=no filtering)
    #[arg(short = 'g', long = "gscore", default_value = "0", help_heading = "Filtering")]
    pub g_score_threshold: f64,
    /// Minimum average read depth per sequenced tag (i.e. species-unique marker coverage,
    /// Sequenced_Reads_Num / Sequenced_Tag_Num). Taxa below this are excluded (0=no filtering).
    #[arg(long = "min-tag-depth", visible_alias = "min-coverage", default_value = "0", help_heading = "Filtering")]
    pub min_tag_depth: f64,

    // ── Options ──
    /// Output per-tag detail files (yes/no)
    #[arg(short = 'v', long = "verbose", default_value = "yes", help_heading = "Options")]
    pub verbose: String,

    // ── Performance ──
    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4", help_heading = "Performance")]
    pub threads: usize,
}

#[derive(Debug, Default)]
struct AbundanceStats {
    theoretical_tag_num: f64,
    sequenced_tag_num: usize,
    sequenced_reads_num: usize,
    sequenced_tag_num_gt1: usize,
}

impl AbundanceStats {
    fn percent(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 { (self.sequenced_tag_num as f64 / self.theoretical_tag_num) * 100.0 } else { 0.0 }
    }
    fn reads_per_theoretical(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 { self.sequenced_reads_num as f64 / self.theoretical_tag_num } else { 0.0 }
    }
    fn reads_per_sequenced(&self) -> f64 {
        if self.sequenced_tag_num > 0 { self.sequenced_reads_num as f64 / self.sequenced_tag_num as f64 } else { 0.0 }
    }
    fn g_score(&self) -> f64 {
        ((self.sequenced_tag_num * self.sequenced_reads_num) as f64).sqrt()
    }
}

// Pre-computed theory stats per taxon (replaces 3-level nested HashMap)
struct TaxonTheory {
    /// Pre-computed: total_unique_tags / gcf_count
    theoretical_tag_num: f64,
    /// Per-GCF unique tag count (for GCF_detected output)
    gcf_unique_tag_count: FxHashMap<Istr, usize>,
}

type TaxonTheoryMap = FxHashMap<Istr, TaxonTheory>;

pub fn run(args: QuantifyArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();
    let verbose = args.verbose.to_lowercase() == "yes";
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme site ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };
    let tax_level = validate_taxonomy_level(&args.taxonomy_level)?;

    tracing::info!("COMMAND: quantify -l {} -d {} -t {} -s {} -o {} -g {} --min-tag-depth {} -v {} -j {}", args.sample_list.display(), args.database_dir.display(), args.taxonomy_level, args.enzyme_site, args.output_dir.display(), args.g_score_threshold, args.min_tag_depth, args.verbose, args.threads);

    std::fs::create_dir_all(&args.output_dir)?;
    let db_file = args.database_dir.join(format!("{}.{}.iibdb", enzyme.name, tax_level));
    let classify_file = args.database_dir.join("abfh_classify_with_speciename.txt.gz");

    if !db_file.exists() { bail!("Database file not found"); }
    if !classify_file.exists() { bail!("Taxonomy file not found"); }

    tracing::info!("### Loading database: {}", db_file.display());
    let (tag_to_gcfs, gcf_to_taxonomy, taxon_theory) = load_database(&db_file, &classify_file, &args.taxonomy_level)?;
    tracing::info!("### Database loaded");

    let samples = read_sample_list(&args.sample_list)?;
    tracing::info!("{} samples to process", samples.len());

    samples.par_iter().for_each(|(sample_name, sample_data)| {
        tracing::info!(">>> ({}) Sample analysis started", sample_name);
        let result = process_sample(sample_name, sample_data, &tag_to_gcfs, &gcf_to_taxonomy, &taxon_theory, enzyme, &args.output_dir, args.g_score_threshold, args.min_tag_depth, verbose);
        match result {
            Ok(_) => tracing::info!("<<< ({}) Sample analysis completed", sample_name),
            Err(e) => tracing::error!("!!! ({}) Error: {}", sample_name, e),
        }
    });
    tracing::info!("\nAll done!");
    Ok(())
}

fn validate_taxonomy_level(level: &str) -> Result<String> {
    let valid = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"];
    if valid.contains(&level) { Ok(level.to_string()) } else { bail!("Invalid taxonomy level"); }
}

fn read_sample_list(list_path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let content = std::fs::read_to_string(list_path)?;
    let mut samples = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 { continue; }
        let sample_path = PathBuf::from(parts[1]);
        if !sample_path.exists() { tracing::warn!("Warning: sample file not found: {}", sample_path.display()); continue; }
        samples.push((parts[0].to_string(), sample_path));
    }
    Ok(samples)
}

fn load_database(
    db_file: &Path,
    classify_file: &Path,
    tax_level: &str,
) -> Result<(FxHashMap<u64, Vec<Istr>>, FxHashMap<Istr, Istr>, TaxonTheoryMap)> {
    let level_index = get_taxonomy_level_index(tax_level);

    // Step 1: Read classify file, intern taxonomy strings
    let content = if classify_file.to_str().unwrap().ends_with(".gz") {
        use flate2::read::GzDecoder; use std::io::Read;
        let file = File::open(classify_file)?;
        let mut decoder = GzDecoder::new(file);
        let mut c = String::new(); decoder.read_to_string(&mut c)?; c
    } else { std::fs::read_to_string(classify_file)? };

    // Intern pool: deduplicates taxonomy strings so identical values share one Arc
    let mut taxonomy_intern: FxHashMap<String, Istr> = FxHashMap::default();
    // Temporary map: gcf_id (String) → taxonomy (Arc<str>)
    let mut classify_map: FxHashMap<String, Istr> = FxHashMap::default();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= level_index { continue; }
        let gcf_id = parts[0].to_string();
        let taxonomy_str = parts[1..=level_index].join("\t");
        let taxonomy_arc = taxonomy_intern.entry(taxonomy_str.clone())
            .or_insert_with(|| Arc::from(taxonomy_str.as_str()))
            .clone();
        classify_map.insert(gcf_id, taxonomy_arc);
    }
    drop(taxonomy_intern); // no longer needed

    // Step 2: Read compact DB, build interned GCF table
    let mut compact_reader = io_utils::open_compact_reader(db_file)?;
    let gcf_intern: Vec<Istr> = compact_reader.gcf_table()
        .iter()
        .map(|s| Arc::from(s.as_str()))
        .collect();

    // Step 3: Build gcf_to_taxonomy with interned keys
    let mut gcf_to_taxonomy: FxHashMap<Istr, Istr> = FxHashMap::default();
    for gcf_arc in &gcf_intern {
        if let Some(taxonomy) = classify_map.get(gcf_arc.as_ref()) {
            gcf_to_taxonomy.insert(gcf_arc.clone(), taxonomy.clone());
        }
    }
    drop(classify_map); // no longer needed

    // Step 4: Load tag records
    let mut tag_to_gcfs: FxHashMap<u64, Vec<Istr>> = FxHashMap::default();
    let mut gcf_tag_counts: FxHashMap<Istr, FxHashMap<Istr, usize>> = FxHashMap::default();
    let mut loaded_count = 0usize;

    while let Some((hash_val, gcf_index)) = compact_reader.next_record()? {
        let gcf_id = &gcf_intern[gcf_index as usize];

        if let Some(taxonomy) = gcf_to_taxonomy.get(gcf_id) {
            // Arc::clone() = atomic refcount increment (cheap, no heap alloc)
            tag_to_gcfs.entry(hash_val).or_default().push(gcf_id.clone());
            *gcf_tag_counts.entry(taxonomy.clone())
                .or_default()
                .entry(gcf_id.clone())
                .or_insert(0) += 1;
            loaded_count += 1;
        }
    }

    // Step 5: Convert to pre-computed TaxonTheory
    let mut taxon_theory: TaxonTheoryMap = FxHashMap::default();
    for (taxonomy, gcf_counts) in gcf_tag_counts {
        let gcf_count = gcf_counts.len();
        let total: usize = gcf_counts.values().sum();
        taxon_theory.insert(taxonomy, TaxonTheory {
            theoretical_tag_num: total as f64 / gcf_count as f64,
            gcf_unique_tag_count: gcf_counts,
        });
    }

    tracing::info!("Successfully loaded {} valid unique tag entries from database", loaded_count);
    Ok((tag_to_gcfs, gcf_to_taxonomy, taxon_theory))
}

fn get_taxonomy_level_index(level: &str) -> usize {
    match level { "kingdom" => 1, "phylum" => 2, "class" => 3, "order" => 4, "family" => 5, "genus" => 6, "species" => 7, "strain" => 8, _ => 7 }
}

fn process_sample(
    sample_name: &str,
    sample_data: &Path,
    tag_to_gcfs: &FxHashMap<u64, Vec<Istr>>,
    gcf_to_taxonomy: &FxHashMap<Istr, Istr>,
    taxon_theory: &TaxonTheoryMap,
    enzyme: &crate::enzymes::Enzyme,
    output_dir: &Path,
    g_score_threshold: f64,
    min_tag_depth: f64,
    verbose: bool,
) -> Result<()> {
    // All .clone() calls on Istr (Arc<str>) are O(1) atomic increments — no heap allocation
    let mut tag_num: FxHashMap<Istr, FxHashMap<u64, usize>> = FxHashMap::default();
    let mut detected_gcf_tag: FxHashMap<Istr, FxHashMap<Istr, fxhash::FxHashSet<u64>>> = FxHashMap::default();
    let mut reader = io_utils::open_binary_reader(sample_data)?;
    let mut ignored_id_buf = String::with_capacity(128);

    while let Some(tag_hash) = reader.next_record_reuse(&mut ignored_id_buf)? {
        if let Some(gcf_list) = tag_to_gcfs.get(&tag_hash) {
            if let Some(first_gcf) = gcf_list.first() {
                if let Some(taxonomy) = gcf_to_taxonomy.get(first_gcf) {
                    *tag_num.entry(taxonomy.clone()).or_default().entry(tag_hash).or_insert(0) += 1;
                    for gcf_id in gcf_list {
                        detected_gcf_tag.entry(taxonomy.clone()).or_default()
                            .entry(gcf_id.clone()).or_default().insert(tag_hash);
                    }
                }
            }
        }
    }

    if tag_num.is_empty() { tracing::warn!("!!! ({}) Warning: no tags detected", sample_name); return Ok(()); }

    let sample_dir = output_dir.join(sample_name);
    std::fs::create_dir_all(&sample_dir)?;
    let gcf_detected_file = sample_dir.join(format!("{}.{}.GCF_detected.xls", sample_name, enzyme.name));
    let mut gcf_writer = BufWriter::new(File::create(&gcf_detected_file)?);

    let mut taxonomy_list: Vec<&Istr> = detected_gcf_tag.keys().collect();
    taxonomy_list.sort_by(|a, b| a.as_ref().cmp(b.as_ref()));

    for taxonomy in taxonomy_list {
        let gcf_map = &detected_gcf_tag[taxonomy];
        let mut gcf_list: Vec<&Istr> = gcf_map.keys().collect();
        gcf_list.sort_by(|a, b| a.as_ref().cmp(b.as_ref()));
        for gcf_id in gcf_list {
            let detected_tags = &gcf_map[gcf_id];
            let detected_tag_num = detected_tags.len();
            let gcf_all_theory_num = if let Some(theory) = taxon_theory.get(taxonomy) {
                theory.gcf_unique_tag_count.get(gcf_id).copied().unwrap_or(0)
            } else { 0 };
            let percent = if gcf_all_theory_num > 0 { detected_tag_num as f64 / gcf_all_theory_num as f64 } else { 0.0 };
            writeln!(gcf_writer, "{}\t{}\t{}\t{}\t{:.4}", taxonomy.as_ref(), gcf_id.as_ref(), gcf_all_theory_num, detected_tag_num, percent)?;
        }
    }

    let output_file = sample_dir.join(format!("{}.{}.xls", sample_name, enzyme.name));
    let mut writer = BufWriter::new(File::create(&output_file)?);
    let level_count = if let Some(first_tax) = gcf_to_taxonomy.values().next() { first_tax.split('\t').count() } else { 7 };
    let tax_levels = get_taxonomy_header_by_level(level_count);
    writeln!(writer, "#{}\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\tSequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\tSequenced_Reads_Num/Sequenced_Tag_Num\tSequenced_Tag_Num(depth>1)\tG_Score", tax_levels)?;

    for (taxonomy, tags) in &tag_num {
        let mut stats = AbundanceStats::default();
        if let Some(theory) = taxon_theory.get(taxonomy) {
            stats.theoretical_tag_num = theory.theoretical_tag_num;
        }
        stats.sequenced_tag_num = tags.len();
        stats.sequenced_reads_num = tags.values().sum();
        stats.sequenced_tag_num_gt1 = tags.values().filter(|&&count| count > 1).count();
        let g_score = stats.g_score();
        let tag_depth = stats.reads_per_sequenced();
        if g_score < g_score_threshold || tag_depth < min_tag_depth { continue; }

        writeln!(writer, "{}\t{:.8}\t{}\t{:.8}%\t{}\t{:.8}\t{:.8}\t{}\t{:.8}", taxonomy.as_ref(), stats.theoretical_tag_num, stats.sequenced_tag_num, stats.percent(), stats.sequenced_reads_num, stats.reads_per_theoretical(), tag_depth, stats.sequenced_tag_num_gt1, g_score)?;

        if verbose {
            let output_name = taxonomy.split('\t').last().unwrap_or("unknown");
            let detail_dir = sample_dir.join(format!("{}.{}", sample_name, enzyme.name));
            std::fs::create_dir_all(&detail_dir)?;
            let detail_file = detail_dir.join(format!("{}.xls", output_name));
            let mut detail_writer = BufWriter::new(File::create(detail_file)?);
            for (tag, &count) in tags { writeln!(detail_writer, "{}\t{}", tag, count)?; }
        }
    }
    Ok(())
}

fn get_taxonomy_header_by_level(level_count: usize) -> String {
    let headers = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"];
    headers[..level_count.min(8)].join("\t")
}
