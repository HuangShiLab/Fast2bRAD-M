use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*; // import rayon
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

/// Select genomes for quantification based on qualitative results
#[derive(Parser, Debug)]
pub struct FindGenomeArgs {
    /// sample list file (TSV format: sample_name<tab>...)
    #[arg(short = 'l', long = "list")]
    pub sample_list: PathBuf,

    /// database directory
    #[arg(short = 'd', long = "database")]
    pub database_dir: PathBuf,

    /// output directory
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// qualitative results directory
    #[arg(long = "qual-dir", alias = "qualdir")]
    pub qual_dir: PathBuf,

    /// G-score threshold (default 5, meaning >5)
    #[arg(long = "gscore", default_value = "5")]
    pub g_score_threshold: i32,

    /// minimum detected tags per GCF (default 1, meaning >1)
    #[arg(long = "gcf", default_value = "1")]
    pub gcf_threshold: i32,

    /// thread count (for parallel sample processing)
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

/// Main function
pub fn run(args: FindGenomeArgs) -> Result<()> {
    // Set up global thread pool
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global();

    tracing::info!(
        "COMMAND: find-genome -l {} -d {} -o {} --qual-dir {} --gscore {} --gcf {} -j {}",
        args.sample_list.display(),
        args.database_dir.display(),
        args.output_dir.display(),
        args.qual_dir.display(),
        args.g_score_threshold,
        args.gcf_threshold,
        args.threads
    );

    // Check database file
    let classify_file = args.database_dir.join("abfh_classify_with_speciename.txt.gz");
    if !classify_file.exists() {
        bail!(
            "Database file not found: {}",
            classify_file.display()
        );
    }

    // Load GCF-to-taxonomy mapping
    let gcf_to_classify = load_gcf_classify(&classify_file)?;
    tracing::info!("Loaded {} genome taxonomy records", gcf_to_classify.len());

    // Create output directory
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("Cannot create output directory: {}", args.output_dir.display()))?;

    // Read sample list
    let samples = read_sample_list(&args.sample_list)?;
    tracing::info!("{} samples to process", samples.len());

    // Process each sample in parallel.
    // We collect errors but let other tasks finish.
    samples.par_iter().for_each(|sample_name| {
        match process_sample(
            sample_name,
            &args.qual_dir,
            &args.database_dir,
            &args.output_dir,
            &gcf_to_classify,
            args.g_score_threshold,
            args.gcf_threshold,
        ) {
            Ok(_) => {},
            Err(e) => tracing::error!("Sample {} processing failed: {}", sample_name, e),
        }
    });

    tracing::info!("\nAll done!");
    Ok(())
}

/// Load GCF-to-taxonomy mapping
fn load_gcf_classify(classify_file: &Path) -> Result<FxHashMap<String, String>> {
    let mut gcf_to_classify = FxHashMap::default();

    let file = File::open(classify_file)
        .with_context(|| format!("Cannot open database file: {}", classify_file.display()))?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.is_empty() {
            continue;
        }

        let gcf_id = parts[0].to_string();
        gcf_to_classify.insert(gcf_id, line.to_string());
    }

    Ok(gcf_to_classify)
}

/// Read sample list
fn read_sample_list(list_file: &Path) -> Result<Vec<String>> {
    let file = File::open(list_file)
        .with_context(|| format!("Cannot open sample list: {}", list_file.display()))?;
    let reader = BufReader::new(file);

    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.is_empty() {
            continue;
        }

        samples.push(parts[0].to_string());
    }

    Ok(samples)
}

/// Process a single sample
fn process_sample(
    sample_name: &str,
    qual_dir: &Path,
    database_dir: &Path,
    output_dir: &Path,
    gcf_to_classify: &FxHashMap<String, String>,
    g_score_threshold: i32,
    gcf_threshold: i32,
) -> Result<()> {
    // Read combine.xls; if absent, fall back to individual enzyme result files
    let combine_file = qual_dir.join(sample_name).join(format!("{}.combine.xls", sample_name));

    let (enzymes, pass_gscore_classes) = if combine_file.exists() {
        // Use combine.xls
        parse_combine_file(&combine_file, g_score_threshold)?
    } else {
        // Fallback: scan individual enzyme result files
        let sample_dir = qual_dir.join(sample_name);
        let mut found_enzymes = Vec::new();
        let mut found_classes = FxHashSet::default();

        // Scan all {sample}.{enzyme}.xls files
        if let Ok(entries) = std::fs::read_dir(&sample_dir) {
            for entry in entries {
                let entry = entry?;
                let path = entry.path();
                if let Some(file_name) = path.file_name().and_then(|n| n.to_str()) {
                    // Match format: {sample}.{enzyme}.xls, excluding GCF_detected.xls
                    if file_name.ends_with(".xls")
                        && !file_name.contains("GCF_detected")
                        && file_name.starts_with(&format!("{}.", sample_name)) {
                        // Extract enzyme name: {sample}.{enzyme}.xls -> {enzyme}
                        let enzyme_part = file_name
                            .strip_prefix(&format!("{}.", sample_name))
                            .and_then(|s| s.strip_suffix(".xls"));

                        if let Some(enzyme) = enzyme_part {
                            found_enzymes.push(enzyme.to_string());
                            // Parse the file and extract taxa passing the G-score threshold
                            let (_, classes) = parse_single_enzyme_file(&path, g_score_threshold)?;
                            found_classes.extend(classes);
                        }
                    }
                }
            }
        }

        if found_enzymes.is_empty() {
            tracing::warn!(
                "!!! {} no qualitative result files (no combine.xls or individual enzyme result files), skipping quantification",
                sample_name
            );
            return Ok(());
        }

        (found_enzymes, found_classes)
    };

    if enzymes.is_empty() {
        tracing::warn!("Warning: {} no enzymes found", sample_name);
        return Ok(());
    }

    tracing::info!("Sample {}: using {} enzymes, {} taxa pass G-score threshold",
             sample_name, enzymes.len(), pass_gscore_classes.len());

    // Create sample output directory
    let sample_output_dir = output_dir.join(sample_name);
    std::fs::create_dir_all(&sample_output_dir)?;

    // Collect all qualifying genomes
    let mut selected_genomes = FxHashSet::default();

    // Iterate over each enzyme
    for enzyme in &enzymes {
        // Check whether the database file exists
        let enzyme_db = database_dir.join(format!("{}.species.iibdb", enzyme));
        if !enzyme_db.exists() {
            tracing::warn!(
                "Warning: Database file not found: {}",
                enzyme_db.display()
            );
            continue;
        }

        // Read GCF_detected.xls file
        let gcf_detected_file = qual_dir
            .join(sample_name)
            .join(format!("{}.{}.GCF_detected.xls", sample_name, enzyme));

        if !gcf_detected_file.exists() {
            tracing::warn!(
                "Warning: {} no GCF_detected.xls for {}: {}",
                sample_name,
                enzyme,
                gcf_detected_file.display()
            );
            continue;
        }

        // Parse GCF_detected.xls
        let gcf_list = parse_gcf_detected_file(
            &gcf_detected_file,
            &pass_gscore_classes,
            gcf_threshold,
        )?;

        for gcf_id in gcf_list {
            selected_genomes.insert(gcf_id);
        }
    }

    // Write sdb.list file (sort and deduplicate)
    let sdb_list_file = sample_output_dir.join("sdb.list");
    let mut writer = BufWriter::new(File::create(&sdb_list_file)?);

    let mut genome_list: Vec<&String> = selected_genomes
        .iter()
        .filter_map(|gcf_id| gcf_to_classify.get(gcf_id))
        .collect();
    genome_list.sort();
    genome_list.dedup();

    let genome_count = genome_list.len();
    for genome_line in genome_list {
        writeln!(writer, "{}", genome_line)?;
    }

    tracing::info!("Sample {}: Selected {} genomes", sample_name, genome_count);
    Ok(())
}

/// Parse combine.xls file, extracting the enzymes used and taxa passing the G-score threshold
fn parse_combine_file(
    combine_file: &Path,
    g_score_threshold: i32,
) -> Result<(Vec<String>, FxHashSet<String>)> {
    let file = File::open(combine_file)
        .with_context(|| format!("Cannot open combine file: {}", combine_file.display()))?;
    let reader = BufReader::new(file);

    let mut enzymes = Vec::new();
    let mut pass_gscore_classes = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        // Skip header line
        if line.to_uppercase().starts_with("#KINGDOM") {
            continue;
        }

        // Parse enzyme info line (e.g. #BcgI CjeI combine)
        if line.starts_with('#') && !line.to_uppercase().starts_with("#KINGDOM") {
            let parts: Vec<&str> = line.trim_start_matches('#').split_whitespace().collect();
            for part in parts {
                let enzyme = part.trim();
                if enzyme != "combine" && !enzyme.is_empty() {
                    enzymes.push(enzyme.to_string());
                }
            }
            continue;
        }

        // Parse data row
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }

        // Last column is G_Score
        if let Ok(g_score) = parts[parts.len() - 1].parse::<f64>() {
            if g_score > g_score_threshold as f64 {
                // Taxonomy info is the first N-8 columns
                let class = parts[0..parts.len() - 8].join("\t");
                pass_gscore_classes.insert(class);
            }
        }
    }

    // Deduplicate enzyme list
    enzymes.sort();
    enzymes.dedup();

    Ok((enzymes, pass_gscore_classes))
}

/// Parse individual enzyme result file (same format as combine.xls but without an enzyme info line)
fn parse_single_enzyme_file(
    enzyme_file: &Path,
    g_score_threshold: i32,
) -> Result<(Vec<String>, FxHashSet<String>)> {
    let file = File::open(enzyme_file)
        .with_context(|| format!("Cannot open enzyme result file: {}", enzyme_file.display()))?;
    let reader = BufReader::new(file);

    let mut pass_gscore_classes = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        // Skip header line
        if line.to_uppercase().starts_with("#KINGDOM") {
            continue;
        }

        // Skip comment line
        if line.starts_with('#') {
            continue;
        }

        // Parse data row
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }

        // Last column is G_Score
        if let Ok(g_score) = parts[parts.len() - 1].parse::<f64>() {
            if g_score > g_score_threshold as f64 {
                // Taxonomy info is the first N-8 columns
                let class = parts[0..parts.len() - 8].join("\t");
                pass_gscore_classes.insert(class);
            }
        }
    }

    // Single-enzyme files contain no enzyme list info; return empty
    Ok((Vec::new(), pass_gscore_classes))
}

/// Parse GCF_detected.xls file and return qualifying GCF IDs
fn parse_gcf_detected_file(
    gcf_detected_file: &Path,
    pass_gscore_classes: &FxHashSet<String>,
    gcf_threshold: i32,
) -> Result<Vec<String>> {
    let file = File::open(gcf_detected_file)
        .with_context(|| format!("Cannot open GCF_detected file: {}", gcf_detected_file.display()))?;
    let reader = BufReader::new(file);

    let mut gcf_list = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }

        // Format: class\tGCF\tGCF_all_theory_num\tdetected_tag_num\tpercent
        // Taxonomy info is the first N-4 columns
        let class = parts[0..parts.len() - 4].join("\t");
        let gcf_id = parts[parts.len() - 4].to_string();
        let detected_tag_num: i32 = parts[parts.len() - 2].parse().unwrap_or(0);

        // Keep entries whose taxon passes the G-score threshold and whose detected tag count > gcf_threshold
        if pass_gscore_classes.contains(&class) && detected_tag_num > gcf_threshold {
            gcf_list.push(gcf_id);
        }
    }

    Ok(gcf_list)
}
