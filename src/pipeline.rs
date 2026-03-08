use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering};
use tracing;

use crate::build_qual_db;
use crate::build_quan_db;
use crate::extract;
use crate::find_genome;
use crate::merge;
use crate::predict;
use crate::quantify;
use crate::enzymes::{enzyme_by_id, enzyme_by_name};

// --- PipelineArgs struct ---

#[derive(Parser, Debug)]
pub struct PipelineArgs {
    /// Run mode: full|db-only|sample-only
    #[arg(long = "mode", default_value = "full")]
    pub mode: String,

    /// Sample list (TSV: sample_name<TAB>path1[<TAB>path2])
    /// Optional in db-only mode; required in all other modes
    #[arg(long = "samples", alias = "list", short = 'l')]
    pub samples: Option<PathBuf>,

    /// Renamed parameter: genome sequence list (TSV: genome_id<TAB>fasta_path)
    /// Used only in db-only mode for batch-extracting 2bRAD features from reference genomes
    #[arg(long = "genome-list")]
    pub genome_list: Option<PathBuf>,

    /// Renamed parameter: genome taxonomy list (TSV: genome_id<TAB>taxonomy...)
    /// Used for database construction
    #[arg(long = "taxonomy")]
    pub taxonomy: Option<PathBuf>,

    /// Enzyme name or ID (e.g. BcgI or 1-16)
    #[arg(long = "site", short = 's')]
    pub site: String,

    /// Taxonomy level (kingdom|phylum|class|order|family|genus|species|strain)
    #[arg(long = "level", short = 't', default_value = "species")]
    pub level: String,

    /// Output directory
    #[arg(long = "outdir", alias = "od")]
    pub outdir: PathBuf,

    /// Global thread count
    #[arg(short = 'j', long = "threads")]
    pub threads: Option<usize>,

    /// G-score threshold (default 5)
    #[arg(long = "gscore", default_value = "5.0")]
    pub g_score: f64,

    /// Minimum detected tags per GCF (default 1)
    #[arg(long = "gcf", default_value = "1")]
    pub gcf_threshold: i32,

    /// Resume from checkpoint
    #[arg(long = "resume", default_value = "no")]
    pub resume: String,

    // --- extract options ---
    #[arg(long = "qc", default_value = "yes")]
    pub quality_control: String,
    #[arg(long = "max-n", default_value = "0.08")]
    pub max_n: f64,
    #[arg(long = "min-qual", default_value = "30")]
    pub min_quality: u8,
    #[arg(long = "min-qual-percent", default_value = "80")]
    pub min_quality_percent: u8,
    #[arg(long = "qual-base", default_value = "33")]
    pub quality_base: u8,
    #[arg(long = "pear-bin")]
    pub pear_bin: Option<String>,

    /// Threads per PEAR process
    #[arg(long = "pc", default_value = "1")]
    pub pear_threads: usize,

    /// Updated help text: whether to use PEAR merging
    #[arg(long = "use-pear", default_value = "no", help = "Whether to use PEAR merging (yes/no, default: no). Choosing yes will significantly slow down analysis.")]
    pub use_pear: String,

    // --- build-db options ---
    #[arg(long = "database")]
    pub database_dir: Option<PathBuf>,
    #[arg(long = "pre-digested-dir")]
    pub pre_digested_dir: Option<PathBuf>,

    // --- merge options ---
    #[arg(long = "prefix", default_value = "Abundance_Stat")]
    pub prefix: String,
    #[arg(long = "mock")]
    pub mock_samples: Option<String>,
    #[arg(long = "control")]
    pub control_samples: Option<String>,

    // --- predict options ---
    /// Species-to-function mapping matrix (TSV: first column = species name, remaining columns = function IDs, values = gene copy numbers)
    /// If provided, pipeline automatically runs functional prediction after merge
    #[arg(long = "ko-mapping")]
    pub ko_mapping: Option<PathBuf>,
}

pub fn run(args: PipelineArgs) -> Result<()> {
    // 1. Explicitly initialize the Rayon global thread pool
    if let Some(n) = args.threads {
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
        {
            tracing::warn!("Note: Rayon thread pool already initialized; --threads may not have fully taken effect (error: {})", e);
        } else {
            tracing::info!("Global parallel thread count set to: {}", n);
        }
    }

    // Validate: modes other than db-only require a sample list
    if args.mode != "db-only" && args.samples.is_none() {
        bail!("Error: run mode '{}' requires a sample list file via --samples or -l.", args.mode);
    }

    // Create directory structure
    let d01 = args.outdir.join("01_extract");
    let d02 = args.outdir.join("02_db_qual");
    let d04 = args.outdir.join("04_quantify");
    let d05 = args.outdir.join("05_merge");
    let d_qual = args.outdir.join("qualitative");
    let d_sdb = args.outdir.join("quantitative_sdb");

    // Output directory for pre-digested files in db-only mode
    let pre_digested_output = args.outdir.join("pre_digested_output");

    if args.mode != "db-only" {
        std::fs::create_dir_all(&d01)?;
        std::fs::create_dir_all(&d04)?;
        std::fs::create_dir_all(&d05)?;
        std::fs::create_dir_all(&d_qual)?;
        std::fs::create_dir_all(&d_sdb)?;
    }
    if args.mode != "sample-only" {
        let _ = std::fs::create_dir_all(&d02);
    }

    let resume = args.resume.eq_ignore_ascii_case("yes");

    // ======================================================
    // Step 0: [db-only] Batch-extract reference genome features (Extract References)
    // ======================================================
    let mut current_pre_digested_dir = args.pre_digested_dir.clone();

    if args.mode == "db-only" {
        if let Some(genome_list_path) = &args.genome_list {
            tracing::info!("\n===== [db-only] Batch-extract reference genome features (Extract References) =====");
            std::fs::create_dir_all(&pre_digested_output)?;

            let extract_done_marker = pre_digested_output.join(".done");
            if !resume || !extract_done_marker.exists() {
                let extract_global_threads = args.threads.unwrap_or(1);

                // Configure ExtractArgs to process reference genomes (Type 1)
                let ext_args = extract::ExtractArgs {
                    genome_list: Some(genome_list_path.clone()), // pass in sequence list
                    input: vec![],
                    input_type: 1, // ReferenceGenome
                    enzyme_site: args.site.clone(),
                    output_dir: pre_digested_output.clone(),
                    output_prefix: vec![],
                    threads: extract_global_threads,
                    // Read-level QC is generally not needed for reference genomes, but keep parameter consistency
                    quality_control: "no".to_string(),
                    max_n: args.max_n,
                    min_quality: args.min_quality,
                    min_quality_percent: args.min_quality_percent,
                    quality_base: args.quality_base,
                    pear_bin:"pear".to_string(),
                    pear_threads: args.pear_threads,
                    use_pear: "no".to_string(),
                };

                extract::run(ext_args).context("db-only batch reference genome extraction failed")?;
                std::fs::write(&extract_done_marker, b"ok")?;
            } else {
                tracing::info!("resume=on: skipping reference genome extraction (already done)");
            }

            // Update the pre-digested directory used by the subsequent build-db step
            current_pre_digested_dir = Some(pre_digested_output);
        }
    }

    // ======================================================
    // Step 1: Extract sample tags (Extract Samples) - (non db-only)
    // ======================================================
    if args.mode != "db-only" {
        let step_num = if args.mode == "sample-only" { "1/5" } else { "1/5" };
        tracing::info!("\n===== [{}] Extract sample tags (extract) =====", step_num);
        tracing::info!("Concurrency: {} samples in parallel, {} PEAR threads each",
            args.threads.unwrap_or(rayon::current_num_threads()),
            args.pear_threads
        );

        let use_pear_bool = args.use_pear.eq_ignore_ascii_case("yes");
        if use_pear_bool {
            tracing::warn!("!!! Warning: PEAR merging enabled (--use-pear yes). This will significantly increase processing time!");
        } else {
            tracing::info!("PEAR merging disabled (default). Paired-end reads will be extracted independently and combined.");
        }

        let extract_done_marker = d01.join(".done");
        if !resume || !extract_done_marker.exists() {
            let extract_global_threads = args.threads.unwrap_or(1);

            let ext_args = extract::ExtractArgs {
                genome_list: args.samples.clone(),
                input: vec![],
                input_type: 2,
                enzyme_site: args.site.clone(),
                output_dir: d01.clone(),
                output_prefix: vec![],
                threads: extract_global_threads,
                quality_control: args.quality_control.clone(),
                max_n: args.max_n,
                min_quality: args.min_quality,
                min_quality_percent: args.min_quality_percent,
                quality_base: args.quality_base,
                pear_bin: args.pear_bin.clone().unwrap_or_else(|| "pear".to_string()),
                pear_threads: args.pear_threads,
                use_pear: args.use_pear.clone(),
            };
            extract::run(ext_args).context("extract step failed")?;
            std::fs::write(&extract_done_marker, b"ok")?;
        } else {
            tracing::info!("resume=on: skipping extract (already done)");
        }
    }

    // ======================================================
    // Step 2: Prepare qualitative database (02_db_qual)
    // ======================================================
    // Logic change: Database construction requires taxonomy information.
    // Prefer --taxonomy; fall back to --genome-list if not provided (assuming the same file contains taxonomy information)
    let effective_taxonomy_list = args.taxonomy.or(args.genome_list.clone());

    let qual_db_dir = if let Some(db) = args.database_dir.as_ref() {
        db.clone()
    } else if let Some(gl) = effective_taxonomy_list {
        tracing::info!("\n===== [2/5] Build qualitative database (build-qual-db) =====");
        let qual_done_marker = d02.join(".done");
        if !resume || !qual_done_marker.exists() {
            let qual_args = build_qual_db::BuildQualDbArgs {
                genome_list: gl.clone(), // pass in taxonomy list
                enzyme_site: args.site.clone(),
                taxonomy_levels: args.level.clone(),
                output_dir: d02.clone(),
                enzyme_file: None,
                // Use the updated pre-digested directory (may be pre_digested_output)
                pre_digested_dir: current_pre_digested_dir.clone(),
                remove_redundant: "yes".to_string(),
                threads: args.threads.unwrap_or(1),
            };
            build_qual_db::run(qual_args).context("build-qual-db step failed")?;
            std::fs::write(&qual_done_marker, b"ok")?;
        } else {
            tracing::info!("resume=on: skipping build-qual-db (already done)");
        }

        ensure_classify_file(&gl, &d02)?;

        d02.clone()
    } else {
        bail!("--taxonomy (or --genome-list) or --database not provided; cannot build database or run analysis");
    };

    if args.mode == "db-only" {
        tracing::info!("\n===== Mode db-only complete (qualitative database built: {}) =====", qual_db_dir.display());
        return Ok(());
    }

    // --------------------------------------------------------------------------
    //  Note: The following steps only run when mode != "db-only".
    //  Already validated that args.samples.is_some() above; safe to unwrap.
    // --------------------------------------------------------------------------
    let samples_path_buf = args.samples.as_ref().unwrap();

    let (step2_prefix, step3_prefix, step4_prefix) = if args.mode == "full" {
        ("2/5", "3/5", "4/5")
    } else {
        ("2a/5", "2b/5", "2c/5")
    };

    // Step 2a/2: Qualitative analysis
    tracing::info!("\n===== [{}] Qualitative analysis (qualitative) =====", step2_prefix);
    let qual_done_marker = d_qual.join(".done");
    if !resume || !qual_done_marker.exists() {
        let sample_list_for_qual = d_qual.join(format!("{}.samples.tsv", args.prefix));

        build_sample_list_for_quantify(samples_path_buf, &d01, &args.site, &sample_list_for_qual)?;

        let q_args = quantify::QuantifyArgs {
            sample_list: sample_list_for_qual,
            database_dir: qual_db_dir.clone(),
            taxonomy_level: args.level.clone(),
            enzyme_site: args.site.clone(),
            output_dir: d_qual.clone(),
            g_score_threshold: 0.0,
            verbose: "yes".to_string(),
            threads: args.threads.unwrap_or(1),
        };
        quantify::run(q_args).context("qualitative analysis step failed")?;
        std::fs::write(&qual_done_marker, b"ok")?;
    } else {
        tracing::info!("resume=on: skipping qualitative analysis (already done)");
    }

    // Step 2b/3: find-genome
    tracing::info!("\n===== [{}] Genome selection (find-genome) =====", step3_prefix);
    let find_genome_done_marker = d_sdb.join(".done");
    if !resume || !find_genome_done_marker.exists() {
        let fg_args = find_genome::FindGenomeArgs {
            sample_list: samples_path_buf.clone(),
            database_dir: qual_db_dir.clone(),
            output_dir: d_sdb.clone(),
            qual_dir: d_qual.clone(),
            g_score_threshold: args.g_score as i32,
            gcf_threshold: args.gcf_threshold,
            threads: args.threads.unwrap_or(1),
        };
        find_genome::run(fg_args).context("find-genome step failed")?;
        std::fs::write(&find_genome_done_marker, b"ok")?;
    } else {
        tracing::info!("resume=on: skipping find-genome (already done)");
    }

    // Step 2c/4: Build per-sample quantitative database + quantification
    tracing::info!("\n===== [{}] Build per-sample quantitative database (parallelized across samples) =====", step4_prefix);

    let samples_vec = read_sample_names(samples_path_buf)?;

    let enzyme_file = qual_db_dir.join(format!("{}.enzyme.iibdb", get_enzyme_name(&args.site)?));
    if !enzyme_file.exists() {
        bail!("Qualitative database enzyme file not found: {}; please verify database integrity", enzyme_file.display());
    }

    let all_quant_finished = AtomicBool::new(true);

    // Parallelized across samples
    samples_vec.par_iter().try_for_each(|sample_name| -> Result<()> {
        let sample_sdb_list = d_sdb.join(sample_name).join("sdb.list");
        if !sample_sdb_list.exists() {
            tracing::warn!("Warning: sample {} has no sdb.list (no genomes may have been selected); skipping quantitative database build and quantification", sample_name);
            all_quant_finished.store(false, Ordering::Relaxed);
            return Ok(());
        }

        let sample_db_dir = d_sdb.join(sample_name).join("database");
        let sample_db_done = sample_db_dir.join(".done");

        // 1. Build quantitative database
        if resume && sample_db_done.exists() {
            tracing::info!("resume=on: skipping quantitative database build for sample {} (already done)", sample_name);
        } else {
            std::fs::create_dir_all(&sample_db_dir)?;

            let classify_file = sample_db_dir.join("abfh_classify_with_speciename.txt.gz");
            if !classify_file.exists() {
                std::fs::copy(&sample_sdb_list, sample_db_dir.join("abfh_classify_with_speciename.txt"))?;
                Command::new("gzip")
                    .arg("-f")
                    .arg(sample_db_dir.join("abfh_classify_with_speciename.txt"))
                    .status()
                    .context("Failed to compress taxonomy file")?;
            }

            let quan_args = build_quan_db::BuildQuanDbArgs {
                genome_list: sample_sdb_list.clone(),
                enzyme_site: args.site.clone(),
                taxonomy_levels: args.level.clone(),
                output_dir: sample_db_dir.clone(),
                enzyme_file: Some(enzyme_file.clone()),
                pre_digested_dir: None,
                remove_redundant: "yes".to_string(),
                threads: args.threads.unwrap_or(1),
            };
            build_quan_db::run(quan_args)
                .with_context(|| format!("Failed to build quantitative database for sample {}", sample_name))?;
            std::fs::write(&sample_db_done, b"ok")?;
        }

        // 2. Quantitative analysis
        let quant_step_num = if args.mode == "full" { "5/5" } else { "4/5" };
        tracing::info!("\n===== [{}] Quantitative analysis (quantify) for {} =====", quant_step_num, sample_name);

        let sample_quant_dir = d04.join(sample_name);
        std::fs::create_dir_all(&sample_quant_dir)?;
        let sample_quant_done = sample_quant_dir.join(".done");

        if resume && sample_quant_done.exists() {
            tracing::info!("resume=on: skipping quantification for sample {} (already done)", sample_name);
        } else {
            let sample_list_file = sample_quant_dir.join(format!("{}.list.tsv", sample_name));
            let sample_iibsp = d01.join(format!("{}.{}.iibsp", sample_name, get_enzyme_name(&args.site)?));
            if !sample_iibsp.exists() {
                tracing::warn!("Warning: extraction output not found for sample {}; skipping quantification", sample_name);
                all_quant_finished.store(false, Ordering::Relaxed);
                return Ok(());
            }
            std::fs::write(&sample_list_file, format!("{}\t{}", sample_name, sample_iibsp.display()))?;

            let q_args = quantify::QuantifyArgs {
                sample_list: sample_list_file,
                database_dir: sample_db_dir,
                taxonomy_level: args.level.clone(),
                enzyme_site: args.site.clone(),
                output_dir: sample_quant_dir.clone(),
                g_score_threshold: args.g_score,
                verbose: "yes".to_string(),
                threads: args.threads.unwrap_or(1),
            };
            quantify::run(q_args)
                .with_context(|| format!("Quantification failed for sample {}", sample_name))?;
            std::fs::write(&sample_quant_done, b"ok")?;
        }
        Ok(())
    })?;

    if !all_quant_finished.load(Ordering::Relaxed) {
         tracing::warn!("Some samples skipped quantification due to missing sdb.list or extraction output.");
    }

    // Step 5/5: merge
    let step_num = "5/5";
    tracing::info!("\n===== [{}] Merge results (merge) =====", step_num);
    let merge_done_marker = d05.join(".done");
    if !resume || !merge_done_marker.exists() {
        let list_path = d05.join(format!("{}.merge_list.tsv", args.prefix));

        build_merge_list_from_sample_quantify(samples_path_buf, &d04, &args.site, &list_path)?;

        let m_args = merge::MergeArgs {
            sample_list: list_path,
            output_dir: d05.clone(),
            prefix: args.prefix.clone(),
            mock_samples: args.mock_samples.clone(),
            control_samples: args.control_samples.clone(),
        };
        merge::run(m_args).context("merge step failed")?;
        std::fs::write(&merge_done_marker, b"ok")?;
    } else {
        tracing::info!("resume=on: skipping merge (already done)");
    }

    // Step 6/6 (optional): Functional abundance prediction
    if let Some(ko_mapping_file) = &args.ko_mapping {
        tracing::info!("\n===== [6/6] Functional abundance prediction (predict) =====");
        let abundance_file = d05.join(format!("{}.all.xls", args.prefix));
        if !abundance_file.exists() {
            tracing::warn!(
                "Warning: species abundance file not found ({}); skipping functional prediction",
                abundance_file.display()
            );
        } else {
            let predict_done_marker = d05.join(".predict.done");
            if !resume || !predict_done_marker.exists() {
                let p_args = predict::PredictArgs {
                    abundance_file,
                    mapping_file: ko_mapping_file.clone(),
                    output_dir: d05.clone(),
                    prefix: args.prefix.clone(),
                };
                predict::run(p_args).context("functional prediction step failed")?;
                std::fs::write(&predict_done_marker, b"ok")?;
            } else {
                tracing::info!("resume=on: skipping functional prediction (already done)");
            }
        }
    }

    tracing::info!("\nPipeline complete: {}", args.outdir.display());
    Ok(())
}

fn ensure_classify_file(genome_list_path: &Path, output_dir: &Path) -> Result<()> {
    let classify_path = output_dir.join("abfh_classify_with_speciename.txt.gz");
    if classify_path.exists() { return Ok(()); }

    tracing::info!("Generating taxonomy mapping file: {}", classify_path.display());

    let file = File::open(genome_list_path)?;
    let reader = BufReader::new(file);

    let out_file = File::create(classify_path)?;
    let mut writer = flate2::write::GzEncoder::new(out_file, flate2::Compression::fast());

    let mut is_gtdb = false;
    let mut first = true;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }

        let parts: Vec<&str> = line.split('\t').collect();
        if first {
             if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                 is_gtdb = true;
                 if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; }
             }
             first = false;
        }

        if is_gtdb {
             if parts.len() < 2 { continue; }
             let gcf = extract_gcf_id(parts[0]);
             let tax = parts[1].replace(';', "\t");
             writeln!(writer, "{}\t{}", gcf, tax)?;
        } else {
             if parts.len() < 9 { continue; }
             let gcf = parts[0];
             let tax = parts[1..9].join("\t");
             writeln!(writer, "{}\t{}", gcf, tax)?;
        }
    }
    Ok(())
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

fn get_enzyme_name(site: &str) -> Result<String> {
     if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme.name.to_string()); }
     if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme.name.to_string()); } }
     bail!("Unknown enzyme: {}", site)
}

fn read_sample_names(list_path: &Path) -> Result<Vec<String>> {
    let file = File::open(list_path)?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if !parts.is_empty() {
            samples.push(parts[0].to_string());
        }
    }
    Ok(samples)
}

fn build_sample_list_for_quantify(original_list: &Path, extract_dir: &Path, site: &str, output_list: &Path) -> Result<()> {
    let file = File::open(original_list)?;
    let reader = BufReader::new(file);
    let mut writer = File::create(output_list)?;

    let enzyme_name = get_enzyme_name(site)?;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        let sample = parts[0];

        let iibsp_path = extract_dir.join(format!("{}.{}.iibsp", sample, enzyme_name));
        if iibsp_path.exists() {
            writeln!(writer, "{}\t{}", sample, iibsp_path.display())?;
        } else {
            tracing::warn!("Warning: extraction result not found for sample {}; skipping in list", sample);
        }
    }
    Ok(())
}

fn build_merge_list_from_sample_quantify(
    original_list: &Path,
    quant_dir: &Path,
    site: &str,
    output_list: &Path
) -> Result<()> {
    let file = File::open(original_list)?;
    let reader = BufReader::new(file);
    let mut writer = File::create(output_list)?;

    let enzyme_name = get_enzyme_name(site)?;

    let mut found_count = 0;
    let mut total_count = 0;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        let sample = parts[0];
        total_count += 1;

        let xls_path = quant_dir.join(sample).join(sample).join(format!("{}.{}.xls", sample, enzyme_name));

        if xls_path.exists() {
            writeln!(writer, "{}\t{}", sample, xls_path.display())?;
            found_count += 1;
        } else {
            tracing::warn!("Warning: quantification result not found for sample {} ({}); skipping in merge list", sample, xls_path.display());
        }
    }

    if found_count == 0 {
        tracing::warn!("Error: no valid quantification results found among {} samples", total_count);
    } else {
        tracing::info!("Merge list generated: found quantification results for {}/{} samples", found_count, total_count);
    }

    Ok(())
}
