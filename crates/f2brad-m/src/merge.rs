use anyhow::{Context, Result, bail};
use clap::Parser;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

/// Merge multi-sample abundance tables
#[derive(Parser, Debug)]
pub struct MergeArgs {
    // ── Input / Output ──
    /// Sample list file (TSV: sample_name<TAB>path_to_{sample}.{enzyme}.xls)
    #[arg(short = 'l', long = "list", help_heading = "Input / Output")]
    pub sample_list: PathBuf,
    /// Output directory
    #[arg(short = 'o', long = "output", help_heading = "Input / Output")]
    pub output_dir: PathBuf,
    /// Output file name prefix
    #[arg(short = 'p', long = "prefix", default_value = "Abundance_Stat", help_heading = "Input / Output")]
    pub prefix: String,

    // ── Filtering (for generating filtered output) ──
    /// Mock community sample names (comma-separated). These samples are excluded from the filtered output
    #[arg(short = 'm', long = "mock", help_heading = "Filtering")]
    pub mock_samples: Option<String>,
    /// Negative control sample names (comma-separated). These samples are excluded, and all taxa detected in controls are removed as contamination
    #[arg(short = 'c', long = "control", help_heading = "Filtering")]
    pub control_samples: Option<String>,
}

/// Species abundance data
#[derive(Debug, Default)]
struct TaxonAbundance {
    /// sample name -> relative abundance
    samples: FxHashMap<String, f64>,
}

pub fn run(args: MergeArgs) -> Result<()> {
    tracing::info!(
        "COMMAND: merge -l {} -o {} -p {}",
        args.sample_list.display(),
        args.output_dir.display(),
        args.prefix
    );

    // creating output directory
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("Cannot create output directory: {}", args.output_dir.display()))?;

    // parsing mock and control sample lists
    let mock_set: FxHashMap<String, ()> = if let Some(ref mock) = args.mock_samples {
        mock.split(',').map(|s| (s.trim().to_string(), ())).collect()
    } else {
        FxHashMap::default()
    };

    let control_set: FxHashMap<String, ()> = if let Some(ref control) = args.control_samples {
        control.split(',').map(|s| (s.trim().to_string(), ())).collect()
    } else {
        FxHashMap::default()
    };

    // reading abundance data for all samples
    let (taxa_abundance, sample_order, header) = read_all_profiles(&args.sample_list)?;

    if sample_order.is_empty() {
        bail!("No valid sample abundance files found");
    }

    tracing::info!("Read {} samples in total", sample_order.len());

    // write merged abundance table (all.xls)
    let all_output = args.output_dir.join(format!("{}.all.xls", args.prefix));
    write_merged_table(&all_output, &taxa_abundance, &sample_order, &header)?;
    tracing::info!("Merged table written: {}", all_output.display());

    // write filtered abundance table (filtered.xls)
    let filtered_output = args.output_dir.join(format!("{}.filtered.xls", args.prefix));
    write_filtered_table(
        &filtered_output,
        &taxa_abundance,
        &sample_order,
        &header,
        &mock_set,
        &control_set,
    )?;
    tracing::info!("Filtered table written: {}", filtered_output.display());

    tracing::info!("\nAll done!");
    Ok(())
}

/// Read abundance data for all samples
fn read_all_profiles(
    list_file: &Path,
) -> Result<(FxHashMap<String, TaxonAbundance>, Vec<String>, String)> {
    let file = File::open(list_file)
        .with_context(|| format!("Cannot open sample list: {}", list_file.display()))?;
    let reader = BufReader::new(file);

    let mut taxa_abundance: FxHashMap<String, TaxonAbundance> = FxHashMap::default();
    let mut sample_order = Vec::new();
    let mut sample_totals: FxHashMap<String, f64> = FxHashMap::default();
    let mut header = String::new();
    let mut classify_col = 0;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        let sample_name = parts[0].to_string();
        let profile_path = Path::new(parts[1]);

        if !profile_path.exists() {
            tracing::warn!("Warning: abundance file for sample {} not found: {}", sample_name, profile_path.display());
            continue;
        }

        sample_order.push(sample_name.clone());

        // read sample abundance file
        let profile_file = File::open(profile_path)
            .with_context(|| format!("Cannot open abundance file: {}", profile_path.display()))?;
        let profile_reader = BufReader::new(profile_file);

        for line in profile_reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.starts_with("#Kingdom") || line.starts_with("Kingdom") {
                // parse header to determine the range of taxonomy columns
                let fields: Vec<&str> = line.trim_start_matches('#').split('\t').collect();
                for (i, field) in fields.iter().enumerate() {
                    if *field == "Theoretical_Tag_Num" {
                        classify_col = i - 1;
                        header = fields[0..=classify_col].join("\t");
                        break;
                    }
                }
                continue;
            }

            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() <= classify_col {
                continue;
            }

            // extract taxonomy ID
            let taxon_id = fields[0..=classify_col].join("\t");

            // extract Sequenced_Reads_Num/Theoretical_Tag_Num (4th column from the end)
            if fields.len() >= 4 {
                let abundance_value: f64 = fields[fields.len() - 4]
                    .parse()
                    .unwrap_or(0.0);

                // record the abundance value of this taxon in this sample
                taxa_abundance
                    .entry(taxon_id.clone())
                    .or_insert_with(TaxonAbundance::default)
                    .samples
                    .insert(sample_name.clone(), abundance_value);

                // accumulate sample total
                *sample_totals.entry(sample_name.clone()).or_insert(0.0) += abundance_value;
            }
        }
    }

    // normalization: compute relative abundance (each taxon as a fraction of sample total)
    for (_taxon_id, abundance) in taxa_abundance.iter_mut() {
        for (sample_name, value) in abundance.samples.iter_mut() {
            if let Some(&total) = sample_totals.get(sample_name) {
                if total > 0.0 {
                    *value = *value / total;
                }
            }
        }
    }

    Ok((taxa_abundance, sample_order, header))
}

/// Write merged abundance table (all.xls)
fn write_merged_table(
    output_path: &Path,
    taxa_abundance: &FxHashMap<String, TaxonAbundance>,
    sample_order: &[String],
    header: &str,
) -> Result<()> {
    let file = File::create(output_path)
        .with_context(|| format!("Cannot create output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // write header
    write!(writer, "{}", header)?;
    for sample in sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // write data (sorted by taxon ID)
    let mut taxa_list: Vec<&String> = taxa_abundance.keys().collect();
    taxa_list.sort();

    for taxon_id in taxa_list {
        let abundance = &taxa_abundance[taxon_id];

        // check whether at least one sample contains this taxon
        let has_data = sample_order.iter().any(|sample| {
            abundance.samples.get(sample).map_or(false, |&v| v > 0.0)
        });

        if !has_data {
            continue;
        }

        write!(writer, "{}", taxon_id)?;
        for sample in sample_order {
            let value = abundance.samples.get(sample).copied().unwrap_or(0.0);
            if value == 0.0 {
                write!(writer, "\t0")?;
            } else {
                write!(writer, "\t{}", value)?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Write filtered abundance table (filtered.xls)
fn write_filtered_table(
    output_path: &Path,
    taxa_abundance: &FxHashMap<String, TaxonAbundance>,
    sample_order: &[String],
    header: &str,
    mock_set: &FxHashMap<String, ()>,
    control_set: &FxHashMap<String, ()>,
) -> Result<()> {
    let file = File::create(output_path)
        .with_context(|| format!("Cannot create output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // filter sample list (remove mock and control)
    let filtered_samples: Vec<String> = sample_order
        .iter()
        .filter(|s| !mock_set.contains_key(*s) && !control_set.contains_key(*s))
        .cloned()
        .collect();

    if filtered_samples.is_empty() {
        tracing::warn!("Warning: no samples remaining after filtering");
        // create empty file
        return Ok(());
    }

    // collect taxa detected in control samples (potential contamination)
    let mut contamination_taxa: FxHashMap<String, ()> = FxHashMap::default();
    for (taxon_id, abundance) in taxa_abundance.iter() {
        for control_sample in control_set.keys() {
            if abundance.samples.get(control_sample).map_or(false, |&v| v > 0.0) {
                contamination_taxa.insert(taxon_id.clone(), ());
                break;
            }
        }
    }

    // write header
    write!(writer, "{}", header)?;
    for sample in &filtered_samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // write data
    let mut taxa_list: Vec<&String> = taxa_abundance.keys().collect();
    taxa_list.sort();

    for taxon_id in taxa_list {
        // skip taxa flagged as potential contamination
        if contamination_taxa.contains_key(taxon_id) {
            continue;
        }

        let abundance = &taxa_abundance[taxon_id];

        // check whether any filtered sample has data
        let has_data = filtered_samples.iter().any(|sample| {
            abundance.samples.get(sample).map_or(false, |&v| v > 0.0)
        });

        if !has_data {
            continue;
        }

        write!(writer, "{}", taxon_id)?;
        for sample in &filtered_samples {
            let value = abundance.samples.get(sample).copied().unwrap_or(0.0);
            if value == 0.0 {
                write!(writer, "\t0")?;
            } else {
                write!(writer, "\t{}", value)?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}
