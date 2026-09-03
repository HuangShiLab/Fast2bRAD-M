use anyhow::{bail, Context, Result};
use clap::Parser;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

/// Functional abundance prediction: t(species abundance table) * species-to-function matrix = functional abundance table
#[derive(Parser, Debug, Clone)]
pub struct PredictArgs {
    // ── Input ──
    /// Species abundance table (.all.xls from merge step)
    #[arg(short = 'a', long = "abundance", help_heading = "Input")]
    pub abundance_file: PathBuf,
    /// Species-to-function mapping matrix (TSV: species name column + function ID columns with gene copy counts)
    #[arg(short = 'm', long = "mapping", help_heading = "Input")]
    pub mapping_file: PathBuf,

    // ── Output ──
    /// Output directory
    #[arg(short = 'o', long = "output", help_heading = "Output")]
    pub output_dir: PathBuf,
    /// Output file name prefix
    #[arg(short = 'p', long = "prefix", default_value = "Abundance_Stat", help_heading = "Output")]
    pub prefix: String,
}

pub fn run(args: PredictArgs) -> Result<()> {
    tracing::info!(
        "COMMAND: predict -a {} -m {} -o {} -p {}",
        args.abundance_file.display(),
        args.mapping_file.display(),
        args.output_dir.display(),
        args.prefix
    );

    std::fs::create_dir_all(&args.output_dir)?;

    // Step 1: Load species abundance matrix
    let (samples, species_abundance) = load_abundance_matrix(&args.abundance_file)?;
    tracing::info!(
        "Species abundance table loaded: {} species, {} samples",
        species_abundance.len(),
        samples.len()
    );
    if species_abundance.is_empty() {
        bail!("No valid data found in species abundance table");
    }

    // Step 2: Load species-to-function mapping matrix
    let (func_names, species_func_map) = load_mapping_matrix(&args.mapping_file)?;
    tracing::info!(
        "Functional mapping matrix loaded: {} functional entries, {} species",
        func_names.len(),
        species_func_map.len()
    );
    if func_names.is_empty() {
        bail!("No functional columns found in mapping matrix");
    }

    // Step 3: Compute functional abundance matrix
    // func_matrix[sample_idx][func_idx] = Σ_species (abundance[species][sample] × mapping[species][func])
    let n_samples = samples.len();
    let n_funcs = func_names.len();
    let mut func_matrix: Vec<Vec<f64>> = vec![vec![0.0; n_funcs]; n_samples];

    let mut matched = 0usize;
    for (species, abundance_vec) in &species_abundance {
        if let Some(func_vec) = species_func_map.get(species.as_str()) {
            matched += 1;
            for (s_idx, &abund) in abundance_vec.iter().enumerate() {
                if abund == 0.0 {
                    continue;
                }
                for (f_idx, &count) in func_vec.iter().enumerate() {
                    func_matrix[s_idx][f_idx] += abund * count;
                }
            }
        }
    }

    tracing::info!(
        "Species matched: {}/{} species found in mapping matrix",
        matched,
        species_abundance.len()
    );

    if matched == 0 {
        tracing::warn!(
            "Warning: no species matched between abundance table and mapping matrix, please check that species names are consistent"
        );
    }

    // Step 3.5: Per-sample normalization — ensure each sample's total functional abundance sums to 1.0
    for s_idx in 0..n_samples {
        let total: f64 = func_matrix[s_idx].iter().sum();
        if total > 0.0 {
            for k_idx in 0..n_funcs {
                func_matrix[s_idx][k_idx] /= total;
            }
        }
    }
    tracing::info!("Functional abundance matrix normalized per-sample (each sample sums to 1.0)");

    // Step 4: Write functional abundance table
    let output_path = args.output_dir.join(format!("{}.func.xls", args.prefix));
    let file = File::create(&output_path)
        .with_context(|| format!("Cannot open/create output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // Header: #Function\tsample1\t...\tsampleN
    write!(writer, "#Function")?;
    for sample in &samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Data rows: one row per functional entry, skip all-zero rows
    let mut written = 0usize;
    for (f_idx, func_name) in func_names.iter().enumerate() {
        let has_nonzero = (0..n_samples).any(|s| func_matrix[s][f_idx] != 0.0);
        if !has_nonzero {
            continue;
        }
        write!(writer, "{}", func_name)?;
        for s_idx in 0..n_samples {
            write!(writer, "\t{:.8}", func_matrix[s_idx][f_idx])?;
        }
        writeln!(writer)?;
        written += 1;
    }

    tracing::info!(
        "Functional abundance table written: {} non-zero functional entries -> {}",
        written,
        output_path.display()
    );
    Ok(())
}

/// Load species abundance matrix (merge output .all.xls format).
///
/// File format (optional '#' prefix in header line):
///   Kingdom\tPhylum\t...\tSpecies\tsample1\tsample2\t...
///   Bacteria\t...\tSpecies_name\t0.66\t0.34\t...
///
/// Returns: (sample name list, species name -> per-sample abundance values)
/// Species name is taken from the last taxonomy column (usually Species or Strain).
fn load_abundance_matrix(path: &Path) -> Result<(Vec<String>, FxHashMap<String, Vec<f64>>)> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open species abundance file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Parse header, find the boundary between taxonomy columns and sample columns
    let header_raw = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Species abundance file is empty: {}", path.display()))??;
    let header_str = header_raw.trim().trim_start_matches('#');
    let headers: Vec<&str> = header_str.split('\t').collect();

    // Known taxonomy level names (case-sensitive, matching merge output)
    const TAX_LEVELS: &[&str] = &[
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain",
    ];

    // Find the last taxonomy column index; samples start at column sample_start
    let mut sample_start = 0usize;
    for (i, h) in headers.iter().enumerate() {
        if TAX_LEVELS.contains(h) {
            sample_start = i + 1;
        } else if sample_start > 0 {
            // Encountered non-taxonomy column after at least one taxonomy column; stop
            break;
        }
    }

    if sample_start == 0 || sample_start >= headers.len() {
        bail!(
            "Cannot parse header of species abundance file, taxonomy or sample columns not found: {}",
            path.display()
        );
    }

    let samples: Vec<String> = headers[sample_start..]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let n_samples = samples.len();

    let mut species_abundance: FxHashMap<String, Vec<f64>> = FxHashMap::default();

    for line_res in lines {
        let line = line_res?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < sample_start + 1 {
            continue;
        }
        // Species name taken from the last taxonomy column
        let species_name = fields[sample_start - 1].trim().to_string();
        let abundances: Vec<f64> = fields[sample_start..]
            .iter()
            .take(n_samples)
            .map(|v| v.parse::<f64>().unwrap_or(0.0))
            .collect();
        species_abundance.insert(species_name, abundances);
    }

    Ok((samples, species_abundance))
}

/// Load species-to-function mapping matrix.
///
/// File format (TSV, optional '#' prefix in header line):
///   #Species\tKO1\tKO2\t...\tKOn
///   Cutibacterium_acnes\t5\t0\t3\t...
///   Escherichia_coli\t2\t8\t0\t...
///
/// Returns: (function name list, species name -> per-function counts)
fn load_mapping_matrix(
    path: &Path,
) -> Result<(Vec<String>, FxHashMap<String, Vec<f64>>)> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open functional mapping matrix file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_raw = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Functional mapping matrix file is empty: {}", path.display()))??;
    let header_str = header_raw.trim().trim_start_matches('#');
    let headers: Vec<&str> = header_str.split('\t').collect();

    if headers.len() < 2 {
        bail!(
            "Functional mapping matrix format error: at least 2 columns required (first column is species name column, remaining columns are function IDs): {}",
            path.display()
        );
    }

    // First column is the species name column; the rest are function IDs
    let func_names: Vec<String> = headers[1..].iter().map(|s| s.to_string()).collect();
    let n_funcs = func_names.len();

    let mut species_map: FxHashMap<String, Vec<f64>> = FxHashMap::default();

    for line_res in lines {
        let line = line_res?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        let species = fields[0].trim().to_string();
        let mut counts: Vec<f64> = fields[1..]
            .iter()
            .map(|v| v.parse::<f64>().unwrap_or(0.0))
            .collect();
        // Pad to function column count (fault-tolerant)
        counts.resize(n_funcs, 0.0);
        species_map.insert(species, counts);
    }

    Ok((func_names, species_map))
}
