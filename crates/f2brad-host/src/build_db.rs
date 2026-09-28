use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Args;
use flate2::read::GzDecoder;
use tracing;

use f2brad_core::enzymes::{enzyme_by_id, enzyme_by_name};
use f2brad_core::extract::Hash;

#[derive(Args, Debug)]
pub struct BuildDbArgs {
    /// Per-locus tag TSV produced by `f2brad-host digest` (sites file)
    #[arg(short = 't', long = "tags", required = true)]
    pub tags: PathBuf,

    /// Optional human mask list from `f2brad-host cross`: one hex hash per line.
    /// Tags whose canonical hash appears in this file are excluded from the DB.
    #[arg(short = 'm', long = "human-mask")]
    pub human_mask: Option<PathBuf>,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    pub enzyme_site: String,

    /// Output path for the host tag database (TSV)
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,
}

pub fn run(args: BuildDbArgs) -> Result<()> {
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };

    tracing::info!("Building host tag DB for {} from {}", enzyme.name, args.tags.display());

    let mask: HashSet<Hash> = if let Some(mask_path) = &args.human_mask {
        tracing::info!("Loading human mask from {}", mask_path.display());
        load_mask(mask_path)?
    } else {
        HashSet::new()
    };
    if !mask.is_empty() {
        tracing::info!("Loaded {} masked human tag hashes", mask.len());
    }

    let file = File::open(&args.tags)
        .with_context(|| format!("Failed to open tags file: {}", args.tags.display()))?;
    let reader: Box<dyn BufRead> = if args.tags.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create output directory: {}", parent.display()))?;
    }
    let out_file = File::create(&args.output)
        .with_context(|| format!("Failed to create output DB: {}", args.output.display()))?;
    let mut writer = BufWriter::new(out_file);

    // Copy the sites header unchanged.
    writeln!(writer, "contig\tpos\tstrand\tseq\tcanonical\thash\tgc_frac\tcpg_count\tcpg_island\tunique")?;

    let mut total = 0usize;
    let mut unique = 0usize;
    let mut kept = 0usize;
    let mut unique_masked = 0usize;
    let mut non_unique = 0usize;

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 10 {
            continue;
        }
        total += 1;

        let is_unique = parts[9] == "1";
        if is_unique {
            unique += 1;
        } else {
            non_unique += 1;
            continue;
        }

        let hash = u64::from_str_radix(parts[5], 16)
            .with_context(|| format!("Invalid hash on line {}: {}", i + 1, parts[5]))?;

        if mask.contains(&hash) {
            unique_masked += 1;
            continue;
        }

        writeln!(writer, "{}", line)?;
        kept += 1;
    }

    writer.flush()?;

    tracing::info!(
        "Host DB summary: {} input loci, {} unique within genome, {} non-unique skipped, {} unique excluded by cross-mask, {} kept",
        total,
        unique,
        non_unique,
        unique_masked,
        kept
    );
    if unique > 0 {
        tracing::info!(
            "Cross-mask removes {:.2}% of unique tags; usable panel = {} tags",
            unique_masked as f64 / unique as f64 * 100.0,
            kept
        );
    }
    tracing::info!("Wrote {}", args.output.display());

    Ok(())
}

fn load_mask(path: &PathBuf) -> Result<HashSet<Hash>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open mask file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut mask = HashSet::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let hash = u64::from_str_radix(line, 16)
            .with_context(|| format!("Invalid mask hash on line {}: {}", i + 1, line))?;
        mask.insert(hash);
    }
    Ok(mask)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tags(path: &PathBuf, lines: &[&str]) {
        let mut f = File::create(path).unwrap();
        writeln!(f, "contig\tpos\tstrand\tseq\tcanonical\thash\tgc_frac\tcpg_count\tcpg_island\tunique").unwrap();
        for line in lines {
            writeln!(f, "{}", line).unwrap();
        }
    }

    #[test]
    fn build_db_filters_non_unique_and_mask() {
        let tmp = std::env::temp_dir();
        let tags_path = tmp.join("f2host_build_db_test_tags.tsv");
        let mask_path = tmp.join("f2host_build_db_test_mask.txt");
        let out_path = tmp.join("f2host_build_db_test_out.tsv");

        write_tags(&tags_path, &[
            "c1\t1\t+\tAAAA\tAAAA\t0000000000000001\t0.0000\t0\t0\t1",
            "c1\t2\t+\tAAAT\tAAAT\t0000000000000002\t0.0000\t0\t0\t0",
            "c1\t3\t+\tACGT\tACGT\t0000000000000003\t0.0000\t0\t0\t1",
            "c2\t10\t-\tTGCA\tTGCA\t0000000000000004\t0.0000\t0\t0\t1",
        ]);
        {
            let mut f = File::create(&mask_path).unwrap();
            writeln!(f, "0000000000000003").unwrap();
        }

        let args = BuildDbArgs {
            tags: tags_path.clone(),
            human_mask: Some(mask_path.clone()),
            enzyme_site: "BcgI".to_string(),
            output: out_path.clone(),
        };
        run(args).unwrap();

        let content = std::fs::read_to_string(&out_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 kept loci
        assert!(lines.iter().any(|l| l.contains("0000000000000001")));
        assert!(lines.iter().any(|l| l.contains("0000000000000004")));
        assert!(!lines.iter().any(|l| l.contains("0000000000000002")));
        assert!(!lines.iter().any(|l| l.contains("0000000000000003")));

        std::fs::remove_file(&tags_path).ok();
        std::fs::remove_file(&mask_path).ok();
        std::fs::remove_file(&out_path).ok();
    }
}
