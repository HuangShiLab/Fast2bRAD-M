use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;

use f2brad_core::io_utils::{open_compact_reader, CompactDatabaseWriter};

/// Convert a qualitative (possibly multi-GCF per tag) compact database into a
/// quantitative database where each tag is kept only if it maps to exactly one
/// GCF. This matches the 2bRAD-M "quan" database used for abundance estimation.
#[derive(Args, Debug)]
pub struct DedupDbArgs {
    /// Input qualitative .iibdb
    #[arg(short = 'i', long = "input", required = true)]
    pub input: PathBuf,

    /// Output quantitative .iibdb
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,
}

pub fn run(args: DedupDbArgs) -> Result<()> {
    // First pass: count how many GCFs carry each hash.
    let mut reader = open_compact_reader(&args.input)
        .with_context(|| format!("Failed to open input: {}", args.input.display()))?;
    let gcf_table = reader.gcf_table().to_vec();

    tracing::info!("Pass 1/2: counting tag occurrences in {}", args.input.display());
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut total_records = 0u64;
    while let Some((hash, _gcf_index)) = reader.next_record()? {
        *counts.entry(hash).or_insert(0) += 1;
        total_records += 1;
        if total_records % 100_000_000 == 0 {
            tracing::info!("  {} records counted", total_records);
        }
    }
    let unique_hashes = counts.values().filter(|&&c| c == 1).count();
    tracing::info!(
        "Pass 1 done. {} total records, {} unique hashes, {} single-GCF hashes",
        total_records,
        counts.len(),
        unique_hashes
    );

    // Second pass: write only single-GCF records.
    tracing::info!("Pass 2/2: writing unique tags to {}", args.output.display());
    let mut reader2 = open_compact_reader(&args.input)
        .with_context(|| format!("Failed to reopen input: {}", args.input.display()))?;
    let gcf_table2 = reader2.gcf_table().to_vec();
    if gcf_table != gcf_table2 {
        anyhow::bail!("GCF table mismatch between passes");
    }

    let out_file = File::create(&args.output)
        .with_context(|| format!("Failed to create output: {}", args.output.display()))?;
    let buf = BufWriter::new(out_file);
    let gcf_refs: Vec<&str> = gcf_table.iter().map(|s| s.as_str()).collect();
    let mut writer = CompactDatabaseWriter::new(buf, &gcf_refs)?;

    let mut written = 0u64;
    while let Some((hash, gcf_index)) = reader2.next_record()? {
        if counts.get(&hash).copied().unwrap_or(0) == 1 {
            writer.write_record(hash, gcf_index)?;
            written += 1;
            if written % 100_000_000 == 0 {
                tracing::info!("  {} unique records written", written);
            }
        }
    }
    writer.finish()?;
    tracing::info!(
        "Wrote {} single-GCF records ({} MB) to {}",
        written,
        std::fs::metadata(&args.output)?.len() / 1_000_000,
        args.output.display()
    );
    Ok(())
}
