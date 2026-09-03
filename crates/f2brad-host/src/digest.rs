use std::path::PathBuf;

use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct DigestArgs {
    /// Reference genome FASTA/FASTQ (may be gzip-compressed)
    #[arg(short = 'i', long = "input", required = true)]
    pub input: PathBuf,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    pub enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    pub output_dir: PathBuf,

    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

pub fn run(args: DigestArgs) -> Result<()> {
    let _ = args;
    todo!("implement in-silico digest analysis")
}
