use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "f2brad-host", version, about = "Host genotyping for fast2bRAD-holo")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// In-silico digest a reference genome and report tag-level statistics.
    Digest(f2brad_host::digest::DigestArgs),
    /// Cross-assignment collision analysis: human tags vs microbial genomes.
    Cross(f2brad_host::cross::CrossArgs),
}

fn main() -> Result<()> {
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();
    match cli.command {
        Commands::Digest(args) => f2brad_host::digest::run(args),
        Commands::Cross(args) => f2brad_host::cross::run(args),
    }
}
