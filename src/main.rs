mod build_qual_db;
mod build_quan_db;
mod enzymes;
mod extract;
mod find_genome;
mod inspect;
mod io_utils;
mod merge;
mod pipeline;
mod predict;
mod quantify;
mod types;

use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing_subscriber;

use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

#[derive(Parser, Debug)]
#[command(
    name = "fast2bRAD-M",
    version,
    about = "Rust rewrite of the 2bRAD-M extract workflow"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// In-silico restriction digestion: extract 2bRAD tags from sequences at enzyme sites
    Extract(extract::ExtractArgs),
    /// Build qualitative database: build a taxon-specific tag database from reference genomes
    BuildQualDb(build_qual_db::BuildQualDbArgs),
    /// Build quantitative database: output unique tags only
    BuildQuanDb(build_quan_db::BuildQuanDbArgs),
    /// Abundance quantification: compute the relative abundance of microbes in a sample
    Quantify(quantify::QuantifyArgs),

    /// Merge multi-sample abundance tables
    Merge(merge::MergeArgs),

    /// Select genomes for quantification based on qualitative results
    FindGenome(find_genome::FindGenomeArgs),

    /// Inspect .iibdb / .iibsp binary files: show format, tag counts and example records
    Inspect(inspect::InspectArgs),

    /// Functional abundance prediction: t(species abundance table) × species function matrix = functional abundance table
    Predict(predict::PredictArgs),

    /// One-command pipeline: extract → build-db → quantify → merge → predict
    Pipeline(pipeline::PipelineArgs),
}

fn main() -> Result<()> {
    // Create a non-blocking writer that outputs to stdout
    let (non_blocking, _guard) = tracing_appender::non_blocking(std::io::stdout());

    tracing_subscriber::fmt()
        .with_writer(non_blocking) // Use async writer
        .with_target(false)
        .with_thread_ids(false)
        .init();

    let cli = Cli::parse();
    match cli.command {
        Commands::Extract(args) => extract::run(args),
        Commands::BuildQualDb(args) => build_qual_db::run(args),
        Commands::BuildQuanDb(args) => build_quan_db::run(args),
        Commands::Quantify(args) => quantify::run(args),
        Commands::Merge(args) => merge::run(args),
        Commands::FindGenome(args) => find_genome::run(args),
        Commands::Inspect(args) => inspect::run(args),
        Commands::Predict(args) => predict::run(args),
        Commands::Pipeline(args) => pipeline::run(args),
    }
}
