mod build_qual_db;
mod build_quan_db;
mod enzymes;
mod extract;
mod io_utils;
mod merge;
mod quantify;
mod types;

use anyhow::Result;
use clap::{Parser, Subcommand};

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
    /// 数字酶切：根据酶位点从序列中提取 2bRAD 标签
    Extract(extract::ExtractArgs),
    /// 构建定性数据库：从参考基因组构建分类特异性标签数据库
    BuildQualDb(build_qual_db::BuildQualDbArgs),
    /// 构建定量数据库：只输出 unique 标签
    BuildQuanDb(build_quan_db::BuildQuanDbArgs),
    /// 丰度计算：计算样品中微生物的相对丰度
    Quantify(quantify::QuantifyArgs),

    /// 合并多样品丰度表
    Merge(merge::MergeArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Extract(args) => extract::run(args),
        Commands::BuildQualDb(args) => build_qual_db::run(args),
        Commands::BuildQuanDb(args) => build_quan_db::run(args),
        Commands::Quantify(args) => quantify::run(args),
        Commands::Merge(args) => merge::run(args),
    }
}
