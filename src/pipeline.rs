use anyhow::{bail, Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::build_qual_db;
use crate::build_quan_db;
use crate::extract;
use crate::merge;
use crate::quantify;

#[derive(Parser, Debug)]
pub struct PipelineArgs {
    /// 运行模式：full|db-only|sample-only
    #[arg(long = "mode", default_value = "full")]
    pub mode: String,

    /// 样品列表（TSV：sample_name<TAB>path1[<TAB>path2]）
    #[arg(long = "samples", alias = "list", short = 'l')]
    pub samples: PathBuf,

    /// 酶名称或编号（如 BcgI 或 1-16）
    #[arg(long = "site", short = 's')]
    pub site: String,

    /// 分类层级（kingdom|phylum|class|order|family|genus|species|strain）
    #[arg(long = "level", short = 't', default_value = "species")]
    pub level: String,

    /// 输出目录
    #[arg(long = "outdir", alias = "od")]
    pub outdir: PathBuf,

    /// 运行前缀（用于合并命名）
    #[arg(long = "prefix", default_value = "run")]
    pub prefix: String,

    /// 线程数（可选）
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    /// 失败时是否跳过已存在产物（yes/no）
    #[arg(long = "resume", default_value = "yes")]
    pub resume: String,

    /// 最低质量分（传给 extract）
    #[arg(long = "min-quality", default_value = "30")]
    pub min_quality: u8,

    /// 最低质量百分比（传给 extract）
    #[arg(long = "min-quality-percent", default_value = "80")]
    pub min_quality_percent: u8,

    /// 最大 N 比例（传给 extract）
    #[arg(long = "max-n", default_value = "0.08")]
    pub max_n: f64,

    /// 质量编码（传给 extract）
    #[arg(long = "quality-base", default_value = "33")]
    pub quality_base: u8,

    /// 预酶切目录（可选，用于构建数据库）
    #[arg(long = "pre-digested-dir")]
    pub pre_digested_dir: Option<PathBuf>,

    /// PEAR 可执行文件（透传给 extract）
    #[arg(long = "pe")]
    pub pear_bin: Option<String>,

    /// PEAR 线程数（透传给 extract）
    #[arg(long = "pc")]
    pub pear_cpu: Option<usize>,

    /// 基因组分类列表（可选，提供则自动构建数据库）
    #[arg(long = "genome-list")]
    pub genome_list: Option<PathBuf>,

    /// 已有数据库目录（可选，若提供则直接用于 quantify）
    #[arg(long = "database", short = 'd')]
    pub database_dir: Option<PathBuf>,

    /// Mock 样品（逗号分隔，可选）
    #[arg(long = "mock")]
    pub mock_samples: Option<String>,

    /// 阴性对照样品（逗号分隔，可选）
    #[arg(long = "control")]
    pub control_samples: Option<String>,

    /// Quantify 的 G-score 阈值
    #[arg(long = "gscore", short = 'g', default_value = "5")]
    pub g_score: f64,
}

pub fn run(args: PipelineArgs) -> Result<()> {
    // 可选：控制 rayon 线程
    if let Some(n) = args.threads {
        std::env::set_var("RAYON_NUM_THREADS", n.to_string());
    }

    // 创建目录结构
    let d01 = args.outdir.join("01_extract");
    let d02 = args.outdir.join("02_db_qual");
    let d03 = args.outdir.join("03_db_quan");
    let d04 = args.outdir.join("04_quantify");
    let d05 = args.outdir.join("05_merge");
    if args.mode != "db-only" {
        std::fs::create_dir_all(&d01)?;
        std::fs::create_dir_all(&d04)?;
        std::fs::create_dir_all(&d05)?;
    }
    if args.mode != "sample-only" {
        // 数据库目录按需创建
        let _ = std::fs::create_dir_all(&d02);
        let _ = std::fs::create_dir_all(&d03);
    }

    let resume = args.resume.eq_ignore_ascii_case("yes");

    // Step 1: extract（batch, type=2）
    if args.mode != "db-only" {
        println!("\n===== [1/4] 提取样品标签 (extract) =====");
        let extract_done_marker = d01.join(".done");
        if !resume || !extract_done_marker.exists() {
            let ext_args = extract::ExtractArgs {
                batch_file: Some(args.samples.clone()),
                input: vec![],                // 由 batch 读取
                input_type: 2,                // Type 2: Shotgun
                enzyme_site: args.site.clone(),
                output_dir: d01.clone(),
                output_prefix: vec![],        // 由 batch 生成
                compress: "yes".to_string(),
                quality_control: "yes".to_string(),
                max_n: args.max_n,
                min_quality: args.min_quality,
                min_quality_percent: args.min_quality_percent,
                quality_base: args.quality_base,
                format: "fa".to_string(),
            pear_bin: args.pear_bin.clone().unwrap_or_else(|| "pear".to_string()),
            pear_cpu: args.pear_cpu.unwrap_or(1),
            };
            extract::run(ext_args).context("extract 阶段失败")?;
            std::fs::write(&extract_done_marker, b"ok")?;
        } else {
            println!("resume=on: 跳过 extract（已存在）");
        }
    }

    // Step 2: 构建数据库（可选）
    let db_dir_in_use = if args.mode != "sample-only" {
        if let Some(gl) = args.genome_list.as_ref() {
            println!("\n===== [2/4] 构建数据库 (build-qual-db / build-quan-db) =====");
            std::fs::create_dir_all(&d02)?;
            std::fs::create_dir_all(&d03)?;

            // build-qual-db
            let qual_done_marker = d02.join(".done");
            if !resume || !qual_done_marker.exists() {
                let qual_args = build_qual_db::BuildQualDbArgs {
                    genome_list: gl.clone(),
                    enzyme_site: args.site.clone(),
                    taxonomy_levels: args.level.clone(),
                    output_dir: d02.clone(),
                    enzyme_file: None,
                    pre_digested_dir: args.pre_digested_dir.clone(),
                    remove_redundant: "yes".to_string(),
                };
                build_qual_db::run(qual_args).context("build-qual-db 阶段失败")?;
                std::fs::write(&qual_done_marker, b"ok")?;
            } else {
                println!("resume=on: 跳过 build-qual-db（已存在）");
            }

            // build-quan-db
            let quan_done_marker = d03.join(".done");
            if !resume || !quan_done_marker.exists() {
                let quan_args = build_quan_db::BuildQuanDbArgs {
                    genome_list: gl.clone(),
                    enzyme_site: args.site.clone(),
                    taxonomy_levels: args.level.clone(),
                    output_dir: d03.clone(),
                    enzyme_file: None,
                    pre_digested_dir: args.pre_digested_dir.clone(),
                    remove_redundant: "yes".to_string(),
                };
                build_quan_db::run(quan_args).context("build-quan-db 阶段失败")?;
                std::fs::write(&quan_done_marker, b"ok")?;
            } else {
                println!("resume=on: 跳过 build-quan-db（已存在）");
            }

            d03.clone()
        } else if let Some(db) = args.database_dir.as_ref() {
            // 直接使用已有数据库
            db.clone()
        } else {
            bail!("未提供 --genome-list 或 --database，无法进行 quantify");
        }
    } else {
        // sample-only 模式必须提供 database
        args.database_dir
            .as_ref()
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("sample-only 模式必须提供 --database"))?
    };

    // 自动生成分类文件（若缺失且 genome-list 为 GTDB）
    let classify_path = db_dir_in_use.join("abfh_classify_with_speciename.txt.gz");
    if !classify_path.exists() {
        if let Some(gl) = args.genome_list.as_ref() {
            if is_gtdb_list(gl).unwrap_or(false) {
                println!("分类文件缺失，检测到 GTDB 清单，尝试自动生成 ...");
                ensure_convert_tool_and_run(gl, &db_dir_in_use)
                    .context("自动生成分类文件失败")?;
            }
        }
    }

    // Step 3: quantify（基于 Step1 产物自动生成样品列表）
    if args.mode != "db-only" {
        println!("\n===== [3/4] 丰度分析 (quantify) =====");
        let quant_done_marker = d04.join(".done");
        if !resume || !quant_done_marker.exists() {
            let sample_list_for_quant = d04.join(format!("{}.samples.tsv", args.prefix));
            build_sample_list_for_quantify(&args.samples, &d01, &args.site, &sample_list_for_quant)?;

            let q_args = quantify::QuantifyArgs {
                sample_list: sample_list_for_quant,
                database_dir: db_dir_in_use.clone(),
                taxonomy_level: args.level.clone(),
                enzyme_site: args.site.clone(),
                output_dir: d04.clone(),
                g_score_threshold: args.g_score,
                verbose: "yes".to_string(),
            };
            quantify::run(q_args).context("quantify 阶段失败")?;
            std::fs::write(&quant_done_marker, b"ok")?;
        } else {
            println!("resume=on: 跳过 quantify（已存在）");
        }
    }

    // Step 4: merge
    if args.mode != "db-only" {
        println!("\n===== [4/4] 合并结果 (merge) =====");
        let merge_done_marker = d05.join(".done");
        if !resume || !merge_done_marker.exists() {
            let list_path = d05.join(format!("{}.merge_list.tsv", args.prefix));
            build_merge_list_from_quantify(&d04, &list_path)?;

            let m_args = merge::MergeArgs {
                sample_list: list_path,
                output_dir: d05.clone(),
                prefix: args.prefix.clone(),
                mock_samples: args.mock_samples.clone(),
                control_samples: args.control_samples.clone(),
            };
            merge::run(m_args).context("merge 阶段失败")?;
            std::fs::write(&merge_done_marker, b"ok")?;
        } else {
            println!("resume=on: 跳过 merge（已存在）");
        }
    }

    println!("\n流水线完成：{}", args.outdir.display());
    Ok(())
}

/// 判断 genome-list 是否为 GTDB 三列表
fn is_gtdb_list(list_path: &Path) -> Result<bool> {
    let f = File::open(list_path)
        .with_context(|| format!("无法打开清单文件：{}", list_path.display()))?;
    let reader = BufReader::new(f);
    for line in reader.lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = t.split('\t').collect();
        if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy" || parts[0] == "accession") {
            return Ok(true);
        } else {
            return Ok(false);
        }
    }
    Ok(false)
}

/// 调用 tools/convert_gtdb_taxonomy.py 生成分类文件到 db_dir
fn ensure_convert_tool_and_run(gtdb_list: &Path, db_dir: &Path) -> Result<()> {
    // 寻找脚本路径
    let mut candidates = Vec::new();
    if let Ok(exe) = std::env::current_exe() {
        if let Some(bin_dir) = exe.parent() {
            candidates.push(bin_dir.join("../../tools/convert_gtdb_taxonomy.py"));
            candidates.push(bin_dir.join("../tools/convert_gtdb_taxonomy.py"));
            candidates.push(bin_dir.join("tools/convert_gtdb_taxonomy.py"));
        }
    }
    candidates.push(PathBuf::from("tools/convert_gtdb_taxonomy.py"));

    let script = candidates.into_iter().find(|p| p.exists()).ok_or_else(|| {
        anyhow::anyhow!("找不到工具脚本 tools/convert_gtdb_taxonomy.py，请确认代码仓目录结构")
    })?;

    std::fs::create_dir_all(db_dir)?;
    let status = Command::new("python3")
        .args([
            script.to_str().unwrap(),
            "--input",
            gtdb_list.to_str().unwrap(),
            "--output-dir",
            db_dir.to_str().unwrap(),
        ])
        .status()
        .context("调用 python3 运行 convert_gtdb_taxonomy.py 失败")?;
    if !status.success() {
        bail!("convert_gtdb_taxonomy.py 运行失败，退出码 {:?}", status.code());
    }
    let out = db_dir.join("abfh_classify_with_speciename.txt.gz");
    if !out.exists() {
        bail!("转换后未找到输出文件：{}", out.display());
    }
    println!("已自动生成分类文件：{}", out.display());
    Ok(())
}

/// 根据原始样品列表与 extract 输出位置，生成 quantify 所需列表
fn build_sample_list_for_quantify(
    original_list: &Path,
    extract_dir: &Path,
    enzyme: &str,
    output_list: &Path,
) -> Result<()> {
    let file = File::open(original_list)
        .with_context(|| format!("无法读取样品列表: {}", original_list.display()))?;
    let reader = BufReader::new(file);

    let mut lines = Vec::new();
    for (ln, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 2 {
            bail!("样品列表第 {} 行格式错误（至少2列）", ln + 1);
        }
        let sample_name = parts[0];
        let sp_path = extract_dir.join(format!("{}.{}.iibsp.gz", sample_name, enzyme));
        if !sp_path.exists() {
            bail!("找不到样品提取产物：{}", sp_path.display());
        }
        lines.push(format!("{}\t{}", sample_name, sp_path.display()));
    }

    let mut w = File::create(output_list)?;
    for l in lines {
        writeln!(w, "{}", l)?;
    }
    Ok(())
}

/// 从 quantify 结果目录生成 merge 所需列表
fn build_merge_list_from_quantify(quant_dir: &Path, output_list: &Path) -> Result<()> {
    let mut entries = Vec::new();
    for entry in std::fs::read_dir(quant_dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }
        let sample = entry.file_name().to_string_lossy().to_string();
        let xls = entry.path().join(format!("{}.BcgI.xls", sample));
        // 允许任意酶名：匹配 *.xls 优先单一文件
        let xls_path = if xls.exists() {
            xls
        } else {
            // 回退扫描一个 .xls
            let mut found = None;
            for sub in std::fs::read_dir(entry.path())? {
                let sub = sub?;
                let p = sub.path();
                if p.extension().map(|e| e == "xls").unwrap_or(false) {
                    found = Some(p);
                    break;
                }
            }
            if let Some(p) = found { p } else { continue }
        };
        entries.push((sample, xls_path));
    }

    if entries.is_empty() {
        bail!("quantify 结果目录为空：{}", quant_dir.display());
    }

    let mut w = File::create(output_list)?;
    for (s, p) in entries {
        writeln!(w, "{}\t{}", s, p.display())?;
    }
    Ok(())
}


