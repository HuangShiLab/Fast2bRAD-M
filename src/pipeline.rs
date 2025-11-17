use anyhow::{bail, Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::build_qual_db;
use crate::build_quan_db;
use crate::extract;
use crate::find_genome;
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

    /// Find-genome 的 GCF 标签数阈值（用于 sample-only 模式）
    #[arg(long = "gcf", default_value = "1")]
    pub gcf_threshold: i32,
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
    let d_qual = args.outdir.join("qualitative");  // 定性分析结果目录
    let d_sdb = args.outdir.join("quantitative_sdb");  // 样品特异性定量数据库目录
    if args.mode != "db-only" {
        std::fs::create_dir_all(&d01)?;
        std::fs::create_dir_all(&d04)?;
        std::fs::create_dir_all(&d05)?;
    }
    if args.mode == "sample-only" {
        // sample-only 模式需要定性分析目录
        std::fs::create_dir_all(&d_qual)?;
        std::fs::create_dir_all(&d_sdb)?;
    }
    if args.mode != "sample-only" {
        // 数据库目录按需创建
        let _ = std::fs::create_dir_all(&d02);
        let _ = std::fs::create_dir_all(&d03);
    }

    let resume = args.resume.eq_ignore_ascii_case("yes");

    // Step 1: extract（batch, type=2）
    if args.mode != "db-only" {
        let step1_num = if args.mode == "sample-only" { "1/5" } else { "1/4" };
        println!("\n===== [{}] 提取样品标签 (extract) =====", step1_num);
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

            // build-quan-db（使用 02_db_qual 中已生成的 enzyme.fa.gz，避免重复生成）
            let quan_done_marker = d03.join(".done");
            if !resume || !quan_done_marker.exists() {
                // 使用 02_db_qual 中已生成的 enzyme.fa.gz
                let qual_enzyme_file = d02.join(format!("{}.enzyme.fa.gz", get_enzyme_name(&args.site)?));
                if !qual_enzyme_file.exists() {
                    bail!("定性数据库酶切文件不存在: {}，请先完成 build-qual-db", qual_enzyme_file.display());
                }
                
                let quan_args = build_quan_db::BuildQuanDbArgs {
                    genome_list: gl.clone(),
                    enzyme_site: args.site.clone(),
                    taxonomy_levels: args.level.clone(),
                    output_dir: d03.clone(),
                    enzyme_file: Some(qual_enzyme_file),
                    pre_digested_dir: None,  // 不再使用 pre_digested_dir，直接使用 02_db_qual 的文件
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
        // sample-only 模式必须提供 database（应该是 02_db_qual 定性数据库）
        let qual_db = args.database_dir
            .as_ref()
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("sample-only 模式必须提供 --database（应为 02_db_qual 定性数据库）"))?;
        
        // 在定性分析之前，检查并自动生成分类文件（若缺失且 genome-list 为 GTDB）
        let classify_path = qual_db.join("abfh_classify_with_speciename.txt.gz");
        if !classify_path.exists() {
            if let Some(gl) = args.genome_list.as_ref() {
                if is_gtdb_list(gl).unwrap_or(false) {
                    println!("分类文件缺失，检测到 GTDB 清单，尝试自动生成 ...");
                    ensure_convert_tool_and_run(gl, &qual_db)
                        .context("自动生成分类文件失败")?;
                }
            }
        }
        
        // Step 2a: 定性分析（使用 02_db_qual）
        println!("\n===== [2a/5] 定性分析 (qualitative) =====");
        let qual_done_marker = d_qual.join(".done");
        if !resume || !qual_done_marker.exists() {
            let sample_list_for_qual = d_qual.join(format!("{}.samples.tsv", args.prefix));
            build_sample_list_for_quantify(&args.samples, &d01, &args.site, &sample_list_for_qual)?;

            let q_args = quantify::QuantifyArgs {
                sample_list: sample_list_for_qual,
                database_dir: qual_db.clone(),
                taxonomy_level: args.level.clone(),
                enzyme_site: args.site.clone(),
                output_dir: d_qual.clone(),
                g_score_threshold: 0.0,  // 定性分析不对 G-score 过滤
                verbose: "yes".to_string(),
            };
            quantify::run(q_args).context("定性分析阶段失败")?;
            std::fs::write(&qual_done_marker, b"ok")?;
        } else {
            println!("resume=on: 跳过定性分析（已存在）");
        }

        // Step 2b: find-genome（根据定性结果筛选基因组）
        println!("\n===== [2b/5] 筛选基因组 (find-genome) =====");
        let find_genome_done_marker = d_sdb.join(".done");
        if !resume || !find_genome_done_marker.exists() {
            let fg_args = find_genome::FindGenomeArgs {
                sample_list: args.samples.clone(),
                database_dir: qual_db.clone(),
                output_dir: d_sdb.clone(),
                qual_dir: d_qual.clone(),
                g_score_threshold: args.g_score as i32,
                gcf_threshold: args.gcf_threshold,
            };
            find_genome::run(fg_args).context("find-genome 阶段失败")?;
            std::fs::write(&find_genome_done_marker, b"ok")?;
        } else {
            println!("resume=on: 跳过 find-genome（已存在）");
        }

        // Step 2c: 为每个样品构建特异性定量数据库
        println!("\n===== [2c/5] 构建样品特异性定量数据库 =====");
        let samples = read_sample_names(&args.samples)?;
        let enzyme_file = qual_db.join(format!("{}.{}.fa.gz", get_enzyme_name(&args.site)?, args.level));
        if !enzyme_file.exists() {
            bail!("定性数据库酶切文件不存在: {}", enzyme_file.display());
        }

        for sample_name in &samples {
            let sample_sdb_list = d_sdb.join(sample_name).join("sdb.list");
            if !sample_sdb_list.exists() {
                eprintln!("警告: 样品 {} 没有 sdb.list，跳过构建定量数据库", sample_name);
                continue;
            }

            let sample_db_dir = d_sdb.join(sample_name).join("database");
            let sample_db_done = sample_db_dir.join(".done");
            if resume && sample_db_done.exists() {
                println!("resume=on: 跳过样品 {} 的定量数据库构建（已存在）", sample_name);
                continue;
            }

            std::fs::create_dir_all(&sample_db_dir)?;
            
            // 复制 sdb.list 为分类文件
            let classify_file = sample_db_dir.join("abfh_classify_with_speciename.txt.gz");
            if !classify_file.exists() {
                std::fs::copy(&sample_sdb_list, sample_db_dir.join("abfh_classify_with_speciename.txt"))?;
                Command::new("gzip")
                    .arg("-f")
                    .arg(sample_db_dir.join("abfh_classify_with_speciename.txt"))
                    .status()
                    .context("压缩分类文件失败")?;
            }

            // 构建定量数据库
            let quan_args = build_quan_db::BuildQuanDbArgs {
                genome_list: sample_sdb_list,
                enzyme_site: args.site.clone(),
                taxonomy_levels: args.level.clone(),
                output_dir: sample_db_dir.clone(),
                enzyme_file: Some(enzyme_file.clone()),
                pre_digested_dir: None,
                remove_redundant: "yes".to_string(),  // 与 full 模式保持一致，去除基因组内冗余标签
            };
            build_quan_db::run(quan_args)
                .with_context(|| format!("构建样品 {} 的定量数据库失败", sample_name))?;
            std::fs::write(&sample_db_done, b"ok")?;
        }

        // 返回第一个样品的数据库目录（用于后续步骤，实际会为每个样品单独处理）
        qual_db
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
        if args.mode == "sample-only" {
            // sample-only 模式：为每个样品使用其特异性定量数据库
            println!("\n===== [3/5] 定量分析 (quantify) =====");
            let samples = read_sample_names(&args.samples)?;
            for sample_name in &samples {
                let sample_db_dir = d_sdb.join(sample_name).join("database");
                if !sample_db_dir.exists() {
                    eprintln!("警告: 样品 {} 没有定量数据库，跳过定量分析", sample_name);
                    continue;
                }

                let sample_quant_dir = d04.join(sample_name);
                std::fs::create_dir_all(&sample_quant_dir)?;
                let sample_quant_done = sample_quant_dir.join(".done");
                if resume && sample_quant_done.exists() {
                    println!("resume=on: 跳过样品 {} 的定量分析（已存在）", sample_name);
                    continue;
                }

                // 为单个样品生成列表
                let sample_list_file = sample_quant_dir.join(format!("{}.list.tsv", sample_name));
                let sample_iibsp = d01.join(format!("{}.{}.iibsp.gz", sample_name, get_enzyme_name(&args.site)?));
                if !sample_iibsp.exists() {
                    eprintln!("警告: 样品 {} 的提取产物不存在，跳过", sample_name);
                    continue;
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
                };
                quantify::run(q_args)
                    .with_context(|| format!("样品 {} 的定量分析失败", sample_name))?;
                std::fs::write(&sample_quant_done, b"ok")?;
            }
        } else {
            // full 模式：使用全局数据库
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
    }

    // Step 4: merge
    if args.mode != "db-only" {
        let step_num = if args.mode == "sample-only" { "4/5" } else { "4/4" };
        println!("\n===== [{}] 合并结果 (merge) =====", step_num);
        let merge_done_marker = d05.join(".done");
        if !resume || !merge_done_marker.exists() {
            let list_path = d05.join(format!("{}.merge_list.tsv", args.prefix));
            if args.mode == "sample-only" {
                // sample-only 模式：从每个样品的子目录中收集结果
                build_merge_list_from_sample_quantify(&d04, &list_path)?;
            } else {
                build_merge_list_from_quantify(&d04, &list_path)?;
            }

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

/// 读取样品列表中的样品名称
fn read_sample_names(list_file: &Path) -> Result<Vec<String>> {
    let file = File::open(list_file)
        .with_context(|| format!("无法读取样品列表: {}", list_file.display()))?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if !parts.is_empty() {
            samples.push(parts[0].to_string());
        }
    }
    Ok(samples)
}

/// 获取酶名称
fn get_enzyme_name(site: &str) -> Result<String> {
    use crate::enzymes::{enzyme_by_id, enzyme_by_name};
    if let Ok(site_num) = site.parse::<u8>() {
        if let Some(enzyme) = enzyme_by_id(site_num) {
            return Ok(enzyme.name.to_string());
        }
    }
    if let Some(enzyme) = enzyme_by_name(site) {
        return Ok(enzyme.name.to_string());
    }
    bail!("无效的酶名称或编号: {}", site);
}

/// 从 sample-only 模式的 quantify 结果目录生成 merge 所需列表
fn build_merge_list_from_sample_quantify(quant_dir: &Path, output_list: &Path) -> Result<()> {
    let mut entries = Vec::new();
    for entry in std::fs::read_dir(quant_dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }
        let sample = entry.file_name().to_string_lossy().to_string();
        // 递归查找样品目录中的 {sample}.{enzyme}.xls 文件（排除 GCF_detected.xls）
        let mut found_xls = None;
        
        // 递归查找函数
        fn find_xls_recursive(dir: &Path, sample: &str) -> Option<PathBuf> {
            if let Ok(entries) = std::fs::read_dir(dir) {
                for entry in entries {
                    if let Ok(entry) = entry {
                        let path = entry.path();
                        if path.is_dir() {
                            // 递归查找子目录
                            if let Some(found) = find_xls_recursive(&path, sample) {
                                return Some(found);
                            }
                        } else if let Some(file_name) = path.file_name().and_then(|n| n.to_str()) {
                            // 查找 {sample}.{enzyme}.xls 格式的文件，排除 GCF_detected.xls
                            if file_name.ends_with(".xls") 
                                && !file_name.contains("GCF_detected")
                                && file_name.starts_with(&format!("{}.", sample)) {
                                return Some(path);
                            }
                        }
                    }
                }
            }
            None
        }
        
        found_xls = find_xls_recursive(&entry.path(), &sample);
        
        if let Some(xls_path) = found_xls {
            entries.push((sample, xls_path));
        }
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

/// 从 quantify 结果目录生成 merge 所需列表（full 模式）
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


