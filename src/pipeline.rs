use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*; 
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering}; 
use tracing;

use crate::build_qual_db;
use crate::build_quan_db;
use crate::extract;
use crate::find_genome;
use crate::merge;
use crate::predict;
use crate::quantify;
use crate::enzymes::{enzyme_by_id, enzyme_by_name};

// --- PipelineArgs 结构体 ---

#[derive(Parser, Debug)]
pub struct PipelineArgs {
    /// 运行模式：full|db-only|sample-only
    #[arg(long = "mode", default_value = "full")]
    pub mode: String,

    /// 样品列表（TSV：sample_name<TAB>path1[<TAB>path2]）
    /// 在 db-only 模式下可选，其他模式必选
    #[arg(long = "samples", alias = "list", short = 'l')]
    pub samples: Option<PathBuf>,

    /// 【参数重命名】基因组序列列表（TSV：genome_id<TAB>fasta_path）
    /// 仅在 db-only 模式下使用，用于批量提取参考基因组的 2bRAD 特征
    #[arg(long = "genome-list")] 
    pub genome_list: Option<PathBuf>,

    /// 【参数重命名】基因组分类信息列表（TSV: genome_id<TAB>taxonomy...）
    /// 用于构建数据库
    #[arg(long = "taxonomy")]
    pub taxonomy: Option<PathBuf>,

    /// 酶名称或编号（如 BcgI 或 1-16）
    #[arg(long = "site", short = 's')]
    pub site: String,

    /// 分类层级（kingdom|phylum|class|order|family|genus|species|strain）
    #[arg(long = "level", short = 't', default_value = "species")]
    pub level: String,

    /// 输出目录
    #[arg(long = "outdir", alias = "od")]
    pub outdir: PathBuf,

    /// 全局线程数
    #[arg(short = 'j', long = "threads")]
    pub threads: Option<usize>,

    /// G-score 阈值（默认 5）
    #[arg(long = "gscore", default_value = "5.0")]
    pub g_score: f64,

    /// GCF 标签数阈值（默认 1）
    #[arg(long = "gcf", default_value = "1")]
    pub gcf_threshold: i32,

    /// 断点续传
    #[arg(long = "resume", default_value = "no")]
    pub resume: String,
    
    // --- extract 相关 ---
    #[arg(long = "qc", default_value = "yes")]
    pub quality_control: String,
    #[arg(long = "max-n", default_value = "0.08")]
    pub max_n: f64,
    #[arg(long = "min-qual", default_value = "30")]
    pub min_quality: u8,
    #[arg(long = "min-qual-percent", default_value = "80")]
    pub min_quality_percent: u8,
    #[arg(long = "qual-base", default_value = "33")]
    pub quality_base: u8,
    #[arg(long = "pear-bin")]
    pub pear_bin: Option<String>,

    /// 每个 PEAR 进程使用的线程数
    #[arg(long = "pc", default_value = "1")]
    pub pear_threads: usize,

    /// 【帮助信息修改】是否使用 PEAR 拼接
    #[arg(long = "use-pear", default_value = "no", help = "是否使用 PEAR 拼接，值为yes/no，默认 no [default: no]，如果选择yes，则会严重降低分析速度")]
    pub use_pear: String,

    // --- build-db 相关 ---
    #[arg(long = "database")]
    pub database_dir: Option<PathBuf>,
    #[arg(long = "pre-digested-dir")]
    pub pre_digested_dir: Option<PathBuf>,

    // --- merge 相关 ---
    #[arg(long = "prefix", default_value = "Abundance_Stat")]
    pub prefix: String,
    #[arg(long = "mock")]
    pub mock_samples: Option<String>,
    #[arg(long = "control")]
    pub control_samples: Option<String>,

    // --- predict 相关 ---
    /// 物种-功能映射矩阵（TSV：第一列=物种名，其余列=功能ID，值=基因拷贝数）
    /// 提供此参数时，pipeline 在 merge 后自动进行功能丰度预测
    #[arg(long = "ko-mapping")]
    pub ko_mapping: Option<PathBuf>,
}

pub fn run(args: PipelineArgs) -> Result<()> {
    // 1. 显式初始化 Rayon 全局线程池
    if let Some(n) = args.threads {
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global() 
        {
            tracing::warn!("注意: Rayon 线程池已被初始化，--threads 参数可能未完全生效 (错误: {})", e);
        } else {
            tracing::info!("已设置全局并行线程数: {}", n);
        }
    }

    // 逻辑校验：非 db-only 模式必须提供样品列表
    if args.mode != "db-only" && args.samples.is_none() {
        bail!("错误: 运行模式为 '{}' 时，必须通过 --samples 或 -l 提供样品列表文件。", args.mode);
    }

    // 创建目录结构
    let d01 = args.outdir.join("01_extract");
    let d02 = args.outdir.join("02_db_qual");
    let d04 = args.outdir.join("04_quantify");
    let d05 = args.outdir.join("05_merge");
    let d_qual = args.outdir.join("qualitative"); 
    let d_sdb = args.outdir.join("quantitative_sdb"); 
    
    // db-only 模式下的预酶切输出目录
    let pre_digested_output = args.outdir.join("pre_digested_output");

    if args.mode != "db-only" {
        std::fs::create_dir_all(&d01)?;
        std::fs::create_dir_all(&d04)?;
        std::fs::create_dir_all(&d05)?;
        std::fs::create_dir_all(&d_qual)?;
        std::fs::create_dir_all(&d_sdb)?;
    }
    if args.mode != "sample-only" {
        let _ = std::fs::create_dir_all(&d02);
    }

    let resume = args.resume.eq_ignore_ascii_case("yes");

    // ======================================================
    // Step 0: [db-only 特有] 批量提取参考基因组特征 (Extract References)
    // ======================================================
    let mut current_pre_digested_dir = args.pre_digested_dir.clone();

    if args.mode == "db-only" {
        if let Some(genome_list_path) = &args.genome_list {
            tracing::info!("\n===== [db-only] 批量提取参考基因组特征 (Extract References) =====");
            std::fs::create_dir_all(&pre_digested_output)?;
            
            let extract_done_marker = pre_digested_output.join(".done");
            if !resume || !extract_done_marker.exists() {
                let extract_global_threads = args.threads.unwrap_or(1);
                
                // 配置 ExtractArgs 处理参考基因组 (Type 1)
                let ext_args = extract::ExtractArgs {
                    genome_list: Some(genome_list_path.clone()), // 传入序列列表
                    input: vec![],
                    input_type: 1, // ReferenceGenome
                    enzyme_site: args.site.clone(),
                    output_dir: pre_digested_output.clone(),
                    output_prefix: vec![],
                    threads: extract_global_threads,
                    // 对于参考基因组，通常不需要 read 级别的 QC，但保持参数一致性
                    quality_control: "no".to_string(), 
                    max_n: args.max_n,
                    min_quality: args.min_quality,
                    min_quality_percent: args.min_quality_percent,
                    quality_base: args.quality_base,
                    pear_bin:"pear".to_string(),
                    pear_threads: args.pear_threads,
                    use_pear: "no".to_string(),
                };
                
                extract::run(ext_args).context("db-only 批量提取参考基因组失败")?;
                std::fs::write(&extract_done_marker, b"ok")?;
            } else {
                tracing::info!("resume=on: 跳过参考基因组提取（已存在）");
            }

            // 更新后续 build-db 使用的目录
            current_pre_digested_dir = Some(pre_digested_output);
        }
    }

    // ======================================================
    // Step 1: 提取样品标签 (Extract Samples) - (非 db-only)
    // ======================================================
    if args.mode != "db-only" {
        let step_num = if args.mode == "sample-only" { "1/5" } else { "1/5" };
        tracing::info!("\n===== [{}] 提取样品标签 (extract) =====", step_num);
        tracing::info!("并发设置: {} 个样品并行, 每个样品使用 {} 个 PEAR 线程", 
            args.threads.unwrap_or(rayon::current_num_threads()), 
            args.pear_threads
        );

        let use_pear_bool = args.use_pear.eq_ignore_ascii_case("yes");
        if use_pear_bool {
            tracing::warn!("!!! 注意：已启用 PEAR 拼接 (--use-pear yes)。这会显著增加处理时间！");
        } else {
            tracing::info!("PEAR 拼接已禁用 (默认)。双端数据将分别提取后合并。");
        }

        let extract_done_marker = d01.join(".done");
        if !resume || !extract_done_marker.exists() {
            let extract_global_threads = args.threads.unwrap_or(1);

            let ext_args = extract::ExtractArgs {
                genome_list: args.samples.clone(),
                input: vec![],          
                input_type: 2,          
                enzyme_site: args.site.clone(),
                output_dir: d01.clone(),
                output_prefix: vec![],  
                threads: extract_global_threads,
                quality_control: args.quality_control.clone(),
                max_n: args.max_n,
                min_quality: args.min_quality,
                min_quality_percent: args.min_quality_percent,
                quality_base: args.quality_base,
                pear_bin: args.pear_bin.clone().unwrap_or_else(|| "pear".to_string()),
                pear_threads: args.pear_threads,
                use_pear: args.use_pear.clone(),
            };
            extract::run(ext_args).context("extract 阶段失败")?;
            std::fs::write(&extract_done_marker, b"ok")?;
        } else {
            tracing::info!("resume=on: 跳过 extract（已存在）");
        }
    }

    // ======================================================
    // Step 2: 准备定性数据库 (02_db_qual)
    // ======================================================
    // 【逻辑修改】建库需要分类信息 (taxonomy)。
    // 优先使用 --taxonomy 参数，如果未提供，尝试使用 --genome-list (假设同一个文件包含了分类信息)
    let effective_taxonomy_list = args.taxonomy.or(args.genome_list.clone());

    let qual_db_dir = if let Some(db) = args.database_dir.as_ref() {
        db.clone()
    } else if let Some(gl) = effective_taxonomy_list {
        tracing::info!("\n===== [2/5] 构建定性数据库 (build-qual-db) =====");
        let qual_done_marker = d02.join(".done");
        if !resume || !qual_done_marker.exists() {
            let qual_args = build_qual_db::BuildQualDbArgs {
                genome_list: gl.clone(), // 传入分类信息列表
                enzyme_site: args.site.clone(),
                taxonomy_levels: args.level.clone(),
                output_dir: d02.clone(),
                enzyme_file: None,
                // 使用更新后的预酶切目录 (可能是 pre_digested_output)
                pre_digested_dir: current_pre_digested_dir.clone(),
                remove_redundant: "yes".to_string(),
                threads: args.threads.unwrap_or(1),
            };
            build_qual_db::run(qual_args).context("build-qual-db 阶段失败")?;
            std::fs::write(&qual_done_marker, b"ok")?;
        } else {
            tracing::info!("resume=on: 跳过 build-qual-db（已存在）");
        }
        
        ensure_classify_file(&gl, &d02)?;

        d02.clone()
    } else {
        bail!("未提供 --taxonomy (或 --genome-list) 或 --database，无法进行建库或分析");
    };

    if args.mode == "db-only" {
        tracing::info!("\n===== 模式 db-only 完成（已构建定性数据库：{}）=====", qual_db_dir.display());
        return Ok(());
    }
    
    // --------------------------------------------------------------------------
    //  注意：以下步骤仅在 mode != "db-only" 时执行。
    //  前面已经校验过 args.samples.is_some()，因此可以安全解包。
    // --------------------------------------------------------------------------
    let samples_path_buf = args.samples.as_ref().unwrap(); 

    let (step2_prefix, step3_prefix, step4_prefix) = if args.mode == "full" { 
        ("2/5", "3/5", "4/5")
    } else { 
        ("2a/5", "2b/5", "2c/5") 
    };

    // Step 2a/2: 定性分析
    tracing::info!("\n===== [{}] 定性分析 (qualitative) =====", step2_prefix);
    let qual_done_marker = d_qual.join(".done");
    if !resume || !qual_done_marker.exists() {
        let sample_list_for_qual = d_qual.join(format!("{}.samples.tsv", args.prefix));
        
        build_sample_list_for_quantify(samples_path_buf, &d01, &args.site, &sample_list_for_qual)?;

        let q_args = quantify::QuantifyArgs {
            sample_list: sample_list_for_qual,
            database_dir: qual_db_dir.clone(), 
            taxonomy_level: args.level.clone(),
            enzyme_site: args.site.clone(),
            output_dir: d_qual.clone(),
            g_score_threshold: 0.0,  
            verbose: "yes".to_string(),
            threads: args.threads.unwrap_or(1),
        };
        quantify::run(q_args).context("定性分析阶段失败")?;
        std::fs::write(&qual_done_marker, b"ok")?;
    } else {
        tracing::info!("resume=on: 跳过定性分析（已存在）");
    }

    // Step 2b/3: find-genome
    tracing::info!("\n===== [{}] 筛选基因组 (find-genome) =====", step3_prefix);
    let find_genome_done_marker = d_sdb.join(".done");
    if !resume || !find_genome_done_marker.exists() {
        let fg_args = find_genome::FindGenomeArgs {
            sample_list: samples_path_buf.clone(), 
            database_dir: qual_db_dir.clone(), 
            output_dir: d_sdb.clone(),
            qual_dir: d_qual.clone(),
            g_score_threshold: args.g_score as i32,
            gcf_threshold: args.gcf_threshold,
            threads: args.threads.unwrap_or(1),
        };
        find_genome::run(fg_args).context("find-genome 阶段失败")?;
        std::fs::write(&find_genome_done_marker, b"ok")?;
    } else {
        tracing::info!("resume=on: 跳过 find-genome（已存在）");
    }

    // Step 2c/4: 构建样品特异性定量数据库 + 定量
    tracing::info!("\n===== [{}] 构建样品特异性定量数据库 (Sample-Parallel) =====", step4_prefix);
    
    let samples_vec = read_sample_names(samples_path_buf)?;
    
    let enzyme_file = qual_db_dir.join(format!("{}.enzyme.iibdb", get_enzyme_name(&args.site)?));
    if !enzyme_file.exists() {
        bail!("定性数据库酶切文件不存在: {}，请确认数据库完整性", enzyme_file.display());
    }

    let all_quant_finished = AtomicBool::new(true);

    // 样品间并行
    samples_vec.par_iter().try_for_each(|sample_name| -> Result<()> {
        let sample_sdb_list = d_sdb.join(sample_name).join("sdb.list");
        if !sample_sdb_list.exists() {
            tracing::warn!("警告: 样品 {} 没有 sdb.list（可能未筛选到基因组），跳过定量数据库构建和定量分析", sample_name);
            all_quant_finished.store(false, Ordering::Relaxed);
            return Ok(());
        }

        let sample_db_dir = d_sdb.join(sample_name).join("database");
        let sample_db_done = sample_db_dir.join(".done");
        
        // 1. 构建定量库
        if resume && sample_db_done.exists() {
            tracing::info!("resume=on: 跳过样品 {} 的定量数据库构建（已存在）", sample_name);
        } else {
            std::fs::create_dir_all(&sample_db_dir)?;
            
            let classify_file = sample_db_dir.join("abfh_classify_with_speciename.txt.gz");
            if !classify_file.exists() {
                std::fs::copy(&sample_sdb_list, sample_db_dir.join("abfh_classify_with_speciename.txt"))?;
                Command::new("gzip")
                    .arg("-f")
                    .arg(sample_db_dir.join("abfh_classify_with_speciename.txt"))
                    .status()
                    .context("压缩分类文件失败")?;
            }

            let quan_args = build_quan_db::BuildQuanDbArgs {
                genome_list: sample_sdb_list.clone(),
                enzyme_site: args.site.clone(),
                taxonomy_levels: args.level.clone(),
                output_dir: sample_db_dir.clone(),
                enzyme_file: Some(enzyme_file.clone()),
                pre_digested_dir: None,
                remove_redundant: "yes".to_string(),
                threads: args.threads.unwrap_or(1),
            };
            build_quan_db::run(quan_args)
                .with_context(|| format!("构建样品 {} 的定量数据库失败", sample_name))?;
            std::fs::write(&sample_db_done, b"ok")?;
        }
        
        // 2. 定量分析
        let quant_step_num = if args.mode == "full" { "5/5" } else { "4/5" };
        tracing::info!("\n===== [{}] 定量分析 (quantify) for {} =====", quant_step_num, sample_name);
        
        let sample_quant_dir = d04.join(sample_name);
        std::fs::create_dir_all(&sample_quant_dir)?;
        let sample_quant_done = sample_quant_dir.join(".done");
        
        if resume && sample_quant_done.exists() {
            tracing::info!("resume=on: 跳过样品 {} 的定量分析（已存在）", sample_name);
        } else {
            let sample_list_file = sample_quant_dir.join(format!("{}.list.tsv", sample_name));
            let sample_iibsp = d01.join(format!("{}.{}.iibsp", sample_name, get_enzyme_name(&args.site)?));
            if !sample_iibsp.exists() {
                tracing::warn!("警告: 样品 {} 的提取产物不存在，跳过定量分析", sample_name);
                all_quant_finished.store(false, Ordering::Relaxed);
                return Ok(());
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
                threads: args.threads.unwrap_or(1),
            };
            quantify::run(q_args)
                .with_context(|| format!("样品 {} 的定量分析失败", sample_name))?;
            std::fs::write(&sample_quant_done, b"ok")?;
        }
        Ok(())
    })?; 
    
    if !all_quant_finished.load(Ordering::Relaxed) {
         tracing::warn!("部分样品由于缺少 sdb.list 或提取产物，定量分析步骤未完成。");
    }

    // Step 5/5: merge
    let step_num = "5/5";
    tracing::info!("\n===== [{}] 合并结果 (merge) =====", step_num);
    let merge_done_marker = d05.join(".done");
    if !resume || !merge_done_marker.exists() {
        let list_path = d05.join(format!("{}.merge_list.tsv", args.prefix));
        
        build_merge_list_from_sample_quantify(samples_path_buf, &d04, &args.site, &list_path)?;

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
        tracing::info!("resume=on: 跳过 merge（已存在）");
    }

    // Step 6/6 (可选): 功能丰度预测
    if let Some(ko_mapping_file) = &args.ko_mapping {
        tracing::info!("\n===== [6/6] 功能丰度预测 (predict) =====");
        let abundance_file = d05.join(format!("{}.all.xls", args.prefix));
        if !abundance_file.exists() {
            tracing::warn!(
                "警告: 物种丰度文件不存在 ({})，跳过功能预测",
                abundance_file.display()
            );
        } else {
            let predict_done_marker = d05.join(".predict.done");
            if !resume || !predict_done_marker.exists() {
                let p_args = predict::PredictArgs {
                    abundance_file,
                    mapping_file: ko_mapping_file.clone(),
                    output_dir: d05.clone(),
                    prefix: args.prefix.clone(),
                };
                predict::run(p_args).context("功能预测阶段失败")?;
                std::fs::write(&predict_done_marker, b"ok")?;
            } else {
                tracing::info!("resume=on: 跳过功能预测（已存在）");
            }
        }
    }

    tracing::info!("\n流水线完成：{}", args.outdir.display());
    Ok(())
}

fn ensure_classify_file(genome_list_path: &Path, output_dir: &Path) -> Result<()> {
    let classify_path = output_dir.join("abfh_classify_with_speciename.txt.gz");
    if classify_path.exists() { return Ok(()); }

    tracing::info!("生成分类映射文件: {}", classify_path.display());
    
    let file = File::open(genome_list_path)?;
    let reader = BufReader::new(file);
    
    let out_file = File::create(classify_path)?;
    let mut writer = flate2::write::GzEncoder::new(out_file, flate2::Compression::fast());
    
    let mut is_gtdb = false;
    let mut first = true;
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if first {
             if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                 is_gtdb = true;
                 if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; }
             }
             first = false;
        }
        
        if is_gtdb {
             if parts.len() < 2 { continue; }
             let gcf = extract_gcf_id(parts[0]);
             let tax = parts[1].replace(';', "\t");
             writeln!(writer, "{}\t{}", gcf, tax)?;
        } else {
             if parts.len() < 9 { continue; }
             let gcf = parts[0];
             let tax = parts[1..9].join("\t");
             writeln!(writer, "{}\t{}", gcf, tax)?;
        }
    }
    Ok(())
}

fn extract_gcf_id(filename: &str) -> String {
    let name = filename.split('/').last().unwrap_or(filename);
    let name_clean = if let Some(pos) = name.find("_genomic") { &name[..pos] } else { name };
    if name_clean.starts_with("GCF_") || name_clean.starts_with("GCA_") {
        let parts: Vec<&str> = name_clean.split('_').collect();
        if parts.len() >= 2 { return format!("{}_{}", parts[0], parts[1]); }
    }
    name_clean.to_string()
}

fn get_enzyme_name(site: &str) -> Result<String> {
     if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme.name.to_string()); }
     if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme.name.to_string()); } }
     bail!("未知的酶：{}", site)
}

fn read_sample_names(list_path: &Path) -> Result<Vec<String>> {
    let file = File::open(list_path)?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if !parts.is_empty() {
            samples.push(parts[0].to_string());
        }
    }
    Ok(samples)
}

fn build_sample_list_for_quantify(original_list: &Path, extract_dir: &Path, site: &str, output_list: &Path) -> Result<()> {
    let file = File::open(original_list)?;
    let reader = BufReader::new(file);
    let mut writer = File::create(output_list)?;
    
    let enzyme_name = get_enzyme_name(site)?;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        let sample = parts[0];
        
        let iibsp_path = extract_dir.join(format!("{}.{}.iibsp", sample, enzyme_name));
        if iibsp_path.exists() {
            writeln!(writer, "{}\t{}", sample, iibsp_path.display())?;
        } else {
            tracing::warn!("警告: 样品 {} 的提取结果未找到，将在列表中跳过", sample);
        }
    }
    Ok(())
}

fn build_merge_list_from_sample_quantify(
    original_list: &Path, 
    quant_dir: &Path, 
    site: &str,
    output_list: &Path
) -> Result<()> {
    let file = File::open(original_list)?;
    let reader = BufReader::new(file);
    let mut writer = File::create(output_list)?;
    
    let enzyme_name = get_enzyme_name(site)?;

    let mut found_count = 0;
    let mut total_count = 0;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        let sample = parts[0];
        total_count += 1;
        
        let xls_path = quant_dir.join(sample).join(sample).join(format!("{}.{}.xls", sample, enzyme_name));
        
        if xls_path.exists() {
            writeln!(writer, "{}\t{}", sample, xls_path.display())?;
            found_count += 1;
        } else {
            tracing::warn!("警告: 样品 {} 的定量结果未找到 ({})，将在合并列表中跳过", sample, xls_path.display());
        }
    }

    if found_count == 0 {
        tracing::warn!("错误：在 {} 个样品中未找到任何有效的定量结果", total_count);
    } else {
        tracing::info!("合并列表生成完毕：找到 {}/{} 个样品的定量结果", found_count, total_count);
    }

    Ok(())
}