use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use fxhash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// 根据定性结果筛选定量基因组
#[derive(Parser, Debug)]
pub struct FindGenomeArgs {
    /// 样品列表文件（TSV格式：sample_name<tab>...）
    #[arg(short = 'l', long = "list")]
    pub sample_list: PathBuf,

    /// 数据库目录
    #[arg(short = 'd', long = "database")]
    pub database_dir: PathBuf,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// 定性结果目录
    #[arg(long = "qual-dir", alias = "qualdir")]
    pub qual_dir: PathBuf,

    /// G-score 阈值（默认 5，表示 >5）
    #[arg(long = "gscore", default_value = "5")]
    pub g_score_threshold: i32,

    /// GCF 标签数阈值（默认 1，表示 >1）
    #[arg(long = "gcf", default_value = "1")]
    pub gcf_threshold: i32,
}

/// 主函数
pub fn run(args: FindGenomeArgs) -> Result<()> {
    println!(
        "COMMAND: find-genome -l {} -d {} -o {} --qual-dir {} --gscore {} --gcf {}",
        args.sample_list.display(),
        args.database_dir.display(),
        args.output_dir.display(),
        args.qual_dir.display(),
        args.g_score_threshold,
        args.gcf_threshold
    );

    // 检查数据库文件
    let classify_file = args.database_dir.join("abfh_classify_with_speciename.txt.gz");
    if !classify_file.exists() {
        bail!(
            "数据库文件不存在: {}",
            classify_file.display()
        );
    }

    // 加载 GCF 到完整分类信息的映射
    let gcf_to_classify = load_gcf_classify(&classify_file)?;
    println!("已加载 {} 个基因组分类信息", gcf_to_classify.len());

    // 创建输出目录
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir.display()))?;

    // 读取样品列表
    let samples = read_sample_list(&args.sample_list)?;
    println!("共 {} 个样品需要处理", samples.len());

    // 处理每个样品
    for sample_name in samples {
        process_sample(
            &sample_name,
            &args.qual_dir,
            &args.database_dir,
            &args.output_dir,
            &gcf_to_classify,
            args.g_score_threshold,
            args.gcf_threshold,
        )?;
    }

    println!("\n全部完成！");
    Ok(())
}

/// 加载 GCF 到完整分类信息的映射
fn load_gcf_classify(classify_file: &Path) -> Result<FxHashMap<String, String>> {
    let mut gcf_to_classify = FxHashMap::default();

    let file = File::open(classify_file)
        .with_context(|| format!("无法打开数据库文件: {}", classify_file.display()))?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.is_empty() {
            continue;
        }

        let gcf_id = parts[0].to_string();
        gcf_to_classify.insert(gcf_id, line.to_string());
    }

    Ok(gcf_to_classify)
}

/// 读取样品列表
fn read_sample_list(list_file: &Path) -> Result<Vec<String>> {
    let file = File::open(list_file)
        .with_context(|| format!("无法打开样品列表: {}", list_file.display()))?;
    let reader = BufReader::new(file);

    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.is_empty() {
            continue;
        }

        samples.push(parts[0].to_string());
    }

    Ok(samples)
}

/// 处理单个样品
fn process_sample(
    sample_name: &str,
    qual_dir: &Path,
    database_dir: &Path,
    output_dir: &Path,
    gcf_to_classify: &FxHashMap<String, String>,
    g_score_threshold: i32,
    gcf_threshold: i32,
) -> Result<()> {
    // 读取 combine.xls 文件
    let combine_file = qual_dir.join(sample_name).join(format!("{}.combine.xls", sample_name));
    
    if !combine_file.exists() {
        eprintln!(
            "!!! {} 没有定性结果文件: {}，跳过定量分析",
            sample_name,
            combine_file.display()
        );
        return Ok(());
    }

    // 解析 combine.xls，提取使用的酶和通过 G-score 阈值的分类
    let (enzymes, pass_gscore_classes) = parse_combine_file(&combine_file, g_score_threshold)?;
    
    if enzymes.is_empty() {
        eprintln!("警告: {} 未找到使用的酶", sample_name);
        return Ok(());
    }

    println!("样品 {}: 使用 {} 个酶，{} 个分类通过 G-score 阈值", 
             sample_name, enzymes.len(), pass_gscore_classes.len());

    // 创建样品输出目录
    let sample_output_dir = output_dir.join(sample_name);
    std::fs::create_dir_all(&sample_output_dir)?;

    // 收集所有符合条件的基因组
    let mut selected_genomes = FxHashSet::default();

    // 遍历每个酶
    for enzyme in &enzymes {
        // 检查数据库文件是否存在
        let enzyme_db = database_dir.join(format!("{}.species.fa.gz", enzyme));
        if !enzyme_db.exists() {
            eprintln!(
                "警告: 数据库文件不存在: {}",
                enzyme_db.display()
            );
            continue;
        }

        // 读取 GCF_detected.xls 文件
        let gcf_detected_file = qual_dir
            .join(sample_name)
            .join(format!("{}.{}.GCF_detected.xls", sample_name, enzyme));

        if !gcf_detected_file.exists() {
            eprintln!(
                "警告: {} 没有 {} 的 GCF_detected.xls 文件: {}",
                sample_name,
                enzyme,
                gcf_detected_file.display()
            );
            continue;
        }

        // 解析 GCF_detected.xls
        let gcf_list = parse_gcf_detected_file(
            &gcf_detected_file,
            &pass_gscore_classes,
            gcf_threshold,
        )?;

        for gcf_id in gcf_list {
            selected_genomes.insert(gcf_id);
        }
    }

    // 输出 sdb.list 文件（排序去重）
    let sdb_list_file = sample_output_dir.join("sdb.list");
    let mut writer = BufWriter::new(File::create(&sdb_list_file)?);

    let mut genome_list: Vec<&String> = selected_genomes
        .iter()
        .filter_map(|gcf_id| gcf_to_classify.get(gcf_id))
        .collect();
    genome_list.sort();
    genome_list.dedup();

    let genome_count = genome_list.len();
    for genome_line in genome_list {
        writeln!(writer, "{}", genome_line)?;
    }

    println!("样品 {}: 筛选出 {} 个基因组", sample_name, genome_count);
    Ok(())
}

/// 解析 combine.xls 文件，提取使用的酶和通过 G-score 阈值的分类
fn parse_combine_file(
    combine_file: &Path,
    g_score_threshold: i32,
) -> Result<(Vec<String>, FxHashSet<String>)> {
    let file = File::open(combine_file)
        .with_context(|| format!("无法打开 combine 文件: {}", combine_file.display()))?;
    let reader = BufReader::new(file);

    let mut enzymes = Vec::new();
    let mut pass_gscore_classes = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        // 跳过表头行
        if line.to_uppercase().starts_with("#KINGDOM") {
            continue;
        }

        // 解析酶信息行（如 #BcgI CjeI combine）
        if line.starts_with('#') && !line.to_uppercase().starts_with("#KINGDOM") {
            let parts: Vec<&str> = line.trim_start_matches('#').split_whitespace().collect();
            for part in parts {
                let enzyme = part.trim();
                if enzyme != "combine" && !enzyme.is_empty() {
                    enzymes.push(enzyme.to_string());
                }
            }
            continue;
        }

        // 解析数据行
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }

        // 最后一列是 G_Score
        if let Ok(g_score) = parts[parts.len() - 1].parse::<f64>() {
            if g_score > g_score_threshold as f64 {
                // 获取分类信息（前 N-8 列）
                let class = parts[0..parts.len() - 8].join("\t");
                pass_gscore_classes.insert(class);
            }
        }
    }

    // 去重酶列表
    enzymes.sort();
    enzymes.dedup();

    Ok((enzymes, pass_gscore_classes))
}

/// 解析 GCF_detected.xls 文件，返回符合条件的 GCF 列表
fn parse_gcf_detected_file(
    gcf_detected_file: &Path,
    pass_gscore_classes: &FxHashSet<String>,
    gcf_threshold: i32,
) -> Result<Vec<String>> {
    let file = File::open(gcf_detected_file)
        .with_context(|| format!("无法打开 GCF_detected 文件: {}", gcf_detected_file.display()))?;
    let reader = BufReader::new(file);

    let mut gcf_list = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }

        // 格式: class\tGCF\tGCF_all_theory_num\tdetected_tag_num\tpercent
        // 分类信息是前 N-4 列
        let class = parts[0..parts.len() - 4].join("\t");
        let gcf_id = parts[parts.len() - 4].to_string();
        let detected_tag_num: i32 = parts[parts.len() - 2].parse().unwrap_or(0);

        // 检查分类是否通过 G-score 阈值，且检测到的标签数 > gcf_threshold
        if pass_gscore_classes.contains(&class) && detected_tag_num > gcf_threshold {
            gcf_list.push(gcf_id);
        }
    }

    Ok(gcf_list)
}

