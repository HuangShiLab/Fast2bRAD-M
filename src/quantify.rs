use anyhow::{bail, Context, Result};
use clap::Parser;
use fxhash::{FxHashMap, FxHashSet};
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

use crate::enzymes::{enzyme_by_id, enzyme_by_name};

/// 丰度计算参数
#[derive(Parser, Debug)]
pub struct QuantifyArgs {
    /// 样品列表文件（TSV格式：sample_name<tab>data_path）
    #[arg(short = 'l', long = "list")]
    pub sample_list: PathBuf,

    /// 数据库目录
    #[arg(short = 'd', long = "database")]
    pub database_dir: PathBuf,

    /// 分类层级
    #[arg(short = 't', long = "taxonomy")]
    pub taxonomy_level: String,

    /// 酶编号（1-16）或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// G-score 阈值（过滤假阳性）
    #[arg(short = 'g', long = "gscore", default_value = "0")]
    pub g_score_threshold: f64,

    /// 是否输出详细信息
    #[arg(short = 'v', long = "verbose", default_value = "yes")]
    pub verbose: String,
}

/// 分类信息
#[derive(Debug, Clone)]
struct TaxonomyInfo {
    gcf_id: String,
    taxonomy: String, // 完整分类路径
}

/// 丰度统计结果
#[derive(Debug, Default)]
struct AbundanceStats {
    theoretical_tag_num: f64,    // 理论标签平均数
    sequenced_tag_num: usize,    // 测到的标签数
    sequenced_reads_num: usize,  // 测到的总reads数
    sequenced_tag_num_gt1: usize, // 深度>1的标签数
}

impl AbundanceStats {
    fn percent(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 {
            (self.sequenced_tag_num as f64 / self.theoretical_tag_num) * 100.0
        } else {
            0.0
        }
    }

    fn reads_per_theoretical(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 {
            self.sequenced_reads_num as f64 / self.theoretical_tag_num
        } else {
            0.0
        }
    }

    fn reads_per_sequenced(&self) -> f64 {
        if self.sequenced_tag_num > 0 {
            self.sequenced_reads_num as f64 / self.sequenced_tag_num as f64
        } else {
            0.0
        }
    }

    fn g_score(&self) -> f64 {
        ((self.sequenced_tag_num * self.sequenced_reads_num) as f64).sqrt()
    }
}

/// 主函数
pub fn run(args: QuantifyArgs) -> Result<()> {
    // 验证参数
    let verbose = args.verbose.to_lowercase() == "yes";

    // 获取酶
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num)
            .ok_or_else(|| anyhow::anyhow!("无效的酶切位点编号: {}", args.enzyme_site))?
    } else {
        enzyme_by_name(&args.enzyme_site)
            .ok_or_else(|| anyhow::anyhow!("无效的酶名称: {}", args.enzyme_site))?
    };

    // 验证分类层级
    let tax_level = validate_taxonomy_level(&args.taxonomy_level)?;

    tracing::info!(
        "COMMAND: quantify -l {} -d {} -t {} -s {} -o {} -g {} -v {}",
        args.sample_list.display(),
        args.database_dir.display(),
        args.taxonomy_level,
        args.enzyme_site,
        args.output_dir.display(),
        args.g_score_threshold,
        args.verbose
    );

    // 创建输出目录
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir.display()))?;

    // 检查数据库文件
    let db_file = args
        .database_dir
        .join(format!("{}.{}.fa.gz", enzyme.name, tax_level));
    let classify_file = args
        .database_dir
        .join("abfh_classify_with_speciename.txt.gz");

    if !db_file.exists() {
        bail!("数据库文件不存在: {}", db_file.display());
    }
    if !classify_file.exists() {
        bail!("分类文件不存在: {}", classify_file.display());
    }

    // 加载数据库
    tracing::info!("### 加载数据库：{}", db_file.display());
    let (tag_to_gcfs, gcf_to_taxonomy, tag_theory_num) =
        load_database(&db_file, &classify_file, &args.taxonomy_level)?;
    tracing::info!("### 数据库加载完成");

    // 处理样品列表
    let samples = read_sample_list(&args.sample_list)?;
    tracing::info!("共 {} 个样品待处理", samples.len());

    for (sample_name, sample_data) in samples {
        tracing::info!("\n### ({}) 样品分析开始", sample_name);

        let result = process_sample(
            &sample_name,
            &sample_data,
            &tag_to_gcfs,
            &gcf_to_taxonomy,
            &tag_theory_num,
            enzyme,
            &args.output_dir,
            args.g_score_threshold,
            verbose,
        );

        match result {
            Ok(_) => tracing::info!("### ({}) 样品分析完成", sample_name),
            Err(e) => tracing::error!("!!! ({}) 错误: {}", sample_name, e),
        }
    }

    tracing::info!("\n全部完成！");
    Ok(())
}

/// 验证分类层级
fn validate_taxonomy_level(level: &str) -> Result<String> {
    let valid_levels = [
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain",
    ];
    if valid_levels.contains(&level) {
        Ok(level.to_string())
    } else {
        bail!("无效的分类层级: {}", level);
    }
}

/// 读取样品列表
fn read_sample_list(list_path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let content = std::fs::read_to_string(list_path)
        .with_context(|| format!("无法读取样品列表: {}", list_path.display()))?;

    let mut samples = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        let sample_name = parts[0].to_string();
        let sample_path = PathBuf::from(parts[1]);

        if !sample_path.exists() {
            tracing::warn!(
                "警告: 样品文件不存在: {} ({})",
                sample_path.display(),
                sample_name
            );
            continue;
        }

        samples.push((sample_name, sample_path));
    }

    Ok(samples)
}

/// 加载数据库
fn load_database(
    db_file: &Path,
    classify_file: &Path,
    tax_level: &str,
) -> Result<(
    FxHashMap<Vec<u8>, Vec<String>>,
    FxHashMap<String, String>,
    FxHashMap<String, FxHashMap<String, FxHashMap<Vec<u8>, usize>>>,
)> {
    let level_index = get_taxonomy_level_index(tax_level);

    // 读取分类信息
    let mut gcf_to_taxonomy: FxHashMap<String, String> = FxHashMap::default();
    let content = if classify_file.to_str().unwrap().ends_with(".gz") {
        use flate2::read::GzDecoder;
        use std::io::Read;
        let file = File::open(classify_file)?;
        let mut decoder = GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;
        content
    } else {
        std::fs::read_to_string(classify_file)?
    };

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= level_index {
            continue;
        }

        let gcf_id = parts[0].to_string();
        let taxonomy = parts[1..=level_index].join("\t");
        gcf_to_taxonomy.insert(gcf_id, taxonomy);
    }

    // 读取数据库标签
    let mut tag_to_gcfs: FxHashMap<Vec<u8>, Vec<String>> = FxHashMap::default();
    let mut tag_theory_num: FxHashMap<String, FxHashMap<String, FxHashMap<Vec<u8>, usize>>> =
        FxHashMap::default();

    let mut reader = parse_fastx_file(db_file)?;

    while let Some(record) = reader.next() {
        let record = record.context("解析数据库文件失败")?;

        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() < 6 {
            continue;
        }

        let gcf_id = parts[0].trim_start_matches('>');
        let unique_flag = parts[5];

        // 只处理 unique 标签
        if unique_flag != "1" {
            continue;
        }

        if let Some(taxonomy) = gcf_to_taxonomy.get(gcf_id) {
            let mut tag_seq = record.seq().to_vec();
            tag_seq.make_ascii_uppercase();

            tag_to_gcfs
                .entry(tag_seq.clone())
                .or_insert_with(Vec::new)
                .push(gcf_id.to_string());

            *tag_theory_num
                .entry(taxonomy.clone())
                .or_insert_with(FxHashMap::default)
                .entry(gcf_id.to_string())
                .or_insert_with(FxHashMap::default)
                .entry(tag_seq)
                .or_insert(0) += 1;
        }
    }

    Ok((tag_to_gcfs, gcf_to_taxonomy, tag_theory_num))
}

/// 获取分类层级索引
fn get_taxonomy_level_index(level: &str) -> usize {
    match level {
        "kingdom" => 1,
        "phylum" => 2,
        "class" => 3,
        "order" => 4,
        "family" => 5,
        "genus" => 6,
        "species" => 7,
        "strain" => 8,
        _ => 7, // 默认 species
    }
}

/// 处理单个样品
fn process_sample(
    sample_name: &str,
    sample_data: &Path,
    tag_to_gcfs: &FxHashMap<Vec<u8>, Vec<String>>,
    gcf_to_taxonomy: &FxHashMap<String, String>,
    tag_theory_num: &FxHashMap<String, FxHashMap<String, FxHashMap<Vec<u8>, usize>>>,
    enzyme: &crate::enzymes::Enzyme,
    output_dir: &Path,
    g_score_threshold: f64,
    verbose: bool,
) -> Result<()> {
    // 统计标签深度
    let mut tag_num: FxHashMap<String, FxHashMap<Vec<u8>, usize>> = FxHashMap::default();
    // 记录每个分类下每个GCF检测到的标签集合 (taxonomy -> GCF -> HashSet<tag>)
    let mut detected_gcf_tag: FxHashMap<String, FxHashMap<String, FxHashSet<Vec<u8>>>> = FxHashMap::default();

    let mut reader = parse_fastx_file(sample_data)?;

    while let Some(record) = reader.next() {
        let record = record.context("解析样品文件失败")?;
        let mut tag_seq = record.seq().to_vec();
        tag_seq.make_ascii_uppercase();

        // 尝试匹配
        if let Some(gcf_list) = tag_to_gcfs.get(&tag_seq) {
            if let Some(first_gcf) = gcf_list.first() {
                if let Some(taxonomy) = gcf_to_taxonomy.get(first_gcf) {
                    *tag_num
                        .entry(taxonomy.clone())
                        .or_insert_with(FxHashMap::default)
                        .entry(tag_seq.clone())
                        .or_insert(0) += 1;

                    for gcf_id in gcf_list {
                        detected_gcf_tag
                            .entry(taxonomy.clone())
                            .or_insert_with(FxHashMap::default)
                            .entry(gcf_id.clone())
                            .or_insert_with(FxHashSet::default)
                            .insert(tag_seq.clone());
                    }
                }
            }
        } else {
            // 尝试反向互补
            let rc_tag = reverse_complement(&tag_seq);
            if let Some(gcf_list) = tag_to_gcfs.get(&rc_tag) {
                if let Some(first_gcf) = gcf_list.first() {
                    if let Some(taxonomy) = gcf_to_taxonomy.get(first_gcf) {
                        *tag_num
                            .entry(taxonomy.clone())
                            .or_insert_with(FxHashMap::default)
                            .entry(rc_tag.clone())
                            .or_insert(0) += 1;

                        for gcf_id in gcf_list {
                            detected_gcf_tag
                                .entry(taxonomy.clone())
                                .or_insert_with(FxHashMap::default)
                                .entry(gcf_id.clone())
                                .or_insert_with(FxHashSet::default)
                                .insert(rc_tag.clone());
                        }
                    }
                }
            }
        }
    }

    if tag_num.is_empty() {
        tracing::warn!(
            "!!! ({}) 警告: {} 级别未检测到任何 2bRAD 标签",
            sample_name,
            enzyme.name
        );
        return Ok(());
    }

    // 创建样品输出目录
    let sample_dir = output_dir.join(sample_name);
    std::fs::create_dir_all(&sample_dir)?;

    // 输出 GCF_detected.xls 文件
    let gcf_detected_file = sample_dir.join(format!("{}.{}.GCF_detected.xls", sample_name, enzyme.name));
    let mut gcf_writer = BufWriter::new(File::create(&gcf_detected_file)?);
    
    // 按分类和GCF排序输出
    let mut taxonomy_list: Vec<&String> = detected_gcf_tag.keys().collect();
    taxonomy_list.sort();
    
    for taxonomy in taxonomy_list {
        let gcf_map = &detected_gcf_tag[taxonomy];
        let mut gcf_list: Vec<&String> = gcf_map.keys().collect();
        gcf_list.sort();
        
        for gcf_id in gcf_list {
            let detected_tags = &gcf_map[gcf_id];
            let detected_tag_num = detected_tags.len();
            
            // 获取理论标签数
            let gcf_all_theory_num = if let Some(gcf_tags) = tag_theory_num.get(taxonomy) {
                if let Some(tags) = gcf_tags.get(gcf_id) {
                    tags.len()
                } else {
                    0
                }
            } else {
                0
            };
            
            let percent = if gcf_all_theory_num > 0 {
                detected_tag_num as f64 / gcf_all_theory_num as f64
            } else {
                0.0
            };
            
            writeln!(
                gcf_writer,
                "{}\t{}\t{}\t{}\t{:.4}",
                taxonomy, gcf_id, gcf_all_theory_num, detected_tag_num, percent
            )?;
        }
    }

    // 输出结果
    let output_file = sample_dir.join(format!("{}.{}.xls", sample_name, enzyme.name));
    let mut writer = BufWriter::new(File::create(&output_file)?);

    // 写入表头
    let level_count = if let Some(first_tax) = gcf_to_taxonomy.values().next() {
        first_tax.split('\t').count()
    } else {
        7 // 默认到 species
    };
    let tax_levels = get_taxonomy_header_by_level(level_count);
    writeln!(
        writer,
        "#{}",
        [
            tax_levels.as_str(),
            "Theoretical_Tag_Num",
            "Sequenced_Tag_Num",
            "Percent",
            "Sequenced_Reads_Num",
            "Sequenced_Reads_Num/Theoretical_Tag_Num",
            "Sequenced_Reads_Num/Sequenced_Tag_Num",
            "Sequenced_Tag_Num(depth>1)",
            "G_Score"
        ]
        .join("\t")
    )?;

    // 计算每个分类的统计信息
    for (taxonomy, tags) in &tag_num {
        let mut stats = AbundanceStats::default();

        // 计算理论标签平均数
        if let Some(gcf_tags) = tag_theory_num.get(taxonomy) {
            let total_theory: usize = gcf_tags
                .values()
                .map(|tags| tags.values().sum::<usize>())
                .sum();
            stats.theoretical_tag_num = total_theory as f64 / gcf_tags.len() as f64;
        }

        // 统计测到的标签
        stats.sequenced_tag_num = tags.len();
        stats.sequenced_reads_num = tags.values().sum();
        stats.sequenced_tag_num_gt1 = tags.values().filter(|&&count| count > 1).count();

        // 计算 G_Score
        let g_score = stats.g_score();

        // 过滤低 G_Score
        if g_score < g_score_threshold {
            continue;
        }

        // 输出
        writeln!(
            writer,
            "{}\t{:.8}\t{}\t{:.8}%\t{}\t{:.8}\t{:.8}\t{}\t{:.8}",
            taxonomy,
            stats.theoretical_tag_num,
            stats.sequenced_tag_num,
            stats.percent(),
            stats.sequenced_reads_num,
            stats.reads_per_theoretical(),
            stats.reads_per_sequenced(),
            stats.sequenced_tag_num_gt1,
            g_score
        )?;

        // 详细输出
        if verbose {
            let output_name = taxonomy.split('\t').last().unwrap_or("unknown");
            let detail_dir = sample_dir.join(format!("{}.{}", sample_name, enzyme.name));
            std::fs::create_dir_all(&detail_dir)?;
            let detail_file = detail_dir.join(format!("{}.xls", output_name));
            let mut detail_writer = BufWriter::new(File::create(detail_file)?);

            for (tag, &count) in tags {
                writeln!(
                    detail_writer,
                    "{}\t{}",
                    std::str::from_utf8(tag).unwrap_or(""),
                    count
                )?;
            }
        }
    }

    Ok(())
}

/// 获取分类表头
fn get_taxonomy_header_by_level(level_count: usize) -> String {
    let headers = [
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain",
    ];
    headers[..level_count.min(8)].join("\t")
}

/// 反向互补
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b,
        })
        .collect()
}

