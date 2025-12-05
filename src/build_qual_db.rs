use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::hash::Hasher;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use flate2::write::GzEncoder;
use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use tracing;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;

// ================== 类型定义与 Hash 工具 ==================

pub type Hash = u64;

// 计算 FxHash
pub fn hash_bytes(bytes: &[u8]) -> Hash {
    let mut hasher = FxHasher::default();
    hasher.write(bytes);
    hasher.finish()
}

// 计算反向互补
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let complement = match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            x => x, 
        };
        rc.push(complement);
    }
    rc
}

// 获取 Canonical 序列（用于去重和哈希）
fn get_canonical_sequence(seq: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(seq);
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

// ================== 参数与结构体 ==================

#[derive(Args, Debug)]
pub struct BuildQualDbArgs {
    /// 基因组分类列表文件（TSV 格式）
    #[arg(short = 'l', long = "list")]
    pub genome_list: PathBuf,

    /// 酶编号（1-16）或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 分类水平（逗号分隔，或 'all'）
    #[arg(short = 't', long = "type")]
    pub taxonomy_levels: String,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// 可选：预酶切文件路径（单个合并文件，跳过酶切步骤）
    #[arg(short = 'e', long = "enzyme-file")]
    pub enzyme_file: Option<PathBuf>,

    /// 可选：预酶切目录（extract 批量输出目录）
    #[arg(long = "pre-digested-dir")]
    pub pre_digested_dir: Option<PathBuf>,

    /// 是否删除基因组内冗余标签
    #[arg(short = 'r', long = "remove-redundant", default_value = "yes")]
    pub remove_redundant: String,
}

/// 分类层级
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TaxonomyLevel {
    Kingdom = 1,
    Phylum = 2,
    Class = 3,
    Order = 4,
    Family = 5,
    Genus = 6,
    Species = 7,
    Strain = 8,
}

impl TaxonomyLevel {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "kingdom" => Some(TaxonomyLevel::Kingdom),
            "phylum" => Some(TaxonomyLevel::Phylum),
            "class" => Some(TaxonomyLevel::Class),
            "order" => Some(TaxonomyLevel::Order),
            "family" => Some(TaxonomyLevel::Family),
            "genus" => Some(TaxonomyLevel::Genus),
            "species" => Some(TaxonomyLevel::Species),
            "strain" => Some(TaxonomyLevel::Strain),
            _ => None,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            TaxonomyLevel::Kingdom => "kingdom",
            TaxonomyLevel::Phylum => "phylum",
            TaxonomyLevel::Class => "class",
            TaxonomyLevel::Order => "order",
            TaxonomyLevel::Family => "family",
            TaxonomyLevel::Genus => "genus",
            TaxonomyLevel::Species => "species",
            TaxonomyLevel::Strain => "strain",
        }
    }

    pub fn all_levels() -> Vec<Self> {
        vec![
            TaxonomyLevel::Kingdom,
            TaxonomyLevel::Phylum,
            TaxonomyLevel::Class,
            TaxonomyLevel::Order,
            TaxonomyLevel::Family,
            TaxonomyLevel::Genus,
            TaxonomyLevel::Species,
            TaxonomyLevel::Strain,
        ]
    }
}

/// 基因组分类记录
#[derive(Debug, Clone)]
pub struct GenomeRecord {
    pub gcf_id: String,
    pub taxonomy: Vec<String>,
    pub genome_path: Option<PathBuf>,
}

pub fn run(args: BuildQualDbArgs) -> Result<()> {
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;
    let remove_redundant = args.remove_redundant.eq_ignore_ascii_case("yes");

    io_utils::ensure_directory(&args.output_dir)?;

    tracing::info!("读取基因组分类列表 ...");
    let genomes = read_genome_list(&args.genome_list)?;
    tracing::info!("共 {} 个基因组", genomes.len());

    // 确定酶切文件（Hash格式）
    let enzyme_file = if let Some(ref file) = args.enzyme_file {
        tracing::info!("使用预酶切文件：{}", file.display());
        file.clone()
    } else if let Some(ref dir) = args.pre_digested_dir {
        tracing::info!("从预酶切目录合并文件：{}", dir.display());
        let output_file = args
            .output_dir
            .join(format!("{}.enzyme.fa.gz", enzyme.name));
        merge_pre_digested_files(&genomes, enzyme, dir, &output_file)?;
        output_file
    } else {
        tracing::info!("开始批量酶切基因组并生成 Hash ...");
        let output_file = args
            .output_dir
            .join(format!("{}.enzyme.fa.gz", enzyme.name));
        digest_genomes(&genomes, enzyme, &output_file)?;
        output_file
    };

    // 构建数据库
    for level in &levels {
        tracing::info!("\n========== 构建 {} 级别数据库 (Hash模式) ==========", level.name());
        
        // 修改这里的调用顺序，以匹配新的函数定义
        build_database_for_level(
            &enzyme_file,       // 1. 酶切文件路径
            enzyme,             // 2. 酶对象
            &args.output_dir,   // 3. 输出目录
            *level,             // 4. 分类层级
            &genomes,           // 5. 基因组列表
            remove_redundant,   // 6. 是否去冗余
        )?;
    }

    tracing::info!("\n全部完成！");
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
    if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme); }
    if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme); } }
    bail!("未知的酶：{}", site)
}

fn parse_taxonomy_levels(levels_str: &str) -> Result<Vec<TaxonomyLevel>> {
    if levels_str.eq_ignore_ascii_case("all") { return Ok(TaxonomyLevel::all_levels()); }
    let mut levels = Vec::new();
    for part in levels_str.split(',') {
        let level = TaxonomyLevel::from_str(part.trim())
            .ok_or_else(|| anyhow::anyhow!("无效的分类水平：{}", part))?;
        levels.push(level);
    }
    if levels.is_empty() { bail!("至少需要指定一个分类水平"); }
    Ok(levels)
}

fn read_genome_list(path: &Path) -> Result<Vec<GenomeRecord>> {
    let file = File::open(path)
        .with_context(|| format!("无法读取基因组列表文件：{}", path.display()))?;
    let reader = BufReader::new(file);
    let mut genomes = Vec::new();
    let mut is_gtdb_format = false;
    let mut first_data_line = true;

    for (_line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }

        let parts: Vec<&str> = trimmed.split('\t').collect();
        
        if first_data_line {
            first_data_line = false;
            if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                is_gtdb_format = true;
                if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; }
            }
        }

        if is_gtdb_format {
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord {
                gcf_id: extract_gcf_id(parts[0]),
                taxonomy: parse_gtdb_taxonomy(parts[1])?,
                genome_path: None,
            });
        } else {
            if parts.len() < 9 { continue; }
            let genome_path = if parts.len() > 9 { Some(PathBuf::from(parts[9])) } else { None };
            genomes.push(GenomeRecord {
                gcf_id: parts[0].to_string(),
                taxonomy: parts[1..9].iter().map(|s| s.to_string()).collect(),
                genome_path,
            });
        }
    }
    if genomes.is_empty() { bail!("基因组列表为空"); }
    Ok(genomes)
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

fn parse_gtdb_taxonomy(gtdb_str: &str) -> Result<Vec<String>> {
    let parts: Vec<&str> = gtdb_str.split(';').collect();
    let mut taxonomy = Vec::new();
    for part in parts.iter() {
        if let Some(pos) = part.find("__") {
            taxonomy.push(part[pos+2..].to_string());
        } else {
            taxonomy.push(part.to_string());
        }
    }
    while taxonomy.len() < 8 {
        if let Some(last) = taxonomy.last() { taxonomy.push(format!("{}_strain", last)); }
        else { taxonomy.push("unknown".to_string()); }
    }
    Ok(taxonomy)
}

// ================== 酶切与合并 (Hash 版本) ==================

fn merge_pre_digested_files(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    pre_digested_dir: &Path,
    output_file: &Path,
) -> Result<()> {

    let file = File::create(output_file)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let total = genomes.len();
    let mut found_count = 0;
    
    // 创建进度条
    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    for genome in genomes {
        // 尝试匹配不同后缀: .iibdb (新版 Extract), .fa.gz (旧版)
        let genome_id = genome.gcf_id.split('.').take(2).collect::<Vec<_>>().join(".");
        
        let patterns = [
            format!("{}.{}.iibdb", genome_id, enzyme.name),
            format!("{}.{}.fa.gz", genome_id, enzyme.name),
            format!("{}.{}.iibdb", genome.gcf_id, enzyme.name),
        ];

        let mut file_found = false;
        for pattern in &patterns {
            let test_path = pre_digested_dir.join(pattern);
            if test_path.exists() {
                // 如果是 .gz 结尾用 Gzip 解压，否则直接读取
                if pattern.ends_with(".gz") {
                     copy_content_with_gcf(&test_path, &genome.gcf_id, &mut writer, true)?;
                } else {
                     copy_content_with_gcf(&test_path, &genome.gcf_id, &mut writer, false)?;
                }
                
                file_found = true;
                found_count += 1;
                break;
            }
        }
        
        if !file_found {
            tracing::warn!("警告：未找到基因组 {} 的预酶切文件", genome.gcf_id);
        }

        pb.inc(1);
    }
    tracing::info!("合并完成，找到 {}/{}", found_count, total);
    pb.finish();
    Ok(())
}

fn copy_content_with_gcf(
    file_path: &Path,
    gcf_id: &str,
    writer: &mut GzEncoder<File>,
    is_gzipped: bool,
) -> Result<()> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut line_iter = reader.lines();
    while let Some(line_res) = line_iter.next() {
        let line = line_res?;
        if line.starts_with('>') {
            let header = &line[1..];
            // 兼容多种 Header 格式
            let (scaffold, pos) = if let Some(idx) = header.rfind(':') {
                (&header[..idx], &header[idx+1..])
            } else if let Some(parts) = header.split_once('-') {
                 (parts.0, parts.1.split('-').next().unwrap_or("0"))
            } else {
                (header, "0")
            };

            // 读取下一行（Hash值）
            if let Some(seq_line_res) = line_iter.next() {
                let seq_line = seq_line_res?;
                if seq_line.trim().parse::<u64>().is_ok() {
                    writeln!(writer, ">{}|0|{}|{}|0|-", gcf_id, scaffold, pos)?;
                    writeln!(writer, "{}", seq_line.trim())?;
                }
            }
        }
    }
    Ok(())
}

fn digest_genomes(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    output_file: &Path,
) -> Result<()> {
    let file = File::create(output_file)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let total = genomes.len();
    
    // 创建进度条
    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    
    for (_i, genome) in genomes.iter().enumerate() {
        let genome_path = genome.genome_path.as_ref()
            .ok_or_else(|| anyhow!("基因组 {} 缺少路径", genome.gcf_id))?;

        if !genome_path.exists() {
             tracing::warn!("警告：基因组文件不存在 {}", genome_path.display());
             pb.inc(1);
             continue;
        }

        let mut reader = parse_fastx_file(genome_path)?;
        let mut tag_index = 0usize;

        while let Some(record) = reader.next() {
            let record = record.context("解析序列失败")?;
            let seq_id = std::str::from_utf8(record.id()).unwrap_or("seq")
                .split_whitespace().next().unwrap_or("seq");
            
            let mut sequence = record.seq().to_vec();
            sequence.make_ascii_uppercase();

            let positions = enzyme.find_all_tags(&sequence);

            for (pos, len) in positions {
                tag_index += 1;
                let tag_seq = &sequence[pos..pos + len];
                let canonical = get_canonical_sequence(tag_seq);
                let hash_val = hash_bytes(&canonical);

                writeln!(writer, ">{}|{}|{}|{}|0|-", genome.gcf_id, tag_index, seq_id, pos + 1)?;
                writeln!(writer, "{}", hash_val)?;
            }
        }

        pb.inc(1);
    }
    pb.finish_with_message("酶切完成");
    Ok(())
}

// ================== 构建数据库核心逻辑 (Hash 版 + 恢复统计) ==================

// ================== 数据库构建 (Hash 版 - 修复分类索引) ==================
// ================== 数据库构建逻辑 ==================

fn build_database_for_level(
    enzyme_file: &Path,        // 对应调用第 1 个参数
    enzyme: &'static Enzyme,   // 对应调用第 2 个参数
    output_dir: &Path,         // 对应调用第 3 个参数
    level: TaxonomyLevel,      // 对应调用第 4 个参数
    genomes: &[GenomeRecord],  // 对应调用第 5 个参数
    remove_redundant: bool,    // 对应调用第 6 个参数
) -> Result<()> {
    
    // 1. 统计阶段
    tracing::info!("第 1 步：统计标签分类信息 ...");

    let mut gcf_to_taxonomy = FxHashMap::default();
    for genome in genomes {
        // [关键修复]：使用 0..level 确保只取到当前层级，不包含更细的层级(Strain)
        // 例如 Species(7)，取 indices 0..7 (即 0,1,2,3,4,5,6)，对应 K,P,C,O,F,G,S
        if level as usize > genome.taxonomy.len() {
             bail!("分类层级超出范围: 需要 {} 但只有 {}", level.name(), genome.taxonomy.len());
        }
        
        // 这里必须是 0..level，不能是 1..=level
        let taxonomy_str = genome.taxonomy[0..level as usize].join("\t");
        gcf_to_taxonomy.insert(genome.gcf_id.clone(), taxonomy_str);
    }

    // Pass 1: Collect
    let (tag_taxonomy, genome_tags) = collect_tag_taxonomy(
        genomes, 
        enzyme_file, 
        level, // 这里主要用于函数签名，实际 taxonomy map 已经在上面构建好了，可能需要微调 collect_tag_taxonomy 内部逻辑，或者直接传 map 进去
        remove_redundant
    )?;

    // 2. 输出阶段
    tracing::info!("第 2 步：识别特异性标签并输出数据库 ...");
    output_database(
        genomes,
        enzyme,
        enzyme_file,
        level,
        &tag_taxonomy,
        &genome_tags,
        remove_redundant,
        output_dir,
    )?;

    Ok(())
}

type TagTaxonomyMap = FxHashMap<Hash, FxHashSet<String>>;
type GenomeTagCountMap = FxHashMap<String, FxHashMap<Hash, usize>>;

fn collect_tag_taxonomy(
    genomes: &[GenomeRecord],
    enzyme_file: &Path,
    level: TaxonomyLevel,
    remove_redundant: bool,
) -> Result<(TagTaxonomyMap, GenomeTagCountMap)>
{
    let gcf_to_taxonomy: HashMap<String, String> = genomes
        .iter()
        .map(|g| {
            let tax = &g.taxonomy[0..level as usize];
            (g.gcf_id.clone(), tax.join("\t"))
        })
        .collect();

    let mut tag_taxonomy: TagTaxonomyMap = FxHashMap::default();
    let mut genome_tags: GenomeTagCountMap = FxHashMap::default();

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut processed_gcfs = FxHashSet::default();

    while let Some(record) = reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.is_empty() { continue; }

        let gcf_id = parts[0].trim_start_matches('>');
        if !gcf_to_taxonomy.contains_key(gcf_id) { continue; }
        
        let taxonomy = gcf_to_taxonomy.get(gcf_id).unwrap();
        processed_gcfs.insert(gcf_id.to_string());

        let seq_bytes = record.seq();
        let hash_str = std::str::from_utf8(&seq_bytes).unwrap_or("0");
        let hash_val: Hash = match hash_str.trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        tag_taxonomy
            .entry(hash_val)
            .or_insert_with(FxHashSet::default)
            .insert(taxonomy.clone());

        if remove_redundant {
            *genome_tags
                .entry(gcf_id.to_string())
                .or_insert_with(FxHashMap::default)
                .entry(hash_val)
                .or_insert(0) += 1;
        }
    }

    let percent = (processed_gcfs.len() * 100) / genomes.len();
    tracing::info!("  处理了 {}/{} 个基因组 ({}%)", processed_gcfs.len(), genomes.len(), percent);
    Ok((tag_taxonomy, genome_tags))
}

fn output_database(
    _genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    enzyme_file: &Path,
    level: TaxonomyLevel,
    tag_taxonomy: &TagTaxonomyMap,
    genome_tags: &GenomeTagCountMap,
    remove_redundant: bool,
    output_dir: &Path,
) -> Result<()> {
    let output_path = output_dir.join(format!("{}.{}.fa.gz", enzyme.name, level.name()));
    let file = File::create(&output_path)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let mut reader = parse_fastx_file(enzyme_file)?;
    
    // 【恢复统计逻辑】
    let mut unique_counts: FxHashMap<String, usize> = FxHashMap::default();
    let mut total_counts: FxHashMap<String, usize> = FxHashMap::default();

    while let Some(record) = reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() < 6 { continue; }

        let gcf_id = parts[0].trim_start_matches('>');
        let tag_index = parts[1];
        let scaffold_id = parts[2];
        let pos = parts[3];
        // Hash 模式下 strand 不再重要，置为 0
        let _original_strand = parts[4]; 

        let seq_bytes = record.seq();
        let hash_str = std::str::from_utf8(&seq_bytes).unwrap_or("0");
        let hash_val: Hash = match hash_str.trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // 【恢复统计计数】：计算该基因组的总标签数
        *total_counts.entry(gcf_id.to_string()).or_insert(0) += 1;

        // 检查唯一性
        let mut is_unique = if let Some(taxonomies) = tag_taxonomy.get(&hash_val) {
            taxonomies.len() == 1
        } else {
            false
        };

        if is_unique && remove_redundant {
            if let Some(counts) = genome_tags.get(gcf_id) {
                if let Some(&count) = counts.get(&hash_val) {
                    if count > 1 { is_unique = false; }
                }
            }
        }

        if is_unique {
            // 【恢复统计计数】：计算该基因组的特异性标签数
            *unique_counts.entry(gcf_id.to_string()).or_insert(0) += 1;
        }

        let unique_flag = if is_unique { "1" } else { "0" };
        
        writeln!(writer, ">{}|{}|{}|{}|0|{}", gcf_id, tag_index, scaffold_id, pos, unique_flag)?;
        writeln!(writer, "{}", hash_val)?;
    }

    tracing::info!("  输出数据库：{}", output_path.display());
    tracing::info!("  包含 {} 个基因组的特异性标签", unique_counts.len());

    Ok(())
}