use anyhow::{bail, Context, Result, anyhow};
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use needletail::parse_fastx_file;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};

// ================== Hash 工具与类型定义 ==================

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
// 逻辑：严格比较序列与其反向互补，取较小者，确保正反链生成同一个 Hash
fn get_canonical_sequence(seq: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(seq);
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

// ================== 参数与结构体 ==================

/// 构建定量数据库参数
#[derive(Parser, Debug)]
pub struct BuildQuanDbArgs {
    /// 基因组分类列表文件（TSV格式）
    #[arg(short = 'l', long = "list")]
    pub genome_list: PathBuf,

    /// 酶编号（1-16）或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 构建数据库的分类层级（逗号分隔，或 'all'）
    #[arg(short = 't', long = "taxonomy")]
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
    #[arg(short = 'r', long = "remove-redundant", default_value = "no")]
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
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "kingdom" => Ok(TaxonomyLevel::Kingdom),
            "phylum" => Ok(TaxonomyLevel::Phylum),
            "class" => Ok(TaxonomyLevel::Class),
            "order" => Ok(TaxonomyLevel::Order),
            "family" => Ok(TaxonomyLevel::Family),
            "genus" => Ok(TaxonomyLevel::Genus),
            "species" => Ok(TaxonomyLevel::Species),
            "strain" => Ok(TaxonomyLevel::Strain),
            _ => bail!("无效的分类层级: {}", s),
        }
    }

    fn as_str(&self) -> &'static str {
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
}

/// 基因组记录
struct GenomeRecord {
    gcf_id: String,
    taxonomy: Vec<String>, 
}

/// 主函数
pub fn run(args: BuildQuanDbArgs) -> Result<()> {
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;

    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num)
            .ok_or_else(|| anyhow!("无效的酶切位点编号: {}", args.enzyme_site))?
    } else {
        enzyme_by_name(&args.enzyme_site)
            .ok_or_else(|| anyhow!("无效的酶名称: {}", args.enzyme_site))?
    };

    println!("读取基因组分类列表 ...");
    let (genome_records, _) = read_genome_list(&args.genome_list, &levels)?;
    println!("共 {} 个基因组", genome_records.len());

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir.display()))?;

    let remove_redundant = args.remove_redundant.to_lowercase() == "yes";

    // 酶切基因组或合并预酶切文件（转换为 Hash 格式）
    let intermediate_enzyme_file = digest_genomes_to_intermediate_file(
        &genome_records,
        enzyme,
        &args.output_dir,
        args.enzyme_file.as_ref(),
        args.pre_digested_dir.as_ref(),
    )?;

    // 为每个分类层级构建数据库
    for level in &levels {
        println!("\n========== 构建 {} 级别数据库 (Hash模式) ==========", level.as_str());
        build_database_for_level(
            &intermediate_enzyme_file,
            enzyme,
            &args.output_dir,
            *level,
            &genome_records,
            remove_redundant,
        )?;
    }

    println!("\n全部完成！");
    Ok(())
}

fn parse_taxonomy_levels(levels_str: &str) -> Result<Vec<TaxonomyLevel>> {
    if levels_str == "all" {
        Ok(vec![
            TaxonomyLevel::Kingdom,
            TaxonomyLevel::Phylum,
            TaxonomyLevel::Class,
            TaxonomyLevel::Order,
            TaxonomyLevel::Family,
            TaxonomyLevel::Genus,
            TaxonomyLevel::Species,
            TaxonomyLevel::Strain,
        ])
    } else {
        levels_str
            .split(',')
            .map(|s| TaxonomyLevel::from_str(s.trim()))
            .collect()
    }
}

fn read_genome_list(
    list_path: &Path,
    levels: &[TaxonomyLevel],
) -> Result<(Vec<GenomeRecord>, FxHashMap<String, usize>)> {
    use std::io::{BufRead, BufReader};
    
    let file = File::open(list_path)
        .with_context(|| format!("无法读取文件: {}", list_path.display()))?;
    let reader = BufReader::new(file);

    let mut genomes = Vec::new();
    let mut taxonomy_levels_map = FxHashMap::default();
    let mut is_gtdb_format = false;
    let mut first_data_line = true;

    for level in levels {
        taxonomy_levels_map.insert(level.as_str().to_string(), *level as usize);
    }

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed_line = line.trim();
        if trimmed_line.is_empty() || trimmed_line.starts_with('#') { continue; }

        let parts: Vec<&str> = trimmed_line.split('\t').collect();
        
        if first_data_line {
            first_data_line = false;
            // 简单的 GTDB 检测
            if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                is_gtdb_format = true;
                if parts[0] == "accession" || parts[0] == "GCF_ID" { continue; }
            }
        }

        if is_gtdb_format {
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord { 
                gcf_id: extract_gcf_id(parts[0].trim()), // Trim ID
                taxonomy: parse_gtdb_taxonomy(parts[1])? 
            });
        } else {
            // 传统格式：GCF, Kingdom, Phylum ...
            if parts.len() < 2 { continue; }
            genomes.push(GenomeRecord { 
                gcf_id: parts[0].trim().to_string(), // Trim ID
                // [FIX]: 跳过 ID (index 0)，从 index 1 开始
                // [FIX]: 增加 .trim() 确保没有隐形空格导致分类不匹配
                taxonomy: parts[1..].iter().map(|s| s.trim().to_string()).collect() 
            });
        }
    }
    
    Ok((genomes, taxonomy_levels_map))
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

// ================== 酶切与合并 (Hash 处理) ==================

fn digest_genomes_to_intermediate_file(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    output_dir: &Path,
    enzyme_file: Option<&PathBuf>,
    pre_digested_dir: Option<&PathBuf>,
) -> Result<PathBuf> {
    if let Some(existing_file) = enzyme_file {
        println!("使用预酶切文件: {}", existing_file.display());
        return Ok(existing_file.clone());
    }

    if let Some(dir) = pre_digested_dir {
        println!("从预酶切目录合并文件：{}", dir.display());
        let output_file = output_dir.join(format!("{}.enzyme.fa.gz", enzyme.name));
        merge_pre_digested_files(genomes, enzyme, dir, &output_file)?;
        return Ok(output_file);
    }

    bail!("请使用 -e 参数提供预酶切文件或 --pre-digested-dir 提供预酶切目录");
}

fn merge_pre_digested_files(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    pre_digested_dir: &Path,
    output_file: &Path,
) -> Result<()> {
    let file = File::create(output_file)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let total = genomes.len();
    let mut processed = 0;
    let mut found = 0;

    for genome in genomes {
        let genome_id = genome.gcf_id.split('.').take(2).collect::<Vec<_>>().join(".");
        
        let patterns = [
            format!("{}.{}.iibdb", genome_id, enzyme.name),
            format!("{}.{}.fa.gz", genome_id, enzyme.name),
            format!("{}.{}.iibdb", genome.gcf_id, enzyme.name),
            format!("{}.{}.fa.gz", genome.gcf_id, enzyme.name),
        ];
        
        let mut file_found = false;
        for pattern in &patterns {
            let test_path = pre_digested_dir.join(pattern);
            if test_path.exists() {
                // 调用处理函数，将序列转换为 Hash (如果是序列) 或保持 Hash (如果是 Hash)
                process_and_write_file(&test_path, &genome.gcf_id, &mut writer)?;
                file_found = true;
                found += 1;
                break;
            }
        }
        
        if !file_found {
            eprintln!("警告：未找到基因组 {} 的预酶切文件", genome.gcf_id);
        }

        processed += 1;
        if processed % (total / 20).max(1) == 0 {
            print!("\r进度: {}/{}", processed, total);
            std::io::stdout().flush()?;
        }
    }

    println!("\n合并完成：{}/{} 个基因组", found, total);
    Ok(())
}

// 读取预酶切文件，转换为统一的 >Header \n Hash 格式写入
fn process_and_write_file(
    path: &Path,
    gcf_id: &str,
    writer: &mut GzEncoder<File>,
) -> Result<()> {
    let mut reader = parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record.context("解析Fastx失败")?;
        let header_bytes = record.id();
        let header_str = std::str::from_utf8(header_bytes).unwrap_or("");
        
        // 尝试解析 header 提取 scaffold, pos 等
        // 兼容 >scaffold-start-end-index 和 >ID:Pos 等格式
        // 目标格式: >GCF|index|scaffold|pos|strand|unique
        
        let mut scaffold = header_str;
        let mut pos = "0";
        let mut tag_index = "0";

        if let Some(parts) = header_str.split_once('|') {
             // 已经是标准格式 >GCF|index|scaffold|...
             // 只需提取需要的信息，或者如果已经是目标格式，我们可能不需要做太多
             // 但为了统一，我们重新组装
             // 这里的处理依赖于上游 Extract 的输出格式
             // 假设 Extract 输出的是 >Scaffold-Pos 或 >Scaffold:Pos
        } 
        
        if let Some(idx) = header_str.rfind(':') {
            scaffold = &header_str[..idx];
            pos = &header_str[idx+1..];
        } else if let Some(parts) = header_str.split_once('-') {
            scaffold = parts.0;
            if let Some(p) = parts.1.split('-').next() {
                pos = p;
            }
             // 尝试提取 index
            if let Some(last_dash) = header_str.rfind('-') {
                 tag_index = &header_str[last_dash+1..];
            }
        }

        let seq_bytes = record.seq();
        let hash_val: Hash;

        // 检查内容是 序列 还是 Hash字符串
        // 简单的检查方法：看是否包含非数字字符 (注意 Hash 字符串是纯数字)
        // 或者是解析 u64 成功
        let seq_str = std::str::from_utf8(&seq_bytes).unwrap_or("");
        
        if let Ok(val) = seq_str.trim().parse::<u64>() {
            // 已经是 Hash
            hash_val = val;
        } else {
            // 是 DNA 序列，计算 Canonical Hash
            let mut seq_vec = seq_bytes.to_vec();
            seq_vec.make_ascii_uppercase();
            let canonical = get_canonical_sequence(&seq_vec);
            hash_val = hash_bytes(&canonical);
        }

        // 写入 Hash
        // Strand 设为 0，Unique 设为 0 (初始状态)
        writeln!(writer, ">{}|{}|{}|{}|0|0", gcf_id, tag_index, scaffold, pos)?;
        writeln!(writer, "{}", hash_val)?;
    }
    Ok(())
}

// ================== 数据库构建 (Hash 版) ==================
fn build_database_for_level(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<()> {
    println!("第 1 步：统计标签分类信息 ...");

    let mut gcf_to_taxonomy = FxHashMap::default();
    
    // Debug: 打印前 3 个基因组的分类字符串，检查是否符合预期
    println!("  [DEBUG] 检查分类字符串格式 (Level={}, Index 0..{}):", level.as_str(), level as usize);
    
    for (i, genome) in genomes.iter().enumerate() {
        // [FIX]: 使用 0..level 确保只取到 Species，不含 Strain
        let end_index = std::cmp::min(level as usize, genome.taxonomy.len());
        let taxonomy_str = genome.taxonomy[0..end_index].join("\t");
        
        println!("  [DEBUG] GCF: {} -> Tax: '{}'", genome.gcf_id, taxonomy_str);

        gcf_to_taxonomy.insert(genome.gcf_id.clone(), taxonomy_str);
    }

    // Pass 1: Collect
    let (tag_taxonomy, genome_tags) = collect_tag_taxonomies(
        enzyme_file,
        &gcf_to_taxonomy,
        genomes,
        remove_redundant,
    )?;

    println!("第 2 步：识别特异性标签并输出数据库 ...");
    identify_and_output_unique_tags(
        enzyme_file,
        enzyme,
        output_dir,
        level,
        &tag_taxonomy,
        &genome_tags,
        remove_redundant,
    )?;

    Ok(())
}

type TagTaxonomyMap = FxHashMap<Hash, FxHashSet<String>>;
type GenomeTagCountMap = FxHashMap<String, FxHashMap<Hash, usize>>;

// ================== 深度调试版核心函数 ==================
fn collect_tag_taxonomies(
    enzyme_file: &Path,
    gcf_to_taxonomy: &FxHashMap<String, String>,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<(TagTaxonomyMap, GenomeTagCountMap)> {
    let mut tag_taxonomy: TagTaxonomyMap = FxHashMap::default();
    let mut genome_tags: GenomeTagCountMap = FxHashMap::default();

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut processed_gcfs = FxHashSet::default();
    let mut total_records = 0;
    let mut valid_hash_records = 0;

    while let Some(record) = reader.next() {
        let record = record.context("解析酶切文件失败")?;
        total_records += 1;
        
        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        // 宽松检查，防止 header 格式轻微差异导致跳过
        if parts.is_empty() { continue; }

        let gcf_id = parts[0].trim_start_matches('>');
        if !gcf_to_taxonomy.contains_key(gcf_id) { continue; }

        let taxonomy = gcf_to_taxonomy.get(gcf_id).unwrap();
        processed_gcfs.insert(gcf_id.to_string());

        // 解析 Hash
        let seq_bytes = record.seq();
        let hash_str = std::str::from_utf8(&seq_bytes).unwrap_or("0");
        
        // [DEBUG]: 打印前几个 Hash 字符串，确认是否真的是数字
        if total_records < 5 {
             println!("  [DEBUG] 读取Hash行: '{}'", hash_str.trim());
        }

        let hash_val: Hash = match hash_str.trim().parse() {
            Ok(v) => {
                valid_hash_records += 1;
                v
            },
            Err(_) => {
                // [DEBUG]: 如果解析失败，打印出来看看是什么怪东西
                if total_records < 10 {
                    println!("  [WARN] Hash解析失败: '{}' 不是有效的 u64", hash_str.trim());
                }
                continue; 
            },
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
    println!("  [DEBUG] 文件读取统计: 总记录={}, 有效Hash={}, 覆盖基因组={}/{} ({}%)", 
        total_records, valid_hash_records, processed_gcfs.len(), genomes.len(), percent);
    println!("  [DEBUG] Map中唯一的Tag(Hash)总数: {}", tag_taxonomy.len());

    // [DEBUG] 随机抽查一个 Tag 看看它的分类情况
    if let Some((hash, taxa)) = tag_taxonomy.iter().next() {
        println!("  [DEBUG] 抽查 Tag {}: 对应 {} 个分类字符串", hash, taxa.len());
        for t in taxa {
            println!("      -> '{}'", t);
        }
    }

    Ok((tag_taxonomy, genome_tags))
}

fn identify_and_output_unique_tags(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    tag_taxonomy: &TagTaxonomyMap,
    genome_tags: &GenomeTagCountMap,
    remove_redundant: bool,
) -> Result<()> {
    let output_path = output_dir.join(format!("{}.{}.fa.gz", enzyme.name, level.as_str()));

    let output_file = File::create(&output_path)?;
    let encoder = GzEncoder::new(output_file, Compression::default());
    let mut writer = BufWriter::new(encoder);

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut unique_counts: FxHashMap<String, usize> = FxHashMap::default();
    
    // [DEBUG] 统计被拒绝的原因
    let mut rejected_not_found = 0;
    let mut rejected_ambiguous = 0;
    let mut rejected_redundant = 0;
    let mut accepted = 0;

    while let Some(record) = reader.next() {
        let record = record.context("解析酶切文件失败")?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        // 宽松检查
        if parts.is_empty() { continue; }

        let gcf_id = parts[0].trim_start_matches('>');
        // 提取 scaffold, pos 等用于输出
        let tag_index = if parts.len() > 1 { parts[1] } else { "0" };
        let scaffold_id = if parts.len() > 2 { parts[2] } else { "scaffold" };
        let pos = if parts.len() > 3 { parts[3] } else { "0" };

        let seq_bytes = record.seq();
        let hash_str = std::str::from_utf8(&seq_bytes).unwrap_or("0");
        let hash_val: Hash = match hash_str.trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // 检查特异性
        let mut is_unique = false;
        if let Some(taxonomies) = tag_taxonomy.get(&hash_val) {
            if taxonomies.len() == 1 {
                is_unique = true;
            } else {
                rejected_ambiguous += 1;
            }
        } else {
            rejected_not_found += 1;
        }

        // 检查基因组内冗余
        if is_unique && remove_redundant {
            if let Some(genome_tag_counts) = genome_tags.get(gcf_id) {
                if let Some(&count) = genome_tag_counts.get(&hash_val) {
                    if count > 1 { 
                        is_unique = false; 
                        rejected_redundant += 1;
                    }
                }
            }
        }

        if is_unique {
            accepted += 1;
            *unique_counts.entry(gcf_id.to_string()).or_insert(0) += 1;
            // 输出格式 >GCF|index|scaffold|pos|strand|unique
            writeln!(writer, ">{}|{}|{}|{}|0|1", gcf_id, tag_index, scaffold_id, pos)?;
            writeln!(writer, "{}", hash_val)?;
        }
    }

    drop(writer);

    println!("  [DEBUG] 筛选统计: 接受={}, 拒绝(未找到)={}, 拒绝(分类模糊/多物种)={}, 拒绝(基因组冗余)={}", 
        accepted, rejected_not_found, rejected_ambiguous, rejected_redundant);
    println!("  输出数据库：{}", output_path.display());
    println!("  包含 {} 个基因组的特异性标签", unique_counts.len());

    Ok(())
}
