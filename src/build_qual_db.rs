use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Args;
use flate2::write::GzEncoder;
use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet};
use needletail::parse_fastx_file;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;

#[derive(Args, Debug)]
pub struct BuildQualDbArgs {
    /// 基因组分类列表文件（TSV 格式）
    /// 支持传统格式（9列）或 GTDB 格式（3列：accession, gtdb_taxonomy, ncbi_taxonomy）
    #[arg(short = 'l', long = "list")]
    pub genome_list: PathBuf,

    /// 酶编号（1-16）或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 分类水平（逗号分隔，或 'all'）
    /// kingdom,phylum,class,order,family,genus,species,strain
    #[arg(short = 't', long = "type")]
    pub taxonomy_levels: String,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// 可选：预酶切文件路径（单个合并文件，跳过酶切步骤）
    #[arg(short = 'e', long = "enzyme-file")]
    pub enzyme_file: Option<PathBuf>,

    /// 可选：预酶切目录（extract 批量输出目录，包含 genome*.fa.gz 文件）
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
    pub taxonomy: Vec<String>, // Kingdom, Phylum, ..., Strain
    pub genome_path: Option<PathBuf>,
}

pub fn run(args: BuildQualDbArgs) -> Result<()> {
    // 解析参数
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;
    let remove_redundant = args.remove_redundant.eq_ignore_ascii_case("yes");

    io_utils::ensure_directory(&args.output_dir)?;

    // 读取基因组分类列表
    println!("读取基因组分类列表 ...");
    let genomes = read_genome_list(&args.genome_list)?;
    println!("共 {} 个基因组", genomes.len());

    // 酶切阶段（或读取预酶切文件/目录）
    let enzyme_file = if let Some(ref file) = args.enzyme_file {
        println!("使用预酶切文件：{}", file.display());
        file.clone()
    } else if let Some(ref dir) = args.pre_digested_dir {
        println!("从预酶切目录合并文件：{}", dir.display());
        let output_file = args
            .output_dir
            .join(format!("{}.enzyme.fa.gz", enzyme.name));
        merge_pre_digested_files(&genomes, enzyme, dir, &output_file)?;
        output_file
    } else {
        println!("开始批量酶切基因组 ...");
        let output_file = args
            .output_dir
            .join(format!("{}.enzyme.fa.gz", enzyme.name));
        digest_genomes(&genomes, enzyme, &output_file)?;
        output_file
    };

    // 为每个分类水平构建数据库
    for level in &levels {
        println!("\n========== 构建 {} 级别数据库 ==========", level.name());
        build_database_for_level(
            &genomes,
            enzyme,
            &enzyme_file,
            *level,
            remove_redundant,
            &args.output_dir,
        )?;
    }

    println!("\n全部完成！");
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
    if let Some(enzyme) = enzyme_by_name(site) {
        return Ok(enzyme);
    }
    if let Ok(id) = site.parse::<u8>() {
        if let Some(enzyme) = enzyme_by_id(id) {
            return Ok(enzyme);
        }
    }
    bail!(
        "未知的酶：{}，支持的酶：{:?}",
        site,
        crate::enzymes::supported_enzyme_names()
    )
}

fn parse_taxonomy_levels(levels_str: &str) -> Result<Vec<TaxonomyLevel>> {
    if levels_str.eq_ignore_ascii_case("all") {
        return Ok(TaxonomyLevel::all_levels());
    }

    let mut levels = Vec::new();
    for part in levels_str.split(',') {
        let level = TaxonomyLevel::from_str(part.trim())
            .ok_or_else(|| anyhow::anyhow!("无效的分类水平：{}", part))?;
        levels.push(level);
    }

    if levels.is_empty() {
        bail!("至少需要指定一个分类水平");
    }

    Ok(levels)
}

/// 读取基因组分类列表
/// 支持两种格式：
/// 1. 传统格式：GCF_ID  Kingdom  Phylum  Class  Order  Family  Genus  Species  Strain  [genome_path]
/// 2. GTDB 格式：accession  gtdb_taxonomy  ncbi_taxonomy （使用 gtdb_taxonomy 列）
fn read_genome_list(path: &Path) -> Result<Vec<GenomeRecord>> {
    let file = File::open(path)
        .with_context(|| format!("无法读取基因组列表文件：{}", path.display()))?;
    let reader = BufReader::new(file);
    let mut genomes = Vec::new();
    let mut is_gtdb_format = false;
    let mut first_data_line = true;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("读取第 {} 行失败", line_no + 1))?;
        let trimmed = line.trim();

        // 跳过注释和空行
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = trimmed.split('\t').collect();
        
        // 自动检测格式（第一行数据）
        if first_data_line {
            first_data_line = false;
            // 检测 GTDB 格式：第二列包含 "d__" 或列名为 "gtdb_taxonomy"
            if parts.len() >= 2 && (parts[1].contains("d__") || parts[1] == "gtdb_taxonomy") {
                is_gtdb_format = true;
                println!("检测到 GTDB 分类格式");
                // 如果是表头，跳过
                if parts[0] == "accession" || parts[0] == "GCF_ID" {
                    continue;
                }
            } else {
                println!("检测到传统分类格式");
            }
        }

        if is_gtdb_format {
            // GTDB 格式：accession  gtdb_taxonomy  [ncbi_taxonomy]
            if parts.len() < 2 {
                bail!(
                    "基因组列表第 {} 行格式错误（GTDB 格式需要至少 2 列）",
                    line_no + 1
                );
            }

            let gcf_id = extract_gcf_id(parts[0]);
            let taxonomy = parse_gtdb_taxonomy(parts[1])?;
            let genome_path = None; // GTDB 格式从 extract 输出目录推断路径

            genomes.push(GenomeRecord {
                gcf_id,
                taxonomy,
                genome_path,
            });
        } else {
            // 传统格式：GCF_ID  Kingdom  Phylum  Class  Order  Family  Genus  Species  Strain  [genome_path]
            if parts.len() < 9 {
                bail!(
                    "基因组列表第 {} 行格式错误，需要至少 9 列（GCFid + 8 级分类）",
                    line_no + 1
                );
            }

            let gcf_id = parts[0].to_string();
            let taxonomy: Vec<String> = parts[1..9].iter().map(|s| s.to_string()).collect();
            let genome_path = if parts.len() > 9 {
                Some(PathBuf::from(parts[9]))
            } else {
                None
            };

            genomes.push(GenomeRecord {
                gcf_id,
                taxonomy,
                genome_path,
            });
        }
    }

    if genomes.is_empty() {
        bail!("基因组列表为空");
    }

    Ok(genomes)
}

/// 从文件名提取 GCF/GCA ID
/// 例如：GCA_000477855.1_Prop_sp_KPL1847_V1_genomic.fna.gz -> GCA_000477855.1
fn extract_gcf_id(filename: &str) -> String {
    let name = filename.split('/').last().unwrap_or(filename);
    
    // 移除 _genomic 及后面的部分
    let name_clean = if let Some(pos) = name.find("_genomic") {
        &name[..pos]
    } else {
        name
    };
    
    // 如果是 GCF/GCA 格式，返回 GCF_XXXXXX.X 格式
    if name_clean.starts_with("GCF_") || name_clean.starts_with("GCA_") {
        let parts: Vec<&str> = name_clean.split('_').collect();
        if parts.len() >= 2 {
            // 返回 GCA_000477855.1 格式（不包含后面的描述）
            return format!("{}_{}", parts[0], parts[1]);
        }
    }
    
    name_clean.to_string()
}

/// 解析 GTDB 分类字符串
/// 例如：d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium_acnes
/// 返回：[Bacteria, Actinobacteriota, Actinomycetia, Propionibacteriales, Propionibacteriaceae, Cutibacterium, Cutibacterium_acnes, Cutibacterium_acnes_strain]
fn parse_gtdb_taxonomy(gtdb_str: &str) -> Result<Vec<String>> {
    let parts: Vec<&str> = gtdb_str.split(';').collect();
    if parts.len() < 7 {
        bail!("GTDB 分类格式错误，需要至少 7 个层级（d__ 到 s__）");
    }

    let mut taxonomy = Vec::new();
    
    for part in parts.iter() {
        // 提取 "d__Bacteria" -> "Bacteria"
        if let Some(pos) = part.find("__") {
            let value = &part[pos+2..];
            taxonomy.push(value.to_string());
        } else {
            taxonomy.push(part.to_string());
        }
    }

    // 补齐到 8 个层级（如果只有 7 个，用 species 复制为 strain）
    while taxonomy.len() < 8 {
        if let Some(last) = taxonomy.last() {
            taxonomy.push(format!("{}_strain", last));
        } else {
            taxonomy.push("unknown".to_string());
        }
    }

    Ok(taxonomy)
}

/// 从 extract 批量输出目录合并预酶切文件
fn merge_pre_digested_files(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    pre_digested_dir: &Path,
    output_file: &Path,
) -> Result<()> {
    use std::io::{BufRead, BufReader};
    use flate2::read::GzDecoder;

    let file = File::create(output_file)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let total = genomes.len();
    let mut processed = 0;
    let mut found = 0;

    for genome in genomes {
        // 构造预酶切文件路径：genome_id.enzyme.fa.gz
        // 去掉扩展名，只保留 GCF/GCA ID
        let genome_id = genome.gcf_id.split('.').take(2).collect::<Vec<_>>().join(".");
        let pre_digested_file = pre_digested_dir.join(format!("{}.{}.fa.gz", genome_id, enzyme.name));

        if !pre_digested_file.exists() {
            // 尝试使用完整的 accession 名称
            let patterns = [
                format!("{}.{}.fa.gz", genome.gcf_id, enzyme.name),
                format!("genome{:02}.{}.fa.gz", processed + 1, enzyme.name),
            ];
            
            let mut file_found = false;
            for pattern in &patterns {
                let test_path = pre_digested_dir.join(pattern);
                if test_path.exists() {
                    copy_gz_content_with_gcf(&test_path, &genome.gcf_id, &mut writer)?;
                    file_found = true;
                    found += 1;
                    break;
                }
            }
            
            if !file_found {
                eprintln!("警告：未找到基因组 {} 的预酶切文件", genome.gcf_id);
            }
        } else {
            copy_gz_content_with_gcf(&pre_digested_file, &genome.gcf_id, &mut writer)?;
            found += 1;
        }

        processed += 1;
        let percent = (processed * 100) / total;
        if processed % (total / 20).max(1) == 0 || processed == total {
            print!("\r进度: {}% ({}/{}, 找到 {})", percent, processed, total, found);
            std::io::stdout().flush()?;
        }
    }

    println!();
    drop(writer);

    println!("合并完成：{}/{} 个基因组", found, total);
    Ok(())
}

/// 复制 gzip 文件内容并转换为标准格式
/// 将 Type 1 格式（>scaffold-start-end-index）转换为数据库格式（>GCF_ID|index|scaffold|pos|strand|unique）
fn copy_gz_content_with_gcf(
    gz_file: &Path,
    gcf_id: &str,
    writer: &mut GzEncoder<File>,
) -> Result<()> {
    use flate2::read::GzDecoder;
    use std::io::{BufRead, BufReader};

    let file = File::open(gz_file)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // 转换格式：>scaffold-start-end-index -> >GCF_ID|index|scaffold|pos|strand|unique
            let header = &line[1..]; // 去掉 '>'
            let parts: Vec<&str> = header.split('-').collect();
            if parts.len() >= 4 {
                let scaffold = parts[0];
                let start = parts[1];
                let index = parts[3];
                // 格式：>GCF_ID|tag_index|scaffold_id|pos|strand|unique
                // strand 和 unique 暂时设为占位值
                writeln!(writer, ">{}|{}|{}|{}|0|0", gcf_id, index, scaffold, start)?;
            } else {
                writeln!(writer, "{}", line)?;
            }
        } else {
            writeln!(writer, "{}", line)?;
        }
    }

    Ok(())
}

/// 复制 gzip 文件内容（原始）
fn copy_gz_content(gz_file: &Path, writer: &mut GzEncoder<File>) -> Result<()> {
    use flate2::read::GzDecoder;
    use std::io::{BufRead, BufReader};

    let file = File::open(gz_file)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    for line in reader.lines() {
        let line = line?;
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

/// 批量酶切基因组
fn digest_genomes(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    output_file: &Path,
) -> Result<()> {
    let file = File::create(output_file)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let total = genomes.len();
    let mut processed = 0;

    for genome in genomes {
        let genome_path = genome
            .genome_path
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("基因组 {} 缺少路径", genome.gcf_id))?;

        if !genome_path.exists() {
            bail!(
                "基因组文件不存在：{}",
                genome_path.display()
            );
        }

        digest_single_genome(&genome.gcf_id, genome_path, enzyme, &mut writer)?;

        processed += 1;
        let percent = (processed * 100) / total;
        if processed % (total / 20).max(1) == 0 || processed == total {
            print!("\r进度: {}% ({}/{})", percent, processed, total);
            std::io::stdout().flush()?;
        }
    }

    println!();
    drop(writer);

    Ok(())
}

/// 酶切单个基因组
fn digest_single_genome(
    gcf_id: &str,
    genome_path: &Path,
    enzyme: &'static Enzyme,
    writer: &mut GzEncoder<File>,
) -> Result<()> {
    let mut reader = parse_fastx_file(genome_path)?;
    let mut tag_index = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;

        let seq_id = std::str::from_utf8(record.id()).unwrap_or("sequence");
        let seq_id = seq_id.split_whitespace().next().unwrap_or("sequence");

        let mut sequence = record.seq().to_vec();
        sequence.make_ascii_uppercase();

        // 查找所有匹配的标签位置
        let positions = enzyme.find_all_tags(&sequence);

        for (pos, len) in positions {
            tag_index += 1;
            let tag_seq = &sequence[pos..pos + len];

            // 格式：>GCFid|tag_index|scaffold_id|pos|strand|unique
            writeln!(
                writer,
                ">{}|{}|{}|{}|0|-",
                gcf_id,
                tag_index,
                seq_id,
                pos + 1
            )?;
            writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
        }
    }

    Ok(())
}

/// 为指定分类水平构建数据库
fn build_database_for_level(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    enzyme_file: &Path,
    level: TaxonomyLevel,
    remove_redundant: bool,
    output_dir: &Path,
) -> Result<()> {
    // 第一步：记录每个标签的分类信息
    println!("第 1 步：统计标签分类信息 ...");
    let (tag_taxonomy, genome_tags) =
        collect_tag_taxonomy(genomes, enzyme_file, level, remove_redundant)?;

    // 第二步：识别特异性标签并输出
    println!("第 2 步：识别特异性标签并输出数据库 ...");
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

/// 收集标签的分类信息
fn collect_tag_taxonomy(
    genomes: &[GenomeRecord],
    enzyme_file: &Path,
    level: TaxonomyLevel,
    remove_redundant: bool,
) -> Result<(FxHashMap<Vec<u8>, FxHashSet<String>>, FxHashMap<String, FxHashMap<Vec<u8>, usize>>)>
{
    let gcf_to_taxonomy: HashMap<String, String> = genomes
        .iter()
        .map(|g| {
            let tax = &g.taxonomy[0..level as usize];
            (g.gcf_id.clone(), tax.join("\t"))
        })
        .collect();

    let mut tag_taxonomy: FxHashMap<Vec<u8>, FxHashSet<String>> = FxHashMap::default();
    let mut genome_tags: FxHashMap<String, FxHashMap<Vec<u8>, usize>> = FxHashMap::default();

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut processed_gcfs = FxHashSet::default();

    while let Some(record) = reader.next() {
        let record = record.context("解析酶切文件失败")?;

        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.is_empty() {
            continue;
        }

        let gcf_id = parts[0].trim_start_matches('>');
        if !gcf_to_taxonomy.contains_key(gcf_id) {
            continue;
        }

        processed_gcfs.insert(gcf_id.to_string());

        let mut tag_seq = record.seq().to_vec();
        tag_seq.make_ascii_uppercase();

        let taxonomy = gcf_to_taxonomy.get(gcf_id).unwrap();

        // Perl 逻辑：如果原始序列已存在就用它，否则用反向互补
        let rc_tag = reverse_complement(&tag_seq);
        let canonical_tag = if tag_taxonomy.contains_key(&tag_seq) {
            tag_seq.clone()
        } else {
            // 如果原始序列不存在，使用反向互补（与Perl版本一致）
            rc_tag.clone()
        };

        tag_taxonomy
            .entry(canonical_tag.clone())
            .or_insert_with(FxHashSet::default)
            .insert(taxonomy.clone());

        if remove_redundant {
            *genome_tags
                .entry(gcf_id.to_string())
                .or_insert_with(FxHashMap::default)
                .entry(canonical_tag)
                .or_insert(0) += 1;
        }
    }

    let percent = (processed_gcfs.len() * 100) / genomes.len();
    println!(
        "  处理了 {}/{} 个基因组 ({}%)",
        processed_gcfs.len(),
        genomes.len(),
        percent
    );

    Ok((tag_taxonomy, genome_tags))
}

/// 输出数据库
fn output_database(
    genomes: &[GenomeRecord],
    enzyme: &'static Enzyme,
    enzyme_file: &Path,
    level: TaxonomyLevel,
    tag_taxonomy: &FxHashMap<Vec<u8>, FxHashSet<String>>,
    genome_tags: &FxHashMap<String, FxHashMap<Vec<u8>, usize>>,
    remove_redundant: bool,
    output_dir: &Path,
) -> Result<()> {
    let output_path = output_dir.join(format!("{}.{}.fa.gz", enzyme.name, level.name()));
    let file = File::create(&output_path)?;
    let mut writer = GzEncoder::new(file, Compression::default());

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut unique_counts: FxHashMap<String, usize> = FxHashMap::default();
    let mut total_counts: FxHashMap<String, usize> = FxHashMap::default();

    while let Some(record) = reader.next() {
        let record = record.context("解析酶切文件失败")?;

        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() < 6 {
            continue;
        }

        let gcf_id = parts[0].trim_start_matches('>');
        let tag_index = parts[1];
        let scaffold_id = parts[2];
        let pos = parts[3];
        let original_strand = parts[4]; // 保持原始 strand

        let mut tag_seq = record.seq().to_vec();
        tag_seq.make_ascii_uppercase();

        *total_counts.entry(gcf_id.to_string()).or_insert(0) += 1;

        // Perl 逻辑：如果原始序列不在 hash 中，用反向互补并翻转 strand
        let rc_tag = reverse_complement(&tag_seq);
        let (final_tag, final_strand) = if tag_taxonomy.contains_key(&tag_seq) {
            (tag_seq.clone(), original_strand.to_string())
        } else {
            // 如果原始序列不在hash中，使用反向互补（与Perl版本一致）
            let flipped_strand = if original_strand == "0" {
                "1".to_string()
            } else {
                "0".to_string()
            };
            (rc_tag.clone(), flipped_strand)
        };

        // 检查是否为特异性标签
        let mut is_unique = if let Some(taxonomies) = tag_taxonomy.get(&final_tag) {
            taxonomies.len() == 1
        } else {
            false
        };

        // 检查基因组内冗余
        if is_unique && remove_redundant {
            if let Some(genome_tag_counts) = genome_tags.get(gcf_id) {
                if let Some(&count) = genome_tag_counts.get(&final_tag) {
                    if count > 1 {
                        is_unique = false; // 标记为非 unique
                    }
                }
            }
        }

        if is_unique {
            *unique_counts.entry(gcf_id.to_string()).or_insert(0) += 1;
        }

        // 输出标签，使用翻转后的 strand，标记 unique 状态
        let unique_flag = if is_unique { "1" } else { "0" };
        writeln!(
            writer,
            ">{}|{}|{}|{}|{}|{}",
            gcf_id, tag_index, scaffold_id, pos, final_strand, unique_flag
        )?;
        writeln!(writer, "{}", std::str::from_utf8(&final_tag).unwrap_or(""))?;
    }

    drop(writer);

    println!("  输出数据库：{}", output_path.display());
    println!(
        "  包含 {} 个基因组的特异性标签",
        unique_counts.len()
    );

    Ok(())
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

