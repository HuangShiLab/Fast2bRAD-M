use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet};
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};

/// 构建定量数据库参数
#[derive(Parser, Debug)]
pub struct BuildQuanDbArgs {
    /// 基因组分类列表文件（TSV格式）
    /// 支持传统格式（9列）或 GTDB 格式（3列：accession, gtdb_taxonomy, ncbi_taxonomy）
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

    /// 可选：预酶切目录（extract 批量输出目录，包含 genome*.fa.gz 文件）
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
    taxonomy: Vec<String>, // 完整分类信息
}

/// 主函数
pub fn run(args: BuildQuanDbArgs) -> Result<()> {
    // 解析分类层级
    let levels = parse_taxonomy_levels(&args.taxonomy_levels)?;

    // 获取酶
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num)
            .ok_or_else(|| anyhow::anyhow!("无效的酶切位点编号: {}", args.enzyme_site))?
    } else {
        enzyme_by_name(&args.enzyme_site)
            .ok_or_else(|| anyhow::anyhow!("无效的酶名称: {}", args.enzyme_site))?
    };

    println!("读取基因组分类列表 ...");
    let (genome_records, _) = read_genome_list(&args.genome_list, &levels)?;
    println!("共 {} 个基因组", genome_records.len());

    // 创建输出目录
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir.display()))?;

    // 判断是否去冗余
    let remove_redundant = args.remove_redundant.to_lowercase() == "yes";

    // 酶切基因组或读取预酶切文件/目录
    let intermediate_enzyme_file = digest_genomes_to_intermediate_file(
        &genome_records,
        enzyme,
        &args.output_dir,
        args.enzyme_file.as_ref(),
        args.pre_digested_dir.as_ref(),
    )?;

    // 为每个分类层级构建数据库
    for level in &levels {
        println!(
            "\n========== 构建 {} 级别数据库 ==========",
            level.as_str()
        );
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

/// 解析分类层级
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

/// 读取基因组分类列表
/// 支持两种格式：
/// 1. 传统格式：GCF_ID  Kingdom  Phylum  Class  Order  Family  Genus  Species  Strain  [genome_path]
/// 2. GTDB 格式：accession  gtdb_taxonomy  ncbi_taxonomy （使用 gtdb_taxonomy 列）
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
                continue;
            }

            let gcf_id = extract_gcf_id(parts[0]);
            let taxonomy = parse_gtdb_taxonomy(parts[1])?;

            genomes.push(GenomeRecord { gcf_id, taxonomy });
        } else {
            // 传统格式
            if parts.len() < 2 {
                continue;
            }

            let gcf_id = parts[0].to_string();
            let taxonomy = parts.iter().map(|s| s.to_string()).collect();

            genomes.push(GenomeRecord { gcf_id, taxonomy });
        }
    }

    Ok((genomes, taxonomy_levels_map))
}

/// 从文件名提取 GCF/GCA ID
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

/// 酶切基因组或读取预酶切文件/目录
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

    println!("开始批量酶切基因组 ...");
    let output_file = output_dir.join(format!("{}.enzyme.fa.gz", enzyme.name));

    // 简化：暂时返回错误，用户需要提供预酶切文件或目录
    bail!("请使用 -e 参数提供预酶切文件或 --pre-digested-dir 提供预酶切目录");
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
        // 构造预酶切文件路径
        let genome_id = genome.gcf_id.split('.').take(2).collect::<Vec<_>>().join(".");
        let pre_digested_file = pre_digested_dir.join(format!("{}.{}.fa.gz", genome_id, enzyme.name));

        if !pre_digested_file.exists() {
            // 尝试其他文件名模式
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
            let header = &line[1..];
            let parts: Vec<&str> = header.split('-').collect();
            if parts.len() >= 4 {
                let scaffold = parts[0];
                let start = parts[1];
                let index = parts[3];
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

/// 为指定分类层级构建数据库
fn build_database_for_level(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<()> {
    println!("第 1 步：统计标签分类信息 ...");

    // 构建 gcf_id -> taxonomy 映射
    let mut gcf_to_taxonomy = FxHashMap::default();
    for genome in genomes {
        let taxonomy_str = genome.taxonomy[1..=level as usize].join("\t");
        gcf_to_taxonomy.insert(genome.gcf_id.clone(), taxonomy_str);
    }

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

/// 统计标签分类信息（第一遍扫描）
fn collect_tag_taxonomies(
    enzyme_file: &Path,
    gcf_to_taxonomy: &FxHashMap<String, String>,
    genomes: &[GenomeRecord],
    remove_redundant: bool,
) -> Result<(
    FxHashMap<Vec<u8>, FxHashSet<String>>,
    FxHashMap<String, FxHashMap<Vec<u8>, usize>>,
)> {
    let mut tag_taxonomy: FxHashMap<Vec<u8>, FxHashSet<String>> = FxHashMap::default();
    let mut genome_tags: FxHashMap<String, FxHashMap<Vec<u8>, usize>> = FxHashMap::default();

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut processed_gcfs = FxHashSet::default();

    while let Some(record) = reader.next() {
        let record = record.context("解析酶切文件失败")?;

        let header = std::str::from_utf8(record.id()).unwrap_or("");
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() < 5 {
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

        // Perl 逻辑：先到先得
        let rc_tag = reverse_complement(&tag_seq);
        let canonical_tag = if tag_taxonomy.contains_key(&tag_seq) {
            tag_seq.clone()
        } else if tag_taxonomy.contains_key(&rc_tag) {
            rc_tag.clone()
        } else {
            tag_seq.clone()
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

/// 识别特异性标签并输出数据库（第二遍扫描）
/// **关键：定量数据库只输出 unique 标签**
fn identify_and_output_unique_tags(
    enzyme_file: &Path,
    enzyme: &'static Enzyme,
    output_dir: &Path,
    level: TaxonomyLevel,
    tag_taxonomy: &FxHashMap<Vec<u8>, FxHashSet<String>>,
    genome_tags: &FxHashMap<String, FxHashMap<Vec<u8>, usize>>,
    remove_redundant: bool,
) -> Result<()> {
    let output_path = output_dir.join(format!("{}.{}.fa.gz", enzyme.name, level.as_str()));

    let output_file = File::create(&output_path)?;
    let encoder = GzEncoder::new(output_file, Compression::default());
    let mut writer = BufWriter::new(encoder);

    let mut reader = parse_fastx_file(enzyme_file)?;
    let mut unique_counts: FxHashMap<String, usize> = FxHashMap::default();

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
        let original_strand = parts[4];

        let mut tag_seq = record.seq().to_vec();
        tag_seq.make_ascii_uppercase();

        // Perl 逻辑：判断使用哪个方向
        let rc_tag = reverse_complement(&tag_seq);
        let (final_tag, final_strand) = if tag_taxonomy.contains_key(&tag_seq) {
            (tag_seq.clone(), original_strand.to_string())
        } else if tag_taxonomy.contains_key(&rc_tag) {
            let flipped_strand = if original_strand == "0" {
                "1".to_string()
            } else {
                "0".to_string()
            };
            (rc_tag.clone(), flipped_strand)
        } else {
            (tag_seq.clone(), original_strand.to_string())
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
                        is_unique = false;
                    }
                }
            }
        }

        // **关键差异：定量数据库只输出 unique 标签**
        if is_unique {
            *unique_counts.entry(gcf_id.to_string()).or_insert(0) += 1;

            writeln!(
                writer,
                ">{}|{}|{}|{}|{}|1",
                gcf_id, tag_index, scaffold_id, pos, final_strand
            )?;
            writeln!(writer, "{}", std::str::from_utf8(&final_tag).unwrap_or(""))?;
        }
    }

    drop(writer);

    println!("  输出数据库：{}", output_path.display());
    println!("  包含 {} 个基因组的特异性标签", unique_counts.len());

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

