use anyhow::{Context, Result, bail};
use clap::Parser;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

/// 合并多样品丰度表
#[derive(Parser, Debug)]
pub struct MergeArgs {
    /// 样品列表文件（TSV格式：sample_name<tab>profile_path）
    #[arg(short = 'l', long = "list")]
    pub sample_list: PathBuf,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// 输出文件前缀
    #[arg(short = 'p', long = "prefix", default_value = "Abundance_Stat")]
    pub prefix: String,

    /// Mock样品名称（逗号分隔，用于过滤）
    #[arg(short = 'm', long = "mock")]
    pub mock_samples: Option<String>,

    /// 阴性对照样品名称（逗号分隔，用于过滤）
    #[arg(short = 'c', long = "control")]
    pub control_samples: Option<String>,
}

/// 物种丰度数据
#[derive(Debug, Default)]
struct TaxonAbundance {
    /// 样品名 -> 相对丰度值
    samples: FxHashMap<String, f64>,
}

pub fn run(args: MergeArgs) -> Result<()> {
    tracing::info!(
        "COMMAND: merge -l {} -o {} -p {}",
        args.sample_list.display(),
        args.output_dir.display(),
        args.prefix
    );

    // 创建输出目录
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir.display()))?;

    // 解析mock和control样品列表
    let mock_set: FxHashMap<String, ()> = if let Some(ref mock) = args.mock_samples {
        mock.split(',').map(|s| (s.trim().to_string(), ())).collect()
    } else {
        FxHashMap::default()
    };

    let control_set: FxHashMap<String, ()> = if let Some(ref control) = args.control_samples {
        control.split(',').map(|s| (s.trim().to_string(), ())).collect()
    } else {
        FxHashMap::default()
    };

    // 读取所有样品的丰度数据
    let (taxa_abundance, sample_order, header) = read_all_profiles(&args.sample_list)?;
    
    if sample_order.is_empty() {
        bail!("未找到有效的样品丰度文件");
    }

    tracing::info!("共读取 {} 个样品", sample_order.len());

    // 输出合并后的丰度表（all.xls）
    let all_output = args.output_dir.join(format!("{}.all.xls", args.prefix));
    write_merged_table(&all_output, &taxa_abundance, &sample_order, &header)?;
    tracing::info!("✅ 合并表已输出：{}", all_output.display());

    // 输出过滤后的丰度表（filtered.xls）
    let filtered_output = args.output_dir.join(format!("{}.filtered.xls", args.prefix));
    write_filtered_table(
        &filtered_output,
        &taxa_abundance,
        &sample_order,
        &header,
        &mock_set,
        &control_set,
    )?;
    tracing::info!("✅ 过滤表已输出：{}", filtered_output.display());

    tracing::info!("\n全部完成！");
    Ok(())
}

/// 读取所有样品的丰度数据
fn read_all_profiles(
    list_file: &Path,
) -> Result<(FxHashMap<String, TaxonAbundance>, Vec<String>, String)> {
    let file = File::open(list_file)
        .with_context(|| format!("无法打开样品列表: {}", list_file.display()))?;
    let reader = BufReader::new(file);

    let mut taxa_abundance: FxHashMap<String, TaxonAbundance> = FxHashMap::default();
    let mut sample_order = Vec::new();
    let mut sample_totals: FxHashMap<String, f64> = FxHashMap::default();
    let mut header = String::new();
    let mut classify_col = 0;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        let sample_name = parts[0].to_string();
        let profile_path = Path::new(parts[1]);

        if !profile_path.exists() {
            tracing::warn!("警告：样品 {} 的丰度文件不存在：{}", sample_name, profile_path.display());
            continue;
        }

        sample_order.push(sample_name.clone());

        // 读取样品丰度文件
        let profile_file = File::open(profile_path)
            .with_context(|| format!("无法打开丰度文件: {}", profile_path.display()))?;
        let profile_reader = BufReader::new(profile_file);

        for line in profile_reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.starts_with("#Kingdom") || line.starts_with("Kingdom") {
                // 解析表头，确定分类列的范围
                let fields: Vec<&str> = line.trim_start_matches('#').split('\t').collect();
                for (i, field) in fields.iter().enumerate() {
                    if *field == "Theoretical_Tag_Num" {
                        classify_col = i - 1;
                        header = fields[0..=classify_col].join("\t");
                        break;
                    }
                }
                continue;
            }

            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() <= classify_col {
                continue;
            }

            // 提取分类ID
            let taxon_id = fields[0..=classify_col].join("\t");
            
            // 提取 Sequenced_Reads_Num/Theoretical_Tag_Num（倒数第4列）
            if fields.len() >= 4 {
                let abundance_value: f64 = fields[fields.len() - 4]
                    .parse()
                    .unwrap_or(0.0);

                // 记录该分类在该样品中的丰度值
                taxa_abundance
                    .entry(taxon_id.clone())
                    .or_insert_with(TaxonAbundance::default)
                    .samples
                    .insert(sample_name.clone(), abundance_value);

                // 累加样品总量
                *sample_totals.entry(sample_name.clone()).or_insert(0.0) += abundance_value;
            }
        }
    }

    // 归一化：计算相对丰度（每个分类占样品总量的百分比）
    for (_taxon_id, abundance) in taxa_abundance.iter_mut() {
        for (sample_name, value) in abundance.samples.iter_mut() {
            if let Some(&total) = sample_totals.get(sample_name) {
                if total > 0.0 {
                    *value = *value / total;
                }
            }
        }
    }

    Ok((taxa_abundance, sample_order, header))
}

/// 输出合并后的丰度表（all.xls）
fn write_merged_table(
    output_path: &Path,
    taxa_abundance: &FxHashMap<String, TaxonAbundance>,
    sample_order: &[String],
    header: &str,
) -> Result<()> {
    let file = File::create(output_path)
        .with_context(|| format!("无法创建输出文件: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // 写入表头
    write!(writer, "{}", header)?;
    for sample in sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // 写入数据（按分类ID排序）
    let mut taxa_list: Vec<&String> = taxa_abundance.keys().collect();
    taxa_list.sort();

    for taxon_id in taxa_list {
        let abundance = &taxa_abundance[taxon_id];
        
        // 检查是否至少有一个样品中存在该分类
        let has_data = sample_order.iter().any(|sample| {
            abundance.samples.get(sample).map_or(false, |&v| v > 0.0)
        });

        if !has_data {
            continue;
        }

        write!(writer, "{}", taxon_id)?;
        for sample in sample_order {
            let value = abundance.samples.get(sample).copied().unwrap_or(0.0);
            if value == 0.0 {
                write!(writer, "\t0")?;
            } else {
                write!(writer, "\t{}", value)?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// 输出过滤后的丰度表（filtered.xls）
fn write_filtered_table(
    output_path: &Path,
    taxa_abundance: &FxHashMap<String, TaxonAbundance>,
    sample_order: &[String],
    header: &str,
    mock_set: &FxHashMap<String, ()>,
    control_set: &FxHashMap<String, ()>,
) -> Result<()> {
    let file = File::create(output_path)
        .with_context(|| format!("无法创建输出文件: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // 过滤样品列表（移除mock和control）
    let filtered_samples: Vec<String> = sample_order
        .iter()
        .filter(|s| !mock_set.contains_key(*s) && !control_set.contains_key(*s))
        .cloned()
        .collect();

    if filtered_samples.is_empty() {
        tracing::warn!("警告：过滤后没有剩余样品");
        // 创建空文件
        return Ok(());
    }

    // 收集在control样品中检测到的分类（潜在污染）
    let mut contamination_taxa: FxHashMap<String, ()> = FxHashMap::default();
    for (taxon_id, abundance) in taxa_abundance.iter() {
        for control_sample in control_set.keys() {
            if abundance.samples.get(control_sample).map_or(false, |&v| v > 0.0) {
                contamination_taxa.insert(taxon_id.clone(), ());
                break;
            }
        }
    }

    // 写入表头
    write!(writer, "{}", header)?;
    for sample in &filtered_samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // 写入数据
    let mut taxa_list: Vec<&String> = taxa_abundance.keys().collect();
    taxa_list.sort();

    for taxon_id in taxa_list {
        // 跳过潜在污染的分类
        if contamination_taxa.contains_key(taxon_id) {
            continue;
        }

        let abundance = &taxa_abundance[taxon_id];
        
        // 检查在过滤后的样品中是否有数据
        let has_data = filtered_samples.iter().any(|sample| {
            abundance.samples.get(sample).map_or(false, |&v| v > 0.0)
        });

        if !has_data {
            continue;
        }

        write!(writer, "{}", taxon_id)?;
        for sample in &filtered_samples {
            let value = abundance.samples.get(sample).copied().unwrap_or(0.0);
            if value == 0.0 {
                write!(writer, "\t0")?;
            } else {
                write!(writer, "\t{}", value)?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

