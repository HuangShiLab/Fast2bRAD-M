use anyhow::{bail, Context, Result};
use clap::Parser;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

/// 功能丰度预测：t(物种丰度表) * 物种功能矩阵 = 功能丰度表
/// 结果会对每个样品的列进行归一化，使得每个样品的功能总丰度之和为 1。
#[derive(Parser, Debug, Clone)]
pub struct PredictArgs {
    /// 物种丰度表（merge 输出的 .all.xls 文件）
    #[arg(short = 'a', long = "abundance")]
    pub abundance_file: PathBuf,

    /// 物种-功能映射矩阵（TSV：第一行为表头，第一列为物种名，其余列为功能ID，值为基因拷贝数）
    #[arg(short = 'm', long = "mapping")]
    pub mapping_file: PathBuf,

    /// 输出目录
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,

    /// 输出文件前缀
    #[arg(short = 'p', long = "prefix", default_value = "Abundance_Stat")]
    pub prefix: String,
}

pub fn run(args: PredictArgs) -> Result<()> {
    tracing::info!(
        "COMMAND: predict -a {} -m {} -o {} -p {}",
        args.abundance_file.display(),
        args.mapping_file.display(),
        args.output_dir.display(),
        args.prefix
    );

    std::fs::create_dir_all(&args.output_dir)?;

    // Step 1: 加载物种丰度矩阵
    let (samples, species_abundance) = load_abundance_matrix(&args.abundance_file)?;
    tracing::info!(
        "加载物种丰度表完成：{} 个物种，{} 个样品",
        species_abundance.len(),
        samples.len()
    );
    if species_abundance.is_empty() {
        bail!("物种丰度表中未读取到有效数据");
    }

    // Step 2: 加载物种-功能映射矩阵
    let (func_names, species_func_map) = load_mapping_matrix(&args.mapping_file)?;
    tracing::info!(
        "加载功能映射矩阵完成：{} 个功能项，{} 个物种",
        func_names.len(),
        species_func_map.len()
    );
    if func_names.is_empty() {
        bail!("功能映射矩阵中未找到功能列");
    }

    // Step 3: 计算功能丰度矩阵
    // func_matrix[sample_idx][func_idx] = Σ_species (abundance[species][sample] × mapping[species][func])
    let n_samples = samples.len();
    let n_funcs = func_names.len();
    let mut func_matrix: Vec<Vec<f64>> = vec![vec![0.0; n_funcs]; n_samples];

    let mut matched = 0usize;
    for (species, abundance_vec) in &species_abundance {
        // 尝试匹配物种名
        if let Some(func_vec) = species_func_map.get(species.as_str()) {
            matched += 1;
            for (s_idx, &abund) in abundance_vec.iter().enumerate() {
                if abund == 0.0 {
                    continue;
                }
                for (f_idx, &count) in func_vec.iter().enumerate() {
                    func_matrix[s_idx][f_idx] += abund * count;
                }
            }
        }
    }

    tracing::info!(
        "物种匹配：{}/{} 个物种在映射矩阵中找到对应记录",
        matched,
        species_abundance.len()
    );

    if matched == 0 {
        tracing::warn!("警告：物种丰度表与映射矩阵无任何物种匹配，请检查物种名称是否一致");
    }

    // --- 核心更新：样品列归一化 (Sum normalization per sample) ---
    tracing::info!("正在对样品功能丰度进行归一化（Sum to 1）...");
    for s_idx in 0..n_samples {
        let sample_total_func_abundance: f64 = func_matrix[s_idx].iter().sum();
        if sample_total_func_abundance > 0.0 {
            for f_idx in 0..n_funcs {
                func_matrix[s_idx][f_idx] /= sample_total_func_abundance;
            }
        }
    }

    // Step 4: 输出功能丰度表
    let output_path = args.output_dir.join(format!("{}.func.xls", args.prefix));
    let file = File::create(&output_path)
        .with_context(|| format!("无法创建输出文件: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // 写入表头：#Function\tsample1\t...\tsampleN
    write!(writer, "#Function")?;
    for sample in &samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // 写入数据行
    let mut written = 0usize;
    for (f_idx, func_name) in func_names.iter().enumerate() {
        // 检查该功能是否在任一样品中非零
        let has_nonzero = (0..n_samples).any(|s| func_matrix[s][f_idx] > 0.0);
        if !has_nonzero {
            continue;
        }
        write!(writer, "{}", func_name)?;
        for s_idx in 0..n_samples {
            write!(writer, "\t{:.8}", func_matrix[s_idx][f_idx])?;
        }
        writeln!(writer)?;
        written += 1;
    }

    tracing::info!(
        "功能丰度表已输出：{} 个有效功能项 -> {}",
        written,
        output_path.display()
    );
    Ok(())
}

fn load_abundance_matrix(path: &Path) -> Result<(Vec<String>, FxHashMap<String, Vec<f64>>)> {
    let file = File::open(path)
        .with_context(|| format!("无法打开物种丰度文件: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_raw = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("物种丰度文件为空: {}", path.display()))??;
    let header_str = header_raw.trim().trim_start_matches('#');
    let headers: Vec<&str> = header_str.split('\t').collect();

    const TAX_LEVELS: &[&str] = &[
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain",
    ];

    let mut sample_start = 0usize;
    for (i, h) in headers.iter().enumerate() {
        if TAX_LEVELS.contains(h) {
            sample_start = i + 1;
        } else if sample_start > 0 {
            break;
        }
    }

    if sample_start == 0 || sample_start >= headers.len() {
        bail!("无法解析物种丰度文件表头，未找到分类列或样品列: {}", path.display());
    }

    let samples: Vec<String> = headers[sample_start..].iter().map(|s| s.to_string()).collect();
    let n_samples = samples.len();
    let mut species_abundance: FxHashMap<String, Vec<f64>> = FxHashMap::default();

    for line_res in lines {
        let line = line_res?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < sample_start + n_samples { continue; }
        
        let species_name = fields[sample_start - 1].trim().to_string();
        let abundances: Vec<f64> = fields[sample_start..sample_start + n_samples]
            .iter()
            .map(|v| v.parse::<f64>().unwrap_or(0.0))
            .collect();
        species_abundance.insert(species_name, abundances);
    }

    Ok((samples, species_abundance))
}

fn load_mapping_matrix(path: &Path) -> Result<(Vec<String>, FxHashMap<String, Vec<f64>>)> {
    let file = File::open(path)
        .with_context(|| format!("无法打开功能映射矩阵文件: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_raw = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("功能映射矩阵文件为空: {}", path.display()))??;
    let header_str = header_raw.trim().trim_start_matches('#');
    let headers: Vec<&str> = header_str.split('\t').collect();

    if headers.len() < 2 {
        bail!("功能映射矩阵格式错误：至少需要 2 列: {}", path.display());
    }

    let func_names: Vec<String> = headers[1..].iter().map(|s| s.to_string()).collect();
    let n_funcs = func_names.len();
    let mut species_map: FxHashMap<String, Vec<f64>> = FxHashMap::default();

    for line_res in lines {
        let line = line_res?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() { continue; }
        
        let species = fields[0].trim().to_string();
        let mut counts: Vec<f64> = fields[1..]
            .iter()
            .map(|v| v.parse::<f64>().unwrap_or(0.0))
            .collect();
        counts.resize(n_funcs, 0.0);
        species_map.insert(species, counts);
    }

    Ok((func_names, species_map))
}