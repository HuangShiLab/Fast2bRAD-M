use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;
use crate::types::{DigestStats, InputType, OutputFormat, QualityControl};

#[derive(Args, Debug)]
pub struct ExtractArgs {
    /// 批量处理：样品列表 TSV（第1列=样品名，第2列=路径，第3列=路径2[可选]）
    #[arg(long = "batch")]
    pub batch_file: Option<PathBuf>,

    /// 输入文件（1-2个，支持 .gz）
    #[arg(short = 'i', long = "input", num_args = 1..=2)]
    pub input: Vec<PathBuf>,

    /// 输入类型：1=参考基因组 2=shotgun 3=单标签 4=5连标签
    #[arg(short = 't', long = "type")]
    pub input_type: u8,

    /// 酶编号（1-16）或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 输出目录
    #[arg(long = "od")]
    pub output_dir: PathBuf,

    /// 输出前缀（样本名，Type 4 需要 5 个）
    #[arg(long = "op", num_args = 1..=5)]
    pub output_prefix: Vec<String>,

    /// 是否压缩输出
    #[arg(long = "gz", default_value = "yes")]
    pub compress: String,

    /// 是否启用质量控制
    #[arg(long = "qc", default_value = "yes")]
    pub quality_control: String,

    /// 最大 N 比例
    #[arg(short = 'n', long, default_value = "0.08")]
    pub max_n: f64,

    /// 最低质量分数
    #[arg(short = 'q', long, default_value = "30")]
    pub min_quality: u8,

    /// 最低质量百分比
    #[arg(short = 'p', long, default_value = "80")]
    pub min_quality_percent: u8,

    /// 质量分数编码
    #[arg(short = 'b', long, default_value = "33")]
    pub quality_base: u8,

    /// 输出格式：fa 或 fq
    #[arg(long = "fm", default_value = "fa")]
    pub format: String,
}

pub fn run(args: ExtractArgs) -> Result<()> {
    // 批量模式：处理样品列表
    if let Some(batch_file) = args.batch_file.clone() {
        return run_batch_mode(args, &batch_file);
    }

    // 单样品模式：原有逻辑
    run_single_sample(args)
}

fn run_single_sample(args: ExtractArgs) -> Result<()> {
    // 解析参数
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let input_type = InputType::from_u8(args.input_type)
        .ok_or_else(|| anyhow!("无效的输入类型：{}", args.input_type))?;
    let output_format = OutputFormat::from_str(&args.format)
        .ok_or_else(|| anyhow!("无效的输出格式：{}", args.format))?;
    let compress = args.compress.eq_ignore_ascii_case("yes");
    let qc = QualityControl {
        enabled: args.quality_control.eq_ignore_ascii_case("yes"),
        max_n: args.max_n,
        min_quality: args.min_quality,
        min_quality_percent: args.min_quality_percent,
        quality_base: args.quality_base,
    };

    // 验证参数组合
    validate_args(&args, input_type)?;

    io_utils::ensure_directory(&args.output_dir)?;

    // 根据输入类型分派处理
    match input_type {
        InputType::ReferenceGenome => {
            extract_reference_genome(&args, enzyme, compress)?;
        }
        InputType::ShotgunMetagenome => {
            extract_shotgun(&args, enzyme, output_format, compress, &qc)?;
        }
        InputType::Single2bRAD => {
            extract_single_tag(&args, enzyme, output_format, compress, &qc)?;
        }
        InputType::Concatenated2bRAD => {
            extract_concatenated_tags(&args, enzyme, output_format, compress, &qc)?;
        }
    }

    Ok(())
}

fn run_batch_mode(base_args: ExtractArgs, batch_file: &Path) -> Result<()> {
    use std::io::{BufRead, BufReader};

    println!("### 批量处理模式：{}", batch_file.display());
    
    let file = File::open(batch_file)
        .with_context(|| format!("无法打开批量文件：{}", batch_file.display()))?;
    let reader = BufReader::new(file);
    
    let mut samples = Vec::new();
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            bail!(
                "批量文件第 {} 行格式错误：需要至少 2 列（样品名、路径）",
                line_num + 1
            );
        }
        
        let sample_name = fields[0].to_string();
        let input1 = PathBuf::from(fields[1]);
        let input2 = if fields.len() > 2 && !fields[2].is_empty() {
            Some(PathBuf::from(fields[2]))
        } else {
            None
        };
        
        samples.push((sample_name, input1, input2));
    }
    
    if samples.is_empty() {
        bail!("批量文件中没有有效样品");
    }
    
    println!("共 {} 个样品待处理", samples.len());
    
    let input_type = InputType::from_u8(base_args.input_type)
        .ok_or_else(|| anyhow!("无效的输入类型：{}", base_args.input_type))?;
    
    // 验证输入类型与文件数是否匹配
    for (idx, (name, _, input2)) in samples.iter().enumerate() {
        let file_count = if input2.is_some() { 2 } else { 1 };
        match input_type {
            InputType::ReferenceGenome | InputType::Single2bRAD => {
                if file_count > 1 {
                    bail!("样品 {} (#{})：Type {} 只支持单个输入文件", name, idx + 1, base_args.input_type);
                }
            }
            InputType::Concatenated2bRAD => {
                if file_count != 2 {
                    bail!("样品 {} (#{})：Type 4 需要 2 个输入文件（R1/R2）", name, idx + 1);
                }
            }
            InputType::ShotgunMetagenome => {
                // Type 2 支持 1-2 个文件，都可以
            }
        }
    }
    
    io_utils::ensure_directory(&base_args.output_dir)?;
    
    // 并行处理样品
    let total_samples = samples.len();
    let completed = Arc::new(AtomicUsize::new(0));
    let failed = Arc::new(AtomicUsize::new(0));
    
    samples.into_par_iter().for_each(|(sample_name, input1, input2)| {
        let completed_count = completed.load(Ordering::Relaxed) + 1;
        println!("\n### [{}/{}] 处理样品：{}", completed_count, total_samples, sample_name);
        
        // 构造每个样品的参数
        let mut sample_args = ExtractArgs {
            batch_file: None,
            input: vec![input1.clone()],
            input_type: base_args.input_type,
            enzyme_site: base_args.enzyme_site.clone(),
            output_dir: base_args.output_dir.clone(),
            output_prefix: vec![],
            compress: base_args.compress.clone(),
            quality_control: base_args.quality_control.clone(),
            max_n: base_args.max_n,
            min_quality: base_args.min_quality,
            min_quality_percent: base_args.min_quality_percent,
            quality_base: base_args.quality_base,
            format: base_args.format.clone(),
        };
        
        if let Some(input2_path) = input2 {
            sample_args.input.push(input2_path);
        }
        
        // Type 4 需要 5 个输出前缀，其他类型只需 1 个
        if input_type == InputType::Concatenated2bRAD {
            sample_args.output_prefix = (0..5)
                .map(|i| format!("{}_tag{}", sample_name, i + 1))
                .collect();
        } else {
            sample_args.output_prefix = vec![sample_name.clone()];
        }
        
        // 处理单个样品
        match run_single_sample(sample_args) {
            Ok(_) => {
                completed.fetch_add(1, Ordering::Relaxed);
                println!("✅ 样品 {} 处理完成", sample_name);
            }
            Err(e) => {
                failed.fetch_add(1, Ordering::Relaxed);
                eprintln!("❌ 样品 {} 处理失败：{}", sample_name, e);
            }
        }
    });
    
    let completed_count = completed.load(Ordering::Relaxed);
    let failed_count = failed.load(Ordering::Relaxed);
    
    println!("\n### 批量处理完成！");
    println!("成功：{} 个，失败：{} 个，总计：{} 个", 
             completed_count, failed_count, total_samples);
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
    // 尝试按名称查找
    if let Some(enzyme) = enzyme_by_name(site) {
        return Ok(enzyme);
    }

    // 尝试按编号查找
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

fn validate_args(args: &ExtractArgs, input_type: InputType) -> Result<()> {
    match input_type {
        InputType::ReferenceGenome | InputType::Single2bRAD => {
            if args.input.len() > 1 {
                bail!("Type {} 只支持单个输入文件", args.input_type);
            }
            if args.output_prefix.len() != 1 {
                bail!("Type {} 需要 1 个输出前缀", args.input_type);
            }
        }
        InputType::ShotgunMetagenome => {
            if args.input.is_empty() || args.input.len() > 2 {
                bail!("Type 2 支持 1-2 个输入文件");
            }
            if args.output_prefix.len() != 1 {
                bail!("Type 2 需要 1 个输出前缀");
            }
        }
        InputType::Concatenated2bRAD => {
            if args.input.len() != 2 {
                bail!("Type 4 需要 2 个输入文件（R1/R2）");
            }
            if args.output_prefix.len() != 5 {
                bail!("Type 4 需要 5 个输出前缀");
            }
        }
    }
    Ok(())
}

// ========== Type 1: 参考基因组 ==========

fn extract_reference_genome(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    compress: bool,
) -> Result<()> {
    let input_path = &args.input[0];
    let output_prefix = &args.output_prefix[0];

    println!("数字酶切参考基因组：{}", input_path.display());

    let output_path = build_output_path(&args.output_dir, output_prefix, enzyme, compress, "fa");
    let stat_path = args
        .output_dir
        .join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));

    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(
            File::create(&output_path)?,
            Compression::default(),
        ))
    } else {
        Box::new(File::create(&output_path)?)
    };

    let mut reader = parse_fastx_file(input_path)?;
    let mut input_sequences = 0usize;
    let mut tag_count = 0usize;
    let mut global_tag_index = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        input_sequences += 1;

        let seq_id = std::str::from_utf8(record.id()).unwrap_or("sequence");
        let seq_id = seq_id.split_whitespace().next().unwrap_or("sequence");

        let mut sequence = record.seq().to_vec();
        sequence.make_ascii_uppercase();

        // 查找所有匹配的标签位置（去重）
        let positions = enzyme.find_all_tags(&sequence);

        for (pos, len) in positions {
            global_tag_index += 1;
            tag_count += 1;
            let tag_seq = &sequence[pos..pos + len];
            let pos_end = pos + len;

            writeln!(
                writer,
                ">{}-{}-{}-{}",
                seq_id, pos + 1, pos_end, global_tag_index
            )?;
            writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
        }
    }

    drop(writer);

    // 写统计
    let stats = DigestStats {
        sample_id: output_prefix.to_string(),
        enzyme: enzyme.name.to_string(),
        input_sequences,
        tag_count,
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    println!(
        "完成：输入 {} 个序列，提取 {} 个标签 ({:.2}%)",
        input_sequences,
        tag_count,
        stats.percent()
    );

    Ok(())
}

// ========== Type 2: Shotgun ==========

fn extract_shotgun(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    output_format: OutputFormat,
    compress: bool,
    qc: &QualityControl,
) -> Result<()> {
    let input_path = &args.input[0];
    let output_prefix = &args.output_prefix[0];

    println!("提取 shotgun 数据标签：{}", input_path.display());

    // 样品文件使用 .iibsp 后缀
    let output_path = build_sample_output_path(&args.output_dir, output_prefix, enzyme, compress);
    let stat_path = args
        .output_dir
        .join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));

    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(
            File::create(&output_path)?,
            Compression::default(),
        ))
    } else {
        Box::new(File::create(&output_path)?)
    };

    let mut reader = parse_fastx_file(input_path)?;
    let mut input_sequences = 0usize;
    let mut tag_count = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        input_sequences += 1;

        let seq_id = record.id();
        let mut sequence = record.seq().to_vec();
        let quality = record.qual();

        // 质量控制
        if !qc.check_n(&sequence) {
            continue;
        }
        if let Some(qual) = quality {
            if !qc.check_quality(qual) {
                continue;
            }
        }

        sequence.make_ascii_uppercase();

        // 在序列中查找所有匹配的标签（去重）
        let positions = enzyme.find_all_tags(&sequence);

        for (pos, len) in positions {
            tag_count += 1;
            let tag_seq = &sequence[pos..pos + len];

            match output_format {
                OutputFormat::Fasta => {
                    let header = format!("@{}-{}", std::str::from_utf8(seq_id).unwrap_or("seq"), pos + 1);
                    writeln!(writer, ">{}", &header[1..])?;
                    writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
                }
                OutputFormat::Fastq => {
                    if let Some(qual) = quality {
                        let tag_qual = &qual[pos..pos + len];
                        writeln!(writer, "@{}-{}", std::str::from_utf8(seq_id).unwrap_or("seq"), pos + 1)?;
                        writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
                        writeln!(writer, "+")?;
                        writeln!(writer, "{}", std::str::from_utf8(tag_qual).unwrap_or(""))?;
                    }
                }
            }
        }
    }

    drop(writer);

    // 写统计
    let stats = DigestStats {
        sample_id: output_prefix.to_string(),
        enzyme: enzyme.name.to_string(),
        input_sequences,
        tag_count,
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    println!(
        "完成：输入 {} 个序列，提取 {} 个标签 ({:.2}%)",
        input_sequences,
        tag_count,
        stats.percent()
    );

    Ok(())
}

// ========== Type 3: 单标签 ==========

fn extract_single_tag(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    output_format: OutputFormat,
    compress: bool,
    qc: &QualityControl,
) -> Result<()> {
    let input_path = &args.input[0];
    let output_prefix = &args.output_prefix[0];

    println!("提取单 2bRAD 标签：{}", input_path.display());

    // 样品文件使用 .iibsp 后缀
    let output_path = build_sample_output_path(&args.output_dir, output_prefix, enzyme, compress);
    let stat_path = args
        .output_dir
        .join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));

    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(
            File::create(&output_path)?,
            Compression::default(),
        ))
    } else {
        Box::new(File::create(&output_path)?)
    };

    let mut reader = parse_fastx_file(input_path)?;
    let mut input_sequences = 0usize;
    let mut enzyme_reads = 0usize;
    let mut qc_passed = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        input_sequences += 1;

        let seq_id = record.id();
        let mut sequence = record.seq().to_vec();
        let quality = record.qual();

        // 截断到前 50bp（PE 数据）
        if sequence.len() > 50 {
            sequence.truncate(50);
        }

        sequence.make_ascii_uppercase();

        // 只查找第一个匹配的标签
        let mut found = false;
        for pattern in enzyme.patterns {
            if sequence.len() < enzyme.tag_length {
                break;
            }

            for offset in 0..sequence.len() {
                if offset + enzyme.tag_length > sequence.len() {
                    break;
                }

                let window = &sequence[offset..offset + enzyme.tag_length];
                if pattern.matches(window) {
                    enzyme_reads += 1;
                    let tag_seq = window;
                    let tag_qual = quality.map(|q| {
                        if offset + enzyme.tag_length <= q.len() {
                            &q[offset..offset + enzyme.tag_length]
                        } else {
                            &[]
                        }
                    });

                    // 质量控制
                    if !qc.check_n(tag_seq) {
                        found = true;
                        break;
                    }
                    if let Some(qual) = tag_qual {
                        if !qual.is_empty() && !qc.check_quality(qual) {
                            found = true;
                            break;
                        }
                    }

                    qc_passed += 1;

                    // 输出标签
                    match output_format {
                        OutputFormat::Fasta => {
                            writeln!(writer, ">{}", std::str::from_utf8(seq_id).unwrap_or("seq"))?;
                            writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
                        }
                        OutputFormat::Fastq => {
                            writeln!(writer, "@{}", std::str::from_utf8(seq_id).unwrap_or("seq"))?;
                            writeln!(writer, "{}", std::str::from_utf8(tag_seq).unwrap_or(""))?;
                            writeln!(writer, "+")?;
                            if let Some(qual) = tag_qual {
                                writeln!(writer, "{}", std::str::from_utf8(qual).unwrap_or(""))?;
                            } else {
                                writeln!(writer)?;
                            }
                        }
                    }

                    found = true;
                    break;
                }
            }

            if found {
                break;
            }
        }
    }

    drop(writer);

    // 写统计
    let mut stat_file = File::create(&stat_path)?;
    writeln!(
        stat_file,
        "sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent"
    )?;
    let percent = if input_sequences > 0 {
        (qc_passed as f64 / input_sequences as f64) * 100.0
    } else {
        0.0
    };
    writeln!(
        stat_file,
        "{}\t{}\t{}\t{}\t{}\t{:.2}",
        output_prefix, enzyme.name, input_sequences, enzyme_reads, qc_passed, percent
    )?;

    println!(
        "完成：输入 {} 个序列，命中 {} 个，质控通过 {} 个 ({:.2}%)",
        input_sequences, enzyme_reads, qc_passed, percent
    );

    Ok(())
}

// ========== Type 4: 5连标签 ==========

fn extract_concatenated_tags(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    output_format: OutputFormat,
    compress: bool,
    qc: &QualityControl,
) -> Result<()> {
    if args.input.len() != 2 {
        bail!("Type 4 需要 R1 和 R2 两个输入文件");
    }
    if args.output_prefix.len() != 5 {
        bail!("Type 4 需要 5 个输出前缀");
    }

    let r1_path = &args.input[0];
    let r2_path = &args.input[1];

    println!("处理 5 连标签数据：R1={}, R2={}", r1_path.display(), r2_path.display());
    println!("注意：Type 4 需要预先用 PEAR 等工具拼接 R1/R2，或直接提供拼接后的 FASTQ");

    // 简化实现：假设用户已经拼接好，只读取 R1（拼接后文件）
    // 完整实现需要调用外部 PEAR 或内置拼接逻辑
    let input_path = r1_path;

    // 统计原始 reads 数
    let mut raw_reads_count = 0usize;
    {
        let mut reader = parse_fastx_file(input_path)?;
        while let Some(_) = reader.next() {
            raw_reads_count += 1;
        }
    }

    // 打开 5 个输出文件
    let mut writers: Vec<Box<dyn Write>> = Vec::with_capacity(5);
    let output_ext = output_format.extension();
    
    for (i, prefix) in args.output_prefix.iter().enumerate() {
        let output_path = build_output_path(&args.output_dir, prefix, enzyme, compress, output_ext);
        let writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(
                File::create(&output_path)?,
                Compression::default(),
            ))
        } else {
            Box::new(File::create(&output_path)?)
        };
        writers.push(writer);
    }

    // 统计信息
    let mut combined_reads = 0usize;
    let mut enzyme_reads = vec![0usize; 5];
    let mut qc_passed = vec![0usize; 5];

    // 处理拼接后的序列
    let mut reader = parse_fastx_file(input_path)?;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        combined_reads += 1;

        let seq_id = record.id();
        let sequence = record.seq();
        let quality = record.qual();

        // 处理 5 个标签位置
        for tag_idx in 0..5 {
            let start = enzyme.concat_starts[tag_idx];
            let end = enzyme.concat_ends[tag_idx];

            if end > sequence.len() {
                continue; // 序列太短，跳过
            }

            let mut tag_seq = sequence[start..=end].to_vec();
            let tag_qual = quality.map(|q| {
                if end < q.len() {
                    &q[start..=end]
                } else {
                    &[]
                }
            });

            tag_seq.make_ascii_uppercase();

            // 检查是否匹配酶切位点（只取第一个匹配）
            let mut matched = false;
            for pattern in enzyme.patterns {
                // 在 tag_seq 中查找第一个匹配
                if tag_seq.len() < enzyme.tag_length {
                    break;
                }

                for offset in 0..=tag_seq.len().saturating_sub(enzyme.tag_length) {
                    let window = &tag_seq[offset..offset + enzyme.tag_length];
                    if pattern.matches(window) {
                        enzyme_reads[tag_idx] += 1;

                        // 提取匹配的标签
                        let final_tag = window;
                        let final_qual = tag_qual.and_then(|q| {
                            if !q.is_empty() && offset + enzyme.tag_length <= q.len() {
                                Some(&q[offset..offset + enzyme.tag_length])
                            } else {
                                None
                            }
                        });

                        // 质量控制
                        if !qc.check_n(final_tag) {
                            matched = true;
                            break;
                        }
                        if let Some(qual) = final_qual {
                            if !qual.is_empty() && !qc.check_quality(qual) {
                                matched = true;
                                break;
                            }
                        }

                        qc_passed[tag_idx] += 1;

                        // 输出标签
                        let writer = &mut writers[tag_idx];
                        let header = format!("{}:{}", std::str::from_utf8(seq_id).unwrap_or("seq"), tag_idx + 1);

                        match output_format {
                            OutputFormat::Fasta => {
                                writeln!(writer, ">{}", header)?;
                                writeln!(writer, "{}", std::str::from_utf8(final_tag).unwrap_or(""))?;
                            }
                            OutputFormat::Fastq => {
                                writeln!(writer, "@{}", header)?;
                                writeln!(writer, "{}", std::str::from_utf8(final_tag).unwrap_or(""))?;
                                writeln!(writer, "+")?;
                                if let Some(qual) = final_qual {
                                    writeln!(writer, "{}", std::str::from_utf8(qual).unwrap_or(""))?;
                                } else {
                                    writeln!(writer)?;
                                }
                            }
                        }

                        matched = true;
                        break;
                    }
                }

                if matched {
                    break;
                }
            }
        }
    }

    // 关闭所有写入器
    drop(writers);

    // 写统计文件
    let stat_name = args.output_prefix.join("-");
    let stat_path = args
        .output_dir
        .join(format!("{}.{}.stat.tsv", stat_name, enzyme.name));

    let mut stat_file = File::create(&stat_path)?;
    writeln!(
        stat_file,
        "sample\tenzyme\tinput_reads_num\tcombine_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent"
    )?;

    for (i, prefix) in args.output_prefix.iter().enumerate() {
        let percent = if raw_reads_count > 0 {
            (qc_passed[i] as f64 / raw_reads_count as f64) * 100.0
        } else {
            0.0
        };
        writeln!(
            stat_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
            prefix,
            enzyme.name,
            raw_reads_count,
            combined_reads,
            enzyme_reads[i],
            qc_passed[i],
            percent
        )?;
    }

    println!(
        "完成：原始 {} reads，拼接后 {} reads，5 个样本分别通过质控：{:?}",
        raw_reads_count, combined_reads, qc_passed
    );

    Ok(())
}

// ========== 辅助函数 ==========

fn build_output_path(
    output_dir: &Path,
    prefix: &str,
    enzyme: &Enzyme,
    compress: bool,
    ext: &str,
) -> PathBuf {
    let filename = if compress {
        format!("{}.{}.{}.gz", prefix, enzyme.name, ext)
    } else {
        format!("{}.{}.{}", prefix, enzyme.name, ext)
    };
    output_dir.join(filename)
}

/// 构建样品输出路径（使用 .iibsp 后缀）
fn build_sample_output_path(
    output_dir: &Path,
    prefix: &str,
    enzyme: &Enzyme,
    compress: bool,
) -> PathBuf {
    let filename = if compress {
        format!("{}.{}.iibsp.gz", prefix, enzyme.name)
    } else {
        format!("{}.{}.iibsp", prefix, enzyme.name)
    };
    output_dir.join(filename)
}
