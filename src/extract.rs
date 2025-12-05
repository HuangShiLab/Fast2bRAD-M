use std::fs::File;
use std::io::{Write, BufReader, BufWriter, BufRead}; // Added BufRead explicitly
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::process::Command;
use std::hash::Hasher;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use fxhash::FxHasher;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;
use crate::types::{DigestStats, InputType, QualityControl};

// ================== 类型定义 ==================

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
    // 比较字典序，取较小者
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

// ================== 参数定义 ==================

#[derive(Args, Debug, Clone)]
pub struct ExtractArgs {
    /// 批量处理：样品列表 TSV
    #[arg(long = "batch")]
    pub batch_file: Option<PathBuf>,

    /// 输入文件
    #[arg(short = 'i', long = "input", num_args = 1..=2)]
    pub input: Vec<PathBuf>,

    /// 输入类型：1=参考基因组 2=shotgun 3=单标签 4=5连标签
    #[arg(short = 't', long = "type")]
    pub input_type: u8,

    /// 酶编号或名称
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,

    /// 输出目录
    #[arg(long = "od")]
    pub output_dir: PathBuf,

    /// 输出前缀
    #[arg(long = "op", num_args = 1..=5)]
    pub output_prefix: Vec<String>,

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

    // 已移除 output format 参数

    /// PEAR 可执行文件路径
    #[arg(long = "pe", default_value = "pear")]
    pub pear_bin: String,

    /// PEAR 线程数
    #[arg(long = "pc", default_value = "1")]
    pub pear_cpu: usize,
}

pub fn run(args: ExtractArgs) -> Result<()> {
    if let Some(batch_file) = args.batch_file.clone() {
        return run_batch_mode(args, &batch_file);
    }
    run_single_sample(args)
}

fn run_single_sample(args: ExtractArgs) -> Result<()> {
    let enzyme = parse_enzyme(&args.enzyme_site)?;
    let input_type = InputType::from_u8(args.input_type)
        .ok_or_else(|| anyhow!("无效的输入类型：{}", args.input_type))?;
    
    let qc = QualityControl {
        enabled: args.quality_control.eq_ignore_ascii_case("yes"),
        max_n: args.max_n,
        min_quality: args.min_quality,
        min_quality_percent: args.min_quality_percent,
        quality_base: args.quality_base,
    };

    validate_args(&args, input_type)?;
    io_utils::ensure_directory(&args.output_dir)?;

    match input_type {
        InputType::ReferenceGenome => {
            // Type 1 -> .iibdb (Text format: >ID\nHash)
            extract_reference_genome(&args, enzyme)?;
        }
        InputType::ShotgunMetagenome => {
            // Type 2 -> .iibsp (Text format: >ID\nHash)
            extract_shotgun(&args, enzyme, &qc)?;
        }
        InputType::Single2bRAD => {
            // Type 3 -> .iibsp (Text format: >ID\nHash)
            extract_single_tag(&args, enzyme, &qc)?;
        }
        InputType::Concatenated2bRAD => {
            // Type 4 -> .iibsp (Text format: >ID\nHash)
            extract_concatenated_tags(&args, enzyme, &qc)?;
        }
    }

    Ok(())
}

fn run_batch_mode(base_args: ExtractArgs, batch_file: &Path) -> Result<()> {
    println!("### 批量处理模式：{}", batch_file.display());
    
    let file = File::open(batch_file)
        .with_context(|| format!("无法打开批量文件：{}", batch_file.display()))?;
    let reader = BufReader::new(file);
    
    let mut samples = Vec::new();
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            eprintln!("Warning: Skipping invalid line {}: {}", line_num + 1, line);
            continue;
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

    let completed = Arc::new(AtomicUsize::new(0));
    let total = samples.len();
    eprintln!("DEBUG: 开始并行处理 {} 个样本", total);

    samples.into_par_iter().for_each(|(sample_name, input1, input2)| {
        let mut sample_args = base_args.clone();
        
        sample_args.batch_file = None;
        sample_args.input = vec![input1];
        if let Some(in2) = input2 {
            sample_args.input.push(in2);
        }
        
        if base_args.input_type == 4 { 
             sample_args.output_prefix = (0..5)
                .map(|i| format!("{}_tag{}", sample_name, i + 1))
                .collect();
        } else {
            sample_args.output_prefix = vec![sample_name.clone()];
        }

        eprintln!("DEBUG: 线程启动处理样本: {}", sample_name);

        match run_single_sample(sample_args) {
            Ok(_) => {
                let count = completed.fetch_add(1, Ordering::Relaxed) + 1;
                println!("✅ [{}/{}] 样品 {} 处理完成", count, total, sample_name);
            }
            Err(e) => {
                eprintln!("❌ 样品 {} 处理失败: {:?}", sample_name, e);
            }
        }
    });
    
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
     if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme); }
     if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme); } }
     bail!("未知的酶：{}", site)
}

fn validate_args(_args: &ExtractArgs, _input_type: InputType) -> Result<()> {
    Ok(())
}

// ========== Type 1: 参考基因组 -> .iibdb (FASTA-like) ==========

fn extract_reference_genome(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
) -> Result<()> {
    let input_path = &args.input[0];
    let output_prefix = &args.output_prefix[0];

    println!("数字酶切参考基因组 (Hash模式)：{}", input_path.display());

    // 生成 .iibdb 文件 (文本格式)
    let output_filename = format!("{}.{}.iibdb", output_prefix, enzyme.name);
    let output_path = args.output_dir.join(output_filename);

    // 打开输出流
    let file = File::create(&output_path).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);

    let mut reader = parse_fastx_file(input_path)
        .context(format!("Failed to open file: {:?}", input_path))?;

    let mut input_sequences = 0usize;
    let mut total_tags = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("读取 Fastx 记录失败")?;
        input_sequences += 1;

        let raw_id = record.id();
        let seq_id = std::str::from_utf8(raw_id)
            .unwrap_or("sequence")
            .split_whitespace()
            .next()
            .unwrap_or("sequence");
        
        let mut sequence = record.seq().to_vec();
        sequence.make_ascii_uppercase();

        let positions_iter = enzyme.find_all_tags(&sequence);
        
        for (pos, len) in positions_iter {
            if pos + len > sequence.len() {
                continue; 
            }

            let tag_seq = &sequence[pos..pos + len];
            let canonical = get_canonical_sequence(tag_seq);
            let hash = hash_bytes(&canonical);
            
            // 写入格式: >ID:Pos\nHash
            writeln!(writer, ">{}:{}\n{}", seq_id, pos, hash)?;
            total_tags += 1;
        }
    }

    // 写统计
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));
    let stats = DigestStats {
        sample_id: output_prefix.to_string(),
        enzyme: enzyme.name.to_string(),
        input_sequences,
        tag_count: total_tags,
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    println!("完成：生成 iibdb 文件 {}，含 {} 个序列，{} 个标签", output_path.display(), input_sequences, total_tags);
    Ok(())
}

// ========== Type 2: Shotgun -> .iibsp (FASTA-like) ==========

fn extract_shotgun(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    qc: &QualityControl,
) -> Result<()> {
    // 若为双端输入（2个文件），优先使用 PEAR 合并
    let mut merged_path: Option<PathBuf> = None;
    if args.input.len() == 2 {
        merged_path = Some(run_pear_and_combine(args, enzyme)?);
    }

    let input_path = merged_path.as_ref().unwrap_or(&args.input[0]);
    let output_prefix = &args.output_prefix[0];

    println!("提取 shotgun 标签 (Hash模式)：{}", input_path.display());

    // 生成 .iibsp 文件 (文本格式)
    let output_filename = format!("{}.{}.iibsp", output_prefix, enzyme.name);
    let output_path = args.output_dir.join(output_filename);

    let file = File::create(&output_path).context("无法创建输出文件")?;
    let mut writer = BufWriter::new(file);

    let mut reader = parse_fastx_file(input_path)?;
    let mut input_sequences = 0usize;
    let mut tag_count = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        input_sequences += 1;

        let seq_id = String::from_utf8_lossy(record.id());
        let mut sequence = record.seq().to_vec();
        let quality = record.qual();

        if !qc.check_n(&sequence) { continue; }
        if let Some(qual) = quality { if !qc.check_quality(qual) { continue; } }

        sequence.make_ascii_uppercase();

        let positions = enzyme.find_all_tags(&sequence);

        for (i, (pos, len)) in positions.iter().enumerate() {
            tag_count += 1;
            let tag_seq = &sequence[*pos..*pos + len];
            let canonical = get_canonical_sequence(tag_seq);
            let tag_hash = hash_bytes(&canonical);

            // 写入格式: >ID_tagIndex\nHash
            writeln!(writer, ">{}_tag{}\n{}", seq_id, i + 1, tag_hash)?;
        }
    }

    // 统计
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));
    let stats = DigestStats {
        sample_id: output_prefix.to_string(),
        enzyme: enzyme.name.to_string(),
        input_sequences,
        tag_count,
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    println!("完成：生成 iibsp 文件，提取 {} 个标签", tag_count);
    Ok(())
}

fn run_pear_and_combine(args: &ExtractArgs, enzyme: &Enzyme) -> Result<PathBuf> {
    let r1 = &args.input[0];
    let r2 = &args.input[1];
    let prefix = &args.output_prefix[0];

    // 输出前缀：<outdir>/<prefix>.<enzyme>
    let base = args.output_dir.join(format!("{}.{}", prefix, enzyme.name));

    // 1) 运行 PEAR
    let status = Command::new(&args.pear_bin)
        .args([
            "-f",
            r1.to_str().unwrap(),
            "-r",
            r2.to_str().unwrap(),
            "-e",
            "-o",
            base.to_str().unwrap(),
            "-j",
            &args.pear_cpu.to_string(),
        ])
        .status()
        .with_context(|| format!("无法执行 PEAR：{}", args.pear_bin))?;
    if !status.success() {
        bail!("PEAR 执行失败，退出码：{}", status.code().unwrap_or(-1));
    }

    // 2) 合并三个输出为 .pear.fastq (不压缩)
    let assembled = args.output_dir.join(format!("{}.{}.assembled.fastq", prefix, enzyme.name));
    let unassembled_f = args.output_dir.join(format!("{}.{}.unassembled.forward.fastq", prefix, enzyme.name));
    let unassembled_r = args.output_dir.join(format!("{}.{}.unassembled.reverse.fastq", prefix, enzyme.name));

    let pear_fastq = args.output_dir.join(format!("{}.{}.pear.fastq", prefix, enzyme.name));
    {
        let mut out = File::create(&pear_fastq)?;
        for p in [&assembled, &unassembled_f, &unassembled_r] {
            if p.exists() {
                let mut f = File::open(p)?;
                std::io::copy(&mut f, &mut out)?;
            }
        }
    }

    // 3) 清理中间文件
    let _ = std::fs::remove_file(assembled);
    let _ = std::fs::remove_file(unassembled_f);
    let _ = std::fs::remove_file(unassembled_r);
    let discarded_path = format!("{}.discarded.fastq", base.to_str().unwrap());
    let _ = std::fs::remove_file(discarded_path);

    Ok(pear_fastq)
}

// ========== Type 3: 单标签 -> .iibsp (FASTA-like) ==========
fn extract_single_tag(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    qc: &QualityControl,
) -> Result<()> {
    let input_path = &args.input[0];
    let output_prefix = &args.output_prefix[0];

    println!("提取单 2bRAD 标签 (Hash模式)：{}", input_path.display());

    let output_filename = format!("{}.{}.iibsp", output_prefix, enzyme.name);
    let output_path = args.output_dir.join(output_filename);
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", output_prefix, enzyme.name));

    let file = File::create(&output_path).context("无法创建输出文件")?;
    let mut writer = BufWriter::new(file);

    let mut reader = parse_fastx_file(input_path)?;
    
    let mut input_sequences = 0usize;
    let mut enzyme_reads = 0usize;
    let mut qc_passed = 0usize;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        input_sequences += 1;

        let seq_id = String::from_utf8_lossy(record.id());
        let mut sequence = record.seq().to_vec();
        let quality = record.qual();

        if sequence.len() > 50 {
            sequence.truncate(50);
        }

        sequence.make_ascii_uppercase();

        let mut found = false;
        for pattern in enzyme.patterns {
            if sequence.len() < enzyme.tag_length { break; }

            for offset in 0..sequence.len() {
                if offset + enzyme.tag_length > sequence.len() { break; }

                let window = &sequence[offset..offset + enzyme.tag_length];
                if pattern.matches(window) {
                    enzyme_reads += 1;
                    let tag_seq = window;
                    
                    let qual_slice = quality.map(|q| {
                        if offset + enzyme.tag_length <= q.len() {
                            &q[offset..offset + enzyme.tag_length]
                        } else {
                            &[]
                        }
                    });

                    if !qc.check_n(tag_seq) { found = true; break; }
                    if let Some(q) = qual_slice {
                        if !q.is_empty() && !qc.check_quality(q) { found = true; break; }
                    }

                    qc_passed += 1;

                    let canonical = get_canonical_sequence(tag_seq);
                    let tag_hash = hash_bytes(&canonical);

                    // 写入格式: >ID\nHash
                    writeln!(writer, ">{}\n{}", seq_id, tag_hash)?;

                    found = true;
                    break;
                }
            }
            if found { break; }
        }
    }

    // 写统计
    let mut stat_file = File::create(&stat_path)?;
    writeln!(stat_file, "sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent")?;
    let percent = if input_sequences > 0 { (qc_passed as f64 / input_sequences as f64) * 100.0 } else { 0.0 };
    writeln!(stat_file, "{}\t{}\t{}\t{}\t{}\t{:.2}%", output_prefix, enzyme.name, input_sequences, enzyme_reads, qc_passed, percent)?;

    println!("完成：输入 {} 个序列，命中 {} 个，质控通过 {} 个 ({:.2}%)", input_sequences, enzyme_reads, qc_passed, percent);
    Ok(())
}

// ========== Type 4: 5连标签 -> .iibsp (FASTA-like) ==========

fn extract_concatenated_tags(
    args: &ExtractArgs,
    enzyme: &'static Enzyme,
    qc: &QualityControl,
) -> Result<()> {
    if args.input.len() != 2 {
        bail!("Type 4 需要 R1 和 R2 两个输入文件");
    }
    if args.output_prefix.len() != 5 {
        bail!("Type 4 需要 5 个输出前缀");
    }

    let r1_path = &args.input[0];
    println!("处理 5 连标签数据 (Hash模式)：R1={}", r1_path.display());
    let input_path = r1_path;

    // 预先打开 5 个输出文件 writers
    let mut writers = Vec::new();
    for i in 0..5 {
        let prefix = &args.output_prefix[i];
        let output_filename = format!("{}.{}.iibsp", prefix, enzyme.name);
        let output_path = args.output_dir.join(output_filename);
        let file = File::create(&output_path).context(format!("无法创建文件 {:?}", output_path))?;
        writers.push(BufWriter::new(file));
    }

    // 统计总 reads 数
    let mut raw_reads_count = 0usize;
    {
        let mut reader = parse_fastx_file(input_path)?;
        while let Some(_) = reader.next() {
            raw_reads_count += 1;
        }
    }

    let mut combined_reads = 0usize;
    let mut enzyme_reads = vec![0usize; 5];
    let mut qc_passed = vec![0usize; 5];

    let mut reader = parse_fastx_file(input_path)?;

    while let Some(record) = reader.next() {
        let record = record.context("解析序列记录失败")?;
        combined_reads += 1;

        let seq_id = String::from_utf8_lossy(record.id());
        let sequence = record.seq();
        let quality = record.qual();

        for tag_idx in 0..5 {
            let start = enzyme.concat_starts[tag_idx];
            let end = enzyme.concat_ends[tag_idx];

            if end > sequence.len() { continue; }

            let mut tag_seq = sequence[start..=end].to_vec();
            tag_seq.make_ascii_uppercase();

            let mut matched = false;
            for pattern in enzyme.patterns {
                if tag_seq.len() < enzyme.tag_length { break; }

                for offset in 0..=tag_seq.len().saturating_sub(enzyme.tag_length) {
                    let window = &tag_seq[offset..offset + enzyme.tag_length];
                    if pattern.matches(window) {
                        enzyme_reads[tag_idx] += 1;
                        
                        let q_start = start + offset;
                        let q_end = q_start + enzyme.tag_length;
                        
                        let qual_slice = quality.and_then(|q| {
                            if q_end <= q.len() { Some(&q[q_start..q_end]) } else { None }
                        });

                        if !qc.check_n(window) { matched = true; break; }
                        if let Some(q) = qual_slice {
                            if !q.is_empty() && !qc.check_quality(q) { matched = true; break; }
                        }

                        qc_passed[tag_idx] += 1;

                        let canonical = get_canonical_sequence(window);
                        let tag_hash = hash_bytes(&canonical);
                        
                        // 写入对应的 writer
                        // 写入格式: >ID:TagIdx\nHash
                        let writer = &mut writers[tag_idx];
                        writeln!(writer, ">{}:{}\n{}", seq_id, tag_idx + 1, tag_hash)?;

                        matched = true;
                        break;
                    }
                }
                if matched { break; }
            }
        }
    }

    // 统计报告
    let stat_name = args.output_prefix.join("-");
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", stat_name, enzyme.name));

    let mut stat_file = File::create(&stat_path)?;
    writeln!(stat_file, "sample\tenzyme\tinput_reads_num\tcombine_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent")?;

    for (i, prefix) in args.output_prefix.iter().enumerate() {
        let percent = if raw_reads_count > 0 { (qc_passed[i] as f64 / raw_reads_count as f64) * 100.0 } else { 0.0 };
        writeln!(stat_file, "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}%", prefix, enzyme.name, raw_reads_count, combined_reads, enzyme_reads[i], qc_passed[i], percent)?;
    }

    println!("完成：原始 {} reads，拼接后 {} reads，5 个样本分别通过质控：{:?}", raw_reads_count, combined_reads, qc_passed);

    Ok(())
}