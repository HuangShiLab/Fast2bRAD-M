use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::process::Command;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord; 
use rayon::prelude::*;
use fxhash::FxHasher;
use std::hash::Hasher;
use tracing;

use crate::enzymes::{Enzyme, enzyme_by_id, enzyme_by_name};
use crate::io_utils;
use crate::types::{DigestStats, InputType, QualityControl};

const BATCH_SIZE: usize = 10000;
const CHANNEL_BUFFER: usize = 16;

pub type Hash = u64;

pub fn hash_bytes(bytes: &[u8]) -> Hash {
    let mut hasher = FxHasher::default();
    hasher.write(bytes);
    hasher.finish()
}

fn get_canonical_sequence(seq: &[u8]) -> Vec<u8> {
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
    
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

// 【优化】RawRecord 结构调整
// 1. id 改为 Vec<u8>，避免 parsing 阶段的 UTF-8 校验和 String 转换开销
// 2. 字段全部保留 Capacity 以供复用
#[derive(Debug, Clone)]
struct RawRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>, // 使用空 Vec 表示 None，避免 Option 的解包开销
}

impl RawRecord {
    fn new() -> Self {
        Self {
            id: Vec::with_capacity(64),
            seq: Vec::with_capacity(150),
            qual: Vec::with_capacity(150),
        }
    }

    // 【核心优化】内存复用逻辑
    // 不分配新内存，直接拷贝数据到已有 buffer
    fn populate_from(&mut self, rec: &SequenceRecord) {
        self.id.clear();
        self.id.extend_from_slice(rec.id());

        self.seq.clear();
        // 【修复】rec.seq() 返回 Cow<[u8]>，extend_from_slice 需要 &[u8]
        // 加上 & 符号借用 Cow，利用 Deref 特性自动转换为 &[u8]
        self.seq.extend_from_slice(&rec.seq());

        self.qual.clear();
        if let Some(q) = rec.qual() {
            self.qual.extend_from_slice(q);
        }
    }
}

struct WriteTask {
    file_index: usize, 
    hash: Hash,
    id_str: String,
}

#[derive(Args, Debug, Clone)]
pub struct ExtractArgs {
    #[arg(long = "genome-list")]
    pub genome_list: Option<PathBuf>,
    #[arg(short = 'i', long = "input", num_args = 1..=2)]
    pub input: Vec<PathBuf>,
    #[arg(short = 't', long = "type")]
    pub input_type: u8,
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,
    #[arg(long = "od")]
    pub output_dir: PathBuf,
    #[arg(long = "op", num_args = 1..=5)]
    pub output_prefix: Vec<String>,
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
    #[arg(long = "qc", default_value = "yes")]
    pub quality_control: String,
    #[arg(short = 'n', long, default_value = "0.08")]
    pub max_n: f64,
    #[arg(short = 'q', long, default_value = "30")]
    pub min_quality: u8,
    #[arg(short = 'p', long, default_value = "80")]
    pub min_quality_percent: u8,
    #[arg(short = 'b', long, default_value = "33")]
    pub quality_base: u8,
    #[arg(long = "pe", default_value = "pear")]
    pub pear_bin: String,
    #[arg(long = "pc", default_value = "1")]
    pub pear_threads: usize,
    
    // 【新增】是否使用 PEAR
    #[arg(long = "use-pear", default_value = "no")]
    pub use_pear: String,
}

pub fn run(args: ExtractArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global();

    if let Some(genome_list) = args.genome_list.clone() {
        return run_batch_mode(args, &genome_list);
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

    io_utils::ensure_directory(&args.output_dir)?;

    match input_type {
        InputType::ReferenceGenome => extract_reference_genome(&args, enzyme)?,
        InputType::ShotgunMetagenome => extract_shotgun(&args, enzyme, &qc)?,
        InputType::Single2bRAD => extract_single_tag(&args, enzyme, &qc)?,
        InputType::Concatenated2bRAD => extract_concatenated_tags(&args, enzyme, &qc)?,
    }
    Ok(())
}

fn run_batch_mode(base_args: ExtractArgs, genome_list: &std::path::Path) -> Result<()> {
    use std::io::BufRead;
    tracing::info!("### 批量处理模式：{}", genome_list.display());
    let file = File::open(genome_list)?;
    let reader = std::io::BufReader::new(file);
    let mut samples = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 { continue; }
        let sample_name = fields[0].to_string();
        let input1 = PathBuf::from(fields[1]);
        let input2 = if fields.len() > 2 && !fields[2].is_empty() { Some(PathBuf::from(fields[2])) } else { None };
        samples.push((sample_name, input1, input2));
    }
    
    let pb = ProgressBar::new(samples.len() as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner} {pos}/{len}").unwrap());

    samples.into_par_iter().for_each(|(sample_name, input1, input2)| {
        let mut sample_args = base_args.clone();
        sample_args.genome_list = None;
        sample_args.input = vec![input1];
        if let Some(in2) = input2 { sample_args.input.push(in2); }
        if base_args.input_type == 4 { 
             sample_args.output_prefix = (0..5).map(|i| format!("{}_tag{}", sample_name, i + 1)).collect();
        } else {
            sample_args.output_prefix = vec![sample_name.clone()];
        }
        
        match run_single_sample(sample_args) {
            Ok(_) => {},
            Err(e) => tracing::error!("样品 {} 处理失败: {}", sample_name, e),
        }
        pb.inc(1);
    });
    pb.finish();
    Ok(())
}

fn parse_enzyme(site: &str) -> Result<&'static Enzyme> {
     if let Some(enzyme) = enzyme_by_name(site) { return Ok(enzyme); }
     if let Ok(id) = site.parse::<u8>() { if let Some(enzyme) = enzyme_by_id(id) { return Ok(enzyme); } }
     bail!("未知的酶：{}", site)
}

// ==========================================
// 【核心优化】通用流水线读取器
// ==========================================
// 启动一个后台线程读取文件，并利用 recycle_rx 接收用过的 Batch 进行复用
type BatchData = (Vec<RawRecord>, usize); // (Buffer容器, 有效数据量)

fn spawn_reader_thread(input_path: PathBuf) -> (
    mpsc::Receiver<Result<BatchData>>, 
    mpsc::Sender<BatchData>, 
    thread::JoinHandle<()>
) {
    // work_tx: 发送填充好的数据给消费者
    let (work_tx, work_rx) = mpsc::sync_channel::<Result<BatchData>>(CHANNEL_BUFFER);
    // recycle_tx: 消费者把用完的容器还给生产者
    let (recycle_tx, recycle_rx) = mpsc::channel::<BatchData>();

    let handle = thread::spawn(move || {
        let mut reader = match parse_fastx_file(&input_path) {
            Ok(r) => r,
            Err(e) => {
                let _ = work_tx.send(Err(anyhow!(e).context("Failed to open input file")));
                return;
            }
        };

        loop {
            // 1. 获取一个 Batch 容器（优先从回收站拿，没有则新建）
            let (mut batch, _) = recycle_rx.try_recv().unwrap_or_else(|_| {
                let mut v = Vec::with_capacity(BATCH_SIZE);
                for _ in 0..BATCH_SIZE { v.push(RawRecord::new()); }
                (v, 0)
            });

            // 2. 填充数据
            let mut count = 0;
            let mut exhausted = false;

            for i in 0..BATCH_SIZE {
                match reader.next() {
                    Some(Ok(rec)) => {
                        // 【复用】这里不会分配新内存，而是复用 batch[i] 里的 Vec
                        batch[i].populate_from(&rec);
                        count += 1;
                    },
                    Some(Err(e)) => {
                        let _ = work_tx.send(Err(anyhow!(e).context("Fastx parse error")));
                        return; 
                    },
                    None => {
                        exhausted = true;
                        break;
                    }
                }
            }

            // 3. 发送数据
            if count > 0 {
                if work_tx.send(Ok((batch, count))).is_err() {
                    break; // 消费者断开
                }
            } else {
                // 如果这次没读到数据且已耗尽，就不发送了
                break;
            }

            if exhausted {
                break;
            }
        }
    });

    (work_rx, recycle_tx, handle)
}

// ========== Type 1: 参考基因组 ==========

fn extract_reference_genome(args: &ExtractArgs, enzyme: &'static Enzyme) -> Result<()> {
    let input_path = args.input[0].clone();
    let output_path = args.output_dir.join(format!("{}.{}.iibdb", args.output_prefix[0], enzyme.name));

    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path).context("Failed to create output file")?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(batch) = write_rx.recv() {
            for task in batch {
                io_utils::write_binary_record(&mut writer, task.hash, &task.id_str)?;
            }
        }
        Ok(())
    });

    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);
    
    let input_sequences = Arc::new(AtomicUsize::new(0));
    let total_tags = Arc::new(AtomicUsize::new(0));

    // 主线程消费循环
    while let Ok(result) = work_rx.recv() {
        let (batch, count) = result?;
        
        // 这里的 &batch[..count] 是 slice，只处理有效数据
        process_genome_batch(&batch[..count], enzyme, &write_tx, &input_sequences, &total_tags)?;

        // 【回收】处理完后，把整个容器还给 Reader
        let _ = recycle_tx.send((batch, 0)); 
    }

    drop(write_tx); 
    let _ = reader_handle.join();
    writer_handle.join().unwrap()?; 
    
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", args.output_prefix[0], enzyme.name));
    let stats = DigestStats {
        sample_id: args.output_prefix[0].clone(),
        enzyme: enzyme.name.to_string(),
        input_sequences: input_sequences.load(Ordering::Relaxed),
        tag_count: total_tags.load(Ordering::Relaxed),
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;

    Ok(())
}

fn process_genome_batch(
    batch: &[RawRecord], // 改为 slice
    enzyme: &Enzyme,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    let results: Vec<WriteTask> = batch.par_iter().flat_map(|record| {
        // record.seq 已经是 Vec<u8>，需要转大写。
        let mut sequence = record.seq.clone(); 
        sequence.make_ascii_uppercase();
        
        let positions_iter = enzyme.find_all_tags(&sequence);
        let mut tasks = Vec::new();
        // ID 处理：RawRecord.id 是 Vec<u8>，这里只在生成 string 时转换，将 UTF8 check 移到了并行线程中
        let id_utf8 = String::from_utf8_lossy(&record.id);
        // Fastx ID 通常包含空格，取第一部分
        let seq_id = id_utf8.split_whitespace().next().unwrap_or("seq");

        for (pos, len) in positions_iter {
            if pos + len > sequence.len() { continue; }
            let tag_seq = &sequence[pos..pos + len];
            let canonical = get_canonical_sequence(tag_seq);
            let hash = hash_bytes(&canonical);
            let id_str = format!("{}:{}", seq_id, pos);
            tasks.push(WriteTask { file_index: 0, hash, id_str });
        }
        tasks
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() {
        tx.send(results).context("Failed to send results to writer")?;
    }
    Ok(())
}

// ========== Type 2: Shotgun ==========

fn extract_shotgun(args: &ExtractArgs, enzyme: &'static Enzyme, qc: &QualityControl) -> Result<()> {
    // 【修改】核心逻辑：根据 use_pear 参数决定是否合并，还是获取所有待处理文件列表
    let inputs_to_process = if args.input.len() == 2 && args.use_pear.eq_ignore_ascii_case("yes") {
        tracing::info!("执行 PEAR 拼接 (use-pear=yes) ...");
        // 合并后产生一个新的文件，只处理这个文件
        vec![run_pear_and_combine(args, enzyme)?]
    } else {
        if args.input.len() == 2 {
            tracing::info!("跳过 PEAR 拼接 (use-pear=no)，依次处理双端文件 ...");
        }
        // 直接处理原始输入文件（可能是1个或2个）
        args.input.clone()
    };
    
    let output_path = args.output_dir.join(format!("{}.{}.iibsp", args.output_prefix[0], enzyme.name));

    // 创建写线程：所有输入文件提取的 tag 都写入这同一个文件
    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path)?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(batch) = write_rx.recv() {
            for task in batch { io_utils::write_binary_record(&mut writer, task.hash, &task.id_str)?; }
        }
        Ok(())
    });

    // 统计数据共享：跨文件累加
    let input_sequences = Arc::new(AtomicUsize::new(0));
    let tag_count = Arc::new(AtomicUsize::new(0));

    // 遍历所有待处理文件（如果是 PEAR 模式则只有1个，否则可能有2个）
    for input_path in inputs_to_process {
        tracing::info!("正在提取文件: {}", input_path.display());
        let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);

        // 处理单个文件的流水线循环
        while let Ok(result) = work_rx.recv() {
            let (batch, count) = result?;
            process_shotgun_batch(&batch[..count], enzyme, qc, &write_tx, &input_sequences, &tag_count)?;
            let _ = recycle_tx.send((batch, 0));
        }
        
        // 等待当前文件的读取线程结束，再处理下一个
        let _ = reader_handle.join();
    }

    // 所有文件处理完毕，关闭写入通道并等待写线程
    drop(write_tx);
    writer_handle.join().unwrap()?;
    
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", args.output_prefix[0], enzyme.name));
    let stats = DigestStats {
        sample_id: args.output_prefix[0].clone(),
        enzyme: enzyme.name.to_string(),
        input_sequences: input_sequences.load(Ordering::Relaxed),
        tag_count: tag_count.load(Ordering::Relaxed),
    };
    io_utils::write_sample_stats(&stat_path, &stats)?;
    
    Ok(())
}

fn process_shotgun_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    count_seq: &AtomicUsize,
    count_tag: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);

    let results: Vec<WriteTask> = batch.par_iter().flat_map(|record| {
        // QC 检查
        if !qc.check_n(&record.seq) { return Vec::new(); }
        if !record.qual.is_empty() {
            if !qc.check_quality(&record.qual) { return Vec::new(); }
        }

        let mut sequence = record.seq.clone(); // Clone for modification (uppercase)
        sequence.make_ascii_uppercase();
        
        let positions = enzyme.find_all_tags(&sequence);
        if positions.is_empty() { return Vec::new(); }

        let id_utf8 = String::from_utf8_lossy(&record.id);
        // 通常 Fastq ID 第一个空格前是 ID
        let seq_id = id_utf8.split_whitespace().next().unwrap_or(&id_utf8);

        let mut tasks = Vec::with_capacity(positions.len());
        for (i, (pos, len)) in positions.iter().enumerate() {
            let tag_seq = &sequence[*pos..*pos + len];
            let canonical = get_canonical_sequence(tag_seq);
            let tag_hash = hash_bytes(&canonical);
            let id_str = format!("{}_tag{}", seq_id, i + 1);
            tasks.push(WriteTask { file_index: 0, hash: tag_hash, id_str });
        }
        tasks
    }).collect();

    count_tag.fetch_add(results.len(), Ordering::Relaxed);
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

fn run_pear_and_combine(args: &ExtractArgs, enzyme: &Enzyme) -> Result<PathBuf> {
    let r1 = &args.input[0];
    let r2 = &args.input[1];
    let prefix = &args.output_prefix[0];
    let base = args.output_dir.join(format!("{}.{}", prefix, enzyme.name));
    let status = Command::new(&args.pear_bin)
        .args(["-f", r1.to_str().unwrap(), "-r", r2.to_str().unwrap(), "-e", "-o", base.to_str().unwrap(), "-j", &args.pear_threads.to_string()])
        .status()?;
    if !status.success() { bail!("PEAR failed"); }
    let pear_fastq = args.output_dir.join(format!("{}.{}.pear.fastq", prefix, enzyme.name));
    {
        let mut out = File::create(&pear_fastq)?;
        for suffix in [".assembled.fastq", ".unassembled.forward.fastq", ".unassembled.reverse.fastq"] {
             let p = args.output_dir.join(format!("{}.{}{}", prefix, enzyme.name, suffix));
             if p.exists() { std::io::copy(&mut File::open(&p)?, &mut out)?; std::fs::remove_file(p)?; }
        }
    }
    let discarded = args.output_dir.join(format!("{}.{}.discarded.fastq", prefix, enzyme.name));
    if discarded.exists() { std::fs::remove_file(discarded)?; }
    Ok(pear_fastq)
}

// ========== Type 3: 单标签 ==========
fn extract_single_tag(args: &ExtractArgs, enzyme: &'static Enzyme, qc: &QualityControl) -> Result<()> {
    let input_path = args.input[0].clone();
    let output_path = args.output_dir.join(format!("{}.{}.iibsp", args.output_prefix[0], enzyme.name));
    
    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&output_path)?;
        let mut writer = BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, file);
        while let Ok(batch) = write_rx.recv() { for task in batch { io_utils::write_binary_record(&mut writer, task.hash, &task.id_str)?; } }
        Ok(())
    });

    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(input_path);

    let input_sequences = Arc::new(AtomicUsize::new(0));
    let enzyme_reads = Arc::new(AtomicUsize::new(0));
    let qc_passed = Arc::new(AtomicUsize::new(0));

    while let Ok(result) = work_rx.recv() {
        let (batch, count) = result?;
        process_single_tag_batch(&batch[..count], enzyme, qc, &write_tx, &input_sequences, &enzyme_reads, &qc_passed)?;
        let _ = recycle_tx.send((batch, 0));
    }

    drop(write_tx);
    let _ = reader_handle.join();
    writer_handle.join().unwrap()?;
    
    let stat_path = args.output_dir.join(format!("{}.{}.stat.tsv", args.output_prefix[0], enzyme.name));
    let seqs = input_sequences.load(Ordering::Relaxed);
    let passed = qc_passed.load(Ordering::Relaxed);
    let mut stat_file = File::create(&stat_path)?;
    writeln!(stat_file, "sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent")?;
    let percent = if seqs > 0 { (passed as f64 / seqs as f64) * 100.0 } else { 0.0 };
    writeln!(stat_file, "{}\t{}\t{}\t{}\t{}\t{:.2}%", 
        args.output_prefix[0], enzyme.name, seqs, enzyme_reads.load(Ordering::Relaxed), passed, percent)?;

    Ok(())
}

fn process_single_tag_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    count_seq: &AtomicUsize,
    count_enz: &AtomicUsize,
    count_qc: &AtomicUsize,
) -> Result<()> {
    count_seq.fetch_add(batch.len(), Ordering::Relaxed);
    let results: Vec<WriteTask> = batch.par_iter().filter_map(|record| {
        // 为了安全起见和逻辑一致，这里需要 clone seq。
        // 但为了性能，如果逻辑允许，可以只引用。
        // 原逻辑中有一个 truncate(50)，这会修改数据。
        let mut sequence = record.seq.clone();
        if sequence.len() > 50 { sequence.truncate(50); }
        sequence.make_ascii_uppercase();

        for pattern in enzyme.patterns {
            if sequence.len() < enzyme.tag_length { break; }
            for offset in 0..sequence.len() {
                if offset + enzyme.tag_length > sequence.len() { break; }
                let window = &sequence[offset..offset + enzyme.tag_length];
                if pattern.matches(window) {
                    let mut pass = true;
                    if qc.check_n(window) { pass = false; }
                    if pass {
                        if !record.qual.is_empty() {
                            if offset + enzyme.tag_length <= record.qual.len() {
                                let qs = &record.qual[offset..offset + enzyme.tag_length];
                                if !qs.is_empty() && !qc.check_quality(qs) { pass = false; }
                            }
                        }
                    }
                    if pass {
                        let canonical = get_canonical_sequence(window);
                        let hash = hash_bytes(&canonical);
                        // 只有在确定通过时才转换 ID string
                        let id_str = String::from_utf8_lossy(&record.id).to_string();
                        return Some((true, WriteTask { file_index: 0, hash, id_str }));
                    } else {
                        // 酶切位点匹配但QC失败
                        return Some((false, WriteTask { file_index: 0, hash: 0, id_str: String::new() })); 
                    }
                }
            }
        }
        None 
    }).map(|(passed, task)| {
        if passed { (1, 1, Some(task)) } else { (1, 0, None) }
    }).collect::<Vec<_>>().into_iter().fold(Vec::new(), |mut acc, (enz, qc_pass, task_opt)| {
        count_enz.fetch_add(enz, Ordering::Relaxed);
        count_qc.fetch_add(qc_pass, Ordering::Relaxed);
        if let Some(t) = task_opt { acc.push(t); }
        acc
    });
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}

// ========== Type 4: 5连标签 ==========
fn extract_concatenated_tags(args: &ExtractArgs, enzyme: &'static Enzyme, qc: &QualityControl) -> Result<()> {
    if args.input.len() != 2 { bail!("Type 4 needs R1 and R2"); }
    let r1_path = args.input[0].clone();
    let output_paths: Vec<PathBuf> = (0..5).map(|i| args.output_dir.join(format!("{}.{}.iibsp", args.output_prefix[i], enzyme.name))).collect();

    let (write_tx, write_rx) = mpsc::sync_channel::<Vec<WriteTask>>(CHANNEL_BUFFER);
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut writers = Vec::new();
        for path in output_paths {
            let f = File::create(&path)?;
            writers.push(BufWriter::with_capacity(io_utils::IO_BUFFER_SIZE, f));
        }
        while let Ok(batch) = write_rx.recv() {
            for task in batch {
                if task.file_index < 5 { io_utils::write_binary_record(&mut writers[task.file_index], task.hash, &task.id_str)?; }
            }
        }
        Ok(())
    });

    let (work_rx, recycle_tx, reader_handle) = spawn_reader_thread(r1_path);
    
    let comb_reads = Arc::new(AtomicUsize::new(0));
    let enz_reads: Vec<Arc<AtomicUsize>> = (0..5).map(|_| Arc::new(AtomicUsize::new(0))).collect();
    let qc_passed: Vec<Arc<AtomicUsize>> = (0..5).map(|_| Arc::new(AtomicUsize::new(0))).collect();

    while let Ok(result) = work_rx.recv() {
        let (batch, count) = result?;
        process_concat_batch(&batch[..count], enzyme, qc, &write_tx, &comb_reads, &enz_reads, &qc_passed)?;
        let _ = recycle_tx.send((batch, 0));
    }
    
    drop(write_tx);
    let _ = reader_handle.join();
    writer_handle.join().unwrap()?;
    // Stats ... (省略)
    Ok(())
}

fn process_concat_batch(
    batch: &[RawRecord],
    enzyme: &Enzyme,
    qc: &QualityControl,
    tx: &mpsc::SyncSender<Vec<WriteTask>>,
    comb_reads: &AtomicUsize,
    enz_counts: &[Arc<AtomicUsize>],
    qc_counts: &[Arc<AtomicUsize>],
) -> Result<()> {
    comb_reads.fetch_add(batch.len(), Ordering::Relaxed);
    let results: Vec<WriteTask> = batch.par_iter().flat_map(|record| {
        let mut tasks = Vec::new();
        let mut work_seq = record.seq.clone();
        work_seq.make_ascii_uppercase();

        let id_utf8 = String::from_utf8_lossy(&record.id);

        for tag_idx in 0..5 {
            let start = enzyme.concat_starts[tag_idx];
            let end = enzyme.concat_ends[tag_idx];
            if end > work_seq.len() { continue; }
            let tag_seq = &work_seq[start..=end]; 
            
            let mut matched = false;
            for pattern in enzyme.patterns {
                if tag_seq.len() < enzyme.tag_length { break; }
                for offset in 0..=tag_seq.len().saturating_sub(enzyme.tag_length) {
                    let window = &tag_seq[offset..offset + enzyme.tag_length];
                    if pattern.matches(window) {
                        enz_counts[tag_idx].fetch_add(1, Ordering::Relaxed);
                        let mut pass = true;
                        if qc.check_n(window) { pass = false; }
                        if pass {
                             if !record.qual.is_empty() {
                                let q_start = start + offset;
                                let q_end = q_start + enzyme.tag_length;
                                if q_end <= record.qual.len() { if !qc.check_quality(&record.qual[q_start..q_end]) { pass = false; } }
                             }
                        }
                        if pass {
                            qc_counts[tag_idx].fetch_add(1, Ordering::Relaxed);
                            let canonical = get_canonical_sequence(window);
                            let hash = hash_bytes(&canonical);
                            tasks.push(WriteTask { file_index: tag_idx, hash, id_str: format!("{}:{}", id_utf8, tag_idx + 1) });
                        }
                        matched = true;
                        break;
                    }
                }
                if matched { break; }
            }
        }
        tasks
    }).collect();
    if !results.is_empty() { tx.send(results)?; }
    Ok(())
}