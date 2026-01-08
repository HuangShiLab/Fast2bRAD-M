use anyhow::{bail, Context, Result, anyhow};
use clap::Parser;
use fxhash::FxHashMap;
use rayon::prelude::*; 
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing;

use crate::enzymes::{enzyme_by_id, enzyme_by_name};
use crate::io_utils;

#[derive(Parser, Debug)]
pub struct QuantifyArgs {
    #[arg(short = 'l', long = "list")]
    pub sample_list: PathBuf,
    #[arg(short = 'd', long = "database")]
    pub database_dir: PathBuf,
    #[arg(short = 't', long = "taxonomy")]
    pub taxonomy_level: String,
    #[arg(short = 's', long = "site")]
    pub enzyme_site: String,
    #[arg(short = 'o', long = "output")]
    pub output_dir: PathBuf,
    #[arg(short = 'g', long = "gscore", default_value = "0")]
    pub g_score_threshold: f64,
    #[arg(short = 'v', long = "verbose", default_value = "yes")]
    pub verbose: String,
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

#[derive(Debug, Default)]
struct AbundanceStats {
    theoretical_tag_num: f64,    
    sequenced_tag_num: usize,    
    sequenced_reads_num: usize,  
    sequenced_tag_num_gt1: usize, 
}

impl AbundanceStats {
    fn percent(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 { (self.sequenced_tag_num as f64 / self.theoretical_tag_num) * 100.0 } else { 0.0 }
    }
    fn reads_per_theoretical(&self) -> f64 {
        if self.theoretical_tag_num > 0.0 { self.sequenced_reads_num as f64 / self.theoretical_tag_num } else { 0.0 }
    }
    fn reads_per_sequenced(&self) -> f64 {
        if self.sequenced_tag_num > 0 { self.sequenced_reads_num as f64 / self.sequenced_tag_num as f64 } else { 0.0 }
    }
    fn g_score(&self) -> f64 {
        ((self.sequenced_tag_num * self.sequenced_reads_num) as f64).sqrt()
    }
}

pub fn run(args: QuantifyArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();
    let verbose = args.verbose.to_lowercase() == "yes";
    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("无效的酶切位点编号"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("无效的酶名称"))?
    };
    let tax_level = validate_taxonomy_level(&args.taxonomy_level)?;

    tracing::info!("COMMAND: quantify -l {} -d {} -t {} -s {} -o {} -g {} -v {} -j {}", args.sample_list.display(), args.database_dir.display(), args.taxonomy_level, args.enzyme_site, args.output_dir.display(), args.g_score_threshold, args.verbose, args.threads);

    std::fs::create_dir_all(&args.output_dir)?;
    // 【修改】读取 .iibdb 后缀
    let db_file = args.database_dir.join(format!("{}.{}.iibdb", enzyme.name, tax_level));
    let classify_file = args.database_dir.join("abfh_classify_with_speciename.txt.gz");

    if !db_file.exists() { bail!("数据库文件不存在"); }
    if !classify_file.exists() { bail!("分类文件不存在"); }

    tracing::info!("### 加载数据库：{}", db_file.display());
    let (tag_to_gcfs, gcf_to_taxonomy, tag_theory_num) = load_database(&db_file, &classify_file, &args.taxonomy_level)?;
    tracing::info!("### 数据库加载完成");

    let samples = read_sample_list(&args.sample_list)?;
    tracing::info!("共 {} 个样品待处理", samples.len());

    samples.par_iter().for_each(|(sample_name, sample_data)| {
        tracing::info!(">>> ({}) 样品分析开始", sample_name);
        let result = process_sample(sample_name, sample_data, &tag_to_gcfs, &gcf_to_taxonomy, &tag_theory_num, enzyme, &args.output_dir, args.g_score_threshold, verbose);
        match result {
            Ok(_) => tracing::info!("<<< ({}) 样品分析完成", sample_name),
            Err(e) => tracing::error!("!!! ({}) 错误: {}", sample_name, e),
        }
    });
    tracing::info!("\n全部完成！");
    Ok(())
}

fn validate_taxonomy_level(level: &str) -> Result<String> {
    let valid = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"];
    if valid.contains(&level) { Ok(level.to_string()) } else { bail!("无效的分类层级"); }
}

fn read_sample_list(list_path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let content = std::fs::read_to_string(list_path)?;
    let mut samples = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 { continue; }
        let sample_path = PathBuf::from(parts[1]);
        if !sample_path.exists() { tracing::warn!("警告: 样品文件不存在: {}", sample_path.display()); continue; }
        samples.push((parts[0].to_string(), sample_path));
    }
    Ok(samples)
}

fn load_database(db_file: &Path, classify_file: &Path, tax_level: &str) -> Result<(FxHashMap<u64, Vec<String>>, FxHashMap<String, String>, FxHashMap<String, FxHashMap<String, FxHashMap<u64, usize>>>)> {
    let level_index = get_taxonomy_level_index(tax_level);
    let mut gcf_to_taxonomy: FxHashMap<String, String> = FxHashMap::default();
    
    let content = if classify_file.to_str().unwrap().ends_with(".gz") {
        use flate2::read::GzDecoder; use std::io::Read;
        let file = File::open(classify_file)?;
        let mut decoder = GzDecoder::new(file);
        let mut c = String::new(); decoder.read_to_string(&mut c)?; c
    } else { std::fs::read_to_string(classify_file)? };

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= level_index { continue; }
        let gcf_id = parts[0].to_string();
        let taxonomy = parts[1..=level_index].join("\t");
        gcf_to_taxonomy.insert(gcf_id, taxonomy);
    }

    let mut tag_to_gcfs: FxHashMap<u64, Vec<String>> = FxHashMap::default();
    let mut tag_theory_num: FxHashMap<String, FxHashMap<String, FxHashMap<u64, usize>>> = FxHashMap::default();
    let mut reader = io_utils::open_binary_reader(db_file)?;
    let mut id_buffer = String::with_capacity(256); // 复用

    while let Some(hash_val) = reader.next_record_reuse(&mut id_buffer)? {
        let mut parts = id_buffer.split('|');
        let gcf_id = parts.next().unwrap_or("");
        // 跳过 index, scaffold, pos, strand (4个)
        let unique_flag = parts.nth(4).unwrap_or("0");
        if unique_flag != "1" { continue; }
        if gcf_id.is_empty() { continue; }

        if let Some(taxonomy) = gcf_to_taxonomy.get(gcf_id) {
            tag_to_gcfs.entry(hash_val).or_insert_with(Vec::new).push(gcf_id.to_string());
            *tag_theory_num.entry(taxonomy.clone()).or_insert_with(FxHashMap::default)
                .entry(gcf_id.to_string()).or_insert_with(FxHashMap::default)
                .entry(hash_val).or_insert(0) += 1;
        }
    }
    Ok((tag_to_gcfs, gcf_to_taxonomy, tag_theory_num))
}

fn get_taxonomy_level_index(level: &str) -> usize {
    match level { "kingdom" => 1, "phylum" => 2, "class" => 3, "order" => 4, "family" => 5, "genus" => 6, "species" => 7, "strain" => 8, _ => 7 }
}

fn process_sample(
    sample_name: &str,
    sample_data: &Path,
    tag_to_gcfs: &FxHashMap<u64, Vec<String>>,
    gcf_to_taxonomy: &FxHashMap<String, String>,
    tag_theory_num: &FxHashMap<String, FxHashMap<String, FxHashMap<u64, usize>>>,
    enzyme: &crate::enzymes::Enzyme,
    output_dir: &Path,
    g_score_threshold: f64,
    verbose: bool,
) -> Result<()> {
    let mut tag_num: FxHashMap<String, FxHashMap<u64, usize>> = FxHashMap::default();
    let mut detected_gcf_tag: FxHashMap<String, FxHashMap<String, fxhash::FxHashSet<u64>>> = FxHashMap::default();
    let mut reader = io_utils::open_binary_reader(sample_data)?;
    let mut ignored_id_buf = String::with_capacity(128); // 复用

    while let Some(tag_hash) = reader.next_record_reuse(&mut ignored_id_buf)? {
        if let Some(gcf_list) = tag_to_gcfs.get(&tag_hash) {
            if let Some(first_gcf) = gcf_list.first() {
                if let Some(taxonomy) = gcf_to_taxonomy.get(first_gcf) {
                    *tag_num.entry(taxonomy.clone()).or_insert_with(FxHashMap::default).entry(tag_hash).or_insert(0) += 1;
                    for gcf_id in gcf_list {
                        detected_gcf_tag.entry(taxonomy.clone()).or_insert_with(FxHashMap::default)
                            .entry(gcf_id.clone()).or_insert_with(fxhash::FxHashSet::default).insert(tag_hash);
                    }
                }
            }
        }
    }

    if tag_num.is_empty() { tracing::warn!("!!! ({}) 警告: 未检测到标签", sample_name); return Ok(()); }

    let sample_dir = output_dir.join(sample_name);
    std::fs::create_dir_all(&sample_dir)?;
    let gcf_detected_file = sample_dir.join(format!("{}.{}.GCF_detected.xls", sample_name, enzyme.name));
    let mut gcf_writer = BufWriter::new(File::create(&gcf_detected_file)?);
    
    let mut taxonomy_list: Vec<&String> = detected_gcf_tag.keys().collect();
    taxonomy_list.sort();
    
    for taxonomy in taxonomy_list {
        let gcf_map = &detected_gcf_tag[taxonomy];
        let mut gcf_list: Vec<&String> = gcf_map.keys().collect();
        gcf_list.sort();
        for gcf_id in gcf_list {
            let detected_tags = &gcf_map[gcf_id];
            let detected_tag_num = detected_tags.len();
            let gcf_all_theory_num = if let Some(gcf_tags) = tag_theory_num.get(taxonomy) {
                if let Some(tags) = gcf_tags.get(gcf_id) { tags.len() } else { 0 }
            } else { 0 };
            let percent = if gcf_all_theory_num > 0 { detected_tag_num as f64 / gcf_all_theory_num as f64 } else { 0.0 };
            writeln!(gcf_writer, "{}\t{}\t{}\t{}\t{:.4}", taxonomy, gcf_id, gcf_all_theory_num, detected_tag_num, percent)?;
        }
    }

    let output_file = sample_dir.join(format!("{}.{}.xls", sample_name, enzyme.name));
    let mut writer = BufWriter::new(File::create(&output_file)?);
    let level_count = if let Some(first_tax) = gcf_to_taxonomy.values().next() { first_tax.split('\t').count() } else { 7 };
    let tax_levels = get_taxonomy_header_by_level(level_count);
    writeln!(writer, "#{}\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\tSequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\tSequenced_Reads_Num/Sequenced_Tag_Num\tSequenced_Tag_Num(depth>1)\tG_Score", tax_levels)?;

    for (taxonomy, tags) in &tag_num {
        let mut stats = AbundanceStats::default();
        if let Some(gcf_tags) = tag_theory_num.get(taxonomy) {
            let total_theory: usize = gcf_tags.values().map(|tags| tags.values().sum::<usize>()).sum();
            stats.theoretical_tag_num = total_theory as f64 / gcf_tags.len() as f64;
        }
        stats.sequenced_tag_num = tags.len();
        stats.sequenced_reads_num = tags.values().sum();
        stats.sequenced_tag_num_gt1 = tags.values().filter(|&&count| count > 1).count();
        let g_score = stats.g_score();
        if g_score < g_score_threshold { continue; }

        writeln!(writer, "{}\t{:.8}\t{}\t{:.8}%\t{}\t{:.8}\t{:.8}\t{}\t{:.8}", taxonomy, stats.theoretical_tag_num, stats.sequenced_tag_num, stats.percent(), stats.sequenced_reads_num, stats.reads_per_theoretical(), stats.reads_per_sequenced(), stats.sequenced_tag_num_gt1, g_score)?;

        if verbose {
            let output_name = taxonomy.split('\t').last().unwrap_or("unknown");
            let detail_dir = sample_dir.join(format!("{}.{}", sample_name, enzyme.name));
            std::fs::create_dir_all(&detail_dir)?;
            let detail_file = detail_dir.join(format!("{}.xls", output_name));
            let mut detail_writer = BufWriter::new(File::create(detail_file)?);
            for (tag, &count) in tags { writeln!(detail_writer, "{}\t{}", tag, count)?; }
        }
    }
    Ok(())
}

fn get_taxonomy_header_by_level(level_count: usize) -> String {
    let headers = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"];
    headers[..level_count.min(8)].join("\t")
}