use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::types::{DigestStats, SampleRecord};

pub fn read_sample_list(path: &Path) -> Result<Vec<SampleRecord>> {
    let file =
        File::open(path).with_context(|| format!("无法读取样本列表文件：{}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("读取样本列表第 {} 行失败", line_no + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 2 {
            bail!("样本列表第 {} 行缺少输入文件路径：{}", line_no + 1, trimmed);
        }
        let sample_id = parts[0].to_string();
        let inputs: Vec<PathBuf> = parts[1..].iter().map(|p| PathBuf::from(p)).collect();
        records.push(SampleRecord { sample_id, inputs });
    }

    if records.is_empty() {
        bail!("样本列表为空：{}", path.display());
    }

    Ok(records)
}

pub fn ensure_directory(path: &Path) -> Result<()> {
    fs::create_dir_all(path).with_context(|| format!("创建输出目录失败：{}", path.display()))
}

pub fn write_sample_stats(path: &Path, stats: &DigestStats) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("写入样本统计失败：{}", path.display()))?;
    writeln!(file, "sample\tenzyme\tinput_sequences\ttag_count\tpercent")?;
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{:.2}",
        stats.sample_id,
        stats.enzyme,
        stats.input_sequences,
        stats.tag_count,
        stats.percent()
    )?;
    Ok(())
}

pub fn write_summary(path: &Path, stats: &[DigestStats]) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("写入汇总统计失败：{}", path.display()))?;
    writeln!(file, "sample\tenzyme\tinput_sequences\ttag_count\tpercent")?;
    for stat in stats {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{:.2}",
            stat.sample_id,
            stat.enzyme,
            stat.input_sequences,
            stat.tag_count,
            stat.percent()
        )?;
    }
    Ok(())
}
