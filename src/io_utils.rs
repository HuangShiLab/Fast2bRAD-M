use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use crate::types::DigestStats;

// 【优化】增大 I/O 缓冲区 (默认 8KB -> 128KB)
// 减少系统调用次数，显著提升大文件读写吞吐量
pub const IO_BUFFER_SIZE: usize = 128 * 1024;

pub fn ensure_directory(path: &Path) -> Result<()> {
    fs::create_dir_all(path).with_context(|| format!("创建输出目录失败：{}", path.display()))
}

pub fn write_sample_stats(path: &Path, stats: &DigestStats) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("写入样本统计失败：{}", path.display()))?;
    writeln!(file, "sample\tenzyme\tinput_sequences\ttag_count\tpercent")?;
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{:.2}%",
        stats.sample_id,
        stats.enzyme,
        stats.input_sequences,
        stats.tag_count,
        stats.percent()
    )?;
    Ok(())
}

// ================== 二进制格式读写工具 ==================

pub fn write_binary_record<W: Write>(writer: &mut W, hash: u64, id: &str) -> io::Result<()> {
    writer.write_all(&hash.to_le_bytes())?;
    let id_bytes = id.as_bytes();
    let id_len = id_bytes.len().min(u16::MAX as usize) as u16;
    writer.write_all(&id_len.to_le_bytes())?;
    writer.write_all(&id_bytes[..id_len as usize])?;
    Ok(())
}

pub struct BinaryRecordReader<R> {
    reader: R,
}

impl<R: Read> BinaryRecordReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    /// 读取下一条记录（旧 API，保留以防万一，但建议用 next_record_reuse）
    pub fn next_record(&mut self) -> Result<Option<(u64, String)>> {
        let mut buffer = String::new();
        if let Some(hash) = self.next_record_reuse(&mut buffer)? {
            Ok(Some((hash, buffer)))
        } else {
            Ok(None)
        }
    }

    /// 【优化】读取下一条记录，复用传入的 String buffer
    /// 这样可以避免数百万次的 String 分配
    pub fn next_record_reuse(&mut self, buffer: &mut String) -> Result<Option<u64>> {
        let mut hash_buf = [0u8; 8];
        if let Err(e) = self.reader.read_exact(&mut hash_buf) {
            if e.kind() == io::ErrorKind::UnexpectedEof {
                return Ok(None);
            }
            return Err(e.into());
        }
        let hash = u64::from_le_bytes(hash_buf);

        let mut len_buf = [0u8; 2];
        self.reader.read_exact(&mut len_buf).context("读取ID长度失败")?;
        let len = u16::from_le_bytes(len_buf) as usize;

        // 复用缓冲区：清空内容但保留容量
        buffer.clear();
        
        // 使用 unsafe 直接写入 vec 以获得最高性能
        // 前提是确信数据源是合法的 UTF-8（通常 fastq id 都是 ascii）
        unsafe {
            let v = buffer.as_mut_vec();
            if v.capacity() < len {
                v.reserve(len - v.len());
            }
            v.set_len(len);
            self.reader.read_exact(v).context("读取ID内容失败")?;
        }

        Ok(Some(hash))
    }
}

pub fn open_binary_reader<P: AsRef<Path>>(
    path: P,
) -> Result<BinaryRecordReader<Box<dyn Read + Send>>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("无法打开文件: {}", path.display()))?;

    // 【优化】应用大缓冲区
    let reader: Box<dyn Read + Send> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::with_capacity(IO_BUFFER_SIZE, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(IO_BUFFER_SIZE, file))
    };

    Ok(BinaryRecordReader::new(reader))
}