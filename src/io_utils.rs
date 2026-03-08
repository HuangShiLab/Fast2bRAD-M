use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use crate::types::DigestStats;

// [Optimization] Increase I/O buffer size (default 8 KB -> 128 KB)
// Reduces the number of system calls and significantly improves throughput for large file I/O
pub const IO_BUFFER_SIZE: usize = 128 * 1024;

pub fn ensure_directory(path: &Path) -> Result<()> {
    fs::create_dir_all(path).with_context(|| format!("Failed to create output directory: {}", path.display()))
}

pub fn write_sample_stats(path: &Path, stats: &DigestStats) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to write sample statistics: {}", path.display()))?;
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

// ================== Binary format read/write utilities ==================

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

    /// Read the next record (legacy API, kept for compatibility; prefer next_record_reuse)
    pub fn next_record(&mut self) -> Result<Option<(u64, String)>> {
        let mut buffer = String::new();
        if let Some(hash) = self.next_record_reuse(&mut buffer)? {
            Ok(Some((hash, buffer)))
        } else {
            Ok(None)
        }
    }

    /// [Optimization] Read the next record, reusing the provided String buffer
    /// This avoids millions of String allocations
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
        self.reader.read_exact(&mut len_buf).context("Failed to read ID length")?;
        let len = u16::from_le_bytes(len_buf) as usize;

        // Reuse the buffer: clear content but retain capacity
        buffer.clear();

        // Use unsafe to write directly into the vec for maximum performance
        // Safe as long as the data source is valid UTF-8 (FASTQ IDs are typically ASCII)
        unsafe {
            let v = buffer.as_mut_vec();
            if v.capacity() < len {
                v.reserve(len - v.len());
            }
            v.set_len(len);
            self.reader.read_exact(v).context("Failed to read ID content")?;
        }

        Ok(Some(hash))
    }
}

pub fn open_binary_reader<P: AsRef<Path>>(
    path: P,
) -> Result<BinaryRecordReader<Box<dyn Read + Send>>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Cannot open file: {}", path.display()))?;

    // [Optimization] Apply large buffer
    let reader: Box<dyn Read + Send> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::with_capacity(IO_BUFFER_SIZE, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(IO_BUFFER_SIZE, file))
    };

    Ok(BinaryRecordReader::new(reader))
}