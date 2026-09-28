use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;


// [Optimization] Increase I/O buffer size (default 8 KB -> 128 KB)
// Reduces the number of system calls and significantly improves throughput for large file I/O
pub const IO_BUFFER_SIZE: usize = 128 * 1024;

pub fn ensure_directory(path: &Path) -> Result<()> {
    fs::create_dir_all(path).with_context(|| format!("Failed to create output directory: {}", path.display()))
}

/// Write the `abfh_classify_with_speciename.txt.gz` taxonomy mapping that
/// quantify/predict/find-genome consume. One gzip-compressed line per genome:
///   GCF_ID<TAB>kingdom<TAB>phylum<TAB>...<TAB>species<TAB>strain
/// Written straight from the parsed `GenomeRecord.taxonomy`, so it always stays
/// consistent with the per-level .iibdb databases (same strain synthesis, etc.).
pub fn write_classify_file<'a, I>(output_dir: &Path, genomes: I) -> Result<()>
where
    I: IntoIterator<Item = (&'a str, &'a [String])>,
{
    let path = output_dir.join("abfh_classify_with_speciename.txt.gz");
    let file = File::create(&path)
        .with_context(|| format!("Failed to create classify file: {}", path.display()))?;
    let buf_writer = BufWriter::with_capacity(IO_BUFFER_SIZE, file);
    let mut encoder = GzEncoder::new(buf_writer, Compression::default());
    for (gcf_id, taxonomy) in genomes {
        encoder.write_all(gcf_id.as_bytes())?;
        for rank in taxonomy {
            encoder.write_all(b"\t")?;
            encoder.write_all(rank.as_bytes())?;
        }
        encoder.write_all(b"\n")?;
    }
    encoder.finish()?;
    Ok(())
}

// ================== Taxonomy rank prefixes ==================

/// Rank prefixes for taxonomy levels 1–7 (kingdom → species). Level 8 (strain)
/// is intentionally left unprefixed (it is typically the accession, or
/// "<species> <accession>").
pub const RANK_PREFIXES: [&str; 7] = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"];

/// Strip a leading rank prefix such as `d__` / `s__` (a short letter code
/// followed by `__`). Anything else is returned unchanged.
fn strip_rank_prefix(s: &str) -> &str {
    match s.find("__") {
        Some(pos) if pos <= 2 => &s[pos + 2..],
        _ => s,
    }
}

/// Normalize a genome's taxonomy to the 2bRAD-M convention: ensure levels 1–7
/// carry `k__/p__/c__/o__/f__/g__/s__` (replacing any existing prefix, e.g.
/// GTDB's `d__`), leaving level 8 (strain) as-is. In place; tolerant of fewer
/// than 8 levels. Bare names get prefixes added; already-prefixed names are
/// normalized to this scheme.
pub fn apply_rank_prefixes(taxonomy: &mut [String]) {
    for (i, prefix) in RANK_PREFIXES.iter().enumerate() {
        if let Some(level) = taxonomy.get_mut(i) {
            let bare = strip_rank_prefix(level).to_string();
            *level = format!("{}{}", prefix, bare);
        }
    }
}

/// Build the normalized 8-rank taxonomy from raw rank strings. `raw_levels` may
/// come from a single semicolon-delimited GTDB column (`d__X;p__Y;..` split on
/// `;`) or from one rank per TSV column; each entry may be bare or already
/// prefixed. Ranks 1–7 are normalized to `k__/p__/.../s__`; rank 8 (strain) is
/// synthesized as "<species> <genome_id>" when only 7 ranks are supplied, so
/// each genome forms its own strain (otherwise every genome of a species would
/// collapse into one strain and the strain-level database would duplicate the
/// species one).
pub fn normalize_taxonomy(raw_levels: &[&str], genome_id: &str) -> Vec<String> {
    let mut tax: Vec<String> = raw_levels
        .iter()
        .map(|s| strip_rank_prefix(s.trim()).to_string())
        .collect();
    while tax.len() < 8 {
        if tax.len() == 7 {
            let species = tax.last().cloned().unwrap_or_else(|| "unknown".to_string());
            // Join species and genome id with '_' (and de-space the species name)
            // so the synthesized strain has no whitespace.
            tax.push(format!("{} {}", species, genome_id).replace(' ', "_"));
        } else if let Some(last) = tax.last().cloned() {
            tax.push(format!("{}_strain", last));
        } else {
            tax.push("unknown".to_string());
        }
    }
    tax.truncate(8);
    apply_rank_prefixes(&mut tax);
    tax
}

/// Heuristic: does this trailing column look like a genome file path rather than
/// a taxonomy rank? Lets the parser peel an optional path column off a
/// tab-separated taxonomy line.
pub fn looks_like_path(s: &str) -> bool {
    let s = s.trim();
    s.contains('/')
        || s.ends_with(".gz")
        || s.ends_with(".fa")
        || s.ends_with(".fna")
        || s.ends_with(".fasta")
        || s.ends_with(".ffn")
        || s.ends_with(".frn")
}

/// Per-genome digest statistics (`-t 1`). Deliberately does *not* carry a
/// "percent": tags/contigs is meaningless (it routinely exceeds 100%, e.g.
/// "5000.00%"), so tag density per Mb is reported instead.
pub fn write_genome_stats(
    path: &Path,
    sample_id: &str,
    enzyme: &str,
    contigs: usize,
    total_bases: usize,
    tag_count: usize,
) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to write sample statistics: {}", path.display()))?;
    writeln!(file, "sample\tenzyme\tcontigs\ttotal_bases\ttag_count\ttags_per_mb")?;
    let tags_per_mb = if total_bases == 0 {
        0.0
    } else {
        tag_count as f64 * 1_000_000.0 / total_bases as f64
    };
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{:.2}",
        sample_id, enzyme, contigs, total_bases, tag_count, tags_per_mb
    )?;
    Ok(())
}

/// Per-sample read statistics (`-t 2`). `percent` is tags per 100 input reads.
pub fn write_read_stats(
    path: &Path,
    sample_id: &str,
    enzyme: &str,
    input_reads: usize,
    tag_count: usize,
) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to write sample statistics: {}", path.display()))?;
    writeln!(file, "sample\tenzyme\tinput_reads_num\ttag_count\tpercent")?;
    let percent = if input_reads == 0 {
        0.0
    } else {
        tag_count as f64 / input_reads as f64 * 100.0
    };
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{:.2}%",
        sample_id, enzyme, input_reads, tag_count, percent
    )?;
    Ok(())
}

// ================== Binary format read/write utilities ==================
//
// Two record streams share this module:
//
//  * **id-carrying stream** (reference genome `*.iibdb`, per-genome digests):
//    `[8] hash | [2] id_len | [id_len] id`. The id is what makes a genome
//    database useful (contig / `gcf|idx|scaffold|pos`), so it is kept.
//
//  * **id-less stream** (sample tags, `*.iibsp`): an 8-byte header
//    (`IIBS` + version) followed by bare `[8] hash` records. Nothing
//    downstream of `extract` ever reads a sample record's id — `quantify`
//    counts hashes — so storing per-read names only burned disk and I/O
//    (~15 bytes/tag, i.e. hundreds of MB on a 10^7-tag sample) and invited
//    duplicate-id confusion between R1/R2. See `write_sample_tag_header`.
//
// `BinaryRecordReader` auto-detects which of the two it is looking at, so
// sample files written by older versions (with ids) still read back fine.

pub fn write_binary_record<W: Write>(writer: &mut W, hash: u64, id: &str) -> io::Result<()> {
    writer.write_all(&hash.to_le_bytes())?;
    let id_bytes = id.as_bytes();
    let id_len = id_bytes.len().min(u16::MAX as usize) as u16;
    writer.write_all(&id_len.to_le_bytes())?;
    writer.write_all(&id_bytes[..id_len as usize])?;
    Ok(())
}

/// Magic marking an id-less sample tag stream. Chosen so it cannot be mistaken
/// for a legacy first record: a legacy stream starts with a tag hash, and the
/// chance of a hash matching these exact 8 bytes is 2^-64.
pub const SAMPLE_TAG_MAGIC: &[u8; 4] = b"IIBS";
pub const SAMPLE_TAG_VERSION: u32 = 1;

/// Write the id-less sample tag stream header. Must be called once, before any
/// `write_sample_tag`.
pub fn write_sample_tag_header<W: Write>(writer: &mut W) -> io::Result<()> {
    writer.write_all(SAMPLE_TAG_MAGIC)?;
    writer.write_all(&SAMPLE_TAG_VERSION.to_le_bytes())
}

/// Write one bare tag hash into an id-less sample tag stream.
#[inline]
pub fn write_sample_tag<W: Write>(writer: &mut W, hash: u64) -> io::Result<()> {
    writer.write_all(&hash.to_le_bytes())
}

/// Fill `buf`, returning how many bytes were actually read. Returns `buf.len()`
/// on success and `0` at a clean end of stream; anything in between means the
/// stream was truncated mid-record. Unlike `read_exact` this tolerates the
/// short reads that decompressors (gzip/zstd) legitimately return.
fn read_full<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

pub struct BinaryRecordReader<R> {
    reader: R,
    /// True for the `IIBS` sample stream: records are bare hashes, no ids.
    id_less: bool,
    /// First 8 bytes of a legacy stream, consumed while sniffing the header.
    pending_hash: Option<[u8; 8]>,
    at_eof: bool,
    /// Scratch for the id bytes of an id-carrying record.
    id_bytes: Vec<u8>,
}

impl<R: Read> BinaryRecordReader<R> {
    /// Sniff the 8-byte header to decide whether this is an id-less sample
    /// stream or a legacy id-carrying one.
    pub fn new(mut reader: R) -> Result<Self> {
        let mut head = [0u8; 8];
        let n = read_full(&mut reader, &mut head).context("Failed to read record stream header")?;
        match n {
            0 => Ok(Self {
                reader,
                id_less: false,
                pending_hash: None,
                at_eof: true,
                id_bytes: Vec::new(),
            }),
            8 if &head[..4] == SAMPLE_TAG_MAGIC => {
                let version = u32::from_le_bytes([head[4], head[5], head[6], head[7]]);
                if version != SAMPLE_TAG_VERSION {
                    anyhow::bail!("Unsupported sample tag stream version: {}", version);
                }
                Ok(Self {
                    reader,
                    id_less: true,
                    pending_hash: None,
                    at_eof: false,
                    id_bytes: Vec::new(),
                })
            }
            8 => Ok(Self {
                reader,
                id_less: false,
                pending_hash: Some(head),
                at_eof: false,
                id_bytes: Vec::with_capacity(128),
            }),
            n => anyhow::bail!("Truncated record stream: only {} of 8 header bytes present", n),
        }
    }

    /// Whether records in this stream carry no id (sample tag stream).
    pub fn is_id_less(&self) -> bool {
        self.id_less
    }

    /// [Optimization] Read the next record, reusing the provided String buffer
    /// This avoids millions of String allocations. For an id-less stream the
    /// buffer is simply left empty.
    pub fn next_record_reuse(&mut self, buffer: &mut String) -> Result<Option<u64>> {
        buffer.clear();

        let hash = match self.pending_hash.take() {
            Some(bytes) => u64::from_le_bytes(bytes),
            None => {
                if self.at_eof {
                    return Ok(None);
                }
                let mut hash_buf = [0u8; 8];
                match read_full(&mut self.reader, &mut hash_buf).context("Failed to read tag hash")? {
                    0 => {
                        self.at_eof = true;
                        return Ok(None);
                    }
                    8 => u64::from_le_bytes(hash_buf),
                    n => anyhow::bail!("Truncated record: {} trailing byte(s) after the last complete record", n),
                }
            }
        };

        if self.id_less {
            return Ok(Some(hash));
        }

        let mut len_buf = [0u8; 2];
        if read_full(&mut self.reader, &mut len_buf).context("Failed to read ID length")? != 2 {
            anyhow::bail!("Truncated record: tag hash is not followed by an ID length");
        }
        let len = u16::from_le_bytes(len_buf) as usize;

        self.id_bytes.clear();
        self.id_bytes.resize(len, 0);
        if read_full(&mut self.reader, &mut self.id_bytes).context("Failed to read ID content")? != len {
            anyhow::bail!("Truncated record: ID is shorter than its declared length ({} bytes)", len);
        }
        // Borrowed (no allocation) for the ASCII ids these files actually
        // contain; only malformed input pays for the lossy conversion.
        buffer.push_str(&String::from_utf8_lossy(&self.id_bytes));

        Ok(Some(hash))
    }
}

/// Zstd frame magic bytes: 0xFD2FB528 (little-endian)
const ZSTD_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];

pub fn open_binary_reader<P: AsRef<Path>>(
    path: P,
) -> Result<BinaryRecordReader<Box<dyn Read + Send>>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Cannot open file: {}", path.display()))?;

    // Gzip: detect by extension
    if path.extension().map_or(false, |ext| ext == "gz") {
        let reader: Box<dyn Read + Send> = Box::new(BufReader::with_capacity(IO_BUFFER_SIZE, GzDecoder::new(file)));
        return BinaryRecordReader::new(reader);
    }

    // Peek first 4 bytes to auto-detect zstd format
    let mut buf_reader = BufReader::with_capacity(IO_BUFFER_SIZE, file);
    let is_zstd = {
        let buf = buf_reader.fill_buf().context("Failed to peek file header")?;
        buf.len() >= 4 && buf[..4] == ZSTD_MAGIC
    };

    let reader: Box<dyn Read + Send> = if is_zstd {
        Box::new(zstd::Decoder::new(buf_reader).context("Failed to create zstd decoder")?)
    } else {
        Box::new(buf_reader)
    };

    BinaryRecordReader::new(reader)
}

// ================== Compact database format ==================
// Optimized for level-specific databases (e.g., BcgI.species.iibdb).
// Stores only hash + GCF index per record (12 bytes vs ~70 bytes in legacy format).
//
// Header (always uncompressed):
//   [4 bytes] magic: b"IIBC"
//   [4 bytes] version: u32 LE
//     - v1: uncompressed records, no record_count
//     - v2: zstd-compressed records, no record_count
//     - v3: zstd-compressed records, record_count present (below)
//   [8 bytes] record_count: u64 LE      (v3 only; total number of records)
//   [4 bytes] gcf_count: u32 LE
//   For each GCF (gcf_count times):
//     [2 bytes] id_len: u16 LE
//     [N bytes] id_bytes (UTF-8)
// Records (repeated until EOF, zstd-compressed in v2/v3):
//   [8 bytes] tag_hash: u64 LE
//   [4 bytes] gcf_index: u32 LE (index into GCF table)

pub const COMPACT_MAGIC: &[u8; 4] = b"IIBC";
/// Current write version: v3 = zstd records + record_count in header.
pub const COMPACT_VERSION: u32 = 3;
/// Byte offset of the v3 `record_count` field (after magic[4] + version[4]).
const RECORD_COUNT_OFFSET: u64 = 8;

// ---- Writer ----

/// Writes compact database files with zstd-compressed records section.
/// Header (magic + version + record_count + GCF table) is always uncompressed.
///
/// The record count is not known until every record has been written, so the
/// header reserves 8 bytes for it (written as 0) and `finish()` seeks back to
/// patch the real value — hence the `Seek` bound. The only callers wrap a
/// `BufWriter<File>`, which is seekable.
pub struct CompactDatabaseWriter<W: Write + Seek> {
    encoder: zstd::Encoder<'static, W>,
    record_count: u64,
}

impl<W: Write + Seek> CompactDatabaseWriter<W> {
    /// Create a new compact database writer. Writes the header immediately.
    pub fn new(mut writer: W, gcf_ids: &[&str]) -> Result<Self> {
        // Write header uncompressed
        writer.write_all(COMPACT_MAGIC)?;
        writer.write_all(&COMPACT_VERSION.to_le_bytes())?;
        // record_count placeholder; patched in finish()
        writer.write_all(&0u64.to_le_bytes())?;
        writer.write_all(&(gcf_ids.len() as u32).to_le_bytes())?;
        for id in gcf_ids {
            let bytes = id.as_bytes();
            let len = bytes.len().min(u16::MAX as usize) as u16;
            writer.write_all(&len.to_le_bytes())?;
            writer.write_all(&bytes[..len as usize])?;
        }
        // Records section: zstd-compressed stream (level 3 = good speed/ratio balance)
        let encoder = zstd::Encoder::new(writer, 3)
            .context("Failed to create zstd encoder")?;
        Ok(Self { encoder, record_count: 0 })
    }

    /// Write a single (hash, gcf_index) record into the compressed stream.
    #[inline]
    pub fn write_record(&mut self, hash: u64, gcf_index: u32) -> io::Result<()> {
        self.encoder.write_all(&hash.to_le_bytes())?;
        self.encoder.write_all(&gcf_index.to_le_bytes())?;
        self.record_count += 1;
        Ok(())
    }

    /// Finalize the zstd stream, patch the record_count into the header, and
    /// flush. Must be called before dropping.
    pub fn finish(self) -> Result<W> {
        let count = self.record_count;
        let mut writer = self.encoder.finish().context("Failed to finalize zstd stream")?;
        // Patch the record_count placeholder now that the total is known.
        writer
            .seek(SeekFrom::Start(RECORD_COUNT_OFFSET))
            .context("Failed to seek to record_count header field")?;
        writer
            .write_all(&count.to_le_bytes())
            .context("Failed to write record_count")?;
        writer.flush().context("Failed to flush compact database")?;
        Ok(writer)
    }
}

// ---- Reader ----

/// Reads compact database files. Supports v1 (uncompressed), v2 (zstd) and v3
/// (zstd + record_count header) formats.
pub struct CompactDatabaseReader {
    reader: Box<dyn Read>,
    gcf_table: Vec<String>,
    record_count: Option<u64>,
}

impl CompactDatabaseReader {
    pub fn gcf_table(&self) -> &[String] {
        &self.gcf_table
    }

    /// Total number of records, read straight from the header (v3+). `None` for
    /// older v1/v2 files, which don't store it — count by iterating instead.
    pub fn record_count(&self) -> Option<u64> {
        self.record_count
    }

    /// Read next record. Returns (hash, gcf_index) or None at EOF.
    pub fn next_record(&mut self) -> Result<Option<(u64, u32)>> {
        let mut hash_buf = [0u8; 8];
        if let Err(e) = self.reader.read_exact(&mut hash_buf) {
            if e.kind() == io::ErrorKind::UnexpectedEof {
                return Ok(None);
            }
            return Err(e.into());
        }
        let hash = u64::from_le_bytes(hash_buf);

        let mut idx_buf = [0u8; 4];
        self.reader.read_exact(&mut idx_buf).context("Failed to read GCF index")?;
        let index = u32::from_le_bytes(idx_buf);

        Ok(Some((hash, index)))
    }
}

pub fn open_compact_reader<P: AsRef<Path>>(path: P) -> Result<CompactDatabaseReader> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Cannot open file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(IO_BUFFER_SIZE, file);

    // Read header (always uncompressed)
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic).context("Failed to read compact DB magic")?;
    if &magic != COMPACT_MAGIC {
        anyhow::bail!("Not a compact database file (invalid magic)");
    }

    let mut ver_buf = [0u8; 4];
    reader.read_exact(&mut ver_buf)?;
    let version = u32::from_le_bytes(ver_buf);
    if version != 1 && version != 2 && version != 3 {
        anyhow::bail!("Unsupported compact database version: {}", version);
    }

    // v3 stores the total record count right after the version.
    let record_count = if version == 3 {
        let mut rc_buf = [0u8; 8];
        reader.read_exact(&mut rc_buf)?;
        Some(u64::from_le_bytes(rc_buf))
    } else {
        None
    };

    let mut count_buf = [0u8; 4];
    reader.read_exact(&mut count_buf)?;
    let gcf_count = u32::from_le_bytes(count_buf) as usize;

    let mut gcf_table = Vec::with_capacity(gcf_count);
    for _ in 0..gcf_count {
        let mut len_buf = [0u8; 2];
        reader.read_exact(&mut len_buf)?;
        let len = u16::from_le_bytes(len_buf) as usize;
        let mut bytes = vec![0u8; len];
        reader.read_exact(&mut bytes)?;
        gcf_table.push(String::from_utf8(bytes).context("Invalid UTF-8 in GCF ID")?);
    }

    // Records section: wrap in zstd decoder for v2/v3, raw for v1
    let records_reader: Box<dyn Read> = if version == 2 || version == 3 {
        Box::new(zstd::Decoder::new(reader).context("Failed to create zstd decoder")?)
    } else {
        Box::new(reader)
    };

    Ok(CompactDatabaseReader { reader: records_reader, gcf_table, record_count })
}