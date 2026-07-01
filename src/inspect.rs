// Inspect 2bRAD binary database / sample files (sylph-style `inspect`).
//
// Fast2bRAD-M writes two on-disk binary formats; this command auto-detects and
// summarizes either:
//
//   1. CompactDatabase  — the taxon-level database, e.g. `BcgI.species.iibdb`.
//      Header `IIBC` + GCF (genome) table + zstd records of (tag_hash, gcf_index).
//
//   2. Binary record stream — per-genome / per-sample digests, e.g.
//      `*.iibdb` (reference) and `*.iibsp` (reads). A stream of (tag_hash, id)
//      records, optionally zstd- or gzip-compressed.
//
// By default only the header / a small preview is read (instant, even for a
// multi-GB database). `--full` does a complete pass to count every tag; add
// `--distinct` to also count distinct tag hashes (memory-heavy on large DBs).

use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Args;
use fxhash::{FxHashMap, FxHashSet};

use crate::enzymes::enzyme_by_name;
use crate::io_utils;

#[derive(Args, Debug)]
pub struct InspectArgs {
    /// One or more .iibdb / .iibsp files to inspect
    #[arg(required = true)]
    pub files: Vec<PathBuf>,

    /// Full scan: count every tag (and per-genome tag counts). Slow on large DBs
    #[arg(short = 'f', long = "full")]
    pub full: bool,

    /// With --full, also count distinct tag hashes (memory-heavy on large DBs)
    #[arg(long = "distinct")]
    pub distinct: bool,

    /// Number of example records to print per file (0 = none)
    #[arg(short = 'n', long = "records", default_value = "5")]
    pub num_records: usize,

    /// Number of genomes to list (top by tag count with --full, otherwise the first N)
    #[arg(short = 't', long = "top", default_value = "10")]
    pub top_genomes: usize,

    /// Write the report to this file instead of stdout
    #[arg(short = 'o', long = "output")]
    pub out_file_name: Option<PathBuf>,
}

pub fn run(args: InspectArgs) -> Result<()> {
    let mut writer: Box<dyn Write> = match &args.out_file_name {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).with_context(|| format!("Failed to create output file: {}", path.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    for file in &args.files {
        if let Err(e) = inspect_one(file, &args, &mut writer) {
            eprintln!("Failed to inspect {}: {:#}", file.display(), e);
        }
    }
    writer.flush()?;
    Ok(())
}

fn inspect_one(path: &Path, args: &InspectArgs, w: &mut Box<dyn Write>) -> Result<()> {
    writeln!(w, "=== File: {} ===", path.display())?;

    // Enzyme / taxon-level hints parsed from the file name (e.g. BcgI.species.iibdb)
    let (enzyme, level) = parse_name_hints(path);
    if let Some(enz) = &enzyme {
        writeln!(w, "Enzyme:        {}", enz)?;
    }
    if let Some(lvl) = &level {
        writeln!(w, "Taxon level:   {}", lvl)?;
    }
    if let Ok(meta) = std::fs::metadata(path) {
        writeln!(w, "File size:     {}", human_bytes(meta.len()))?;
    }

    if is_compact(path)? {
        inspect_compact(path, args, w)?;
    } else {
        inspect_stream(path, args, w)?;
    }
    writeln!(w)?;
    Ok(())
}

/// CompactDatabase header magic is written uncompressed, so a 4-byte peek tells
/// the two formats apart. zstd/gzip streams start with their own magic instead.
fn is_compact(path: &Path) -> Result<bool> {
    let mut f = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut magic = [0u8; 4];
    match f.read_exact(&mut magic) {
        Ok(()) => Ok(&magic == io_utils::COMPACT_MAGIC),
        Err(_) => Ok(false), // shorter than 4 bytes -> treat as (empty) stream
    }
}

fn inspect_compact(path: &Path, args: &InspectArgs, w: &mut Box<dyn Write>) -> Result<()> {
    let mut reader = io_utils::open_compact_reader(path)?;
    let gcf_table: Vec<String> = reader.gcf_table().to_vec();
    let header_count = reader.record_count();

    writeln!(w, "Format:        CompactDatabase (taxon-level database)")?;
    writeln!(w, "Genomes:       {}", gcf_table.len())?;
    if let Some(rc) = header_count {
        writeln!(w, "Total tags:    {} (from header)", rc)?;
    }

    if args.full {
        // Full pass: count tags per genome (cheap counters) and optionally distinct.
        let mut per_genome = vec![0usize; gcf_table.len()];
        let mut distinct: FxHashSet<u64> = FxHashSet::default();
        let mut total = 0usize;
        let mut examples: Vec<(u64, u32)> = Vec::new();
        while let Some((hash, idx)) = reader.next_record()? {
            total += 1;
            if (idx as usize) < per_genome.len() {
                per_genome[idx as usize] += 1;
            }
            if args.distinct {
                distinct.insert(hash);
            }
            if examples.len() < args.num_records {
                examples.push((hash, idx));
            }
        }
        // Only print the scanned total if the header didn't already provide one.
        match header_count {
            None => writeln!(w, "Total tags:    {}", total)?,
            Some(rc) if rc != total as u64 => {
                writeln!(w, "Total tags:    {} (scanned; header said {})", total, rc)?
            }
            Some(_) => {}
        }
        if args.distinct {
            writeln!(w, "Distinct tags: {}", distinct.len())?;
        }
        write_top_genomes_by_count(w, &gcf_table, &per_genome, args.top_genomes)?;
        write_compact_examples(w, &gcf_table, &examples)?;
    } else {
        // Fast: list the first genomes from the header + a few example records.
        if args.top_genomes > 0 {
            let show = args.top_genomes.min(gcf_table.len());
            writeln!(w, "First {} genomes:", show)?;
            for id in gcf_table.iter().take(show) {
                writeln!(w, "  {}", id)?;
            }
        }
        let mut examples: Vec<(u64, u32)> = Vec::new();
        while examples.len() < args.num_records {
            match reader.next_record()? {
                Some(rec) => examples.push(rec),
                None => break,
            }
        }
        write_compact_examples(w, &gcf_table, &examples)?;
        if header_count.is_some() {
            writeln!(w, "(run with --full --distinct for per-genome and unique-tag counts)")?;
        } else {
            writeln!(w, "(run with --full to count all tags; add --distinct for unique-tag counts)")?;
        }
    }
    Ok(())
}

fn write_top_genomes_by_count(
    w: &mut Box<dyn Write>,
    gcf_table: &[String],
    per_genome: &[usize],
    top: usize,
) -> Result<()> {
    if top == 0 || gcf_table.is_empty() {
        return Ok(());
    }
    let mut idx: Vec<usize> = (0..gcf_table.len()).collect();
    idx.sort_by(|&a, &b| per_genome[b].cmp(&per_genome[a]));
    let show = top.min(idx.len());
    writeln!(w, "Top {} genomes by tag count:", show)?;
    for &i in idx.iter().take(show) {
        writeln!(w, "  {:<24} {}", gcf_table[i], per_genome[i])?;
    }
    Ok(())
}

fn write_compact_examples(w: &mut Box<dyn Write>, gcf_table: &[String], examples: &[(u64, u32)]) -> Result<()> {
    if examples.is_empty() {
        return Ok(());
    }
    writeln!(w, "Example records (tag_hash -> genome):")?;
    for (hash, idx) in examples {
        let name = gcf_table.get(*idx as usize).map(|s| s.as_str()).unwrap_or("<out-of-range>");
        writeln!(w, "  {:016x}  {}", hash, name)?;
    }
    Ok(())
}

/// What a binary record stream represents, inferred from the file name and the
/// record id structure. Determines how records are grouped and labelled.
#[derive(Clone, Copy, PartialEq)]
enum StreamKind {
    /// `*.iibsp`: per-read sample tags (ids are read names) — not grouped.
    SampleReads,
    /// extract default `*.iibdb`: id = contig.
    GenomeContig,
    /// extract `--record-pos` `*.pos.iibdb`: id = `contig|offset`.
    GenomeContigPos,
    /// build-qual-db `*.enzyme.iibdb` digest: id = `gcf|idx|scaffold|pos|..`.
    GenomeDigest,
}

impl StreamKind {
    fn label(self) -> &'static str {
        match self {
            StreamKind::GenomeDigest => "Genomes",
            _ => "Contigs",
        }
    }
    fn describe(self) -> &'static str {
        match self {
            StreamKind::SampleReads => "sample reads (per-read 2bRAD tags)",
            StreamKind::GenomeContig => "reference genome tags, grouped by contig (no positions)",
            StreamKind::GenomeContigPos => "reference genome tags with positions (id = contig|offset)",
            StreamKind::GenomeDigest => "per-genome digest (id = gcf|idx|scaffold|pos|..)",
        }
    }
    /// Grouping key for one record id, or None when records aren't grouped.
    fn group_key<'a>(self, id: &'a str) -> Option<&'a str> {
        match self {
            StreamKind::SampleReads => None,
            StreamKind::GenomeContig => Some(id),
            StreamKind::GenomeContigPos | StreamKind::GenomeDigest => {
                Some(id.split('|').next().unwrap_or(id))
            }
        }
    }
}

/// Classify from the file name. Generic `*.iibdb` returns `None` and is
/// disambiguated (contig vs digest) from the first record id.
fn classify_by_name(path: &Path) -> Option<StreamKind> {
    let name = path.to_string_lossy();
    if name.ends_with(".iibsp") {
        Some(StreamKind::SampleReads)
    } else if name.ends_with(".pos.iibdb") {
        Some(StreamKind::GenomeContigPos)
    } else {
        None
    }
}

fn inspect_stream(path: &Path, args: &InspectArgs, w: &mut Box<dyn Write>) -> Result<()> {
    let mut reader = io_utils::open_binary_reader(path)?;
    // Generic `.iibdb` (not `.pos.iibdb`/`.iibsp`): contig ids have no `|`,
    // digest ids do — decide on the first record.
    let decide = |id: &str| {
        if id.contains('|') { StreamKind::GenomeDigest } else { StreamKind::GenomeContig }
    };

    let mut kind = classify_by_name(path);
    let mut examples: Vec<(u64, String)> = Vec::new();
    let mut buf = String::with_capacity(128);

    if args.full {
        let mut distinct: FxHashSet<u64> = FxHashSet::default();
        let mut groups: FxHashMap<String, usize> = FxHashMap::default();
        let mut total = 0usize;
        while let Some(hash) = reader.next_record_reuse(&mut buf)? {
            total += 1;
            if kind.is_none() {
                kind = Some(decide(&buf));
            }
            if args.distinct {
                distinct.insert(hash);
            }
            if let Some(key) = kind.unwrap().group_key(&buf) {
                *groups.entry(key.to_string()).or_insert(0) += 1;
            }
            if examples.len() < args.num_records {
                examples.push((hash, buf.clone()));
            }
        }
        let kind = kind.unwrap_or(StreamKind::GenomeContig);
        writeln!(w, "Format:        Binary record stream — {}", kind.describe())?;
        writeln!(w, "{:<15}{}", "Total tags:", total)?;
        if args.distinct {
            writeln!(w, "{:<15}{}", "Distinct tags:", distinct.len())?;
        }
        if !groups.is_empty() && args.top_genomes > 0 {
            writeln!(w, "{:<15}{}", format!("{}:", kind.label()), groups.len())?;
            let mut v: Vec<(&String, &usize)> = groups.iter().collect();
            v.sort_by(|a, b| b.1.cmp(a.1));
            let show = args.top_genomes.min(v.len());
            writeln!(w, "Top {} {} by tag count:", show, kind.label().to_lowercase())?;
            for (key, count) in v.into_iter().take(show) {
                writeln!(w, "  {:<24} {}", key, count)?;
            }
        }
    } else {
        while examples.len() < args.num_records {
            match reader.next_record_reuse(&mut buf)? {
                Some(hash) => {
                    if kind.is_none() {
                        kind = Some(decide(&buf));
                    }
                    examples.push((hash, buf.clone()));
                }
                None => break,
            }
        }
        let kind = kind.unwrap_or(StreamKind::GenomeContig);
        writeln!(w, "Format:        Binary record stream — {}", kind.describe())?;
    }

    if !examples.is_empty() {
        writeln!(w, "Example records (tag_hash : id):")?;
        for (hash, id) in &examples {
            writeln!(w, "  {:016x}  {}", hash, id)?;
        }
    }
    if !args.full {
        writeln!(w, "(run with --full to count all tags; add --distinct for unique-tag counts)")?;
    }
    Ok(())
}

/// Pull an enzyme name and/or taxonomy level out of a file name like
/// `BcgI.species.iibdb` or `GCF_000001.1.BcgI.iibdb`.
fn parse_name_hints(path: &Path) -> (Option<String>, Option<String>) {
    const LEVELS: [&str; 8] = [
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain",
    ];
    let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
    let mut enzyme = None;
    let mut level = None;
    for token in name.split('.') {
        if enzyme.is_none() {
            if let Some(enz) = enzyme_by_name(token) {
                enzyme = Some(enz.name.to_string());
                continue;
            }
        }
        if level.is_none() && LEVELS.contains(&token) {
            level = Some(token.to_string());
        }
    }
    (enzyme, level)
}

fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KB", "MB", "GB", "TB"];
    let mut v = n as f64;
    let mut u = 0;
    while v >= 1024.0 && u < UNITS.len() - 1 {
        v /= 1024.0;
        u += 1;
    }
    if u == 0 {
        format!("{} {}", n, UNITS[u])
    } else {
        format!("{:.2} {}", v, UNITS[u])
    }
}
