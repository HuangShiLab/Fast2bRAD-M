use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow, bail};
use clap::Args;
use flate2::read::GzDecoder;
use needletail::parse_fastx_file;
use rayon;
use tracing;

use f2brad_core::enzymes::{enzyme_by_id, enzyme_by_name};

#[derive(Args, Debug)]
pub struct GenotypeArgs {
    /// Host tag database TSV produced by `f2brad-host build-db`
    #[arg(short = 'd', long = "db", required = true)]
    pub db: PathBuf,

    /// Read 1 FASTQ/FASTA (may be gzip-compressed)
    #[arg(short = '1', long = "r1", required = true)]
    pub r1: PathBuf,

    /// Optional read 2 FASTQ/FASTA (paired-end)
    #[arg(short = '2', long = "r2")]
    pub r2: Option<PathBuf>,

    /// Enzyme name (e.g. BcgI, BsaXI, AlfI) or numeric ID (1–16)
    #[arg(short = 's', long = "site", required = true)]
    pub enzyme_site: String,

    /// Output directory
    #[arg(short = 'o', long = "output", required = true)]
    pub output_dir: PathBuf,

    /// Minimum base quality (Phred) for a base to contribute to a genotype
    #[arg(short = 'q', long = "min-qual", default_value = "20")]
    pub min_qual: u8,

    /// Minimum mapping quality proxy: maximum Hamming distance to a reference tag
    #[arg(long = "max-mismatch", default_value = "2")]
    pub max_mismatch: usize,

    /// Minimum per-locus depth to emit a genotype record
    #[arg(long = "min-depth", default_value = "4")]
    pub min_depth: usize,

    /// Number of parallel threads
    #[arg(short = 'j', long = "threads", default_value = "4")]
    pub threads: usize,
}

#[derive(Debug, Clone)]
pub struct Locus {
    pub contig: String,
    pub pos: usize,
    pub seq: Vec<u8>,
    pub canonical: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct HostDb {
    pub loci: Vec<Locus>,
    pub index: HostTagIndex,
}

#[derive(Debug, Clone)]
pub struct HostTagIndex {
    tag_len: usize,
    part_size: usize,
    buckets: Vec<HashMap<Vec<u8>, Vec<usize>>>,
    canonical: Vec<Vec<u8>>,
}

impl HostTagIndex {
    pub fn new(loci: &[Locus], max_mismatch: usize) -> Self {
        if loci.is_empty() {
            return Self {
                tag_len: 0,
                part_size: 0,
                buckets: Vec::new(),
                canonical: Vec::new(),
            };
        }
        let tag_len = loci[0].canonical.len();
        let parts = max_mismatch + 1;
        let part_size = (tag_len + parts - 1) / parts;
        let mut buckets: Vec<HashMap<Vec<u8>, Vec<usize>>> = vec![HashMap::new(); parts];
        for (i, locus) in loci.iter().enumerate() {
            for p in 0..parts {
                let start = p * part_size;
                let end = ((p + 1) * part_size).min(tag_len);
                if start >= tag_len {
                    continue;
                }
                let part = locus.canonical[start..end].to_vec();
                buckets[p].entry(part).or_default().push(i);
            }
        }
        Self {
            tag_len,
            part_size,
            buckets,
            canonical: loci.iter().map(|l| l.canonical.clone()).collect(),
        }
    }

    /// Find the closest reference tag within `max_mismatch` of `query`. Returns
    /// (locus index, Hamming distance) for the closest match.
    pub fn find(&self, query: &[u8], max_mismatch: usize) -> Option<(usize, usize)> {
        let parts = max_mismatch + 1;
        let mut candidates: Vec<usize> = Vec::new();
        for p in 0..parts {
            let start = p * self.part_size;
            let end = ((p + 1) * self.part_size).min(self.tag_len);
            if start >= self.tag_len {
                continue;
            }
            let part = &query[start..end];
            if let Some(idxs) = self.buckets[p].get(part) {
                candidates.extend(idxs.iter().copied());
            }
        }
        candidates.sort_unstable();
        candidates.dedup();

        let mut best: Option<(usize, usize)> = None;
        for &idx in &candidates {
            let dist = hamming_distance(query, &self.canonical[idx]);
            if dist <= max_mismatch {
                if best.is_none() || dist < best.unwrap().1 {
                    best = Some((idx, dist));
                    if dist == 0 {
                        break;
                    }
                }
            }
        }
        best
    }
}

pub fn run(args: GenotypeArgs) -> Result<()> {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();

    let enzyme = if let Ok(site_num) = args.enzyme_site.parse::<u8>() {
        enzyme_by_id(site_num).ok_or_else(|| anyhow!("Invalid enzyme ID"))?
    } else {
        enzyme_by_name(&args.enzyme_site).ok_or_else(|| anyhow!("Invalid enzyme name"))?
    };

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("Failed to create output directory: {}", args.output_dir.display()))?;

    tracing::info!("Loading host DB from {}", args.db.display());
    let db = load_host_db(&args.db)?;
    tracing::info!("Loaded {} host tags", db.loci.len());

    let max_mismatch = args.max_mismatch.max(1).min(3);
    tracing::info!("Building host tag index for Hamming distance <= {}", max_mismatch);
    let db = HostDb {
        index: HostTagIndex::new(&db.loci, max_mismatch),
        loci: db.loci,
    };

    tracing::info!("Genotyping reads from {}", args.r1.display());
    if let Some(r2) = &args.r2 {
        tracing::info!("Paired-end with {}", r2.display());
    }

    let mut pileup = Pileup::new(db.loci.len(), db.loci[0].canonical.len());

    // Process reads. For paired-end input, process mates in lockstep so a tag
    // observed on both R1 and R2 is counted once per fragment, not twice.
    if let Some(r2) = &args.r2 {
        process_paired(&args.r1, r2, enzyme, &db, max_mismatch, args.min_qual, &mut pileup)?;
    } else {
        process_file(&args.r1, enzyme, &db, max_mismatch, args.min_qual, &mut pileup)?;
    }

    tracing::info!(
        "Pileup complete: {} loci covered",
        pileup.depth.iter().filter(|d| d.iter().any(|&x| x > 0)).count()
    );

    // Call genotypes and write VCF.
    let vcf_path = args.output_dir.join("genotypes.vcf");
    write_vcf(&vcf_path, &db.loci, &pileup, args.min_depth)?;
    tracing::info!("Wrote {}", vcf_path.display());

    // Write GEMMA BIMBAM mean-genotype file.
    let bimbam_path = args.output_dir.join("dosages.bimbam");
    write_bimbam(&bimbam_path, &db.loci, &pileup, args.min_depth)?;
    tracing::info!("Wrote {}", bimbam_path.display());

    Ok(())
}

pub fn load_host_db(path: &PathBuf) -> Result<HostDb> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open host DB: {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut loci = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 10 {
            continue;
        }
        let contig = parts[0].to_string();
        let pos: usize = parts[1].parse().with_context(|| format!("Invalid pos on line {}: {}", i + 1, parts[1]))?;
        let seq = parts[3].as_bytes().to_vec();
        let canonical = parts[4].as_bytes().to_vec();

        loci.push(Locus { contig, pos, seq, canonical });
    }

    if loci.is_empty() {
        bail!("Host DB contains no loci");
    }

    Ok(HostDb { loci, index: HostTagIndex::new(&[], 2) })
}

/// One tag-to-locus assignment from a single read.
#[derive(Debug, Clone)]
pub struct ReadMatch {
    locus_idx: usize,
    dist: usize,
    tag_seq: Vec<u8>,
    qual: Vec<u8>,
    same_strand: bool,
}

pub fn extract_matches(
    seq_bytes: &[u8],
    qual: Option<&[u8]>,
    enzyme: &f2brad_core::enzymes::Enzyme,
    db: &HostDb,
    max_mismatch: usize,
) -> Vec<ReadMatch> {
    let mut matches = Vec::new();
    let tags = enzyme.find_all_tags(seq_bytes);
    for (offset, len) in tags {
        let tag_seq = &seq_bytes[offset..offset + len];
        let canonical = canonicalize(tag_seq);
        if let Some((idx, dist)) = db.index.find(&canonical, max_mismatch) {
            let locus = &db.loci[idx];
            let ref_tag = &locus.seq;
            let ref_rc = reverse_complement(ref_tag);
            // Determine orientation from the raw read tag, not the canonical form.
            let same_strand = hamming_distance(tag_seq, ref_tag) <= hamming_distance(tag_seq, &ref_rc);
            let qual_window = qual.map(|q| q[offset..offset + len].to_vec()).unwrap_or_default();
            matches.push(ReadMatch {
                locus_idx: idx,
                dist,
                tag_seq: tag_seq.to_vec(),
                qual: qual_window,
                same_strand,
            });
        }
    }
    matches
}

fn add_match_to_pileup(m: &ReadMatch, min_qual: u8, db: &HostDb, pileup: &mut Pileup) {
    let locus = &db.loci[m.locus_idx];
    pileup.add_read(&locus.seq, &m.tag_seq, &m.qual, m.same_strand, min_qual, m.locus_idx, m.dist);
}

fn process_file(
    path: &PathBuf,
    enzyme: &f2brad_core::enzymes::Enzyme,
    db: &HostDb,
    max_mismatch: usize,
    min_qual: u8,
    pileup: &mut Pileup,
) -> Result<()> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to open reads: {}", path.display()))?;

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("Failed to read record from {}", path.display()))?;
        let seq = record.seq();
        let seq_bytes = seq.as_ref();
        let qual = record.qual();
        for m in extract_matches(seq_bytes, qual, enzyme, db, max_mismatch) {
            add_match_to_pileup(&m, min_qual, db, pileup);
        }
    }

    Ok(())
}

fn process_paired(
    path1: &PathBuf,
    path2: &PathBuf,
    enzyme: &f2brad_core::enzymes::Enzyme,
    db: &HostDb,
    max_mismatch: usize,
    min_qual: u8,
    pileup: &mut Pileup,
) -> Result<()> {
    let mut reader1 = parse_fastx_file(path1)
        .with_context(|| format!("Failed to open reads: {}", path1.display()))?;
    let mut reader2 = parse_fastx_file(path2)
        .with_context(|| format!("Failed to open reads: {}", path2.display()))?;

    let mut n_pairs = 0usize;
    loop {
        let rec1 = reader1.next();
        let rec2 = reader2.next();
        match (rec1, rec2) {
            (None, None) => break,
            (Some(r1), Some(r2)) => {
                let r1 = r1.with_context(|| format!("Failed to read record from {}", path1.display()))?;
                let r2 = r2.with_context(|| format!("Failed to read record from {}", path2.display()))?;
                let seq1 = r1.seq();
                let qual1 = r1.qual();
                let seq2 = r2.seq();
                let qual2 = r2.qual();

                let matches1 = extract_matches(seq1.as_ref(), qual1, enzyme, db, max_mismatch);
                let matches2 = extract_matches(seq2.as_ref(), qual2, enzyme, db, max_mismatch);

                // Merge matches by locus, keeping the observation with the
                // fewest mismatches. This prevents a tag seen on both mates of
                // the same fragment from being counted twice.
                let mut merged: HashMap<usize, ReadMatch> = HashMap::new();
                for m in matches1.into_iter().chain(matches2.into_iter()) {
                    merged
                        .entry(m.locus_idx)
                        .and_modify(|existing| {
                            if m.dist < existing.dist {
                                *existing = m.clone();
                            }
                        })
                        .or_insert(m);
                }
                for m in merged.values() {
                    add_match_to_pileup(m, min_qual, db, pileup);
                }
                n_pairs += 1;
            }
            _ => bail!(
                "Paired input files have different numbers of reads: {} vs {}",
                path1.display(),
                path2.display()
            ),
        }
    }

    tracing::info!("Processed {} paired-end fragments", n_pairs);
    Ok(())
}

pub fn canonicalize(seq: &[u8]) -> Vec<u8> {
    let fwd: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
    let rc = reverse_complement(&fwd);
    if fwd <= rc { fwd } else { rc }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            x => x,
        })
        .collect()
}

pub fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Per-locus, per-position pileup.
struct Pileup {
    /// depth[locus][position] = number of contributing reads
    depth: Vec<Vec<usize>>,
    /// counts[locus][position][base] = number of observations
    counts: Vec<Vec<[usize; 4]>>,
    /// qual_sums[locus][position][base] = sum of quality scores
    qual_sums: Vec<Vec<[usize; 4]>>,
    /// mismatch_counts[locus] = how many reads matched with >0 mismatches
    mismatch_reads: Vec<usize>,
}

const BASE_TO_IDX: [usize; 256] = {
    let mut arr = [4usize; 256];
    arr[b'A' as usize] = 0;
    arr[b'C' as usize] = 1;
    arr[b'G' as usize] = 2;
    arr[b'T' as usize] = 3;
    arr
};

const IDX_TO_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl Pileup {
    fn new(n_loci: usize, tag_len: usize) -> Self {
        Self {
            depth: vec![vec![0; tag_len]; n_loci],
            counts: vec![vec![[0; 4]; tag_len]; n_loci],
            qual_sums: vec![vec![[0; 4]; tag_len]; n_loci],
            mismatch_reads: vec![0; n_loci],
        }
    }

    fn add_read(
        &mut self,
        ref_tag: &[u8],
        read_tag: &[u8],
        qual: &[u8],
        same_strand: bool,
        min_qual: u8,
        locus_idx: usize,
        dist: usize,
    ) {
        if dist > 0 {
            self.mismatch_reads[locus_idx] += 1;
        }
        let tag_len = ref_tag.len();
        for i in 0..tag_len {
            let _ref_base = ref_tag[i];
            let (read_base, q) = if same_strand {
                (read_tag[i], qual.get(i).copied().unwrap_or(0))
            } else {
                // Read is reverse-complement; position i in reference corresponds
                // to position tag_len-1-i in read, and base is RC.
                let ri = tag_len - 1 - i;
                (complement(read_tag[ri]), qual.get(ri).copied().unwrap_or(0))
            };
            if q < min_qual || q == 0 {
                continue;
            }
            if read_base == b'N' || read_base == b'n' {
                continue;
            }
            // Skip if read base does not match reference and we are requiring
            // exact matching at this position. (We already allow up to max_mismatch
            // across the tag, so no extra skip here.)
            let bidx = BASE_TO_IDX[read_base as usize];
            if bidx > 3 {
                continue;
            }
            self.depth[locus_idx][i] += 1;
            self.counts[locus_idx][i][bidx] += 1;
            self.qual_sums[locus_idx][i][bidx] += (q - 33).min(93) as usize;
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        x => x,
    }
}

fn write_vcf(
    path: &PathBuf,
    loci: &[Locus],
    pileup: &Pileup,
    min_depth: usize,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create VCF: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=f2brad-host genotype")?;
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(writer, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods (log10)\">")?;
    writeln!(writer, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">")?;
    writeln!(writer, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">")?;
    writeln!(writer, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">")?;
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")?;

    // Collect callable records and sort by contig/position.
    let mut records: Vec<(usize, usize, usize)> = Vec::new(); // (locus_idx, ref_pos, tag_pos)
    for (locus_idx, locus) in loci.iter().enumerate() {
        for (tag_pos, &depth) in pileup.depth[locus_idx].iter().enumerate() {
            if depth >= min_depth {
                records.push((locus_idx, locus.pos + tag_pos, tag_pos));
            }
        }
    }
    records.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
    // Stable grouping by contig is implicit because locus.pos is monotonic within
    // a contig in the input DB, but contig order is alphabetical. We keep the
    // DB order (records are already grouped by locus).

    let mut prev_locus = usize::MAX;
    for (locus_idx, genomic_pos, tag_pos) in records {
        if locus_idx != prev_locus {
            // Nothing special needed; contig is carried in the locus.
            prev_locus = locus_idx;
        }
        let locus = &loci[locus_idx];
        let ref_base = locus.seq[tag_pos];
        let counts = pileup.counts[locus_idx][tag_pos];
        let depth = pileup.depth[locus_idx][tag_pos];

        // Pick the strongest non-REF allele as ALT.
        let ref_idx = BASE_TO_IDX[ref_base as usize];
        if ref_idx > 3 {
            continue; // skip positions with non-ACGT reference
        }
        let mut alt_idx = None;
        let mut alt_count = 0usize;
        for i in 0..4 {
            if i == ref_idx {
                continue;
            }
            if counts[i] > alt_count {
                alt_count = counts[i];
                alt_idx = Some(i);
            }
        }

        let (alt_base, alt_str) = match alt_idx {
            Some(idx) => (IDX_TO_BASE[idx], String::from_utf8_lossy(&[IDX_TO_BASE[idx]]).to_string()),
            None => (0, ".".to_string()),
        };

        // Genotype likelihoods under diploid model.
        let gl = genotype_likelihoods(ref_base, alt_base, &counts, &pileup.qual_sums[locus_idx][tag_pos]);
        let pl = phred_scale(&gl);
        let best_g = best_genotype(&gl);

        let gt = match best_g {
            0 => "0/0",
            1 => "0/1",
            2 => "1/1",
            _ => "./.",
        };

        let ad_ref = counts[ref_idx];
        let ad_alt = alt_idx.map(|i| counts[i]).unwrap_or(0);

        let gl_str = gl.iter().map(|x| format!("{:.3}", x)).collect::<Vec<_>>().join(",");
        let pl_str = pl.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");

        writeln!(
            writer,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:GL:PL:DP:AD\t{}:{}:{}:{}:{}",
            locus.contig,
            genomic_pos + 1, // VCF is 1-based
            ref_base as char,
            alt_str,
            gt,
            gl_str,
            pl_str,
            depth,
            format!("{},{}", ad_ref, ad_alt)
        )?;
    }

    writer.flush()?;
    Ok(())
}

/// Write a GEMMA BIMBAM-format mean-genotype file.
/// Columns: SNP_id A1_allele A2_allele dosage
fn write_bimbam(
    path: &PathBuf,
    loci: &[Locus],
    pileup: &Pileup,
    min_depth: usize,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create BIMBAM file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    for (locus_idx, locus) in loci.iter().enumerate() {
        for (tag_pos, &depth) in pileup.depth[locus_idx].iter().enumerate() {
            if depth < min_depth {
                continue;
            }
            let ref_base = locus.seq[tag_pos];
            let ref_idx = BASE_TO_IDX[ref_base as usize];
            if ref_idx > 3 {
                continue;
            }
            let counts = pileup.counts[locus_idx][tag_pos];

            // Pick the strongest non-REF allele as A2.
            let mut alt_idx = None;
            let mut alt_count = 0usize;
            for i in 0..4 {
                if i == ref_idx {
                    continue;
                }
                if counts[i] > alt_count {
                    alt_count = counts[i];
                    alt_idx = Some(i);
                }
            }
            let alt_base = alt_idx.map(|i| IDX_TO_BASE[i]).unwrap_or(b'.');

            let gl = genotype_likelihoods(ref_base, alt_base, &counts, &pileup.qual_sums[locus_idx][tag_pos]);
            let dosage = genotype_dosage(&gl);

            let snp_id = format!("{}_{}", locus.contig, locus.pos + tag_pos + 1);
            writeln!(
                writer,
                "{}\t{}\t{}\t{:.4}",
                snp_id,
                ref_base as char,
                alt_base as char,
                dosage
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Compute expected dosage = P(het) + 2*P(hom-alt) from log10 likelihoods.
fn genotype_dosage(gl: &[f64; 3]) -> f64 {
    // Convert log10 likelihoods to probabilities.
    let p: Vec<f64> = gl.iter().map(|&x| 10f64.powf(x)).collect();
    let sum: f64 = p.iter().sum();
    if sum == 0.0 {
        return 0.0;
    }
    let _p0 = p[0] / sum;
    let p1 = p[1] / sum;
    let p2 = p[2] / sum;
    p1 + 2.0 * p2
}

/// Compute log10 genotype likelihoods for genotypes 0/0, 0/1, 1/1.
/// If alt_base == 0 (no alt allele observed), likelihoods are still computed
/// against an arbitrary alt but only 0/0 is emitted.
fn genotype_likelihoods(ref_base: u8, alt_base: u8, counts: &[usize; 4], qual_sums: &[usize; 4]) -> [f64; 3] {
    let ref_idx = BASE_TO_IDX[ref_base as usize];
    let alt_idx = if alt_base == 0 { 3 } else { BASE_TO_IDX[alt_base as usize] };

    let mut log_likes = [0.0f64; 3];
    for i in 0..4 {
        if counts[i] == 0 {
            continue;
        }
        let n = counts[i];
        let avg_q = qual_sums[i] as f64 / n as f64;
        let e = 10f64.powf(-avg_q / 10.0).clamp(1e-6, 0.75);
        let correct = 1.0 - e;
        let wrong = e / 3.0;

        let p_obs_given_ref = if i == ref_idx { correct } else { wrong };
        let p_obs_given_alt = if i == alt_idx { correct } else { wrong };

        // Diploid genotype probabilities: P(obs | AA) = P(obs | A)^2,
        // P(obs | AB) = P(obs | A)*P(obs | B), P(obs | BB) = P(obs | B)^2.
        // We use the per-read likelihood averaged over the two haplotypes.
        let like_hom_ref = p_obs_given_ref;
        let like_het = (p_obs_given_ref + p_obs_given_alt) / 2.0;
        let like_hom_alt = p_obs_given_alt;

        log_likes[0] += (like_hom_ref.ln() / std::f64::consts::LN_10) * n as f64;
        log_likes[1] += (like_het.ln() / std::f64::consts::LN_10) * n as f64;
        log_likes[2] += (like_hom_alt.ln() / std::f64::consts::LN_10) * n as f64;
    }

    log_likes
}

fn phred_scale(gl: &[f64; 3]) -> [i32; 3] {
    let max = gl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0i32; 3];
    for i in 0..3 {
        let diff = (gl[i] - max) * -10.0;
        pl[i] = diff.clamp(0.0, 999.0) as i32;
    }
    pl
}

fn best_genotype(gl: &[f64; 3]) -> usize {
    let mut best = 0usize;
    let mut best_val = gl[0];
    for i in 1..3 {
        if gl[i] > best_val {
            best_val = gl[i];
            best = i;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileup_records_same_strand_read() {
        let mut pileup = Pileup::new(1, 4);
        let ref_tag = b"ACGT";
        let read_tag = b"ACGT";
        let qual = b"IIII"; // Phred+33 Q40
        pileup.add_read(ref_tag, read_tag, qual, true, 20, 0, 0);
        assert_eq!(pileup.depth[0], vec![1, 1, 1, 1]);
        assert_eq!(pileup.counts[0][0], [1, 0, 0, 0]);
        assert_eq!(pileup.counts[0][3], [0, 0, 0, 1]);
    }

    #[test]
    fn pileup_rc_read_matches_reference_orientation() {
        let mut pileup = Pileup::new(1, 4);
        // Use a non-palindrome: ref_tag = "AAGG", canonical = "AAGG", RC = "CCTT".
        let ref_tag = b"AAGG";
        let read_tag = b"CCTT"; // raw read is RC of ref
        let qual = b"IIII";
        pileup.add_read(ref_tag, read_tag, qual, false, 20, 0, 0);
        // After reverse-complementing the read, bases should align to ref AAGG.
        assert_eq!(pileup.depth[0], vec![1, 1, 1, 1]);
        assert_eq!(pileup.counts[0][0], [1, 0, 0, 0]); // A at ref position 0
        assert_eq!(pileup.counts[0][1], [1, 0, 0, 0]); // A at ref position 1
        assert_eq!(pileup.counts[0][2], [0, 0, 1, 0]); // G at ref position 2
        assert_eq!(pileup.counts[0][3], [0, 0, 1, 0]); // G at ref position 3
    }

    #[test]
    fn genotype_likelihood_prefers_homozygote_when_all_match() {
        let counts = [10, 0, 0, 0]; // 10 A observations
        let qual_sums = [330, 0, 0, 0]; // avg Q 33 (after -33)
        let gl = genotype_likelihoods(b'A', b'C', &counts, &qual_sums);
        assert!(gl[0] > gl[1]);
        assert!(gl[0] > gl[2]);
    }

    #[test]
    fn phred_scale_zeros_best_genotype() {
        let gl = [-10.0, -5.0, -20.0];
        let pl = phred_scale(&gl);
        assert_eq!(pl[1], 0);
        assert!(pl[0] > 0);
        assert!(pl[2] > 0);
    }
}
