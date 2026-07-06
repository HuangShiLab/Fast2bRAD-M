use fxhash::FxHashSet;
use regex::bytes::Regex;
use std::sync::OnceLock;

#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub offset: usize,
    pub motif: &'static [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct Pattern {
    pub anchors: &'static [Anchor],
}

impl Pattern {
    pub fn matches(&self, window: &[u8]) -> bool {
        self.anchors.iter().all(|anchor| {
            let start = anchor.offset;
            let end = start + anchor.motif.len();
            end <= window.len() && &window[start..end] == anchor.motif
        })
    }
}

/// How an enzyme's recognition site is matched against a candidate window.
///
/// Most Type IIB enzymes have recognition sites made up entirely of fixed
/// literal bases, so their patterns are checked with cheap byte-slice
/// equality at fixed offsets (`Exact`). Three enzymes (BaeI, HaeIV, Hin4I)
/// have IUPAC-degenerate positions in their recognition sequence (e.g. a
/// pyrimidine-only or purine-only position) that cannot be expressed as a
/// literal byte string; these are matched with a regular expression instead
/// (`Degenerate`), mirroring the `[CT]`/`[AG]`/`[GAC]` character classes
/// used in the original `2bRADExtraction.pl`. The regexes are intentionally
/// left unanchored (no `^...$`) and matched with `find`/`find_at` rather
/// than `is_match` on a fixed window, so the regex engine's literal
/// prefiltering can skip ahead between candidate positions instead of
/// re-checking the full pattern at every offset — see `find_all_tags`.
#[derive(Debug, Clone, Copy)]
pub enum MatchRule {
    Exact(&'static [Pattern]),
    /// Function returning the lazily-compiled, process-wide-cached list of
    /// regexes (one per strand orientation) for this enzyme.
    Degenerate(fn() -> &'static [Regex]),
}

#[derive(Debug, Clone, Copy)]
pub struct Enzyme {
    pub name: &'static str,
    pub id: u8,
    pub tag_length: usize,
    pub rule: MatchRule,
}

// Create a static lookup table of size 256 (the range of u8)
// Only positions for ATCGatcg are true; all others (including N) are false
const ATCG_TABLE: [bool; 256] = {
    let mut table = [false; 256];
    table[b'A' as usize] = true; table[b'a' as usize] = true;
    table[b'T' as usize] = true; table[b't' as usize] = true;
    table[b'C' as usize] = true; table[b'c' as usize] = true;
    table[b'G' as usize] = true; table[b'g' as usize] = true;
    table
};

#[inline]
fn is_pure_atcg(window: &[u8]) -> bool {
    window.iter().all(|&b| ATCG_TABLE[b as usize])
}

impl Enzyme {
    /// Find all matching tag positions and lengths in the sequence
    /// (deduplicated). Used for reference-genome digestion (Type 1) and
    /// shotgun-read scanning (Type 2), both of which look for *every*
    /// occurrence of the site, including overlapping ones — this mirrors
    /// `2bRADExtraction.pl`'s `Electronic_enzyme`/`fastq` subroutines, which
    /// rewind the regex cursor to `match_start + 1` after each hit instead
    /// of continuing from the end of the match.
    pub fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let mut positions = FxHashSet::default();
        let len = self.tag_length;
        if sequence.len() < len {
            return Vec::new();
        }

        match self.rule {
            MatchRule::Exact(patterns) => {
                for pattern in patterns {
                    for offset in 0..=(sequence.len() - len) {
                        let window = &sequence[offset..offset + len];
                        if pattern.matches(window) && is_pure_atcg(window) {
                            positions.insert((offset, len));
                        }
                    }
                }
            }
            MatchRule::Degenerate(get_regexes) => {
                // Uses Regex::find_at + rewind-to-(match_start + 1), mirroring
                // 2bRADExtraction.pl's overlap-scanning technique, instead of
                // testing `is_match` at every single offset: that naive
                // per-offset approach defeats the regex engine's literal
                // prefiltering and was measured ~15-17x slower on realistic
                // genome/read-scale inputs for no behavioral benefit — the
                // patterns use only fixed-count `{n}` repetitions, so a match
                // found anywhere by find_at already has length == tag_length.
                for re in get_regexes() {
                    let mut start = 0usize;
                    while start <= sequence.len() {
                        match re.find_at(sequence, start) {
                            Some(m) => {
                                positions.insert((m.start(), m.end() - m.start()));
                                start = m.start() + 1;
                            }
                            None => break,
                        }
                    }
                }
            }
        }

        let mut result: Vec<_> = positions.into_iter().collect();
        result.sort_unstable();
        result
    }

    /// Find only the left-most matching window, trying each pattern in
    /// declaration order and exhausting all offsets for one pattern before
    /// moving to the next. This mirrors `Single_Lable`'s per-site, leftmost,
    /// non-greedy match-then-`last` behaviour used for Type 3 (single-tag)
    /// input, where only the first tag found in a read is kept.
    pub fn find_first_tag(&self, sequence: &[u8]) -> Option<(usize, usize)> {
        let len = self.tag_length;
        if sequence.len() < len {
            return None;
        }

        match self.rule {
            MatchRule::Exact(patterns) => {
                for pattern in patterns {
                    for offset in 0..=(sequence.len() - len) {
                        let window = &sequence[offset..offset + len];
                        if pattern.matches(window) && is_pure_atcg(window) {
                            return Some((offset, len));
                        }
                    }
                }
                None
            }
            MatchRule::Degenerate(get_regexes) => {
                // A single Regex::find gives the left-most match directly —
                // no offset loop needed (see find_all_tags for why find_at
                // rather than a per-offset is_match scan).
                for re in get_regexes() {
                    if let Some(m) = re.find(sequence) {
                        return Some((m.start(), m.end() - m.start()));
                    }
                }
                None
            }
        }
    }
}

// ========== Definitions for 16 enzymes ==========
// Tag lengths and anchor positions below are derived directly from the
// `@site` regex patterns in `2bRADExtraction.pl` (shihuang047/2bRAD-M), so
// that the extracted tag sequence/length matches the original Perl
// implementation exactly for every enzyme, not just the default BcgI.

// 1. CspCI (tag_length=33): [AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}
const CSPCI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 11, motif: b"CAA" },
    Anchor { offset: 19, motif: b"GTGG" },
];
const CSPCI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 10, motif: b"CCAC" },
    Anchor { offset: 19, motif: b"TTG" },
];
const CSPCI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &CSPCI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &CSPCI_PATTERN_REVERSE_ANCHORS },
];
pub const CSPCI: Enzyme = Enzyme {
    name: "CspCI",
    id: 1,
    tag_length: 33,
    rule: MatchRule::Exact(&CSPCI_PATTERNS),
};

// 2. AloI (tag_length=27): [AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}
const ALOI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GAAC" },
    Anchor { offset: 17, motif: b"TCC" },
];
const ALOI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GGA" },
    Anchor { offset: 16, motif: b"GTTC" },
];
const ALOI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &ALOI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &ALOI_PATTERN_REVERSE_ANCHORS },
];
pub const ALOI: Enzyme = Enzyme {
    name: "AloI",
    id: 2,
    tag_length: 27,
    rule: MatchRule::Exact(&ALOI_PATTERNS),
};

// 3. BsaXI (tag_length=27): [AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}
const BSAXI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 9, motif: b"AC" },
    Anchor { offset: 16, motif: b"CTCC" },
];
const BSAXI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GGAG" },
    Anchor { offset: 16, motif: b"GT" },
];
const BSAXI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &BSAXI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &BSAXI_PATTERN_REVERSE_ANCHORS },
];
pub const BSAXI: Enzyme = Enzyme {
    name: "BsaXI",
    id: 3,
    tag_length: 27,
    rule: MatchRule::Exact(&BSAXI_PATTERNS),
};

// 4. BaeI (tag_length=28, degenerate): recognition site contains IUPAC
// pyrimidine (Y=[CT]) and purine (R=[AG]) positions, so it is matched with
// an unanchored regex instead of fixed-byte anchors.
//   fwd: [AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}
//   rev: [AGCT]{7}G[AG]TAC[AGCT]{4}GT[AGCT]{10}
static BAEI_REGEXES: OnceLock<Vec<Regex>> = OnceLock::new();
fn baei_regexes() -> &'static [Regex] {
    BAEI_REGEXES.get_or_init(|| {
        vec![
            Regex::new(r"[ACGT]{10}AC[ACGT]{4}GTA[CT]C[ACGT]{7}")
                .expect("valid BaeI forward regex"),
            Regex::new(r"[ACGT]{7}G[AG]TAC[ACGT]{4}GT[ACGT]{10}")
                .expect("valid BaeI reverse regex"),
        ]
    })
}
pub const BAEI: Enzyme = Enzyme {
    name: "BaeI",
    id: 4,
    tag_length: 28,
    rule: MatchRule::Degenerate(baei_regexes),
};

// 5. BcgI (tag_length=32): [AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}
const BCGI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 10, motif: b"CGA" },
    Anchor { offset: 19, motif: b"TGC" },
];
const BCGI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 10, motif: b"GCA" },
    Anchor { offset: 19, motif: b"TCG" },
];
const BCGI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &BCGI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &BCGI_PATTERN_REVERSE_ANCHORS },
];
pub const BCGI: Enzyme = Enzyme {
    name: "BcgI",
    id: 5,
    tag_length: 32,
    rule: MatchRule::Exact(&BCGI_PATTERNS),
};

// 6. CjeI (tag_length=28): [AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}
const CJEI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"CCA" },
    Anchor { offset: 17, motif: b"GT" },
];
const CJEI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 9, motif: b"AC" },
    Anchor { offset: 17, motif: b"TGG" },
];
const CJEI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &CJEI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &CJEI_PATTERN_REVERSE_ANCHORS },
];
pub const CJEI: Enzyme = Enzyme {
    name: "CjeI",
    id: 6,
    tag_length: 28,
    rule: MatchRule::Exact(&CJEI_PATTERNS),
};

// 7. PpiI (tag_length=27): [AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}
// NOTE: forward anchor2 offset corrected from 17 -> 16 (the previous value
// implied a 6bp spacer; the Perl pattern requires exactly 5bp).
const PPII_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GAAC" },
    Anchor { offset: 16, motif: b"CTC" },
];
const PPII_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"GAG" },
    Anchor { offset: 16, motif: b"GTTC" },
];
const PPII_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &PPII_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &PPII_PATTERN_REVERSE_ANCHORS },
];
pub const PPII: Enzyme = Enzyme {
    name: "PpiI",
    id: 7,
    tag_length: 27,
    rule: MatchRule::Exact(&PPII_PATTERNS),
};

// 8. PsrI (tag_length=27): [AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}
const PSRI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GAAC" },
    Anchor { offset: 17, motif: b"TAC" },
];
const PSRI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"GTA" },
    Anchor { offset: 16, motif: b"GTTC" },
];
const PSRI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &PSRI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &PSRI_PATTERN_REVERSE_ANCHORS },
];
pub const PSRI: Enzyme = Enzyme {
    name: "PsrI",
    id: 8,
    tag_length: 27,
    rule: MatchRule::Exact(&PSRI_PATTERNS),
};

// 9. BplI (tag_length=27, palindrome): [AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}
const BPLI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"GAG" },
    Anchor { offset: 16, motif: b"CTC" },
];
const BPLI_PATTERNS: [Pattern; 1] = [Pattern { anchors: &BPLI_PATTERN_ANCHORS }];
pub const BPLI: Enzyme = Enzyme {
    name: "BplI",
    id: 9,
    tag_length: 27,
    rule: MatchRule::Exact(&BPLI_PATTERNS),
};

// 10. FalI (tag_length=27, palindrome): [AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}
const FALI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"AAG" },
    Anchor { offset: 16, motif: b"CTT" },
];
const FALI_PATTERNS: [Pattern; 1] = [Pattern { anchors: &FALI_PATTERN_ANCHORS }];
pub const FALI: Enzyme = Enzyme {
    name: "FalI",
    id: 10,
    tag_length: 27,
    rule: MatchRule::Exact(&FALI_PATTERNS),
};

// 11. Bsp24I (tag_length=27): [AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}
const BSP24I_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"GAC" },
    Anchor { offset: 17, motif: b"TGG" },
];
const BSP24I_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"CCA" },
    Anchor { offset: 16, motif: b"GTC" },
];
const BSP24I_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &BSP24I_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &BSP24I_PATTERN_REVERSE_ANCHORS },
];
pub const BSP24I: Enzyme = Enzyme {
    name: "Bsp24I",
    id: 11,
    tag_length: 27,
    rule: MatchRule::Exact(&BSP24I_PATTERNS),
};

// 12. HaeIV (tag_length=27, degenerate):
//   fwd: [AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}
//   rev: [AGCT]{9}GA[CT][AGCT]{5}[AG]TC[AGCT]{7}
static HAEIV_REGEXES: OnceLock<Vec<Regex>> = OnceLock::new();
fn haeiv_regexes() -> &'static [Regex] {
    HAEIV_REGEXES.get_or_init(|| {
        vec![
            Regex::new(r"[ACGT]{7}GA[CT][ACGT]{5}[AG]TC[ACGT]{9}")
                .expect("valid HaeIV forward regex"),
            Regex::new(r"[ACGT]{9}GA[CT][ACGT]{5}[AG]TC[ACGT]{7}")
                .expect("valid HaeIV reverse regex"),
        ]
    })
}
pub const HAEIV: Enzyme = Enzyme {
    name: "HaeIV",
    id: 12,
    tag_length: 27,
    rule: MatchRule::Degenerate(haeiv_regexes),
};

// 13. CjePI (tag_length=27): [AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}
const CJEPI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 7, motif: b"CCA" },
    Anchor { offset: 17, motif: b"TC" },
];
const CJEPI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 8, motif: b"GA" },
    Anchor { offset: 17, motif: b"TGG" },
];
const CJEPI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &CJEPI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &CJEPI_PATTERN_REVERSE_ANCHORS },
];
pub const CJEPI: Enzyme = Enzyme {
    name: "CjePI",
    id: 13,
    tag_length: 27,
    rule: MatchRule::Exact(&CJEPI_PATTERNS),
};

// 14. Hin4I (tag_length=27, degenerate):
//   fwd: [AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}
//   rev: [AGCT]{8}GA[CTG][AGCT]{5}[AG]TC[AGCT]{8}
static HIN4I_REGEXES: OnceLock<Vec<Regex>> = OnceLock::new();
fn hin4i_regexes() -> &'static [Regex] {
    HIN4I_REGEXES.get_or_init(|| {
        vec![
            Regex::new(r"[ACGT]{8}GA[CT][ACGT]{5}[GAC]TC[ACGT]{8}")
                .expect("valid Hin4I forward regex"),
            Regex::new(r"[ACGT]{8}GA[CTG][ACGT]{5}[AG]TC[ACGT]{8}")
                .expect("valid Hin4I reverse regex"),
        ]
    })
}
pub const HIN4I: Enzyme = Enzyme {
    name: "Hin4I",
    id: 14,
    tag_length: 27,
    rule: MatchRule::Degenerate(hin4i_regexes),
};

// 15. AlfI (tag_length=32, palindrome): [AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}
const ALFI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor { offset: 10, motif: b"GCA" },
    Anchor { offset: 19, motif: b"TGC" },
];
const ALFI_PATTERNS: [Pattern; 1] = [Pattern { anchors: &ALFI_PATTERN_ANCHORS }];
pub const ALFI: Enzyme = Enzyme {
    name: "AlfI",
    id: 15,
    tag_length: 32,
    rule: MatchRule::Exact(&ALFI_PATTERNS),
};

// 16. BslFI (tag_length=25): [AGCT]{6}GGGAC[AGCT]{14}
const BSLFI_PATTERN_FORWARD_ANCHORS: [Anchor; 1] = [Anchor { offset: 6, motif: b"GGGAC" }];
const BSLFI_PATTERN_REVERSE_ANCHORS: [Anchor; 1] = [Anchor { offset: 14, motif: b"GTCCC" }];
const BSLFI_PATTERNS: [Pattern; 2] = [
    Pattern { anchors: &BSLFI_PATTERN_FORWARD_ANCHORS },
    Pattern { anchors: &BSLFI_PATTERN_REVERSE_ANCHORS },
];
pub const BSLFI: Enzyme = Enzyme {
    name: "BslFI",
    id: 16,
    tag_length: 25,
    rule: MatchRule::Exact(&BSLFI_PATTERNS),
};

// ========== Enzyme lookup functions ==========

static ENZYMES: &[&Enzyme] = &[
    &CSPCI, &ALOI, &BSAXI, &BAEI, &BCGI, &CJEI, &PPII, &PSRI, &BPLI, &FALI, &BSP24I, &HAEIV,
    &CJEPI, &HIN4I, &ALFI, &BSLFI,
];

pub fn enzyme_by_name(name: &str) -> Option<&'static Enzyme> {
    ENZYMES
        .iter()
        .copied()
        .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
}

pub fn enzyme_by_id(id: u8) -> Option<&'static Enzyme> {
    ENZYMES.iter().copied().find(|enzyme| enzyme.id == id)
}
