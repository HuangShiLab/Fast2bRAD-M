use fxhash::FxHashSet;

#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub offset: usize,
    pub motif: &'static [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct Pattern {
    pub anchors: &'static [Anchor],
}

#[derive(Debug, Clone, Copy)]
pub struct Enzyme {
    pub name: &'static str,
    pub id: u8,
    pub tag_length: usize,
    pub patterns: &'static [Pattern],
    // Type 4 (5连标签) 的位置参数（暂未使用）
    #[allow(dead_code)]
    pub concat_starts: &'static [usize],
    #[allow(dead_code)]
    pub concat_ends: &'static [usize],
    #[allow(dead_code)]
    pub min_pear: Option<usize>,
    #[allow(dead_code)]
    pub max_pear: Option<usize>,
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

// 创建一个静态查找表，大小为 256 (u8 的范围)
// 只有 ATCGatcg 的位置是 true，其他（包括 N）都是 false
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
    // 这种查表法是 O(1) 的指令，通常比分支判断更快
    window.iter().all(|&b| ATCG_TABLE[b as usize])
}

impl Enzyme {

    /// 在序列中查找所有匹配的标签位置和长度（去重）
    pub fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let mut positions = FxHashSet::default();

        for pattern in self.patterns {
            for offset in 0..sequence.len() {
                if offset + self.tag_length > sequence.len() {
                    break;
                }
                let window = &sequence[offset..offset + self.tag_length];
                if pattern.matches(window) {
                    if is_pure_atcg(window) {
                        positions.insert((offset, self.tag_length));
                    }
                }
            }
        }

        let mut result: Vec<_> = positions.into_iter().collect();
        result.sort_unstable();
        result
    }
}

// ========== 16 种酶的定义 ==========

// 1. CspCI (tag_length=36)
const CSPCI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 11,
        motif: b"CAA",
    },
    Anchor {
        offset: 19,
        motif: b"GTGG",
    },
];
const CSPCI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 10,
        motif: b"CCAC",
    },
    Anchor {
        offset: 19,
        motif: b"TTG",
    },
];
const CSPCI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &CSPCI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &CSPCI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const CSPCI: Enzyme = Enzyme {
    name: "CspCI",
    id: 1,
    tag_length: 36,
    patterns: &CSPCI_PATTERNS,
    concat_starts: &[0, 37, 78, 119, 160],
    concat_ends: &[41, 82, 123, 164, 205],
    min_pear: None,
    max_pear: None,
};

// 2. AloI (tag_length=37)
const ALOI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GAAC",
    },
    Anchor {
        offset: 17,
        motif: b"TCC",
    },
];
const ALOI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GGA",
    },
    Anchor {
        offset: 16,
        motif: b"GTTC",
    },
];
const ALOI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &ALOI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &ALOI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const ALOI: Enzyme = Enzyme {
    name: "AloI",
    id: 2,
    tag_length: 37,
    patterns: &ALOI_PATTERNS,
    concat_starts: &[0, 38, 80, 122, 164],
    concat_ends: &[42, 84, 126, 168, 210],
    min_pear: None,
    max_pear: None,
};

// 3. BsaXI (tag_length=32)
const BSAXI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 9,
        motif: b"AC",
    },
    Anchor {
        offset: 16,
        motif: b"CTCC",
    },
];
const BSAXI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GGAG",
    },
    Anchor {
        offset: 16,
        motif: b"GT",
    },
];
const BSAXI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &BSAXI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &BSAXI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const BSAXI: Enzyme = Enzyme {
    name: "BsaXI",
    id: 3,
    tag_length: 32,
    patterns: &BSAXI_PATTERNS,
    concat_starts: &[0, 33, 69, 105, 141],
    concat_ends: &[35, 71, 107, 143, 180],
    min_pear: Some(173),
    max_pear: Some(181),
};

// 4. BaeI (tag_length=36)
const BAEI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 10,
        motif: b"AC",
    },
    Anchor {
        offset: 16,
        motif: b"GTA",
    },
];
const BAEI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"G",
    },
    Anchor {
        offset: 9,
        motif: b"TAC",
    },
];
const BAEI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &BAEI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &BAEI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const BAEI: Enzyme = Enzyme {
    name: "BaeI",
    id: 4,
    tag_length: 36,
    patterns: &BAEI_PATTERNS,
    concat_starts: &[0, 38, 79, 120, 161],
    concat_ends: &[40, 81, 122, 163, 205],
    min_pear: Some(198),
    max_pear: Some(206),
};

// 5. BcgI (tag_length=32)
const BCGI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 10,
        motif: b"CGA",
    },
    Anchor {
        offset: 19,
        motif: b"TGC",
    },
];
const BCGI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 10,
        motif: b"GCA",
    },
    Anchor {
        offset: 19,
        motif: b"TCG",
    },
];
const BCGI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &BCGI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &BCGI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const BCGI: Enzyme = Enzyme {
    name: "BcgI",
    id: 5,
    tag_length: 32,
    patterns: &BCGI_PATTERNS,
    concat_starts: &[0, 36, 75, 114, 153],
    concat_ends: &[38, 77, 116, 155, 195],
    min_pear: Some(188),
    max_pear: Some(196),
};

// 6. CjeI (tag_length=37)
const CJEI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"CCA",
    },
    Anchor {
        offset: 17,
        motif: b"GT",
    },
];
const CJEI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 9,
        motif: b"AC",
    },
    Anchor {
        offset: 17,
        motif: b"TGG",
    },
];
const CJEI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &CJEI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &CJEI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const CJEI: Enzyme = Enzyme {
    name: "CjeI",
    id: 6,
    tag_length: 37,
    patterns: &CJEI_PATTERNS,
    concat_starts: &[0, 40, 83, 126, 169],
    concat_ends: &[42, 85, 128, 171, 214],
    min_pear: None,
    max_pear: None,
};

// 7. PpiI (tag_length=35)
const PPII_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GAAC",
    },
    Anchor {
        offset: 17,
        motif: b"CTC",
    },
];
const PPII_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"GAG",
    },
    Anchor {
        offset: 16,
        motif: b"GTTC",
    },
];
const PPII_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &PPII_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &PPII_PATTERN_REVERSE_ANCHORS,
    },
];
pub const PPII: Enzyme = Enzyme {
    name: "PpiI",
    id: 7,
    tag_length: 35,
    patterns: &PPII_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 199],
    min_pear: None,
    max_pear: None,
};

// 8. PsrI (tag_length=35)
const PSRI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GAAC",
    },
    Anchor {
        offset: 17,
        motif: b"TAC",
    },
];
const PSRI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"GTA",
    },
    Anchor {
        offset: 16,
        motif: b"GTTC",
    },
];
const PSRI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &PSRI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &PSRI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const PSRI: Enzyme = Enzyme {
    name: "PsrI",
    id: 8,
    tag_length: 35,
    patterns: &PSRI_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 199],
    min_pear: None,
    max_pear: None,
};

// 9. BplI (tag_length=35, palindrome)
const BPLI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"GAG",
    },
    Anchor {
        offset: 16,
        motif: b"CTC",
    },
];
const BPLI_PATTERNS: [Pattern; 1] = [Pattern {
    anchors: &BPLI_PATTERN_ANCHORS,
}];
pub const BPLI: Enzyme = Enzyme {
    name: "BplI",
    id: 9,
    tag_length: 35,
    patterns: &BPLI_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 199],
    min_pear: None,
    max_pear: None,
};

// 10. FalI (tag_length=36, palindrome)
const FALI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"AAG",
    },
    Anchor {
        offset: 16,
        motif: b"CTT",
    },
];
const FALI_PATTERNS: [Pattern; 1] = [Pattern {
    anchors: &FALI_PATTERN_ANCHORS,
}];
pub const FALI: Enzyme = Enzyme {
    name: "FalI",
    id: 10,
    tag_length: 36,
    patterns: &FALI_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 200],
    min_pear: Some(193),
    max_pear: Some(201),
};

// 11. Bsp24I (tag_length=36)
const BSP24I_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"GAC",
    },
    Anchor {
        offset: 17,
        motif: b"TGG",
    },
];
const BSP24I_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"CCA",
    },
    Anchor {
        offset: 16,
        motif: b"GTC",
    },
];
const BSP24I_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &BSP24I_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &BSP24I_PATTERN_REVERSE_ANCHORS,
    },
];
pub const BSP24I: Enzyme = Enzyme {
    name: "Bsp24I",
    id: 11,
    tag_length: 36,
    patterns: &BSP24I_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 200],
    min_pear: None,
    max_pear: None,
};

// 12. HaeIV (tag_length=37)
const HAEIV_PATTERN_FORWARD_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 7,
    motif: b"GA",
}];
const HAEIV_PATTERN_REVERSE_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 9,
    motif: b"GA",
}];
const HAEIV_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &HAEIV_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &HAEIV_PATTERN_REVERSE_ANCHORS,
    },
];
pub const HAEIV: Enzyme = Enzyme {
    name: "HaeIV",
    id: 12,
    tag_length: 37,
    patterns: &HAEIV_PATTERNS,
    concat_starts: &[0, 38, 79, 120, 161],
    concat_ends: &[40, 81, 122, 163, 204],
    min_pear: None,
    max_pear: None,
};

// 13. CjePI (tag_length=38)
const CJEPI_PATTERN_FORWARD_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 7,
        motif: b"CCA",
    },
    Anchor {
        offset: 17,
        motif: b"TC",
    },
];
const CJEPI_PATTERN_REVERSE_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 8,
        motif: b"GA",
    },
    Anchor {
        offset: 17,
        motif: b"TGG",
    },
];
const CJEPI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &CJEPI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &CJEPI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const CJEPI: Enzyme = Enzyme {
    name: "CjePI",
    id: 13,
    tag_length: 38,
    patterns: &CJEPI_PATTERNS,
    concat_starts: &[0, 39, 81, 123, 165],
    concat_ends: &[41, 83, 125, 167, 209],
    min_pear: None,
    max_pear: None,
};

// 14. Hin4I (tag_length=35)
const HIN4I_PATTERN_FORWARD_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 8,
    motif: b"GA",
}];
const HIN4I_PATTERN_REVERSE_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 8,
    motif: b"GA",
}];
const HIN4I_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &HIN4I_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &HIN4I_PATTERN_REVERSE_ANCHORS,
    },
];
pub const HIN4I: Enzyme = Enzyme {
    name: "Hin4I",
    id: 14,
    tag_length: 35,
    patterns: &HIN4I_PATTERNS,
    concat_starts: &[0, 37, 77, 117, 157],
    concat_ends: &[39, 79, 119, 159, 199],
    min_pear: None,
    max_pear: None,
};

// 15. AlfI (tag_length=33, palindrome)
const ALFI_PATTERN_ANCHORS: [Anchor; 2] = [
    Anchor {
        offset: 10,
        motif: b"GCA",
    },
    Anchor {
        offset: 19,
        motif: b"TGC",
    },
];
const ALFI_PATTERNS: [Pattern; 1] = [Pattern {
    anchors: &ALFI_PATTERN_ANCHORS,
}];
pub const ALFI: Enzyme = Enzyme {
    name: "AlfI",
    id: 15,
    tag_length: 33,
    patterns: &ALFI_PATTERNS,
    concat_starts: &[0, 36, 75, 114, 153],
    concat_ends: &[38, 77, 116, 155, 194],
    min_pear: None,
    max_pear: None,
};

// 16. BslFI (tag_length=33)
const BSLFI_PATTERN_FORWARD_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 6,
    motif: b"GGGAC",
}];
const BSLFI_PATTERN_REVERSE_ANCHORS: [Anchor; 1] = [Anchor {
    offset: 14,
    motif: b"GTCCC",
}];
const BSLFI_PATTERNS: [Pattern; 2] = [
    Pattern {
        anchors: &BSLFI_PATTERN_FORWARD_ANCHORS,
    },
    Pattern {
        anchors: &BSLFI_PATTERN_REVERSE_ANCHORS,
    },
];
pub const BSLFI: Enzyme = Enzyme {
    name: "BslFI",
    id: 16,
    tag_length: 33,
    patterns: &BSLFI_PATTERNS,
    concat_starts: &[0, 34, 72, 110, 148],
    concat_ends: &[38, 76, 114, 152, 190],
    min_pear: None,
    max_pear: None,
};

// ========== 酶查找函数 ==========

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
