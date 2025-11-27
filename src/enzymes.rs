use fxhash::FxHashSet;
use regex::Regex;

// ========== 备份：锚点匹配实现（保留以备后用） ==========

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Anchor {
    pub offset: usize,
    pub motif: &'static [u8],
}

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct AnchorPattern {
    pub anchors: &'static [Anchor],
}

impl AnchorPattern {
    #[allow(dead_code)]
    pub fn matches(&self, window: &[u8]) -> bool {
        self.anchors.iter().all(|anchor| {
            let start = anchor.offset;
            let end = start + anchor.motif.len();
            end <= window.len() && &window[start..end] == anchor.motif
        })
    }
}

// ========== 正则表达式匹配实现（与Perl版本一致） ==========

/// 正则表达式模式（编译后的Regex）
pub struct RegexPattern {
    pub regex: Regex,
}

impl RegexPattern {
    fn new(pattern: &str) -> Result<Self, regex::Error> {
        Ok(RegexPattern {
            regex: Regex::new(pattern)?,
        })
    }

    /// 在序列中查找所有匹配的位置和长度
    fn find_all(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let seq_str = match std::str::from_utf8(sequence) {
            Ok(s) => s,
            Err(_) => return Vec::new(), // 无效UTF-8，跳过
        };

        let mut positions = Vec::new();
        for mat in self.regex.find_iter(seq_str) {
            positions.push((mat.start(), mat.end() - mat.start()));
        }
        positions
    }
}

/// 预编译的正则表达式模式
struct CompiledRegexPatterns {
    patterns: Vec<RegexPattern>,
}

impl CompiledRegexPatterns {
    fn new(patterns: &[&'static str]) -> Result<Self, regex::Error> {
        let compiled: Result<Vec<_>, _> = patterns.iter().map(|p| RegexPattern::new(p)).collect();
        Ok(CompiledRegexPatterns {
            patterns: compiled?,
        })
    }

    fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let mut all_positions = FxHashSet::default();
        for pattern in &self.patterns {
            for (pos, len) in pattern.find_all(sequence) {
                all_positions.insert((pos, len));
            }
        }
        let mut result: Vec<_> = all_positions.into_iter().collect();
        result.sort_unstable();
        result
    }
}

// 为CompiledRegexPatterns添加访问patterns的方法
impl CompiledRegexPatterns {
    fn patterns(&self) -> &[RegexPattern] {
        &self.patterns
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Enzyme {
    pub name: &'static str,
    pub id: u8,
    pub tag_length: usize,
    /// 备份：锚点模式（保留以备后用）
    #[allow(dead_code)]
    pub patterns: &'static [AnchorPattern],
    /// 正则表达式模式字符串（与Perl版本一致）
    pub regex_patterns: &'static [&'static str],
    // Type 4 (5连标签) 的位置参数
    #[allow(dead_code)]
    pub concat_starts: &'static [usize],
    #[allow(dead_code)]
    pub concat_ends: &'static [usize],
    #[allow(dead_code)]
    pub min_pear: Option<usize>,
    #[allow(dead_code)]
    pub max_pear: Option<usize>,
}

// 使用线程局部存储缓存编译后的正则表达式，避免重复编译
thread_local! {
    static PATTERN_CACHE: std::cell::RefCell<fxhash::FxHashMap<&'static str, CompiledRegexPatterns>> = 
        std::cell::RefCell::new(fxhash::FxHashMap::default());
}

impl Enzyme {
    /// 在序列中查找所有匹配的标签位置和长度（去重）
    /// 使用正则表达式匹配，与Perl版本一致，确保[AGCT]不匹配N
    pub fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        PATTERN_CACHE.with(|cache| {
            let mut cache = cache.borrow_mut();
            let compiled = cache.entry(self.name).or_insert_with(|| {
                CompiledRegexPatterns::new(self.regex_patterns)
                    .unwrap_or_else(|e| {
                        panic!("Failed to compile regex patterns for {}: {}", self.name, e)
                    })
            });
            compiled.find_all_tags(sequence)
        })
    }

    /// 在序列中查找第一个匹配的标签位置和长度（用于Type 3和Type 4）
    /// 使用正则表达式匹配，与Perl版本一致
    pub fn find_first_tag(&self, sequence: &[u8]) -> Option<(usize, usize)> {
        let seq_str = match std::str::from_utf8(sequence) {
            Ok(s) => s,
            Err(_) => return None,
        };

        PATTERN_CACHE.with(|cache| {
            let mut cache = cache.borrow_mut();
            let compiled = cache.entry(self.name).or_insert_with(|| {
                CompiledRegexPatterns::new(self.regex_patterns)
                    .unwrap_or_else(|e| {
                        panic!("Failed to compile regex patterns for {}: {}", self.name, e)
                    })
            });

            // 按顺序检查每个正则表达式模式，返回第一个匹配
            for pattern in compiled.patterns() {
                if let Some(mat) = pattern.regex.find(seq_str) {
                    return Some((mat.start(), mat.end() - mat.start()));
                }
            }
            None
        })
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
const CSPCI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &CSPCI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &CSPCI_PATTERN_REVERSE_ANCHORS,
    },
];
const CSPCI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}",
    r"[AGCT]{10}CCAC[AGCT]{5}TTG[AGCT]{11}",
];
pub const CSPCI: Enzyme = Enzyme {
    name: "CspCI",
    id: 1,
    tag_length: 36,
    patterns: &CSPCI_PATTERNS,
    regex_patterns: CSPCI_REGEX_PATTERNS,
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
const ALOI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &ALOI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &ALOI_PATTERN_REVERSE_ANCHORS,
    },
];
const ALOI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}",
    r"[AGCT]{7}GGA[AGCT]{6}GTTC[AGCT]{7}",
];
pub const ALOI: Enzyme = Enzyme {
    name: "AloI",
    id: 2,
    tag_length: 37,
    patterns: &ALOI_PATTERNS,
    regex_patterns: ALOI_REGEX_PATTERNS,
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
const BSAXI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &BSAXI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &BSAXI_PATTERN_REVERSE_ANCHORS,
    },
];
const BSAXI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}",
    r"[AGCT]{7}GGAG[AGCT]{5}GT[AGCT]{9}",
];
pub const BSAXI: Enzyme = Enzyme {
    name: "BsaXI",
    id: 3,
    tag_length: 32,
    patterns: &BSAXI_PATTERNS,
    regex_patterns: BSAXI_REGEX_PATTERNS,
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
const BAEI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &BAEI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &BAEI_PATTERN_REVERSE_ANCHORS,
    },
];
const BAEI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}",
    r"[AGCT]{7}G[AG]TAC[AGCT]{4}GT[AGCT]{10}",
];
pub const BAEI: Enzyme = Enzyme {
    name: "BaeI",
    id: 4,
    tag_length: 36,
    patterns: &BAEI_PATTERNS,
    regex_patterns: BAEI_REGEX_PATTERNS,
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
const BCGI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &BCGI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &BCGI_PATTERN_REVERSE_ANCHORS,
    },
];
const BCGI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}",
    r"[AGCT]{10}GCA[AGCT]{6}TCG[AGCT]{10}",
];
pub const BCGI: Enzyme = Enzyme {
    name: "BcgI",
    id: 5,
    tag_length: 32,
    patterns: &BCGI_PATTERNS,
    regex_patterns: BCGI_REGEX_PATTERNS,
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
const CJEI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &CJEI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &CJEI_PATTERN_REVERSE_ANCHORS,
    },
];
const CJEI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}",
    r"[AGCT]{9}AC[AGCT]{6}TGG[AGCT]{8}",
];
pub const CJEI: Enzyme = Enzyme {
    name: "CjeI",
    id: 6,
    tag_length: 37,
    patterns: &CJEI_PATTERNS,
    regex_patterns: CJEI_REGEX_PATTERNS,
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
const PPII_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &PPII_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &PPII_PATTERN_REVERSE_ANCHORS,
    },
];
const PPII_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}",
    r"[AGCT]{8}GAG[AGCT]{5}GTTC[AGCT]{7}",
];
pub const PPII: Enzyme = Enzyme {
    name: "PpiI",
    id: 7,
    tag_length: 35,
    patterns: &PPII_PATTERNS,
    regex_patterns: PPII_REGEX_PATTERNS,
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
const PSRI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &PSRI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &PSRI_PATTERN_REVERSE_ANCHORS,
    },
];
const PSRI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}",
    r"[AGCT]{7}GTA[AGCT]{6}GTTC[AGCT]{7}",
];
pub const PSRI: Enzyme = Enzyme {
    name: "PsrI",
    id: 8,
    tag_length: 35,
    patterns: &PSRI_PATTERNS,
    regex_patterns: PSRI_REGEX_PATTERNS,
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
const BPLI_PATTERNS: [AnchorPattern; 1] = [AnchorPattern {
    anchors: &BPLI_PATTERN_ANCHORS,
}];
const BPLI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}",
];
pub const BPLI: Enzyme = Enzyme {
    name: "BplI",
    id: 9,
    tag_length: 35,
    patterns: &BPLI_PATTERNS,
    regex_patterns: BPLI_REGEX_PATTERNS,
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
const FALI_PATTERNS: [AnchorPattern; 1] = [AnchorPattern {
    anchors: &FALI_PATTERN_ANCHORS,
}];
const FALI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}",
];
pub const FALI: Enzyme = Enzyme {
    name: "FalI",
    id: 10,
    tag_length: 36,
    patterns: &FALI_PATTERNS,
    regex_patterns: FALI_REGEX_PATTERNS,
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
const BSP24I_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &BSP24I_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &BSP24I_PATTERN_REVERSE_ANCHORS,
    },
];
const BSP24I_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}",
    r"[AGCT]{7}CCA[AGCT]{6}GTC[AGCT]{8}",
];
pub const BSP24I: Enzyme = Enzyme {
    name: "Bsp24I",
    id: 11,
    tag_length: 36,
    patterns: &BSP24I_PATTERNS,
    regex_patterns: BSP24I_REGEX_PATTERNS,
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
const HAEIV_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &HAEIV_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &HAEIV_PATTERN_REVERSE_ANCHORS,
    },
];
const HAEIV_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}",
    r"[AGCT]{9}GA[CT][AGCT]{5}[AG]TC[AGCT]{7}",
];
pub const HAEIV: Enzyme = Enzyme {
    name: "HaeIV",
    id: 12,
    tag_length: 37,
    patterns: &HAEIV_PATTERNS,
    regex_patterns: HAEIV_REGEX_PATTERNS,
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
const CJEPI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &CJEPI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &CJEPI_PATTERN_REVERSE_ANCHORS,
    },
];
const CJEPI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}",
    r"[AGCT]{8}GA[AGCT]{7}TGG[AGCT]{7}",
];
pub const CJEPI: Enzyme = Enzyme {
    name: "CjePI",
    id: 13,
    tag_length: 38,
    patterns: &CJEPI_PATTERNS,
    regex_patterns: CJEPI_REGEX_PATTERNS,
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
const HIN4I_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &HIN4I_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &HIN4I_PATTERN_REVERSE_ANCHORS,
    },
];
const HIN4I_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}",
    r"[AGCT]{8}GA[CTG][AGCT]{5}[AG]TC[AGCT]{8}",
];
pub const HIN4I: Enzyme = Enzyme {
    name: "Hin4I",
    id: 14,
    tag_length: 35,
    patterns: &HIN4I_PATTERNS,
    regex_patterns: HIN4I_REGEX_PATTERNS,
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
const ALFI_PATTERNS: [AnchorPattern; 1] = [AnchorPattern {
    anchors: &ALFI_PATTERN_ANCHORS,
}];
const ALFI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}",
];
pub const ALFI: Enzyme = Enzyme {
    name: "AlfI",
    id: 15,
    tag_length: 33,
    patterns: &ALFI_PATTERNS,
    regex_patterns: ALFI_REGEX_PATTERNS,
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
const BSLFI_PATTERNS: [AnchorPattern; 2] = [
    AnchorPattern {
        anchors: &BSLFI_PATTERN_FORWARD_ANCHORS,
    },
    AnchorPattern {
        anchors: &BSLFI_PATTERN_REVERSE_ANCHORS,
    },
];
const BSLFI_REGEX_PATTERNS: &[&str] = &[
    r"[AGCT]{6}GGGAC[AGCT]{14}",
    r"[AGCT]{14}GTCCC[AGCT]{6}",
];
pub const BSLFI: Enzyme = Enzyme {
    name: "BslFI",
    id: 16,
    tag_length: 33,
    patterns: &BSLFI_PATTERNS,
    regex_patterns: BSLFI_REGEX_PATTERNS,
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

pub fn supported_enzyme_names() -> Vec<&'static str> {
    ENZYMES.iter().map(|e| e.name).collect()
}
