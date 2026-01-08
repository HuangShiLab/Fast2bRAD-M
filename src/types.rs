#[derive(Debug, Clone)]
pub struct DigestStats {
    pub sample_id: String,
    pub enzyme: String,
    pub input_sequences: usize,
    pub tag_count: usize,
}

impl DigestStats {
    pub fn percent(&self) -> f64 {
        if self.input_sequences == 0 {
            0.0
        } else {
            (self.tag_count as f64 / self.input_sequences as f64) * 100.0
        }
    }
}

/// 输入数据类型
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputType {
    /// Type 1: 参考基因组 FASTA
    ReferenceGenome = 1,
    /// Type 2: Shotgun 测序数据（SE/PE）
    ShotgunMetagenome = 2,
    /// Type 3: 单条 2bRAD 标签（SE/PE，只取第一个匹配）
    Single2bRAD = 3,
    /// Type 4: 5 个连接的 2bRAD 标签（PE，按固定位置切分）
    Concatenated2bRAD = 4,
}

impl InputType {
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            1 => Some(InputType::ReferenceGenome),
            2 => Some(InputType::ShotgunMetagenome),
            3 => Some(InputType::Single2bRAD),
            4 => Some(InputType::Concatenated2bRAD),
            _ => None,
        }
    }
}

/// 质量控制配置
#[derive(Debug, Clone)]
pub struct QualityControl {
    /// 是否启用质量控制
    pub enabled: bool,
    /// 最大 N 比例（0.0-1.0 表示百分比，>=1.0 表示绝对数量）
    pub max_n: f64,
    /// 最低质量分数
    pub min_quality: u8,
    /// 最低质量百分比（0-100）
    pub min_quality_percent: u8,
    /// 质量分数编码类型（通常是 33 或 64）
    pub quality_base: u8,
}

impl Default for QualityControl {
    fn default() -> Self {
        Self {
            enabled: true,
            max_n: 0.08,
            min_quality: 30,
            min_quality_percent: 80,
            quality_base: 33,
        }
    }
}

impl QualityControl {
    /// 检查序列中 N 的比例是否满足要求
    pub fn check_n(&self, sequence: &[u8]) -> bool {
        if !self.enabled {
            return true;
        }

        let n_count = sequence.iter().filter(|&&b| b == b'N').count();

        if self.max_n > 0.0 && self.max_n < 1.0 {
            // 百分比模式
            let ratio = n_count as f64 / sequence.len() as f64;
            ratio <= self.max_n
        } else {
            // 绝对数量模式
            (n_count as f64) <= self.max_n
        }
    }

    /// 检查质量分数是否满足要求
    pub fn check_quality(&self, quality: &[u8]) -> bool {
        if !self.enabled {
            return true;
        }

        let min_phred = self.min_quality + self.quality_base;
        let passed_count = quality.iter().filter(|&&q| q >= min_phred).count();
        let passed_percent = (passed_count * 100) / quality.len();

        passed_percent >= self.min_quality_percent as usize
    }
}
