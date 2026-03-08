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

/// Input data type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputType {
    /// Type 1: Reference genome FASTA
    ReferenceGenome = 1,
    /// Type 2: Shotgun sequencing data (SE/PE)
    ShotgunMetagenome = 2,
    /// Type 3: Single 2bRAD tag (SE/PE, take only the first match)
    Single2bRAD = 3,
}

impl InputType {
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            1 => Some(InputType::ReferenceGenome),
            2 => Some(InputType::ShotgunMetagenome),
            3 => Some(InputType::Single2bRAD),
            _ => None,
        }
    }
}

/// Quality control configuration
#[derive(Debug, Clone)]
pub struct QualityControl {
    /// Whether quality control is enabled
    pub enabled: bool,
    /// Maximum N ratio (0.0–1.0 for fraction, >=1.0 for absolute count)
    pub max_n: f64,
    /// Minimum quality score
    pub min_quality: u8,
    /// Minimum quality percentage (0–100)
    pub min_quality_percent: u8,
    /// Quality score encoding base (typically 33 or 64)
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
    /// Check whether the N ratio in a sequence satisfies the threshold
    pub fn check_n(&self, sequence: &[u8]) -> bool {
        if !self.enabled {
            return true;
        }

        let n_count = sequence.iter().filter(|&&b| b == b'N').count();

        if self.max_n > 0.0 && self.max_n < 1.0 {
            // Fraction mode
            let ratio = n_count as f64 / sequence.len() as f64;
            ratio <= self.max_n
        } else {
            // Absolute count mode
            (n_count as f64) <= self.max_n
        }
    }

    /// Check whether the quality scores satisfy the threshold
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
