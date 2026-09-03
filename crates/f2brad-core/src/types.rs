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
    /// Check whether the N ratio in a sequence satisfies the threshold.
    /// Accepts both cases: sequences reach here straight from the reader in
    /// some paths, and a lower-case `n` is just as unusable as an upper-case
    /// one. An empty slice has no N and passes.
    pub fn check_n(&self, sequence: &[u8]) -> bool {
        if !self.enabled {
            return true;
        }
        if sequence.is_empty() {
            return true;
        }

        let n_count = sequence.iter().filter(|&&b| b == b'N' || b == b'n').count();

        if self.max_n > 0.0 && self.max_n < 1.0 {
            // Fraction mode
            let ratio = n_count as f64 / sequence.len() as f64;
            ratio <= self.max_n
        } else {
            // Absolute count mode
            (n_count as f64) <= self.max_n
        }
    }

    /// Check whether the quality scores satisfy the threshold. An empty slice
    /// has no failing base and passes (guarding the division below, which the
    /// callers used to have to protect against themselves).
    pub fn check_quality(&self, quality: &[u8]) -> bool {
        if !self.enabled {
            return true;
        }
        if quality.is_empty() {
            return true;
        }

        let min_phred = self.min_quality.saturating_add(self.quality_base);
        let passed_count = quality.iter().filter(|&&q| q >= min_phred).count();
        let passed_percent = (passed_count * 100) / quality.len();

        passed_percent >= self.min_quality_percent as usize
    }
}
