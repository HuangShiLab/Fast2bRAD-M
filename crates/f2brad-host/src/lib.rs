//! Host genotyping primitives for fast2bRAD-holo.
//!
//! Input: a reference genome + a Type IIB restriction enzyme.
//! Output: a locus-keyed tag database, per-position genotype likelihoods, and
//! VCF/BIMBAM files.

pub mod digest;
