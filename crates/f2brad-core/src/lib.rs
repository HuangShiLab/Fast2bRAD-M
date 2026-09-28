//! Core primitives shared by the Fast2bRAD-M microbial pipeline and the
//! fast2bRAD-holo host-genotyping arm.
//!
//! This crate intentionally contains no CLI entry points — only the enzyme
//! definitions, tag-extraction logic, binary formats, and shared types.

pub mod enzymes;
pub mod extract;
pub mod io_utils;
pub mod types;
