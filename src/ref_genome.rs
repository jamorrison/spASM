use std::{
    path::Path,
};
use rust_htslib::faidx;

/// Open reference genome
///
/// Inputs
///     path: path to genome (must be coercible to std::path::Path)
/// Returns
///     faidx::Reader
pub fn open_ref<P: AsRef<Path>>(path: P) -> faidx::Reader {
    faidx::Reader::from_path(path).ok().unwrap()
}

/// Retrieve base from FASTA
///
/// Inputs
///     r: opened FASTA
///     pos: 0-based position in 
pub fn get_base(r: &faidx::Reader, chr: &str, pos: usize) -> Option<char> {
    let seq = r.fetch_seq(chr, pos, pos).unwrap();
    let base: Vec<char> = seq.iter().map(|b| *b as char).collect();

    base.into_iter().nth(0)
}
