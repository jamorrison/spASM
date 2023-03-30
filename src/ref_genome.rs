use std::{
    path::Path,
    collections::HashMap,
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

/// Create lookup table and its reverse to translate chromosomes into IDs
///
/// Inputs
///     path: path to reference (must be FAI indexed)
/// Returns
///     ( HashMap::<String, u32>, HashMap::<u32, String> )
pub fn chromosome_lookup_tables<P: AsRef<Path>>(path: P) -> (HashMap::<String, u32>, HashMap::<u32, String>) {
    let mut k_chr: HashMap::<String, u32> = HashMap::new();
    let mut k_int: HashMap::<u32, String> = HashMap::new();
    let genome = open_ref(path);

    let nseq = genome.n_seqs();
    for i in 0..nseq {
        let seq = genome.seq_name(i as i32).unwrap();
        k_chr.insert(seq.clone(), i as u32);
        k_int.insert(i as u32, seq);
    }

    (k_chr, k_int)
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
