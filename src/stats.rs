use std::{
    str::FromStr,
    fmt,
    cmp::Ordering,
    collections::HashMap,
};

use thiserror::Error;

use crate::utils;

use crate::fishers_exact_test;
use crate::constants::{
    CLOSE_TO_ZERO,
    N_METH_STATES,
    Base,
    CpgType,
};

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
pub struct PValMetadata {
    /// SNP1
    snp1: Base,
    /// SNP2
    snp2: Base,
    /// CpG1
    cpg1: CpgType,
    /// CpG2
    cpg2: CpgType,
    /// Contingency table
    mat: (u32, u32, u32, u32),
}

impl PValMetadata {
    pub fn new(s1: Base, s2: Base, c1: CpgType, c2: CpgType, m: (u32, u32, u32, u32)) -> PValMetadata {
        PValMetadata { snp1: s1, snp2: s2, cpg1: c1, cpg2: c2, mat: m }
    }
}

/// SNP-CpG pair with P-value
#[derive(Debug, PartialEq, Clone)]
pub struct SnpCpgData {
    /// Chromosome
    chr: u32,
    /// SNP location (0-based)
    snp_pos: u32,
    /// CpG location (0-based)
    cpg_pos: u32,
    /// p-value
    p: f64,
    /// Metadata on the p-value
    pdata: PValMetadata,
}

impl SnpCpgData {
    pub fn new(chr: u32, s_pos: u32, c_pos: u32, p: f64, pdata: PValMetadata) -> SnpCpgData {
        SnpCpgData { chr: chr, snp_pos: s_pos, cpg_pos: c_pos, p: p, pdata: pdata }
    }

    pub fn get_p(&self) -> &f64 {
        &self.p
    }

    pub fn get_snp1(&self) -> &Base {
        &self.pdata.snp1
    }

    pub fn get_snp2(&self) -> &Base {
        &self.pdata.snp2
    }

    /// Write in BEDPE format
    /// SNP chr, SNP start, SNP end, CpG chrom, CpG start, CpG end, name, p-value, SNP strand, CpG
    /// strand, SNP1, SNP2, CpG1, CpG2, m11, m12, m21, m22
    pub fn to_bedpe(&self, k_int: &HashMap::<u32, String>, cutoff: f64) -> String {
        let name: &str = if self.p < cutoff {
            "candidate"
        } else {
            "non_candidate"
        };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.10e}\t.\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            k_int.get(&self.chr).unwrap(),
            self.snp_pos,
            self.snp_pos+1,
            k_int.get(&self.chr).unwrap(),
            self.cpg_pos,
            self.cpg_pos+2,
            name,
            self.p,
            self.pdata.snp1,
            self.pdata.snp2,
            self.pdata.cpg1,
            self.pdata.cpg2,
            self.pdata.mat.0,
            self.pdata.mat.1,
            self.pdata.mat.2,
            self.pdata.mat.3,
        )
    }

    /// Write in BISCUIT ASM format
    /// chr, snp position, cpg position, SNP1, SNP2, CPG1, CPG2, m11, m12, m21, m22, p-value
    pub fn to_biscuit_asm(&self, k_int: &HashMap::<u32, String>) -> String {
        format!(
            "{}\t{}\t{}\t{}/{}\t{}/{}\t{}\t{}\t{}\t{}\t{:.10e}\t.\n",
            k_int.get(&self.chr).unwrap(),
            self.snp_pos,
            self.cpg_pos,
            self.pdata.snp1,
            self.pdata.snp2,
            self.pdata.cpg1,
            self.pdata.cpg2,
            self.pdata.mat.0,
            self.pdata.mat.1,
            self.pdata.mat.2,
            self.pdata.mat.3,
            self.p,
        )
    }
}

impl fmt::Display for SnpCpgData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "chr: {}, SNP pos: {}, CpG pos: {}, SNP1: {} SNP2: {}, CpG1: {}, CpG2: {}, contingency table: (m11: {}, m12: {}, m21: {}, m22: {}) p-value: {}",
            self.chr,
            self.snp_pos,
            self.cpg_pos,
            self.pdata.snp1,
            self.pdata.snp2,
            self.pdata.cpg1,
            self.pdata.cpg2,
            self.pdata.mat.0,
            self.pdata.mat.1,
            self.pdata.mat.2,
            self.pdata.mat.3,
            self.p
        );

        write!(f, "{}", s)
    }
}

impl PartialOrd for SnpCpgData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if (self.p - &other.p).abs() < CLOSE_TO_ZERO {
            if self.chr == other.chr {
                if self.snp_pos == other.snp_pos {
                    self.cpg_pos.partial_cmp(&other.cpg_pos)
                } else {
                    self.snp_pos.partial_cmp(&other.snp_pos)
                }
            } else {
                self.chr.partial_cmp(&other.chr)
            }
        } else {
            self.p.partial_cmp(&other.p)
        }
    }
}

pub fn calculate_p_values(fm: &[u32], nrow: usize, ncol: usize) -> Option<(f64, PValMetadata)> {
    if nrow < 2 || ncol < 2 {
        eprintln!("nrow ({}) or ncol ({}) not large enough for Fisher's exact test", nrow, ncol);
        quit::with_code(1);
    }

    // Find row and column sums
    let rs: Vec<u32> = utils::row_sums(fm, nrow, ncol);
    let cs: Vec<u32> = utils::col_sums(fm, nrow, ncol);

    // Find top two values in rows and columns
    let (row_one, row_two) = utils::top_two(&rs, nrow);
    let (col_one, col_two) = utils::top_two(&cs, ncol);

    if rs[row_one] > 0 && rs[row_two] > 0 && cs[col_one] > 0 && cs[col_two] > 0 {
        let p_value = fishers_exact_test::fishers_exact_test(
            fm[N_METH_STATES*row_one+col_one] as i64,
            fm[N_METH_STATES*row_one+col_two] as i64,
            fm[N_METH_STATES*row_two+col_one] as i64,
            fm[N_METH_STATES*row_two+col_two] as i64
        );

        return Some(
            (p_value,
             PValMetadata::new(
                 Base::from_usize(row_one).unwrap(),
                 Base::from_usize(row_two).unwrap(),
                 CpgType::from_usize(col_one).unwrap(),
                 CpgType::from_usize(col_two).unwrap(),
                 (fm[N_METH_STATES*row_one+col_one], fm[N_METH_STATES*row_one+col_two],
                  fm[N_METH_STATES*row_two+col_one], fm[N_METH_STATES*row_two+col_two]),
             )
            )
        );
    }

    None
}

/// Available types of FDR corrections
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
enum FdrType {
    Bh,
    Holm,
    Hochberg,
    Bonferroni,
    By,
    No,
}

/// FdrType errors
#[derive(Error, Debug, PartialEq)]
pub enum FdrTypeError {
    /// unknown type
    #[error("unknown FDR type. try spasm --help for options")]
    UnknownType,
}


impl fmt::Display for FdrType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = match self {
            FdrType::Bh => "Benjamini-Hochberg".to_string(),
            FdrType::Holm => "Holm".to_string(),
            FdrType::Hochberg => "Hochberg".to_string(),
            FdrType::Bonferroni => "Bonferroni".to_string(),
            FdrType::By => "Benjamini-Yekutieli".to_string(),
            FdrType::No => "No Correction".to_string(),
        };

        write!(f, "{}", s)
    }
}

impl FromStr for FdrType {
    type Err = FdrTypeError;

    fn from_str(fdr: &str) -> Result<Self, Self::Err> {
        let out: FdrType = match fdr {
            "BH" => FdrType::Bh,
            "HOLM" => FdrType::Holm,
            "HOCHBERG" => FdrType::Hochberg,
            "BONFERRONI" => FdrType::Bonferroni,
            "BY" => FdrType::By,
            "NO" => FdrType::No,
            _ => return Err(FdrTypeError::UnknownType),
        };

        Ok(out)
    }
}

pub fn validate_fdr_type(fdr: &str) -> Result<String, FdrTypeError> {
    match FdrType::from_str(fdr.to_uppercase().as_str()) {
        Ok(_) => Ok(String::from(fdr)),
        Err(e) => Err(e),
    }
}

/// Perform Bonferroni p-value adjustment
fn bonferroni(p: &mut Vec<SnpCpgData>, n: usize) {
    for v in p.iter_mut() {
        let tmp: f64 = if v.p*n as f64 > 1.0 {
            1.0 
        } else { 
            v.p * n as f64
        };

        v.p = tmp;
    }
}

/// Perform BH p-value adjustment, sorting done in-function
fn benjamini_hochberg(p: &mut Vec<SnpCpgData>, n: usize) {
    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Sort inputs in descending order
    utils::sort_with_floats(p, true);

    let mut curr_min = f64::MAX;
    for (i, v) in p.iter_mut().enumerate() {
        let tmp = (n as f64 / rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        v.p = curr_min.min(1.0);
    }
}

/// Performs BY p-value adjustment, sorting done in-function
fn benjamini_yekutieli(p: &mut Vec<SnpCpgData>, n: usize) {
    // Euler-Mascheroni constant
    const GAMMA: f64 = 0.577215664901532;

    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Sort inputs in descending order
    utils::sort_with_floats(p, true);

    // Calculate harmonic number
    // Starting at 125, the percent difference between the approximation of the harmonic constant
    // and the actual value is <0.0001%, so at this point, switch over to using the approximation
    // and leave off calculating the sum
    let cm: f64 = if n >= 125 {
        (n as f64).ln() + GAMMA + 1.0/(2.0*n as f64)
    } else {
        (1..n+1).map(|x| 1.0/(x as f64)).sum()
    };

    let mut curr_min = f64::MAX;
    let mut tmp;
    for (i, v) in p.iter_mut().enumerate() {
        tmp = (cm * n as f64 / rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        v.p = curr_min.min(1.0);
    }
}

/// Perform Hochberg p-value adjustment, sorting done in-function
fn hochberg(p: &mut Vec<SnpCpgData>, n: usize) {
    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Sort inputs in descending order
    utils::sort_with_floats(p, true);

    let mut curr_min = f64::MAX;
    let mut tmp;
    for (i, v) in p.iter_mut().enumerate() {
        tmp = (n as f64 + 1.0 - rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        v.p = curr_min.min(1.0);
    }
}

/// Perform Holm p-value adjustment, sorting done in-function
fn holm(p: &mut Vec<SnpCpgData>, n: usize) {
    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).collect();

    // Sort inputs in descending order
    utils::sort_with_floats(p, false);

    let mut curr_max = f64::MIN;
    let mut tmp;
    for (i, v) in p.iter_mut().enumerate() {
        tmp = (n as f64 + 1.0 - rank[i] as f64) * v.p;
        if tmp > curr_max {
            curr_max = tmp;
        }

        v.p = curr_max.min(1.0);
    }
}

pub fn false_discovery_correction(p: &mut Vec<SnpCpgData>, typ: &String, n: usize) {
    // No need to correct things if nothing there or only one entry
    if p.len() <= 1 {
        return;
    }

    // Get FDR type
    let t = match FdrType::from_str(typ.to_uppercase().as_str()) {
        Ok(f) => f,
        Err(err) => {
            eprintln!("{}", err);
            quit::with_code(1);
        },
    };

    match t {
        FdrType::Bh => {
            benjamini_hochberg(p, n);
        },
        FdrType::Holm => {
            holm(p, n);
        },
        FdrType::Hochberg => {
            hochberg(p, n);
        },
        FdrType::Bonferroni => {
            bonferroni(p, n);
        },
        FdrType::By => {
            benjamini_yekutieli(p, n);
        },
        FdrType::No => {
            return;
        },
    }
}
