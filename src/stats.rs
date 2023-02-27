use std::{
    str::FromStr,
    fmt,
    cmp::Ordering,
    collections::HashMap,
};

use thiserror::Error;

use crate::utils::{
    row_sums,
    col_sums,
    top_two,
};

use crate::constants::{
    CLOSE_TO_ZERO,
    N_METH_STATES,
    Base,
    CpgType,
};

use crate::sorting::{
    unsort_sorted_indexes,
    sort_with_floats,
    resort,
};

/// calculate Fisher's exact test
/// Method from (fastFishersExactTest):
///     https://genome.sph.umich.edu/w/images/b/b3/Bios615-fa12-lec03-presentation.pdf
pub fn log_hypergeometric_dist(lf: &Vec<f64>, m11: i64, m12: i64, m21: i64, m22: i64) -> f64 {
    lf[(m11+m12) as usize] +
    lf[(m21+m22) as usize] +
    lf[(m11+m21) as usize] +
    lf[(m12+m22) as usize] -
    lf[(m11) as usize] -
    lf[(m12) as usize] -
    lf[(m21) as usize] -
    lf[(m22) as usize] -
    lf[(m11+m12+m21+m22) as usize]
}

pub fn fishers_exact(m11: i64, m12: i64, m21: i64, m22: i64) -> f64 {
    let n: i64 = m11 + m12 + m21 + m22;

    // Pre-store factorial values as ln(N!) to handle potentially large values
    let mut log_factorials: Vec<f64> = vec![0.0; (n+1) as usize];
    for i in 1..n+1 {
        log_factorials[i as usize] = log_factorials[(i-1) as usize] + (i as f64).ln();
    }

    let p_cutoff = log_hypergeometric_dist(&log_factorials, m11, m12, m21, m22);

    let mut p: f64 = 0.0;
    for i in 0..n+1 {
        if m11+m12-i >= 0 && m11+m21-i >= 0 && m22-m11+i >= 0 {
            let l = log_hypergeometric_dist(&log_factorials, i, m11+m12-i, m11+m21-i, m22-m11+i);
            if l <= p_cutoff {
                p += (l as f64 - p_cutoff).exp();
            }
        }
    }

    let log_p = p_cutoff + p.ln();

    let mut p = log_p.exp();
    if p > 1.0 {
        p = 1.0;
    }

    p
}

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
    let rs: Vec<u32> = row_sums(fm, nrow, ncol);
    let cs: Vec<u32> = col_sums(fm, nrow, ncol);

    // Find top two values in rows and columns
    let (row_one, row_two) = top_two(&rs, nrow);
    let (col_one, col_two) = top_two(&cs, ncol);

    if rs[row_one] > 0 && rs[row_two] > 0 && cs[col_one] > 0 && cs[col_two] > 0 {
        let p_value = fishers_exact(
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
enum FdrTypeError {
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

/// Perform Bonferroni p-value adjustment
fn bonferroni(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();
    for v in p {
        let tmp: f64 = if v.p*n as f64 > 1.0 {
            1.0 
        } else { 
            v.p * n as f64
        };

        out.push( SnpCpgData { chr: v.chr, snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: tmp, pdata: v.pdata } );
    }

    out
}

/// Perform BH p-value adjustment, sorting done in-function
fn benjamini_hochberg(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();

    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Indexes for returning to original order
    let original_order = unsort_sorted_indexes(p, true);

    // Sort inputs in descending order
    let sorted = sort_with_floats(p, true);

    let mut curr_min = f64::MAX;
    let mut tmp;
    for (i, v) in sorted.iter().enumerate() {
        tmp = (n as f64 / rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        out.push( SnpCpgData { chr: v.chr, snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
    }

    resort(&out, &original_order)
}

/// Performs BY p-value adjustment, sorting done in-function
fn benjamini_yekutieli(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    // Euler-Mascheroni constant
    const GAMMA: f64 = 0.577215664901532;

    let mut out: Vec<SnpCpgData> = Vec::new();

    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Indexes for returning to original order
    let original_order = unsort_sorted_indexes(p, true);

    // Sort inputs in descending order
    let sorted = sort_with_floats(p, true);

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
    for (i, v) in sorted.iter().enumerate() {
        tmp = (cm * n as f64 / rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        out.push( SnpCpgData { chr: v.chr, snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
    }

    resort(&out, &original_order)
}

/// Perform Hochberg p-value adjustment, sorting done in-function
fn hochberg(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();

    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).rev().collect();

    // Indexes for returning to original order
    let original_order = unsort_sorted_indexes(p, true);

    // Sort inputs in descending order
    let sorted = sort_with_floats(p, true);

    let mut curr_min = f64::MAX;
    let mut tmp;
    for (i, v) in sorted.iter().enumerate() {
        tmp = (n as f64 + 1.0 - rank[i] as f64) * v.p;
        if tmp < curr_min {
            curr_min = tmp;
        }

        out.push( SnpCpgData { chr: v.chr, snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
    }

    resort(&out, &original_order)
}

/// Perform Holm p-value adjustment, sorting done in-function
fn holm(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();

    // Rank of p-values
    let rank: Vec<usize> = (1..p.len()+1).collect();

    // Indexes for returning to original order
    let original_order = unsort_sorted_indexes(p, false);

    // Sort inputs in descending order
    let sorted = sort_with_floats(p, false);

    let mut curr_max = f64::MIN;
    let mut tmp;
    for (i, v) in sorted.iter().enumerate() {
        tmp = (n as f64 + 1.0 - rank[i] as f64) * v.p;
        if tmp > curr_max {
            curr_max = tmp;
        }

        out.push( SnpCpgData { chr: v.chr, snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_max.min(1.0), pdata: v.pdata } );
    }

    resort(&out, &original_order)
}

pub fn false_discovery_correction(p: &Vec<SnpCpgData>, typ: &str, n: usize) -> Vec<SnpCpgData> {
    // No need to correct things if nothing there or only one entry
    if p.len() <= 1 {
        return p.to_vec();
    }

    // Get FDR type
    let t = match FdrType::from_str(typ) {
        Ok(f) => f,
        Err(err) => {
            eprintln!("{}", err);
            quit::with_code(1);
        },
    };

    let out: Vec<SnpCpgData>;
    match t {
        FdrType::Bh => {
            out = benjamini_hochberg(p, n);
        },
        FdrType::Holm => {
            out = holm(p, n);
        },
        FdrType::Hochberg => {
            out = hochberg(p, n);
        },
        FdrType::Bonferroni => {
            out = bonferroni(p, n);
        },
        FdrType::By => {
            out = benjamini_yekutieli(p, n);
        },
        FdrType::No => {
            out = p.clone();
        },
    };

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants;

    fn float_equality(x: f64, y: f64) -> bool {
        (x - y).abs() < constants::CLOSE_TO_ZERO
    }

    #[test]
    fn test_fishers_exact_sig() {
        let test = fishers_exact(15, 0, 0, 15);
        assert!(float_equality(test, 1.289e-8));
    }

    #[test]
    fn test_fishers_exact_not_sig() {
        let test = fishers_exact(15, 0, 15, 0);
        assert!(float_equality(test, 1.0));
    }
}
