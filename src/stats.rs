use std::{
    str::FromStr,
    fmt,
    cmp::Ordering,
};

use crate::utils::{
    row_sums,
    col_sums,
    top_two,
};

use crate::constants::{
    CLOSE_TO_ZERO,
    N_METH_STATES,
    Base,
};

use crate::pairs::{
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

fn kinda_round(f: f64, n: u32) -> f64 {
    let exp = 10_i32.pow(n) as f64;
    (f * exp).round() / exp
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
    mat: (i64, i64, i64, i64),
}

impl PValMetadata {
    pub fn new(s1: Base, s2: Base, c1: CpgType, c2: CpgType, m: (i64, i64, i64, i64)) -> PValMetadata {
        PValMetadata { snp1: s1, snp2: s2, cpg1: c1, cpg2: c2, mat: m }
    }
}

/// SNP-CpG pair with P-value
#[derive(Debug, PartialEq, Clone)]
pub struct SnpCpgData {
    /// Chromosome
    chr: String,
    /// SNP location (0-based)
    snp_pos: u64,
    /// CpG location (0-based)
    cpg_pos: u64,
    /// p-value
    p: f64,
    /// Metadata on the p-value
    pdata: PValMetadata,
}

impl SnpCpgData {
    pub fn new(chr: String, s_pos: u64, c_pos: u64, p: f64, pdata: PValMetadata) -> SnpCpgData {
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
    pub fn to_bedpe(&self, cutoff: f64) -> String {
        let name: &str = if self.p < cutoff {
            "candidate"
        } else {
            "non_candidate"
        };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            //"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:10.16e}\t.\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.chr,
            self.snp_pos,
            self.snp_pos+1,
            self.chr,
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
    pub fn to_biscuit_asm(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            //"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:10.6e}\n",
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

pub fn calculate_p_values(fm: &Vec<i64>, nrow: usize, ncol: usize, verbose: &usize) -> Option<(f64, PValMetadata)> {
    if nrow < 2 || ncol < 2 {
        eprintln!("nrow ({}) or ncol ({}) not large enough for Fisher's exact test", nrow, ncol);
        quit::with_code(1);
    }

    // Find row and column sums
    let rs: Vec<i64> = row_sums(fm, nrow, ncol);
    let cs: Vec<i64> = col_sums(fm, nrow, ncol);

    // Find top two values in rows and columns
    let (row_one, row_two) = top_two(&rs, nrow);
    let (col_one, col_two) = top_two(&cs, ncol);

    if rs[row_one] > 0 && rs[row_two] > 0 && cs[col_one] > 0 && cs[col_two] > 0 {
        let p_value = fishers_exact(
            fm[N_METH_STATES*row_one+col_one],
            fm[N_METH_STATES*row_one+col_two],
            fm[N_METH_STATES*row_two+col_one],
            fm[N_METH_STATES*row_two+col_two]
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
    type Err = std::string::ParseError;

    fn from_str(fdr: &str) -> Result<Self, Self::Err> {
        let tmp = fdr.to_uppercase();
        let out: FdrType = if tmp == "BH".to_string()  {
            FdrType::Bh
        } else if tmp == "HOLM" {
            FdrType::Holm
        } else if tmp == "HOCHBERG" {
            FdrType::Hochberg
        } else if tmp == "BONFERRONI" {
            FdrType::Bonferroni
        } else if tmp == "BY" {
            FdrType::By
        } else if tmp == "NO" {
            FdrType::No
        } else {
            FdrType::No
        };

        Ok(out)
    }
}

/// Perform Bonferonni p-value adjustment
fn bonferroni(p: &Vec<SnpCpgData>, n: usize) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();
    for v in p {
        let tmp: f64 = if v.p*n as f64 > 1.0 {
            1.0 
        } else { 
            v.p * n as f64
        };

        out.push( SnpCpgData { chr: v.chr.clone(), snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: tmp, pdata: v.pdata } );
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

        out.push( SnpCpgData { chr: v.chr.clone(), snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
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

        out.push( SnpCpgData { chr: v.chr.clone(), snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
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

        out.push( SnpCpgData { chr: v.chr.clone(), snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_min.min(1.0), pdata: v.pdata } );
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

        out.push( SnpCpgData { chr: v.chr.clone(), snp_pos: v.snp_pos, cpg_pos: v.cpg_pos, p: curr_max.min(1.0), pdata: v.pdata } );
    }

    resort(&out, &original_order)
}

pub fn false_discovery_correction(p: &Vec<SnpCpgData>, typ: &str, n: usize) -> Vec<SnpCpgData> {
    // No need to correct things if nothing there or only one entry
    if p.len() <= 1 {
        return p.to_vec();
    }

    // Get FDR type
    let t = FdrType::from_str(typ).unwrap();

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

    //#[test]
    //fn test_ln_gamma() {
    //    let x: f64 = ln_gamma(3);

    //    assert!(float_equality(x, 0.693147181));
    //}

    //#[test]
    //fn test_ln_binomial() {
    //    let x: f64 = ln_binomial(3, 2);

    //    assert!(float_equality(x, 1.098612289));
    //}

    //#[test]
    //fn test_hypergeometric() {
    //    let x: f64 = hypergeometric(2, 3, 4, 6);

    //    assert!(float_equality(x, 0.6));
    //}

    //#[test]
    //fn test_fisher_exact_left() {
    //    let (l, _, _) = fisher_exact(2, 1, 2, 1);

    //    assert!(float_equality(l, 0.8));
    //}

    //#[test]
    //fn test_fisher_exact_right() {
    //    let (_, r, _) = fisher_exact(2, 1, 2, 1);

    //    assert!(float_equality(r, 0.8));
    //}

    //#[test]
    //fn test_fisher_exact_two() {
    //    let (_, _, t) = fisher_exact(2, 1, 2, 1);

    //    assert!(float_equality(t, 1.0));
    //}

    #[test]
    fn test_other_implementation_1() {
        let t = other_implementation(6, 2, 0, 1);
        println!("{}", t-0.3333);

        assert!(float_equality(t, 0.3333));
    }

    #[test]
    fn test_other_implementation_2() {
        let t = other_implementation(12345, 67890, 54321, 9876);
        println!("{:.50}", t);

        assert!(float_equality(t, 0.02));
        //assert!(float_equality(t, 0.00000000000000002));
    }
}
