use std::{
    str::FromStr,
    fmt,
    cmp::Ordering,
};

use snafu::prelude::*;

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

/// calculate fisher's exact test p-value
/// m11 m12 | r1
/// m21 m22 | r2
/// --------+----
/// c1  c2  | n
pub struct HyperGeoAcc {
    m11: i64,
    r1 : i64,
    c1 : i64,
    n  : i64,
    p  : f64,
}

impl HyperGeoAcc {
    fn new() -> HyperGeoAcc {
        HyperGeoAcc {
            m11: 0,
            r1:  0,
            c1:  0,
            n:   0,
            p:   0.0
        }
    }
}

#[derive(Debug, Snafu)]
#[snafu(display("Invalid value given to ln_gamma. Z must be >= 0 (z = {z})"))]
pub struct GammaError{ z: i64 }

/// Reference: "Lanczos, C. 'A precision approximation of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
/// Transposition of  Alan Miller's FORTRAN-implementation. See:  http://lib.stat.cmu.edu/apstat/245
pub fn ln_gamma(z: i64) -> Result<f64, GammaError> {
    const COEFS: &[f64] = &[
        0.9999999999995183,
        676.5203681218835,
        -1259.139216722289,
        771.3234287757674,
        -176.6150291498386,
        12.50734324009056,
        -0.1385710331296526,
        0.9934937113930748e-05,
        0.1659470187408462e-06
    ];
    const LN_SQRT_2PI: f64 = 0.9189385332046727;

    if z < 0 {
        return Err(GammaError{ z });
    }

    let     z        = z as f64;
    let mut out: f64 = 0.0;
    let mut tmp: f64 = 7.0 + z;

    for i in (1..9).rev() {
        out += COEFS[i] / tmp;
        tmp -= 1.0;
    }
    out += COEFS[0];

    Ok(out.ln() + LN_SQRT_2PI - (6.5 + z) + (z - 0.5)*(z + 6.5).ln())
}

pub fn ln_binomial(n: i64, k: i64) -> Result<f64, GammaError> {
    if k == 0 || n == 0 || n-k <= 0 {
        return Ok(0.0);
    }

    let n1  = match ln_gamma(n+1) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let k1  = match ln_gamma(k+1) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let nk1 = match ln_gamma(n-k+1) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };

    Ok(n1 - k1 - nk1)
}

pub fn hypergeometric(m11: i64, r1: i64, c1: i64, n: i64) -> Result<f64, GammaError> {
    let a = match ln_binomial(r1, m11) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let b = match ln_binomial(n-r1, c1-m11) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let c = match ln_binomial(n, c1) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };

    Ok((a + b - c).exp())
}

pub fn hypergeometric_accumulator(m11: i64, r1: i64, c1: i64, n: i64, acc: &mut HyperGeoAcc) -> Result<f64, GammaError> {
    if r1 > 0 || c1 > 0 || n > 0 {
        acc.m11 = m11;
        acc.r1  = r1;
        acc.c1  = c1;
        acc.n   = n;
    } else {
        // only m11 changed, the others are fixed
        if m11 % 11 > 0 && m11+acc.n-acc.r1-acc.c1 > 0 {
            if m11 == acc.m11+1 { // incremental
                acc.p *= ((acc.r1 - acc.m11) / m11) as f64 * ((acc.c1 - acc.m11) / (m11 + acc.n - acc.r1 - acc.c1)) as f64;
                acc.m11 = m11;

                return Ok(acc.p);
            }
            if m11 == acc.m11-1 { // incremental
                acc.p *= (acc.m11 / (acc.r1 - m11)) as f64 * ((acc.m11 + acc.n - acc.r1 - acc.c1) / (acc.c1 - m11)) as f64;
                acc.m11 = m11;

                return Ok(acc.p);
            }
        }

        acc.m11 = m11;
    }

    acc.p = match hypergeometric(acc.m11, acc.r1, acc.c1, acc.n) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };

    return Ok(acc.p);
}

pub fn logHypergeometricProb(logFacs: &Vec<f64>, m11: i64, m12: i64, m21: i64, m22: i64) -> f64 {
    logFacs[(m11+m12) as usize] +
    logFacs[(m21+m22) as usize] +
    logFacs[(m11+m21) as usize] +
    logFacs[(m12+m22) as usize] -
    logFacs[(m11) as usize] -
    logFacs[(m12) as usize] -
    logFacs[(m21) as usize] -
    logFacs[(m22) as usize] -
    logFacs[(m11+m12+m21+m22) as usize]
}

pub fn other_implementation(m11: i64, m12: i64, m21: i64, m22: i64) -> f64 {
    let n: i64 = m11 + m12 + m21 + m22;
    let mut logFacs: Vec<f64> = vec![0.0; (n+1) as usize];
    for i in 1..n+1 {
        logFacs[i as usize] = logFacs[(i-1) as usize] + (i as f64).ln();
    }

    let logpCutoff = logHypergeometricProb(&logFacs, m11, m12, m21, m22);

    let mut pFraction: f64 = 0.0;
    for i in 0..n+1 {
        if m11+m12-i >= 0 && m11+m21-i >= 0 && m22-m11+i >= 0 {
            let l = logHypergeometricProb(&logFacs, i, m11+m12-i, m11+m21-i, m22-m11+i);
            if l <= logpCutoff {
                pFraction += (l as f64 - logpCutoff).exp();
            }
        }
    }

    let logpValues = logpCutoff + pFraction.ln();

    logpValues.exp()
}

fn kinda_round(f: f64, n: u32) -> f64 {
    let exp = 10_i32.pow(n) as f64;
    (f * exp).round() / exp
}

pub fn fisher_exact(m11: i64, m12: i64, m21: i64, m22: i64) -> Result<(f64, f64, f64), GammaError> {
    let r1: i64 = m11 + m12; // sum of first row
    let c1: i64 = m11 + m21; // sum of first column
    let n : i64 = m11 + m12 + m21 + m22; // sum of all values

    let max: i64 = if c1 < r1 { c1 } else { r1 }; // for right tail
    let min: i64 = if r1 + c1 - n < 0 { 0 } else { r1 + c1 - n }; // for left tail

    if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
        println!("r1: {} c1: {} n: {} min: {} max: {}", &r1, &c1, &n, &min, &max);
    }

    // No need to do test
    if min == max {
        return Ok((1.0, 1.0, 1.0));
    }

    // Accumulator of hypergeometric test values
    let mut acc = HyperGeoAcc::new();

    // Probability of current table
    let q = match hypergeometric_accumulator(m11, r1, c1, n, &mut acc) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };

    // Left tail
    let mut left: f64 = 0.0;
    let mut p: f64 = match hypergeometric_accumulator(min, 0, 0, 0, &mut acc) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let mut i: i64 = min + 1;

    if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
        println!("i-while before - i: {} p: {} q: {} <q: {} >q: {}", &i, &p, &q, 0.99999999*q, 1.00000001*q);
    }
    while i <= max && p < 0.99999999*q {
        left += p;
        p = match hypergeometric_accumulator(i, 0, 0, 0, &mut acc) {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        i += 1;

        if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
            println!("i-while - p: {} q: {} left: {}", &p, &q, &left);
        }
    }

    i -= 1;
    if p < 1.00000001*q {
        left += p;
    } else {
        i -= 1;
    }

    if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
        println!("left - p: {} q: {} left: {}", &p, &q, &left);
    }

    // Right tail
    p = match hypergeometric_accumulator(max, 0, 0, 0, &mut acc) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let mut right: f64 = 0.0;
    let mut j: i64 = max - 1;
    while j >= 0 && p < 0.99999999*q {
        right += p;
        p = match hypergeometric_accumulator(j, 0, 0, 0, &mut acc) {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        j -= 1;
    }

    j += 1;
    if p < 1.00000001*q {
        right += p;
    } else {
        j += 1;
    }

    // Two tail
    let mut two = left + right;
    if two > 1.0 {
        two = 1.0;
    }

    if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
        println!("{} {} {}", &left, &right, &two);
    }

    // Adjust left and right tails
    if (i-m11).abs() < (j-m11).abs() {
        right = 1.0 - left + q;
    } else {
        left = 1.0 - right + q;
    }

    let p_other = other_implementation(m11, m12, m21, m22);
    if kinda_round(p_other, 8) - kinda_round(two, 8) > CLOSE_TO_ZERO {
        println!("{} {} {} {} {} {}", m11, m12, m21, m22, two, p_other);
    }
    if m11 == 15 && m12 == 13 && m21 == 0 && m22 == 1 {
        println!("max val: {}", f64::MAX / 100000000.0);
        println!("other implementation: {}", p_other);
    }

    Ok((left, right, two))
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
            //"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:10.16e}\t.\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
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
            //"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:10.6e}\n",
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

pub fn calculate_p_values(fm: &Vec<i64>, nrow: usize, ncol: usize, verbose: &usize) -> Option<(f64, f64, f64, PValMetadata)> {
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
        let (left, right, two) = match fisher_exact(
            fm[N_METH_STATES*row_one+col_one],
            fm[N_METH_STATES*row_one+col_two],
            fm[N_METH_STATES*row_two+col_one],
            fm[N_METH_STATES*row_two+col_two]
        ) {
            Ok((l, r, t)) => (l, r, t),
            Err(e) => {
                if *verbose >= 2 {
                    eprintln!("NOTE: {}. Skipping.", e);
                }
                return None;
            },
        };

        return Some(
            (left, right, two,
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
