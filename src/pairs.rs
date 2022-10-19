use std::{
    fmt,
    cmp::Ordering,
};

use crate::stats::PValMetadata;

const CLOSE_TO_ZERO: f64 = 0.000001;

pub const N_METH_STATES: usize = 2;

/// Type of SNP (A/C/G/T/N)
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
pub enum SnpType {
    SnpA = 0,
    SnpC = 1,
    SnpG = 2,
    SnpT = 3,
    SnpN = 4,
}

impl fmt::Display for SnpType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = match self {
            SnpType::SnpA => "A".to_string(),
            SnpType::SnpC => "C".to_string(),
            SnpType::SnpG => "G".to_string(),
            SnpType::SnpT => "T".to_string(),
            SnpType::SnpN => "N".to_string(),
        };

        write!(f, "{}", s)
    }
}

impl Into<usize> for SnpType {
    fn into(self) -> usize {
        self as usize
    }
}

impl SnpType {
    /// SnpType::from('A') == Ok(SnpType::SnpA)
    pub fn from(c: char) -> Result<SnpType, ()> {
        match c {
            'A' => Ok(SnpType::SnpA),
            'C' => Ok(SnpType::SnpC),
            'G' => Ok(SnpType::SnpG),
            'T' => Ok(SnpType::SnpT),
            'N' => Ok(SnpType::SnpN),
            _   => Err(()),
        }
    }

    /// SnpType::from_usize(0) == Ok(SnpType::SnpA)
    pub fn from_usize(i: usize) -> Result<SnpType, ()> {
        match i {
            0 => Ok(SnpType::SnpA),
            1 => Ok(SnpType::SnpC),
            2 => Ok(SnpType::SnpG),
            3 => Ok(SnpType::SnpT),
            4 => Ok(SnpType::SnpN),
            _   => Err(()),
        }
    }
}

/// Location and type of SNP
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Snp {
    /// chromosome
    pub chr: String,
    /// 0-based location
    pub pos: u64,
    /// SNP type
    pub typ: SnpType,
}

impl fmt::Display for Snp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "Chromosome: {} Position: {} Type: {:?}",
            self.chr,
            self.pos,
            self.typ,
        );

        write!(f, "{}", s)
    }
}

/// Methylation status
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
pub enum CpgType {
    NotMeth = 0,
    Meth = 1,
}

impl fmt::Display for CpgType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = match self {
            CpgType::NotMeth => "U".to_string(),
            CpgType::Meth => "M".to_string(),
        };

        write!(f, "{}", s)
    }
}

impl Into<usize> for CpgType {
    fn into(self) -> usize {
        self as usize
    }
}

impl CpgType {
    /// CpgType::from('M') == Ok(CpgType::Meth)
    pub fn from(c: char) -> Result<CpgType, ()> {
        match c {
            'U' => Ok(CpgType::NotMeth),
            'M' => Ok(CpgType::Meth),
            _   => Err(()),
        }
    }

    /// CpgType::from_usize(1) == Ok(CpgType::Meth)
    pub fn from_usize(i: usize) -> Result<CpgType, ()> {
        match i {
            0 => Ok(CpgType::NotMeth),
            1 => Ok(CpgType::Meth),
            _   => Err(()),
        }
    }
}

/// Location and status of CpG
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Cpg {
    /// chromosome
    pub chr: String,
    /// 0-based location
    pub pos: u64,
    /// SNP type
    pub typ: CpgType,
}

impl fmt::Display for Cpg {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "Chromosome: {} Position: {} Type: {:?}",
            self.chr,
            self.pos,
            self.typ,
        );

        write!(f, "{}", s)
    }
}

/// Store SNP-CpG pair information
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Pair {
    /// SNP
    pub snp: Snp,
    /// CpG
    pub cpg: Cpg,
}

impl fmt::Display for Pair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "SNP: {} CpG: {}",
            self.snp,
            self.cpg,
        );

        write!(f, "{}", s)
    }
}

impl Pair {
    pub fn get_index(&self) -> usize {
        N_METH_STATES*<SnpType as Into<usize>>::into(self.snp.typ) + <CpgType as Into<usize>>::into(self.cpg.typ)
    }
}

/// SNP-CpG pair with P-value
// TODO: This might be a confusing name as the Pair and PairP structs are related, but PairP isn't
//       a superset of Pair. Maybe look at renaming this?
#[derive(Debug, PartialEq, Clone)]
pub struct PairP {
    /// Chromosome
    pub chr: String,
    /// SNP location (0-based)
    pub snp_pos: u64,
    /// CpG location (0-based)
    pub cpg_pos: u64,
    /// p-value
    pub p: f64,
    /// Metadata on the p-value
    pub pdata: PValMetadata,
}

impl PairP {
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

impl fmt::Display for PairP {
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

impl PartialOrd for PairP {
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
