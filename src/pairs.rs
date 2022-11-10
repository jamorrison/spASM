use std::fmt;

use crate::constants::N_METH_STATES;

/// Type of SNP (A/C/G/T/N)
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
pub enum SnpType {
    SnpA = 0,
    SnpC = 1,
    SnpG = 2,
    SnpT = 3,
    SnpN = 4,
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

/// Location and type of SNP
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Snp {
    /// chromosome
    chr: String,
    /// 0-based location
    pos: u64,
    /// SNP type
    typ: SnpType,
}

impl Snp {
    pub fn new(chr: String, pos: u64, typ: char) -> Snp {
        Snp { chr: chr, pos: pos, typ: SnpType::from(typ).unwrap() }
    }

    pub fn get_chr(&self) -> &String {
        &self.chr
    }

    pub fn get_pos(&self) -> &u64 {
        &self.pos
    }
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

/// Location and status of CpG
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Cpg {
    /// chromosome
    chr: String,
    /// 0-based location
    pos: u64,
    /// CpG type
    typ: CpgType,
}

impl Cpg {
    pub fn new(chr: String, pos: u64, typ: char) -> Cpg {
        Cpg { chr: chr, pos: pos, typ: CpgType::from(typ).unwrap() }
    }

    pub fn get_pos(&self) -> &u64 {
        &self.pos
    }
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
    snp: Snp,
    /// CpG
    cpg: Cpg,
}

impl Pair {
    pub fn new(s: Snp, c: Cpg) -> Pair {
        Pair { snp: s, cpg: c }
    }

    pub fn get_snp(&self) -> &Snp {
        &self.snp
    }

    pub fn get_cpg(&self) -> &Cpg {
        &self.cpg
    }

    pub fn get_index(&self) -> usize {
        N_METH_STATES*<SnpType as Into<usize>>::into(self.snp.typ) + <CpgType as Into<usize>>::into(self.cpg.typ)
    }
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
