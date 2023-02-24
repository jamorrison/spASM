use std::{
    fmt,
    cmp::Ordering,
};

use crate::constants::{
    N_METH_STATES,
    Base,
};

/// Location and type of SNP
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
pub struct Snp {
    /// chromosome
    chr: u32,
    /// 0-based location
    pos: u32,
    /// SNP type
    typ: Base,
}

impl Snp {
    pub fn new(chr: u32, pos: u32, typ: char) -> Snp {
        Snp { chr: chr, pos: pos, typ: Base::from(typ).unwrap() }
    }

    pub fn get_chr(&self) -> &u32 {
        &self.chr
    }

    pub fn get_pos(&self) -> &u32 {
        &self.pos
    }
}

impl fmt::Display for Snp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "Chromosome: {} Position: {} Base: {}",
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
    chr: u32,
    /// 0-based location
    pos: u32,
    /// CpG type
    typ: CpgType,
}

impl Cpg {
    pub fn new(chr: u32, pos: u32, typ: char) -> Cpg {
        Cpg { chr: chr, pos: pos, typ: CpgType::from(typ).unwrap() }
    }

    pub fn get_pos(&self) -> &u32 {
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
#[derive(Debug, Eq, PartialEq, Ord, Clone)]
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
        N_METH_STATES*<Base as Into<usize>>::into(self.snp.typ) + <CpgType as Into<usize>>::into(self.cpg.typ)
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

impl PartialOrd for Pair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // snp chromosome -> snp position -> cpg chromosome -> cpg position
        if self.snp.chr == other.snp.chr {
            if self.snp.pos == other.snp.pos {
                if self.cpg.chr == other.cpg.chr {
                    self.cpg.pos.partial_cmp(&other.cpg.pos)
                } else {
                    self.cpg.chr.partial_cmp(&other.cpg.chr)
                }
            } else {
                self.snp.pos.partial_cmp(&other.snp.pos)
            }
        } else {
            self.snp.chr.partial_cmp(&other.snp.chr)
        }
    }
}
