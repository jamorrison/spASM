use std::fmt;

/// Small number for comparing floats
pub const CLOSE_TO_ZERO: f64 = 0.00000000001;

/// Number of methylation states (methylated or unmethylated)
pub const N_METH_STATES: usize = 2;

/// Number of SNPs available (A/T/C/G/N)
pub const N_BASES: usize = 7;

/// Flat matrix length for number of bases and methylation states
pub const LENGTH: usize = N_METH_STATES * N_BASES;

/// Base support
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone, Copy)]
pub enum Base {
    BaseR = 0,
    BaseY = 1,
    BaseA = 2,
    BaseC = 3,
    BaseG = 4,
    BaseT = 5,
    BaseN = 6,
}

impl Base {
    pub fn from(c: char) -> Result<Base, ()> {
        match c {
            'A' => Ok(Base::BaseA),
            'C' => Ok(Base::BaseC),
            'G' => Ok(Base::BaseG),
            'T' => Ok(Base::BaseT),
            'R' => Ok(Base::BaseR),
            'Y' => Ok(Base::BaseY),
            'N' => Ok(Base::BaseN),
            _   => Err(()),
        }
    }

    pub fn from_usize(i: usize) -> Result<Base, ()> {
        match i {
            0 => Ok(Base::BaseR),
            1 => Ok(Base::BaseY),
            2 => Ok(Base::BaseA),
            3 => Ok(Base::BaseC),
            4 => Ok(Base::BaseG),
            5 => Ok(Base::BaseT),
            6 => Ok(Base::BaseN),
            _ => Err(()),
        }
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = match self {
            Base::BaseA => "A".to_string(),
            Base::BaseC => "C".to_string(),
            Base::BaseG => "G".to_string(),
            Base::BaseT => "T".to_string(),
            Base::BaseR => "R".to_string(),
            Base::BaseY => "Y".to_string(),
            Base::BaseN => "N".to_string(),
        };

        write!(f, "{}", s)
    }
}

impl Into<usize> for Base {
    fn into(self) -> usize{
        self as usize
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
