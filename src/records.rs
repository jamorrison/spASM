use std::{
    str::FromStr,
    fmt,
};

use crate::utils;

/// Errors for parsing record line
#[derive(Debug, PartialEq)]
pub enum RecordParseError {
    /// Record is empty
    Empty,
    /// Record has an incorrect number of entries
    BadLen,
}

/// Store the values from each line in the epiBED
#[derive(Debug, PartialEq, Clone)]
pub struct Record {
    /// chromosome
    chr: String,
    /// 0-based start position
    start: u64,
    /// 1-based, non-inclusive end position
    end: u64,
    /// read name
    name: String,
    /// read number in paired-end sequencing (0 = fragment, 1 = read 1, 2 = read 2)
    read_number: u8,
    /// BS strand (+ for OT/CTOT, - for OB/CTOB)
    bs_strand: char,
    /// string for CpG methylation/SNPs
    cpg: String,
    /// string for GpC methylation/SNPs (None if not included)
    gpc: Option<String>,
}

impl Record {
    pub fn get_chr(&self) -> &String {
        &self.chr
    }

    pub fn get_start(&self) -> u64 {
        self.start
    }

    pub fn get_end(&self) -> u64 {
        self.end
    }

    pub fn get_name(&self) -> &String {
        &self.name
    }

    pub fn get_read_number(&self) -> u8 {
        self.read_number
    }

    pub fn get_bs_strand(&self) -> &char {
        &self.bs_strand
    }

    pub fn get_cpg(&self) -> &String {
        &self.cpg
    }

    pub fn get_gpc(&self) -> &Option<String> {
        &self.gpc
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s: String = format!(
            "Chromosome: {}\nStart: {}\nEnd: {}\nName: {}\nRead Number: {}\nBS Strand: {}\nCpG: {}",
            self.chr,
            self.start,
            self.end,
            self.name,
            self.read_number,
            self.bs_strand,
            self.cpg,
        );

        let gpc = match &self.gpc {
            Some(rle) => rle.clone(),
            None => String::from("NA"),
        };
        s = format!("{}\nGpC: {}\n", s, gpc);

        write!(f, "{}", s)
    }
}

impl FromStr for Record {
    type Err = RecordParseError;

    fn from_str(record: &str) -> Result<Self, Self::Err> {
        if record.len() == 0 {
            return Err(RecordParseError::Empty);
        }

        let vec: Vec<&str> = record.split("\t").collect();
        if vec.len() < 7 || vec.len() > 8 {
            return Err(RecordParseError::BadLen);
        }

        let tmp: Vec<char> = vec[5].chars().collect();

        let chr: String     = String::from(vec[0]);
        let start: u64      = vec[1].parse().unwrap();
        let end: u64        = vec[2].parse().unwrap();
        let name: String    = String::from(vec[3]);
        let read_number: u8 = vec[4].parse().unwrap();
        let bs_strand: char = tmp[0];
        let cpg: String     = utils::decode_rle(&String::from(vec[6]));

        let gpc = if vec.len() == 8 {
            Some(utils::decode_rle(&String::from(vec[7])))
        } else {
            None
        };

        Ok (
            Record {
                chr: chr,
                start: start,
                end: end,
                name: name,
                read_number: read_number,
                bs_strand: bs_strand,
                cpg: cpg,
                gpc: gpc,
            }
        )
    }
}
