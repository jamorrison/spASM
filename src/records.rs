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
    start: u32,
    /// 1-based, non-inclusive end position
    end: u32,
    /// read name
    name: String,
    /// read number in paired-end sequencing (0 = fragment, 1 = read 1, 2 = read 2)
    read_number: u8,
    /// BS strand (true for OT/CTOT, false for OB/CTOB)
    bs_strand: bool,
    /// string for CpG methylation
    cpg: String,
    /// string for SNPs and other variants
    snp: String,
}

impl Record {
    pub fn new(chr: String, start: u32, end: u32, name: String, rn: u8, bss: bool, cpg: String, snp: String) -> Record {
        Record {
            chr: chr,
            start: start,
            end: end,
            name: name,
            read_number: rn,
            bs_strand: bss,
            cpg: cpg,
            snp: snp
        }
    }

    pub fn get_chr(&self) -> &String {
        &self.chr
    }

    pub fn get_start(&self) -> &u32 {
        &self.start
    }

    pub fn get_end(&self) -> &u32 {
        &self.end
    }

    pub fn get_name(&self) -> &String {
        &self.name
    }

    pub fn get_read_number(&self) -> &u8 {
        &self.read_number
    }

    pub fn set_read_number(&mut self, n: u8) {
        self.read_number = n;
    }

    pub fn get_bs_strand(&self) -> &bool {
        &self.bs_strand
    }

    pub fn get_cpg(&self) -> &String {
        &self.cpg
    }

    pub fn get_snp(&self) -> &String {
        &self.snp
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: String = format!(
            "Chromosome: {}\nStart: {}\nEnd: {}\nName: {}\nRead Number: {}\nBS Strand: {}\nCpG: {}\nSNP: {}",
            self.chr,
            self.start,
            self.end,
            self.name,
            self.read_number,
            self.bs_strand,
            self.cpg,
            self.snp,
        );

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
        if vec.len() != 9 {
            return Err(RecordParseError::BadLen);
        }

        // TODO: Remove insert values (acgtni) from cpg and snp
        // let s = s.replace(&['a', 'c', 'g', 't', 'n', 'i'][..], "");
        let tmp: Vec<char> = vec[5].chars().collect();

        let chr: String     = String::from(vec[0]);
        let start: u32      = vec[1].parse().unwrap();
        let end: u32        = vec[2].parse().unwrap();
        let name: String    = String::from(vec[3]);
        let read_number: u8 = vec[4].parse().unwrap();
        let bs_strand: bool = if tmp[0] == '+' { true } else { false };
        let cpg: String     = utils::decode_rle(&String::from(vec[6])).replace(&['a', 'c', 'g', 't', 'n', 'i'][..], "");
        let snp: String     = utils::decode_rle(&String::from(vec[8])).replace(&['a', 'c', 'g', 't', 'n', 'i'][..], "");

        Ok (
            Record {
                chr: chr,
                start: start,
                end: end,
                name: name,
                read_number: read_number,
                bs_strand: bs_strand,
                cpg: cpg,
                snp: snp,
            }
        )
    }
}
