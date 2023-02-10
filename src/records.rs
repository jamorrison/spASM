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
    /// string for CpG methylation
    cpg: String,
    /// string for GpC methylation (None if not included)
    gpc: Option<String>,
    /// string for SNPs and other variants
    snp: String,
}

impl Record {
    pub fn new(chr: String, start: u64, end: u64, name: String, rn: u8, bss: char, cpg: String, gpc: Option<String>, snp: String) -> Record {
        Record {
            chr: chr,
            start: start,
            end: end,
            name: name,
            read_number: rn,
            bs_strand: bss,
            cpg: cpg,
            gpc: gpc,
            snp: snp
        }
    }

    pub fn get_chr(&self) -> &String {
        &self.chr
    }

    pub fn get_start(&self) -> &u64 {
        &self.start
    }

    pub fn get_end(&self) -> &u64 {
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

    pub fn get_bs_strand(&self) -> &char {
        &self.bs_strand
    }

    pub fn get_cpg(&self) -> &String {
        &self.cpg
    }

    pub fn get_gpc(&self) -> &Option<String> {
        &self.gpc
    }

    pub fn get_snp(&self) -> &String {
        &self.snp
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
        s = format!("{}\nGpC: {}", s, gpc);

        s = format!("{}\nSNP: {}\n", s, self.snp);

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

        // TODO: Remove insert values (acgtni) from cpg, gpc, and snp
        // let s = s.replace(&['a', 'c', 'g', 't', 'n', 'i'][..], "");
        let tmp: Vec<char> = vec[5].chars().collect();

        let chr: String     = String::from(vec[0]);
        let start: u64      = vec[1].parse().unwrap();
        let end: u64        = vec[2].parse().unwrap();
        let name: String    = String::from(vec[3]);
        let read_number: u8 = vec[4].parse().unwrap();
        let bs_strand: char = tmp[0];
        let cpg: String     = utils::decode_rle(&String::from(vec[6]));
        let snp: String     = utils::decode_rle(&String::from(vec[8]));

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
                snp: snp,
            }
        )
    }
}
