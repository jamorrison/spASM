use std::{
    io::{BufRead, BufReader, Write},
    str::FromStr,
    collections::HashMap,
};

use anyhow::Result;
use clap::Parser;
use rust_htslib::bgzf;

mod utils;

mod records;
use crate::records::Record;

#[derive(Parser, Debug)]
#[clap(name = "cis_611_project")]
#[clap(author = "J. Morrison")]
#[clap(version = "version 1.1")]
struct Cli {
    /// path to epiBED file
    #[clap(parse(from_os_str))]
    path: std::path::PathBuf,

    /// region to extract (chr:start-end or chr)
    #[clap(short, long, default_value_t = String::from("all"))]
    region: String,

    /// collapse reads to fragment-level
    #[clap(short, long, action)]
    fragment: bool,

    /// assume epiBED is sorted (lexicographically by chr, then start, then end)
    #[clap(short, long, action)]
    sorted: bool,

    /// type of false discovery rate correction to perform
    /// possibilities:
    ///     BH (Benjamini-Hochberg),
    ///     BY (Benjamini-Yekutieli),
    ///     Bonferroni,
    ///     Hochberg,
    ///     Holm,
    ///     No (do not apply false discovery correction)
    #[clap(long, default_value_t = String::from("BH"))]
    fdr: String,

    /// p-value significance cutoff
    #[clap(short, long, default_value_t = 0.05)]
    pcutoff: f64,

    /// output file name, compression level based on file name (stdout if nothing given)
    #[clap(short, long)]
    output: Option<String>,

    /// write in BISCUIT ASM output format
    #[clap(short, long, action)]
    biscuit: bool,
}

/// Read file and put results in a HashMap keyed on the read name
fn process_file(fh: &mut BufReader<bgzf::Reader>, args: &Cli, chr: &String, start: &u64, end: &u64) -> Result<HashMap::<String, Vec<Record>>> {
    let mut line    = String::new();
    let mut records = HashMap::<String, Vec<Record>>::new();

    loop {
        match fh.read_line(&mut line) {
            Ok(bytes_read) => {
                if bytes_read == 0 {
                    break;
                }

                // Strip off newline
                let len = line.len() - 1;
                line.truncate(len);

                // Parse record
                let r = match Record::from_str(&line) {
                    Ok(rec) => rec,
                    Err(_) => {
                        eprintln!("Error found when reading record!");
                        quit::with_code(1);
                    },
                };

                if !r.get_gpc().is_none() {
                    eprintln!("This tool doesn't work with NOMe-seq data. Sorry!");
                    quit::with_code(1);
                }

                // TODO: Could make this tabix, but then would have to figure out how to make that
                // work, and that's not the easiest thing to do
                // If input is sorted and we've made it past the end of the region, then we can
                // break out of out loop
                if args.sorted && r.get_chr() == chr && r.get_start() > *end {
                    line.clear();
                    break;
                }

                // Skip records that don't fall in region (if requested)
                if chr != "all" && (r.get_chr() != chr || r.get_end() <= *start || r.get_start() > *end) {
                    line.clear();
                    continue;
                }

                let name: String = r.get_name().clone();
                if !records.contains_key(&name) {
                    records.insert(name, vec!(r));
                } else {
                    records.get_mut(&name).unwrap().push(r);
                }

                line.clear();
            }

            Err(err) => {
                return Err(anyhow::Error::from(err));
            }
        };
    }

    Ok(records)
}

fn main() {
    // Command line arguments
    let args = Cli::parse();

    // Chromosome, start, and end of region of interest
    let (r_chr, r_start, r_end) = utils::parse_region(&args.region);

    // Check for file and create reader
    let file = match bgzf::Reader::from_path(&args.path) {
        Ok(f) => f,
        Err(err) => {
            eprintln!("Unable to find file: {:?}", err);
            quit::with_code(1);
        },
    };
    let mut reader = BufReader::new(file);

    let file_records = process_file(&mut reader, &args, &r_chr, &r_start, &r_end).expect("Error parsing file.");

    for fr in &file_records {
        println!("{:?}", fr);
    }
}
