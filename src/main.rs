use std::{
    io::{BufRead, BufReader, Write},
    str::FromStr,
    collections::HashMap,
};

use anyhow::Result;
use clap::Parser;
use rust_htslib::bgzf;

mod utils;
mod sorting;
mod collapse;
mod constants;

mod records;
use crate::records::Record;

mod pairs;
use crate::pairs::{
    Snp,
    Cpg,
    Pair,
};

mod stats;
use crate::stats::{
    SnpCpgData,
};

#[derive(Parser, Debug)]
#[clap(name = "spasm")]
#[clap(author = "Jacob Morrison <jacob.morrison@vai.org>")]
#[clap(version = "version 1.0.0")]
struct Args {
    /// path to indexed reference FASTA
    genome: std::path::PathBuf,

    /// path to epiBED file
    path: std::path::PathBuf,

    /// region to extract (chr:start-end or chr)
    #[clap(short = 'g', long, default_value_t = String::from("all"))]
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
    #[clap(short = 'c', long, default_value_t = String::from("BH"))]
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
fn process_file(fh: &mut BufReader<bgzf::Reader>, args: &Args, chr: &String, start: &u64, end: &u64) -> Result<HashMap::<String, Vec<Record>>> {
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
                if args.sorted && r.get_chr() == chr && r.get_start() > end {
                    line.clear();
                    break;
                }

                // Skip records that don't fall in region (if requested)
                if chr != "all" && (r.get_chr() != chr || r.get_end() <= start || r.get_start() > end) {
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

/// Loop over records, collapose to fragment if requested, then loop over each vector entry to pull
/// out CpGs and SNPs to form input to ASM calculation
fn create_snp_cpg_pairs(fr: &HashMap::<String, Vec<Record>>, args: &Args, chr: &String, start: &u64, end: &u64) -> Vec<Pair> {
    let mut out: Vec<Pair> = Vec::new();

    for val in fr.values() {
        // Only want primary reads, but it's hard to tell if a read is primary or secondary without
        // the samtools FLAG, which we don't have in the epiBED format
        if val.len() > 2 {
            eprintln!("Read name with secondary reads found ({}). Skipping!", val[0].get_name());
            continue;
        }

        // Collapse to fragment if desired and able
        let process: Vec<Record> = if args.fragment && val.len() != 1 {
            collapse::collapse_to_fragment(val)
        } else {
            val.to_vec()
        };

        let mut ins_count: u64;
        let mut pos: u64;
        let mut snps: Vec<Snp> = Vec::new();
        let mut cpgs: Vec<Cpg> = Vec::new();
        for v in process {
            // Clear values for next read
            ins_count = 0;
            snps.clear();
            cpgs.clear();

            // Loop over characters in CpG decoded RLE string
            for (i, c) in v.get_cpg().chars().enumerate() {
                pos = v.get_start() + i as u64 - ins_count;

                if (chr != "all") && (v.get_chr() != chr || pos < *start || pos >= *end) {
                    continue;
                }

                match c {
                    'F' | 'x' | 'D' | 'P' => {
                        continue;
                    }
                    'a' | 'c' | 'g' | 't' | 'n' => {
                        ins_count += 1;
                        continue;
                    }
                    'A' | 'C' | 'G' | 'T' | 'N' => {
                        snps.push(Snp::new(v.get_chr().clone(), pos, c));
                    }
                    'M' | 'U' => {
                        cpgs.push(Cpg::new(v.get_chr().clone(), pos, c));
                    }
                    _ => {
                        continue;
                    }
                }
            }

            // Ignore reads with no SNPs in them
            if snps.len() == 0 {
                continue;
            }

            // Create SNP-CpG pairs
            for s in &snps {
                for c in &cpgs {
                    out.push(Pair::new(s.clone(), c.clone()));
                }
            }
        }
    }

    // Need the pairs to be sorted to correctly find ASM
    out.sort();

    out
}

fn find_p_values(locs: &Vec<Pair>) -> Vec<SnpCpgData> {
    let length = constants::N_METH_STATES * constants::N_SNP_STATES;

    let mut chrm = String::from("");
    let mut snp_prev: u64 = u64::MAX;
    let mut cpg_prev: u64 = u64::MAX;
    let mut snp_curr: u64 = u64::MAX;
    let mut cpg_curr: u64 = u64::MAX;
    let mut flat_matrix: Vec<i64> = vec![0; length];
    let mut index: usize;

    let mut out: Vec<SnpCpgData> = Vec::new();
    for l in locs {
        snp_curr = *l.get_snp().get_pos();
        cpg_curr = *l.get_cpg().get_pos();

        if chrm == "" || cpg_curr != cpg_prev || snp_curr != snp_prev || chrm != *l.get_snp().get_chr() {
            if chrm != "" {
                let p_vals = stats::calculate_p_values(&flat_matrix, constants::N_SNP_STATES, constants::N_METH_STATES);

                if !p_vals.is_none() {
                    out.push(
                        SnpCpgData::new(
                            chrm.clone(),
                            snp_prev,
                            cpg_prev,
                            p_vals.unwrap().2,
                            p_vals.unwrap().3
                        )
                    );
                }
            }

            chrm = l.get_snp().get_chr().clone();
            cpg_prev = cpg_curr;
            snp_prev = snp_curr;

            flat_matrix = vec![0; length];
        }

        index = l.get_index();
        flat_matrix[index] += 1;
    }

    // Catch remaining values
    if chrm != "" {
        let p_vals = stats::calculate_p_values(&flat_matrix, constants::N_SNP_STATES, constants::N_METH_STATES);

        if !p_vals.is_none() {
            out.push(
                SnpCpgData::new(
                    chrm.clone(),
                    snp_curr,
                    cpg_curr,
                    p_vals.unwrap().2,
                    p_vals.unwrap().3
                )
            );
        }
    }

    out
}

fn setup_output(fname: &Option<String>) -> bgzf::Writer {
    if fname.is_none() {
        return bgzf::Writer::from_stdout_with_compression(bgzf::CompressionLevel::NoCompression).unwrap();
    }

    let len  = fname.as_ref().unwrap().len();
    let is_gz = &fname.as_ref().unwrap()[len-3..] == ".gz";

    if is_gz {
        bgzf::Writer::from_path(&fname.as_ref().unwrap()).unwrap()
    } else {
        bgzf::Writer::from_path_with_level(&fname.as_ref().unwrap(), bgzf::CompressionLevel::NoCompression).unwrap()
    }
}

fn write_data(fh: &mut bgzf::Writer, data: &Vec<SnpCpgData>, is_biscuit: bool, cutoff: f64) -> () {
    for p in data.iter() {
        // Format string as requested
        let tmp: String = if is_biscuit {
            p.to_biscuit_asm()
        } else {
            p.to_bedpe(cutoff)
        };

        // Write to buffer
        fh.write(tmp.as_bytes()).unwrap();
    }
}

fn main() {
    // Command line arguments
    let args = Args::parse();

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

    // Read epiBED and put into records for processing
    let file_records = process_file(&mut reader, &args, &r_chr, &r_start, &r_end).expect("Error parsing file.");

    // Pull out matched SNP-CpG pairs from reads/fragments
    let locations = create_snp_cpg_pairs(&file_records, &args, &r_chr, &r_start, &r_end);

    // Find p-values from inputs
    let mut p_vals = find_p_values(&locations);

    // Perform p-value false discovery rate correction
    let n = p_vals.len();
    let p_corrected = stats::false_discovery_correction(&mut p_vals, &args.fdr, n);

    // Write data to output
    let mut writer = setup_output(&args.output);

    write_data(&mut writer, &p_corrected, args.biscuit, args.pcutoff);
}
