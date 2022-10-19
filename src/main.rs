use std::{
    io::{BufRead, BufReader, Write},
    str::FromStr,
    collections::HashMap,
};

use anyhow::Result;
use clap::Parser;
use rust_htslib::bgzf;

mod utils;

mod stats;

mod collapse;

mod records;
use crate::records::Record;

mod pairs;
use crate::pairs::{
    SnpType,
    Snp,
    CpgType,
    Cpg,
    Pair,
    PairP,
};

mod sorting;

#[derive(Parser)]
#[clap(name = "cis_611_project")]
#[clap(author = "J. Morrison")]
#[clap(version = "version 1.0")]
struct Cli {
    /// path to epiBED file
    #[clap(parse(from_os_str))]
    path: std::path::PathBuf,

    /// region to extract
    #[clap(short, long, default_value_t = String::from("all"))]
    region: String,

    /// collapse reads to fragment-level
    #[clap(short, long, action)]
    fragment: bool,

    /// assume epiBED is sorted (lexicographically by chr, then start, then end)
    #[clap(short, long, action)]
    sorted: bool,

    /// type of false discovery rate correction to perform
    #[clap(long, default_value_t = String::from("BY"))]
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

/// Read file and put results into a HashMap based on read name
fn process_file(fh: &mut BufReader<bgzf::Reader>, args: &Cli, r_chr: &String, r_start: &u64, r_end: &u64) -> Result<HashMap::<String, Vec<Record>>> {
    let mut line    = String::new();
    let mut records = HashMap::<String, Vec<Record>>::new();

    loop {
        match fh.read_line(&mut line) {
            Ok(bytes_read) => {
                if bytes_read == 0 {
                    break;
                }

                // Strip off newline
                let len = line.len()-1;
                line.truncate(len);

                // Parse record
                let record = match Record::from_str(&line) {
                    Ok(rec) => rec,
                    Err(_) => {
                        eprintln!("Error found when reading record!");
                        quit::with_code(1);
                    },
                };

                if !record.gpc.is_none() {
                    eprintln!("This tool doesn't work with NOMe-seq data. Sorry!");
                    quit::with_code(1);
                }

                // TODO: Could make this tabix, but then would have to figure out how to make that
                // work, and that's not the easiest thing to do
                // If input is sorted and we've made it past the end of the region, then we can
                // break out of out loop
                if args.sorted && record.chr == *r_chr && record.start > *r_end {
                    line.clear();
                    break;
                }

                // Skip records that don't fall in region (if requested)
                if (r_chr != "all") && (record.chr != *r_chr || record.end <= *r_start || record.start > *r_end) {
                    line.clear();
                    continue;
                }

                let name: String = record.name.clone();
                if !records.contains_key(&record.name) {
                    records.insert(name, vec!(record));
                } else {
                    records.get_mut(&name).unwrap().push(record);
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

/// Loop over file_records, collapse to fragment if requested, then loop over each vector
/// entry to pull out CpGs and SNPs to form input to ASM calculation
fn create_snp_cpg_pairs(fr: &HashMap::<String, Vec<Record>>, args: &Cli, r_chr: &String, r_start: &u64, r_end: &u64) -> Vec<Pair> {
    let mut out: Vec<Pair> = Vec::new();

    for val in fr.values() {
        // We only want primary reads, but it's hard to tell if a read is secondary or primary
        // without the flag, which we don't have in the epiBED format
        if val.len() > 2 {
            eprintln!("Read name with secondary reads found ({}). Skipping!", val[0].name);
            continue;
        }

        // Collapse to fragment if desired and able
        let process = if args.fragment && val.len() != 1 {
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
            for (i, c) in v.cpg.chars().enumerate() {
                pos = v.start + i as u64 - ins_count;

                // Skip bases that don't fall in requested region
                if (r_chr != "all") && (v.chr != *r_chr || pos < *r_start || pos >= *r_end) {
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
                        snps.push(
                            Snp {
                                chr: v.chr.clone(),
                                pos: pos,
                                typ: SnpType::from(c).unwrap(),
                            }
                        );
                    }
                    'M' | 'U' => {
                        cpgs.push(
                            Cpg {
                                chr: v.chr.clone(),
                                pos: pos,
                                typ: CpgType::from(c).unwrap(),
                            }
                        );
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
                    out.push( Pair { snp: s.clone(), cpg: c.clone() } );
                }
            }
        }
    }

    // Need the pairs to be sorted to correctly find ASM
    out.sort();

    out
}

fn find_p_values(locs: &Vec<Pair>) -> Vec<PairP> {
    let mut chrm = String::from("");
    let mut snp_prev: u64 = u64::MAX;
    let mut cpg_prev: u64 = u64::MAX;
    let mut snp_curr: u64 = u64::MAX;
    let mut cpg_curr: u64 = u64::MAX;
    let mut flat_matrix: Vec<i64> = vec![0; 10];
    let mut index: usize;

    let mut out: Vec<PairP> = Vec::new();
    for l in locs {
        snp_curr = l.snp.pos;
        cpg_curr = l.cpg.pos;

        if chrm == "" || cpg_curr != cpg_prev || snp_curr != snp_prev || chrm != l.snp.chr {
            if chrm != "" {
                let p_vals = stats::calculate_p_values(&flat_matrix, 5, 2);

                if !p_vals.is_none() {
                    out.push(
                        PairP {
                            chr: chrm.clone(),
                            snp_pos: snp_prev,
                            cpg_pos: cpg_prev,
                            p: p_vals.unwrap().2,
                            pdata: p_vals.unwrap().3
                        }
                    );
                }
            }

            chrm = l.snp.chr.clone();
            cpg_prev = cpg_curr;
            snp_prev = snp_curr;

            flat_matrix = vec![0; 10];
        }

        index = l.get_index();
        flat_matrix[index] += 1;
    }
    
    // Catch remaining values
    if chrm != "" {
        let p_vals = stats::calculate_p_values(&flat_matrix, 5, 2);

        if !p_vals.is_none() {
            out.push(
                PairP {
                    chr: chrm.clone(),
                    snp_pos: snp_curr,
                    cpg_pos: cpg_curr,
                    p: p_vals.unwrap().2,
                    pdata: p_vals.unwrap().3
                }
            );
        }
    }

    out
}

fn setup_output(fname: &Option<String>) -> bgzf::Writer {
    if fname.is_none() {
        return bgzf::Writer::from_stdout_with_compression(bgzf::CompressionLevel::NoCompression).unwrap();
    }

    let filename: String = fname.clone().unwrap();

    let len  = filename.len();
    let is_gz = &filename[len-3..] == ".gz";

    if is_gz {
        bgzf::Writer::from_path(&filename).unwrap()
    } else {
        bgzf::Writer::from_path_with_level(&filename, bgzf::CompressionLevel::NoCompression).unwrap()
    }
}

fn write_data(fh: &mut bgzf::Writer, data: &Vec<PairP>, is_biscuit: bool, cutoff: f64) -> () {
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
    let args = Cli::parse();

    // Chromosome, start, and end of region of interest
    let (r_chr, r_start, r_end) = utils::parse_region(&args.region);

    // Prep file to be read
    let     file   = bgzf::Reader::from_path(&args.path).unwrap();
    let mut reader = BufReader::new(file);

    // Read epiBED file and put into records for processing
    let file_records = process_file(&mut reader, &args, &r_chr, &r_start, &r_end).expect("Error parsing file");

    // Pull out matched SNP-CpG pairs from reads/fragments
    let locations = create_snp_cpg_pairs(&file_records, &args, &r_chr, &r_start, &r_end);

    // Find p-values from inputs
    let mut p_vals = find_p_values(&locations);

    // Perform p-value false discovery rate correction
    let n = p_vals.len();
    let mut fdr = stats::FdrType::from_str(&args.fdr).unwrap();
    let p_corrected = stats::false_discovery_correction(&mut p_vals, &mut fdr, n);

    // Write data to output
    let mut writer = setup_output(&args.output);

    write_data(&mut writer, &p_corrected, args.biscuit, args.pcutoff);
}
