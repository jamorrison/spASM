use std::{
    io::{BufRead, BufReader, Write},
    collections::HashMap,
    path::PathBuf,
};

use anyhow::Result;
use clap::Parser;
use rust_htslib::{
    bgzf,
    tbx::{self, Read},
};

mod utils;
mod collapse;
mod constants;
use crate::constants::{
    Base,
    CpgType,
};
mod ref_genome;
mod snp;
mod fishers_exact_test;

mod records;
use crate::records::Record;

mod stats;
use crate::stats::{
    SnpCpgData,
};

#[derive(Parser, Debug)]
#[clap(name = "spasm")]
#[clap(author = "Jacob Morrison <jacob.morrison@vai.org>")]
#[clap(version = "version 1.1.0")]
struct Args {
    /// path to indexed reference FASTA
    genome: std::path::PathBuf,

    /// path to epiBED file
    path: std::path::PathBuf,

    /// region to extract (chr:start-end or chr)
    #[clap(short = 'g', long, default_value_t = String::from("all"))]
    region: String,

    /// do not merge mate reads together into a single DNA fragment
    #[clap(short, long, action)]
    no_mate_merging: bool,

    /// type of false discovery rate correction to perform
    /// possibilities:
    ///     BH (Benjamini-Hochberg),
    ///     BY (Benjamini-Yekutieli),
    ///     Bonferroni,
    ///     Hochberg,
    ///     Holm,
    ///     No (do not apply false discovery correction)
    #[clap(short = 'c', long, default_value_t = String::from("BH"), value_parser = stats::validate_fdr_type)]
    fdr: String,

    /// p-value significance cutoff
    #[clap(short, long, default_value_t = 0.05)]
    pcutoff: f64,

    /// output file name, compression level based on file name [stdout]
    #[clap(short, long)]
    output: Option<String>,

    /// write only candidate locations, FDR-corrected p-value based on all locations probed
    #[clap(short = 'O', long, action)]
    candidate: bool,

    /// write only locations with no ambiguous SNPs
    #[clap(short = 'N', long, action)]
    no_ambiguous: bool,

    /// write in BISCUIT ASM output format
    #[clap(short, long, action)]
    biscuit: bool,

    /// verbosity level (0: ERRORS ONLY | 1: WARNINGS + ERRORS | 2+: ALL)
    #[clap(short, long, default_value_t = 1)]
    verbose: usize,
}

fn process_file(fname: &PathBuf, genome: &PathBuf, k_chr: &HashMap::<String, u32>, k_int: &HashMap::<u32, String>, chr: &str, start: &u32, end: &u32, verbose: &usize) -> Result<(HashMap::<String, Vec<Record>>, HashMap::<String, [Option<char>; 2]>)> {
    let mut line    = String::new();
    let mut records = HashMap::<String, Vec<Record>>::new();
    let mut support = HashMap::<String, Vec<u16>>::new();

    let mut count: usize = 0;
    if chr == "all" {
        if *verbose >= 2 {
            eprintln!("NOTE: You requested that all reads be read. This may take a long time for large epiBEDs");
        }

        let file = match bgzf::Reader::from_path(fname) {
            Ok(f) => f,
            Err(_) => {
                eprintln!("Unable to open file: {}", fname.display());
                quit::with_code(1);
            },
        };
        let mut reader = BufReader::new(file);

        loop {
            match reader.read_line(&mut line) {
                Ok(bytes_read) => {
                    if bytes_read == 0 {
                        break;
                    }

                    // Strip newline
                    let len = line.len() - 1;
                    line.truncate(len);

                    // Parse record
                    let r = match Record::create(&line, &k_chr) {
                        Ok(rec) => rec,
                        Err(err) => return Err(anyhow::Error::from(err)),
                    };

                    snp::snp_support(&r, &mut support, k_int);

                    let name: String = r.get_name().clone();
                    if !records.contains_key(&name) {
                        records.insert(name, vec!(r));
                    } else {
                        records.get_mut(&name).unwrap().push(r);
                    }

                    line.clear();

                    count += 1;
                }

                Err(err) => {
                    return Err(anyhow::Error::from(err));
                }
            };
        }
    } else {
        let mut tbx_reader = match tbx::Reader::from_path(fname) {
            Ok(t) => t,
            Err(_) => {
                eprintln!("Unable to open file: {}", fname.display());
                quit::with_code(1);
            },
        };

        let tid = match tbx_reader.tid(&chr) {
            Ok(tid) => tid,
            Err(_) => {
                eprintln!("Could not find '{}' in {}", &chr, fname.display());
                quit::with_code(1);
            },
        };

        tbx_reader.fetch(tid, *start as u64, *end as u64).expect(&format!("Could not seek to {}:{}-{}", chr, start, end));

        for record in tbx_reader.records() {
            let rec = record.ok().unwrap();
            line = rec.iter().map(|b| *b as char).collect();

            let r = match Record::create(&line, &k_chr) {
                Ok(rec) => rec,
                Err(err) => return Err(anyhow::Error::from(err)),
            };

            snp::snp_support(&r, &mut support, k_int);

            let name: String = r.get_name().clone();
            if !records.contains_key(&name) {
                records.insert(name, vec!(r));
            } else {
                records.get_mut(&name).unwrap().push(r);
            }

            line.clear();

            count += 1;
        }
    }

    let dist_ambig = snp::redistribute_ambiguous_calls(&support, &genome);

    if *verbose >= 2 {
        eprintln!("NOTE: {} reads processed", count);
        eprintln!("NOTE: {} unique read names processed", records.len());
    }

    Ok((records, dist_ambig))
}

/// Loop over records, collapse mates if requested, then loop over each vector entry to pull out
/// CpGs and SNPs to form input to ASM calculation
fn create_snp_cpg_pairs(fr: HashMap::<String, Vec<Record>>, redist: HashMap::<String, [Option<char>; 2]>, merge: &bool, chr_id: &Option<&u32>, start: &u32, end: &u32, verbose: &usize) -> HashMap<(u32, u32, u32), [u32; constants::LENGTH]> {
    let mut out: HashMap<(u32, u32, u32), [u32; constants::LENGTH]> = HashMap::new();

    for val in fr.values() {
        // Only want primary reads, but it's hard to tell if a read is primary or secondary without
        // the samtools FLAG, which we don't have in the epiBED format
        if val.len() > 2 {
            if *verbose >= 1 {
                eprintln!("WARNING: Read name with secondary reads found ({}). Skipping!", val[0].get_name());
            }
            continue;
        }

        // Collapse to single fragment if desired and able
        let process: Vec<Record> = if *merge && val.len() != 1 {
            collapse::collapse_to_fragment(val)
        } else {
            val.to_vec()
        };

        let mut i: u32;
        let mut pos: u32;
        let mut chr: u32;
        let mut snps: Vec<(u32, Base)> = Vec::new();
        let mut cpgs: Vec<(u32, CpgType)> = Vec::new();
        for v in process {
            // Clear values for next read
            i = 0;
            snps.clear();
            cpgs.clear();
            chr = *v.get_chr_id();

            // Loop over characters in CpG and variant RLE strings
            for it in v.get_cpg().chars().zip(v.get_snp().chars()) {
                let (cg, vr) = it;
                pos = v.get_start() + i;

                i += 1;

                if !chr_id.is_none() && (&Some(v.get_chr_id()) != chr_id || pos < *start || pos >= *end) {
                    continue;
                }

                match cg {
                    'M' | 'U' => {
                        let cg_pos = if *v.get_bs_strand() { pos } else { pos-1 };
                        cpgs.push( (cg_pos, CpgType::from(cg).unwrap()) );
                    },
                    _ => {},
                }

                match vr {
                    'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'N' => {
                        let mut s: Base = Base::from(vr).unwrap();
                        match redist.get(&format!("{}:{}", v.get_chr_id(), pos)) {
                            Some(vec) => {
                                if vr == 'R' || vr == 'Y' {
                                    match vec[Base::from(vr).unwrap() as usize] {
                                        Some(c) => { s = Base::from(c).unwrap(); },
                                        None => {},
                                    };
                                }
                            },
                            None => {},
                        };

                        snps.push( (pos, s) );
                    },
                    _ => {},
                }
            }

            if cpgs.len() == 0 || snps.len() == 0 {
                continue;
            }

            for s in &snps {
                for c in &cpgs {
                    let tmp = (chr, s.0, c.0);
                    let idx = utils::get_index(s.1, c.1);

                    out.entry(tmp).and_modify(|fm| fm[idx] += 1).or_insert_with(|| { let mut fm = [0; constants::LENGTH]; fm[idx] += 1; fm } );
                }
            }
        }
    }

    out
}

fn find_p_values(locs: HashMap<(u32, u32, u32), [u32; constants::LENGTH]>, fdr: &String) -> Vec<SnpCpgData> {
    let mut out: Vec<SnpCpgData> = Vec::new();

    for (k, v) in locs {
        let p_vals = stats::calculate_p_values(&v, constants::N_BASES, constants::N_METH_STATES);

        if !p_vals.is_none() {
            out.push(
                SnpCpgData::new(
                    k.0,
                    k.1,
                    k.2,
                    p_vals.unwrap().0,
                    p_vals.unwrap().1
                )
            );
        }
    }

    // Perform p-value false discovery rate correction
    let n = out.len();
    stats::false_discovery_correction(&mut out, fdr, n);

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

fn write_data(fh: &mut bgzf::Writer, data: Vec<SnpCpgData>, k_int: HashMap::<u32, String>, is_biscuit: bool, cutoff: f64, candidate: bool, no_ambiguous: bool) -> () {
    for p in data.iter() {
        // Format string as requested
        let tmp: String = if is_biscuit {
            p.to_biscuit_asm(&k_int)
        } else {
            p.to_bedpe(&k_int, cutoff)
        };

        let print_val: bool = if no_ambiguous {
            !(*p.get_snp1() == Base::BaseR ||
              *p.get_snp1() == Base::BaseY ||
              *p.get_snp2() == Base::BaseR ||
              *p.get_snp2() == Base::BaseY)
        } else {
            true
        };

        // Write to buffer
        if print_val && (!candidate || (candidate && *p.get_p() < cutoff)) {
            fh.write(tmp.as_bytes()).unwrap();
        }
    }
}

#[quit::main]
fn main() {
    // Command line arguments
    // TODO: Add an option to read a BED file of desired locations
    let args = Args::parse();

    // Lookup tables for chromosome names and IDs
    let (k_chr, k_int) = ref_genome::chromosome_lookup_tables(&args.genome);

    // Chromosome, start, and end of region of interest
    let (r_chr, r_start, r_end) = utils::parse_region(&args.region);
    let r_chr_id = k_chr.get(&r_chr);

    // Read epiBED and put into records for processing
    let (file_records, redist) = match process_file(&args.path, &args.genome, &k_chr, &k_int, &r_chr, &r_start, &r_end, &args.verbose) {
        Ok((f, r)) => (f, r),
        Err(err) => {
            eprintln!("Error parsing file: {}", err);
            quit::with_code(1);
        },
    };

    if file_records.is_empty() {
        eprintln!("No records to process from epiBED file");
        quit::with_code(0);
    }

    // Pull out matched SNP-CpG pairs from reads/fragments
    // default is to merge mates, so when --no-mate-merging is given, the value is set to true
    // therefore, we want to take the opposite of what no_mate_merging to get the code to do what
    // we want it to
    let locations = create_snp_cpg_pairs(file_records, redist, &(!args.no_mate_merging), &r_chr_id, &r_start, &r_end, &args.verbose);

    if locations.is_empty() {
        eprintln!("No CpG-SNP pairs to process from epiBED file");
        quit::with_code(0);
    }

    // Find p-values and perform p-value false discovery rate correction from inputs
    let p_vals = find_p_values(locations, &args.fdr);

    // Write data to output
    let mut writer = setup_output(&args.output);

    write_data(&mut writer, p_vals, k_int, args.biscuit, args.pcutoff, args.candidate, args.no_ambiguous);
}
