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
mod sorting;
mod collapse;
mod constants;
use crate::constants::{
    Base,
};
mod ref_genome;
mod snp;

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
#[clap(version = "version 1.1.0")]
struct Args {
    /// path to indexed reference FASTA
    genome: std::path::PathBuf,

    /// path to epiBED file
    path: std::path::PathBuf,

    /// region to extract (chr:start-end or chr)
    #[clap(short = 'g', long, default_value_t = String::from("all"))]
    region: String,

    /// read mates come from the same original DNA fragment, therefore they represent the same
    /// "epi-haplotype". To recover this correlation, mate reads are merged into a single fragment.
    /// However, there may be cases where mate reads need to be treated individually. This option
    /// provides the ability to avoid merging mate reads.
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
    #[clap(short = 'c', long, default_value_t = String::from("BH"))]
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

        // TODO: Count number of reads and print the number under DEBUG verbose
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
fn create_snp_cpg_pairs(fr: HashMap::<String, Vec<Record>>, redist: HashMap::<String, [Option<char>; 2]>, merge: &bool, chr_id: &Option<&u32>, start: &u32, end: &u32, verbose: &usize) -> Vec<Pair> {
    let mut out: Vec<Pair> = Vec::new();

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
        let mut snps: Vec<Snp> = Vec::new();
        let mut cpgs: Vec<Cpg> = Vec::new();
        for v in process {
            // Clear values for next read
            i = 0;
            snps.clear();
            cpgs.clear();

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
                        cpgs.push(Cpg::new(*v.get_chr_id(), cg_pos, cg));
                    },
                    _ => {}
                }

                match vr {
                    'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'N' => {
                        let mut s: char = vr;
                        match redist.get(&format!("{}:{}", v.get_chr_id(), pos)) {
                            Some(vec) => {
                                if vr == 'R' || vr == 'Y' {
                                    match vec[Base::from(vr).unwrap() as usize] {
                                        Some(c) => { s = c; },
                                        None => {},
                                    };
                                }
                            },
                            None => {},
                        };

                        snps.push(Snp::new(*v.get_chr_id(), pos, s));
                    },
                    _ => {
                        continue;
                    }
                }
            }

            if cpgs.len() == 0 || snps.len() == 0 {
                continue;
            }

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

fn find_p_values(locs: Vec<Pair>) -> Vec<SnpCpgData> {
    const LENGTH: usize = constants::N_METH_STATES * constants::N_BASES;

    let mut chrm: Option<u32> = None;
    let mut snp_prev: u32 = u32::MAX;
    let mut cpg_prev: u32 = u32::MAX;
    let mut snp_curr: u32 = u32::MAX;
    let mut cpg_curr: u32 = u32::MAX;
    let mut flat_matrix = [0; LENGTH];
    let mut index: usize;

    let mut out: Vec<SnpCpgData> = Vec::new();
    for l in locs {
        snp_curr = *l.get_snp().get_pos();
        cpg_curr = *l.get_cpg().get_pos();

        if chrm.is_none() || cpg_curr != cpg_prev || snp_curr != snp_prev || chrm != Some(*l.get_snp().get_chr()) {
            if !chrm.is_none() {
                let p_vals = stats::calculate_p_values(&flat_matrix, constants::N_BASES, constants::N_METH_STATES);

                if !p_vals.is_none() {
                    out.push(
                        SnpCpgData::new(
                            chrm.unwrap(),
                            snp_prev,
                            cpg_prev,
                            p_vals.unwrap().0,
                            p_vals.unwrap().1
                        )
                    );
                }
            }

            chrm = Some(*l.get_snp().get_chr());
            cpg_prev = cpg_curr;
            snp_prev = snp_curr;

            flat_matrix = [0; LENGTH];
        }

        index = l.get_index();
        flat_matrix[index] += 1;
    }

    // Catch remaining values
    if !chrm.is_none() {
        let p_vals = stats::calculate_p_values(&flat_matrix, constants::N_BASES, constants::N_METH_STATES);

        if !p_vals.is_none() {
            out.push(
                SnpCpgData::new(
                    chrm.unwrap(),
                    snp_curr,
                    cpg_curr,
                    p_vals.unwrap().0,
                    p_vals.unwrap().1
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

fn write_data(fh: &mut bgzf::Writer, data: &Vec<SnpCpgData>, k_int: &HashMap::<u32, String>, is_biscuit: bool, cutoff: f64, candidate: bool, no_ambiguous: bool) -> () {
    for p in data.iter() {
        // Format string as requested
        let tmp: String = if is_biscuit {
            p.to_biscuit_asm(k_int)
        } else {
            p.to_bedpe(k_int, cutoff)
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

    // Pull out matched SNP-CpG pairs from reads/fragments
    // default is to merge mates, so when --no-mate-merging is given, the value is set to true
    // therefore, we want to take the opposite of what no_mate_merging to get the code to do what
    // we want it to
    let locations = create_snp_cpg_pairs(file_records, redist, &(!args.no_mate_merging), &r_chr_id, &r_start, &r_end, &args.verbose);

    // Find p-values from inputs
    let mut p_vals = find_p_values(locations);

    // Perform p-value false discovery rate correction
    let n = p_vals.len();
    // TODO: I think I can do this p-correction in place, which would save on some memory
    let p_corrected = stats::false_discovery_correction(&mut p_vals, &args.fdr, n);

    // Write data to output
    let mut writer = setup_output(&args.output);

    write_data(&mut writer, &p_corrected, &k_int, args.biscuit, args.pcutoff, args.candidate, args.no_ambiguous);
}
