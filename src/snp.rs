use std::{
    collections::HashMap,
    path::PathBuf,
};

use crate::utils;
use crate::ref_genome;
use crate::constants::{
    Base,
    N_BASES,
};
use crate::records::Record;

/// count the number of reads with a given base at each snp location
pub fn snp_support(r: &Record, support: &mut HashMap::<String, Vec<u16>>) {
    let mut pos: u32;
    for (i, c) in r.get_snp().chars().enumerate() {
        pos = r.get_start() + i as u32;
        match c {
            'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'N' => {
                let name = format!("{}:{}-{}", r.get_chr(), pos, pos);
                match support.get_mut(&name) {
                    Some(v) => {
                        // snp already in hashmap, so just increment
                        v[Base::from(c).unwrap() as usize] += 1;
                    },
                    None => {
                        // snp not found, so init, increment appropriate base, and push to hashmap
                        let mut vector: Vec<u16> = vec![0; N_BASES];
                        vector[Base::from(c).unwrap() as usize] += 1;

                        support.insert(name, vector);
                    },
                }
            },
            _ => {
                continue;
            }
        }
    }
}

/// find if ambiguous bases (R/Y) can be redistributed to their respective bases
///     (Y -> C/T, R -> G/A)
pub fn redistribute_ambiguous_calls(support: &HashMap<String, Vec<u16>>, ref_fn: &PathBuf) -> HashMap::<String, [Option<char>; 2]> {
    let mut out: HashMap::<String, [Option<char>; 2]> = HashMap::new();

    // Open reference
    let genome = ref_genome::open_ref(ref_fn);

    // loop over snps
    for (k, v) in support {
        let (chr, start, _end) = utils::parse_region(k);

        // if a snp exists, then we need to have a vector to look up
        //     None       -> do not change base
        //     Some(char) -> change ambiguous bases to char
        let name = format!("{}:{}", chr, start);
        out.insert(name.clone(), [None; 2]);

        // reference base, quit nicely if location can't be found
        let rb = match ref_genome::get_base(&genome, &chr, start as usize) {
            Some(c) => c,
            None => {
                eprintln!("Could not find base at location: {}:{}", chr, start);
                quit::with_code(1);
            },
        };

        // Only check for Y redistribution if Y's count is non-zero
        if v[Base::BaseY as usize] > 0 {
            if (rb == 'T' || v[Base::BaseT as usize] > 0) && v[Base::BaseC as usize] == 0 && rb != 'C' {
                let tmp = out.get_mut(&name).unwrap();
                tmp[Base::BaseY as usize] = Some('T');
            } else if (rb == 'C' || v[Base::BaseC as usize] > 0) && v[Base::BaseT as usize] == 0 && rb != 'T' {
                let tmp = out.get_mut(&name).unwrap();
                tmp[Base::BaseY as usize] = Some('C');
            }
        }

        // Only check for R redistribution if R's count is non-zero
        if v[Base::BaseR as usize] > 0 {
            if (rb == 'A' || v[Base::BaseA as usize] > 0) && v[Base::BaseG as usize] == 0 && rb != 'G' {
                let tmp = out.get_mut(&name).unwrap();
                tmp[Base::BaseR as usize] = Some('A');
            } else if (rb == 'G' || v[Base::BaseG as usize] > 0) && v[Base::BaseA as usize] == 0 && rb != 'A' {
                let tmp = out.get_mut(&name).unwrap();
                tmp[Base::BaseR as usize] = Some('G');
            }
        }
    }

    out
}
