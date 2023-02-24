use crate::records::Record;

/// Collapse paired-reads to individual fragments
pub fn collapse_to_fragment(reads: &Vec<Record>) -> Vec<Record> {
    // Figure out which read in vector is read 1 and read 2
    let (idx1, idx2) = if reads[0].get_read_number() == &1 { (0, 1) } else { (1, 0) };

    let mut r1: Record = reads[idx1].clone();
    let     r2: Record = reads[idx2].clone();

    // If reads are on separate chromosomes, then return without processing
    if r1.get_chr() != r2.get_chr() {
        return Vec::from([r1, r2]);
    }

    // Figure out where reads are relative to one another for overlapping
    if r1.get_start() > r2.get_start() { // dovetail only
        return vec!(collapse_dovetail(r1, r2));
    } else if r1.get_start() < r2.get_start() {
        if r1.get_end() >= r2.get_end() { // read 1 entirely overlaps read 2
            r1.set_read_number(0);
            return vec!(r1);
        } else { // canonical overlap
            return vec!(collapse_canonical_proper_pair(r1, r2));
        }
    } else {
        if r1.get_end() >= r2.get_end() { // read 1 entirely overlaps read 2
            r1.set_read_number(0);
            return vec!(r1);
        } else { // canonical overlap
            return vec!(collapse_canonical_proper_pair(r1, r2));
        }
    }
}

/// Collapse dovetail reads
fn collapse_dovetail(r1: Record, r2: Record) -> Record {
    let     new_start: u32 = *r2.get_start();
    let mut new_end: u32 = *r1.get_end();
    let mut new_cpg;
    let mut new_snp;

    if r2.get_end() > r1.get_start() {
        // Difference in start locations
        let diff: usize = (r1.get_start() - r2.get_start()).try_into().unwrap();

        // Pull out substring of read 2 to tack on to read 1
        let r2_cpg = r2.get_cpg()[..diff].to_string();
        let r2_snp = r2.get_snp()[..diff].to_string();
        new_cpg    = format!("{}{}", r2_cpg, r1.get_cpg());
        new_snp    = format!("{}{}", r2_snp, r1.get_snp());

        // Handle case where read 2 starts before read 1 and ends after it
        if r2.get_end() > r1.get_end() {
            new_end = *r2.get_end();

            let tmp_start: usize = diff + r1.get_cpg().len();
            let r2_cpg           = r2.get_cpg()[tmp_start..].to_string();
            let r2_snp           = r2.get_snp()[tmp_start..].to_string();
            new_cpg              = format!("{}{}", new_cpg, r2_cpg);
            new_snp              = format!("{}{}", new_snp, r2_snp);
        }
    } else {
        // Difference in start locations
        let diff: usize = (r1.get_start() - r2.get_end()).try_into().unwrap();

        // Padding added between end of read 2 and read 1
        let pad: String = std::iter::repeat("x").take(diff).collect();

        new_cpg = format!("{}{}{}", r2.get_cpg(), pad, r1.get_cpg());
        new_snp = format!("{}{}{}", r2.get_snp(), pad, r1.get_snp());
    }

    if new_end - new_start != new_cpg.len() as u32 {
        eprintln!("Malformed collapsed fragment. (dovetail)",);
        eprintln!("Read 1: {}", r1);
        eprintln!("Read 2: {}", r2);
        eprintln!("Fragment length: {}, end-start: {}", new_cpg.len(), new_end-new_start);
        quit::with_code(1);
    }

    Record::new(
        r1.get_chr().to_string(),
        *r1.get_chr_id(),
        new_start,
        new_end,
        r1.get_name().to_string(),
        0,
        *r1.get_bs_strand(),
        new_cpg,
        new_snp,
    )
}

/// Collapse canonically-paired reads
fn collapse_canonical_proper_pair(r1: Record, r2: Record) -> Record {
    let     new_start: u32 = *r1.get_start();
    let     new_end: u32 = *r2.get_end();
    let     new_cpg;
    let     new_snp;

    if r2.get_start() > r1.get_end() {
        // Difference between end of read 1 and start of read 2
        let diff: usize = (r2.get_start() - r1.get_end()).try_into().unwrap();

        // Padding added between read 1 and read 2
        let pad: String = std::iter::repeat("x").take(diff).collect();

        new_cpg = format!("{}{}{}", r1.get_cpg(), pad, r2.get_cpg());
        new_snp = format!("{}{}{}", r1.get_snp(), pad, r2.get_snp());
    } else {
        let diff: usize = (r1.get_end() - r2.get_start()).try_into().unwrap();

        let r2_cpg = r2.get_cpg()[diff..].to_string();
        let r2_snp = r2.get_snp()[diff..].to_string();
        new_cpg    = format!("{}{}", r1.get_cpg(), r2_cpg);
        new_snp    = format!("{}{}", r1.get_snp(), r2_snp);
    }

    if new_end - new_start != new_cpg.len() as u32 {
        eprintln!("Malformed collapsed fragment (proper pair).",);
        eprintln!("Read 1: {}", r1);
        eprintln!("Read 2: {}", r2);
        eprintln!("Fragment length: {}, end-start: {}", new_cpg.len(), new_end-new_start);
        quit::with_code(1);
    }

    Record::new(
        r1.get_chr().to_string(),
        *r1.get_chr_id(),
        new_start,
        new_end,
        r1.get_name().to_string(),
        0,
        *r1.get_bs_strand(),
        new_cpg,
        new_snp,
    )
}
