use crate::records::Record;

/// Collapse paired-reads to individual fragments
pub fn collapse_to_fragment(reads: &Vec<Record>) -> Vec<Record> {
    // Figure out which read in vector is read 1 and read 2
    let idx1 = if reads[0].read_number == 1 { 0 } else { 1 };
    let idx2 = if reads[0].read_number == 2 { 0 } else { 1 };

    let mut r1: Record = reads[idx1].clone();
    let     r2: Record = reads[idx2].clone();

    // If reads are on separate chromosomes, then return without processing
    if r1.chr != r2.chr {
        return Vec::from([r1, r2]);
    }

    // Figure out where reads are relative to one another for overlapping
    if r1.start > r2.start { // dovetail only
        return vec!(collapse_dovetail(r1, r2));
    } else if r1.start < r2.start {
        if r1.end >= r2.end { // read 1 entirely overlaps read 2
            r1.read_number = 0;
            return vec!(r1);
        } else { // canonical overlap
            return vec!(collapse_canonical_proper_pair(r1, r2));
        }
    } else {
        if r1.end >= r2.end { // read 1 entirely overlaps read 2
            r1.read_number = 0;
            return vec!(r1);
        } else { // canonical overlap
            return vec!(collapse_canonical_proper_pair(r1, r2));
        }
    }
}

/// Collapse dovetail reads
fn collapse_dovetail(r1: Record, r2: Record) -> Record {
    let     new_start: u64 = r2.start;
    let mut new_end: u64 = r1.end;
    let mut new_cpg;
    let mut new_gpc: Option<String> = None;

    if r2.end > r1.start {
        // Difference in start locations
        let diff: usize = (r1.start - r2.start).try_into().unwrap();

        // Pull out substring of read 2 to tack on to read 1
        let r2_cpg = &r2.cpg[..diff].to_string();
        new_cpg    = format!("{}{}", r2_cpg, r1.cpg);

        // Set GpC string if it exists
        if !r1.gpc.is_none() {
            let r1_gpc = r1.gpc.as_ref().unwrap();
            let r2_gpc = r2.gpc.as_ref().unwrap();

            let tmp = &r2_gpc[..diff].to_string();
            new_gpc = Some(format!("{}{}", tmp, r1_gpc));
        }

        // Handle case where read 2 starts before read 1 and ends after it
        if r2.end > r1.end {
            new_end = r2.end;

            let tmp_start: usize = diff + r1.cpg.len();
            let r2_cpg           = &r2.cpg[tmp_start..].to_string();
            new_cpg              = format!("{}{}", new_cpg, r2_cpg);

            if !r1.gpc.is_none() {
                let r2_gpc = r2.gpc.as_ref().unwrap();
                let tmp    = &r2_gpc[tmp_start..].to_string();
                new_gpc = Some(format!("{}{}", new_gpc.unwrap(), tmp));
            }
        }
    } else {
        // Difference in start locations
        let diff: usize = (r1.start - r2.end).try_into().unwrap();

        // Padding added between end of read 2 and read 1
        let pad: String = std::iter::repeat("x").take(diff).collect();

        new_cpg = format!("{}{}{}", r2.cpg, pad, r1.cpg);

        // Handle GpC if it exists
        if !r1.gpc.is_none() {
            let r1_gpc = r1.gpc.as_ref().unwrap();
            let r2_gpc = r2.gpc.as_ref().unwrap();
            new_gpc = Some(format!("{}{}{}", r2_gpc, pad, r1_gpc));
        }
    }

    if new_end - new_start != new_cpg.len() as u64 {
        eprintln!("Malformed collapsed fragment.",);
        eprintln!("Read 1: {}", r1);
        eprintln!("Read 2: {}", r2);
        quit::with_code(1);
    }

    Record {
        chr: r1.chr,
        start: new_start,
        end: new_end,
        name: r1.name,
        read_number: 0,
        bs_strand: r1.bs_strand,
        cpg: new_cpg,
        gpc: new_gpc,
    }
}

/// Collapse canonically-paired reads
fn collapse_canonical_proper_pair(r1: Record, r2: Record) -> Record {
    let     new_start: u64 = r1.start;
    let     new_end: u64 = r2.end;
    let     new_cpg;
    let mut new_gpc: Option<String> = None;

    if r2.start > r1.end {
        // Difference between end of read 1 and start of read 2
        let diff: usize = (r2.start - r1.end).try_into().unwrap();

        // Padding added between read 1 and read 2
        let pad: String = std::iter::repeat("x").take(diff).collect();

        new_cpg = format!("{}{}{}", r1.cpg, pad, r2.cpg);

        // Handle GpC if it exists
        if !r1.gpc.is_none() {
            let r1_gpc = r1.gpc.as_ref().unwrap();
            let r2_gpc = r2.gpc.as_ref().unwrap();
            new_gpc = Some(format!("{}{}{}", r1_gpc, pad, r2_gpc));
        }
    } else {
        let diff: usize = (r1.end - r2.start).try_into().unwrap();

        let r2_cpg = &r2.cpg[diff..].to_string();
        new_cpg    = format!("{}{}", r1.cpg, r2_cpg);

        // Handle GpC if it exists
        if !r1.gpc.is_none() {
            let r1_gpc = r1.gpc.as_ref().unwrap();
            let r2_gpc = r2.gpc.as_ref().unwrap();

            let tmp = &r2_gpc[diff..].to_string();
            new_gpc = Some(format!("{}{}", r1_gpc, tmp));
        }
    }

    if new_end - new_start != new_cpg.len() as u64 {
        eprintln!("Malformed collapsed fragment.",);
        eprintln!("Read 1: {}", r1);
        eprintln!("Read 2: {}", r2);
        quit::with_code(1);
    }

    Record {
        chr: r1.chr,
        start: new_start,
        end: new_end,
        name: r1.name,
        read_number: 0,
        bs_strand: r1.bs_strand,
        cpg: new_cpg,
        gpc: new_gpc,
    }
}

