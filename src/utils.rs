use crate::constants::{
    Base,
    CpgType,
    N_METH_STATES,
};

/// Parse region string
pub fn parse_region(r: &String) -> (String, u32, u32) {
    if r == "all" {
        return ("all".to_owned(), 0, u32::MAX)
    }

    if !r.contains(":") && !r.contains("-") {
        return (r.clone(), 0, u32::MAX);
    }

    let r_1: Vec<&str> = r.split(":").collect();
    if r_1.len() != 2 {
        eprintln!("Invalid region string ({}) given.", r);
        quit::with_code(1);
    }
    let chr: String = String::from(r_1[0]);

    let r_2: Vec<&str> = r_1[1].split("-").collect();
    if r_2.len() != 2 {
        eprintln!("Invalid region string ({}) given.", r);
        quit::with_code(1);
    }

    let start: u32 = r_2[0].parse::<u32>().unwrap();
    let end: u32   = r_2[1].parse::<u32>().unwrap();

    (chr, start, end)
}

/// Decode RLE string with missing 1's and starts with letter
pub fn decode_rle(s: &str) -> String {
    // Pull out lengths from string
    let l_tmp: Vec<&str> = s.split(char::is_alphabetic).collect();
    let lengths: Vec<&str> = l_tmp[1..]
        .iter()
        .map(|x| {
            if x.is_empty() {
                "1"
            } else {
                x
            }
        })
        .collect();

    // Pull out characters from string
    let tmp: String = s.chars().filter(|c| !c.is_digit(10)).collect();
    let letters: Vec<&str> = tmp.split("").filter(|c| !c.is_empty()).collect();

    // Form output string
    let mut out = String::new();
    for (i, l) in letters.iter().enumerate() {
        out = format!("{}{}", out, l.to_string().repeat(lengths[i].parse::<usize>().unwrap()));
    }

    out
}

/// Calculate index based on Base and CpgType
pub fn get_index(snp: Base, cpg: CpgType) -> usize {
    N_METH_STATES*<Base as Into<usize>>::into(snp) + <CpgType as Into<usize>>::into(cpg)
}

/// find sum of rows in flattened matrix
pub fn row_sums(fm: &[u32], nrow: usize, ncol: usize) -> Vec<u32> {
    let mut sums: Vec<u32> = vec![0; nrow];

    for r in 0..nrow {
        for c in 0..ncol {
            sums[r] += fm[ncol*r + c];
        }
    }

    sums
}

/// find sum of columns in flattened matrix
pub fn col_sums(fm: &[u32], nrow: usize, ncol: usize) -> Vec<u32> {
    let mut sums: Vec<u32> = vec![0; ncol];

    for r in 0..nrow {
        for c in 0..ncol {
            sums[c] += fm[ncol*r + c];
        }
    }

    sums
}

/// find indexes with the largest two values - takes the first occurrences of any matching largest
/// values
pub fn top_two(v: &Vec<u32>, n: usize) -> (usize, usize) {
    let (mut winner, mut first_loser): (usize, usize) = if v[0] > v[1] { (0, 1) } else { (1, 0) };

    if n > 2 {
        for i in 2..n {
            if v[i] > v[winner] {
                first_loser = winner;
                winner = i;
            } else if v[i] > v[first_loser] {
                first_loser = i;
            }
        }
    }

    (winner, first_loser)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_region_chr() {
        let test: String = "chr1".to_string();
        let p = parse_region(&test);

        assert_eq!(p, ("chr1".to_string(), 0, u32::MAX));
    }

    #[test]
    fn test_parse_region_full() {
        let test: String = "chr1:1-1000".to_string();
        let p = parse_region(&test);

        assert_eq!(p, ("chr1".to_string(), 1, 1000));
    }

    #[test]
    fn test_decode_rle() {
        let test = "A3BC3";
        let d = decode_rle(test);

        assert_eq!(d, "AAABCCC".to_string());
    }

    #[test]
    fn test_row_sums() {
        let matrix: Vec<u32> = vec![1, 2, 3, 4, 5, 0, 1, 2, 3, 4];
        let rs = row_sums(&matrix, 2, 5);

        assert_eq!(rs, [15, 10]);
    }

    #[test]
    fn test_col_sums() {
        let matrix: Vec<u32> = vec![1, 2, 3, 4, 5, 0, 1, 2, 3, 4];
        let cs = col_sums(&matrix, 2, 5);

        assert_eq!(cs, [1, 3, 5, 7, 9]);
    }

    #[test]
    fn test_top_two_input_2() {
        let x: Vec<u32> = vec![15, 10];
        let (one, two) = top_two(&x, 2);

        assert_eq!((one, two), (0, 1))
    }

    #[test]
    fn test_top_two_input_5() {
        let x: Vec<u32> = vec![1, 3, 5, 7, 9];
        let (one, two) = top_two(&x, 5);

        assert_eq!((one, two), (4, 3))
    }
}
