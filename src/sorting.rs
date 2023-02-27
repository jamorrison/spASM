/// Returns indexes of original input in sorted order
/// reverse == true returns descending order
/// reverse == false returns ascending order
pub fn sorted_indexes<T: PartialOrd>(p: &Vec<T>, reverse: bool) -> Vec<usize> {
    let mut idxs = (0..p.len()).collect::<Vec<usize>>();
    idxs.sort_by(|a, b| p[*a].partial_cmp(&p[*b]).unwrap());

    if reverse {
        idxs.reverse();
    }

    idxs
}

/// Rank of sorted vector (reverses sorting back to original order)
pub fn unsort_sorted_indexes<T: PartialOrd>(p: &Vec<T>, reverse: bool) -> Vec<usize> {
    let idxs = sorted_indexes(&p, reverse);

    sorted_indexes(&idxs, false)
}

/// Sort vector of type generic that contains a float in it
pub fn sort_with_floats<T: PartialOrd + Clone>(v: &mut Vec<T>, reverse: bool) {
    if reverse {
        v.sort_by(|a, b| b.partial_cmp(a).unwrap())
    } else {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap())
    }
}

/// Resorts a vector based on provided ranks
pub fn resort<T: Clone>(v: &[T], ranks: &[usize]) -> Vec<T> {
    ranks.iter().map(|x| v[*x].clone()).collect()
}
