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
pub fn sort_with_floats<T: PartialOrd + Clone>(v: &Vec<T>, reverse: bool) -> Vec<T> {
    let mut out = v.to_vec();
    if reverse {
        out.sort_by(|a, b| b.partial_cmp(a).unwrap());
    } else {
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    out
}

/// Resorts a vector based on provided ranks
pub fn resort<T: Clone>(v: &[T], ranks: &[usize]) -> Vec<T> {
    ranks.iter().map(|x| v[*x].clone()).collect()
}

//#[cfg(test)]
//mod tests {
//    use super::*;
//    use crate::pairs::{
//        SnpType,
//        CpgType,
//    };
//    use crate::stats::{
//        SnpCpgData,
//        PValMetadata,
//    };
//
//    #[test]
//    fn test_sorted_indexes_forward() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let sort: Vec<usize> = sorted_indexes(&test, false);
//
//        assert_eq!(sort, vec![0, 3, 2, 1]);
//    }
//
//    #[test]
//    fn test_sorted_indexes_reverse() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let sort: Vec<usize> = sorted_indexes(&test, true);
//
//        assert_eq!(sort, vec![1, 2, 3, 0]);
//    }
//
//    #[test]
//    fn test_sorted_ranks_forward() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let sort: Vec<usize> = unsort_sorted_indexes(&test, false);
//
//        assert_eq!(sort, vec![0, 3, 2, 1]);
//    }
//
//    #[test]
//    fn test_sorted_ranks_reverse() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let sort: Vec<usize> = unsort_sorted_indexes(&test, true);
//
//        assert_eq!(sort, vec![3, 0, 1, 2]);
//    }
//
//    #[test]
//    fn test_sort_with_floats_forward() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let compare: Vec<SnpCpgData> = vec![
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//        ];
//
//        let sort: Vec<SnpCpgData> = sort_with_floats(&test, false);
//
//        assert_eq!(sort, compare);
//    }
//
//    #[test]
//    fn test_sort_with_floats_reverse() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let compare: Vec<SnpCpgData> = vec![
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//        ];
//
//        let sort: Vec<SnpCpgData> = sort_with_floats(&test, true);
//
//        assert_eq!(sort, compare);
//    }
//
//    #[test]
//    fn test_resort() {
//        let mut test: Vec<SnpCpgData> = Vec::new();
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//        test.push(
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            )
//        );
//
//        let compare: Vec<SnpCpgData> = vec![
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1010,
//                1.0,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1000,
//                1010,
//                0.5,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr2".to_string(),
//                1000,
//                1010,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//            SnpCpgData::new(
//                "chr1".to_string(),
//                1100,
//                1110,
//                0.4,
//                PValMetadata::new(
//                    SnpType::SnpA,
//                    SnpType::SnpT,
//                    CpgType::Meth,
//                    CpgType::NotMeth,
//                    (1, 0, 0, 1)
//                )
//            ),
//        ];
//
//        let rank: Vec<usize> = vec![1, 2, 3, 0];
//
//        let check = resort(&test, &rank);
//
//        assert_eq!(check, compare);
//    }
//}
