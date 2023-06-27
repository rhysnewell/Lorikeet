#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::pair_hmm::pair_hmm_model::PairHMMModel;
use lorikeet_genome::utils::quality_utils::QualityUtils;
use ndarray::Array2;

lazy_static! {
    static ref INS_QUALS: Vec<u8> = vec![30, 45, 20, 10, 5, 60, 123];
    static ref DEL_QUALS: Vec<u8> = vec![30, 45, 20, 10, 5, 60, 123];
    static ref GAP_QUALS: Vec<u8> = vec![10, 20, 5];
    static ref READ_LENGTHS: Vec<usize> = vec![0, 1, 5, 20, 100, 250, 1000, 10000];
}

static TOLERANCE: f64 = 1e-9;

fn test_qual_to_probs(ins_qual: u8, del_qual: u8, gcp: u8, expected: Vec<f64>) {
    let mut actual = vec![0.0; PairHMMModel::TRANS_PROB_ARRAY_LENGTH];
    PairHMMModel::new().qual_to_trans_probs_with_vec(&mut actual, ins_qual, del_qual, gcp);
    assert_eq!(actual.len(), PairHMMModel::TRANS_PROB_ARRAY_LENGTH);
    assert_equals_double_array(&actual, &expected, TOLERANCE);
}

fn test_qual_to_probs_log10(ins_qual: u8, del_qual: u8, gcp: u8, expected: Vec<f64>) {
    let log_expected = expected
        .into_iter()
        .map(|v| v.log10())
        .collect::<Vec<f64>>();

    let mut actual = vec![0.0; PairHMMModel::TRANS_PROB_ARRAY_LENGTH];
    PairHMMModel::new().qual_to_trans_probs_log10_with_vec(&mut actual, ins_qual, del_qual, gcp);
    assert_eq!(actual.len(), PairHMMModel::TRANS_PROB_ARRAY_LENGTH);
    assert_equals_double_array(&actual, &log_expected, TOLERANCE);
}

fn test_quals_to_trans_probs(
    ins_quals: Vec<u8>,
    del_quals: Vec<u8>,
    gap_quals: Vec<u8>,
    expected: Array2<f64>,
) {
    let actual =
        PairHMMModel::new().qual_to_trans_probs_return_array(&ins_quals, &del_quals, &gap_quals);
    assert_eq!(actual.nrows(), expected.nrows());
    assert_eq!(actual.row(0).len(), expected.row(0).len());
    for i in 0..actual.nrows() {
        assert_equals_double_array(
            &actual.row(i).to_vec(),
            &expected.row(i).to_vec(),
            TOLERANCE,
        );
    }
}

fn assert_equals_double_array(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "Arrays of differing lengths actual {} -> expected {}",
        actual.len(),
        expected.len()
    );
    for (i, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            relative_eq!(a, e, epsilon = tolerance),
            "Actual {}, Expected {} at index {}",
            a,
            e,
            i
        );
    }
}

#[test]
fn make_test_qual_to_probs() {
    let qual_to_probs_data_provider = QualsToProbsDataProvider::new();
    for (ins_qual, del_qual, gcp, expected) in qual_to_probs_data_provider.into_iter() {
        test_qual_to_probs(ins_qual, del_qual, gcp, expected);
    }
}

#[test]
fn make_test_qual_to_probs_log10() {
    let qual_to_probs_data_provider = QualsToProbsDataProvider::new();
    for (ins_qual, del_qual, gcp, expected) in qual_to_probs_data_provider.into_iter() {
        test_qual_to_probs_log10(ins_qual, del_qual, gcp, expected);
    }
}

#[test]
fn make_test_quals_to_trans_probs() {
    let qual_to_trans_probs_data_provider = QualsToTransProbsDataProvider::new();
    for (ins_quals, del_quals, gap_quals, expected) in qual_to_trans_probs_data_provider.into_iter()
    {
        test_quals_to_trans_probs(ins_quals, del_quals, gap_quals, expected)
    }
}

// fn qual_to_probs_data_provider() {
//     let total_count = INS_QUALS.len() * DEL_QUALS.len() * GAP_QUALS.len();
//     for i in 0..total_count {
//         let gap = i % GAP_QUALS.len();
//         let indel_group = i % DEL_QUALS.len();
//         let del = indel_group % DEL_QUALS.len();
//         let ins = indel_group % DEL_QUALS.len();
//
//         let ins_qual = INS_QUALS[ins];
//         let del_qual = DEL_QUALS[del];
//         let gap_qual = GAP_QUALS[gap];
//
//         let trans = quals_to_probs(ins_qual, del_qual, gap_qual);
//     }
// }
//
fn quals_to_probs(ins_qual: u8, del_qual: u8, gap_qual: u8) -> Vec<f64> {
    let mut trans = vec![0.0; PairHMMModel::TRANS_PROB_ARRAY_LENGTH];
    let match_to_match = PairHMMModel::match_to_match_prob_static(ins_qual, del_qual);
    let match_to_insert = QualityUtils::qual_to_error_prob(ins_qual);
    let match_to_deletion = QualityUtils::qual_to_error_prob(del_qual);
    let indel_to_match = QualityUtils::qual_to_prob(gap_qual);
    let indel_to_indel = QualityUtils::qual_to_error_prob(gap_qual);

    trans[PairHMMModel::match_to_match] = match_to_match;
    trans[PairHMMModel::match_to_insertion] = match_to_insert;
    trans[PairHMMModel::match_to_deletion] = match_to_deletion;
    trans[PairHMMModel::indel_to_match] = indel_to_match;
    trans[PairHMMModel::deletion_to_deletion] = indel_to_indel;
    trans[PairHMMModel::insertion_to_insertion] = indel_to_indel;

    trans
}

struct QualsToTransProbsDataProvider {
    read_length_iterator: Box<dyn Iterator<Item = usize>>,
    quals_iterator: QualIterator,
}

impl QualsToTransProbsDataProvider {
    fn new() -> Self {
        Self {
            read_length_iterator: Box::new(READ_LENGTHS.clone().into_iter()),
            quals_iterator: QualIterator::new(),
        }
    }
}

impl Iterator for QualsToTransProbsDataProvider {
    type Item = (Vec<u8>, Vec<u8>, Vec<u8>, Array2<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        let read_length = self.read_length_iterator.next();
        match read_length {
            Some(read_length) => {
                let mut matrix =
                    Array2::zeros((read_length + 1, PairHMMModel::TRANS_PROB_ARRAY_LENGTH));
                let mut ins_quals = vec![0; read_length];
                let mut del_quals = vec![0; read_length];
                let mut gap_quals = vec![0; read_length];

                for i in 0..read_length {
                    let quals = match self.quals_iterator.next() {
                        Some(quals) => quals,
                        None => {
                            self.quals_iterator = QualIterator::new();
                            self.quals_iterator.next().unwrap()
                        }
                    };
                    let ins_qual = quals.0;
                    let del_qual = quals.1;
                    let gap_qual = quals.2;

                    let trans = quals_to_probs(ins_qual, del_qual, gap_qual);
                    let mut row = matrix.row_mut(i + 1);
                    row.iter_mut().zip(trans.iter()).for_each(|(m, t)| *m = *t);

                    ins_quals[i] = ins_qual;
                    del_quals[i] = del_qual;
                    gap_quals[i] = gap_qual;
                }

                Some((ins_quals, del_quals, gap_quals, matrix))
            }
            None => None,
        }
    }
}

struct QualsToProbsDataProvider {
    quals_iterator: QualIterator,
}

impl QualsToProbsDataProvider {
    fn new() -> Self {
        Self {
            quals_iterator: QualIterator::new(),
        }
    }
}

impl Iterator for QualsToProbsDataProvider {
    type Item = (u8, u8, u8, Vec<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        let quals = self.quals_iterator.next();
        match quals {
            None => None,
            Some(quals) => {
                let ins_qual = quals.0;
                let del_qual = quals.1;
                let gap_qual = quals.2;

                let trans = quals_to_probs(ins_qual, del_qual, gap_qual);

                Some((ins_qual, del_qual, gap_qual, trans))
            }
        }
    }
}

struct QualIterator {
    total_count: usize,
    i: usize,
}

impl QualIterator {
    fn new() -> Self {
        Self {
            total_count: INS_QUALS.len() * DEL_QUALS.len() * GAP_QUALS.len(),
            i: 0,
        }
    }
}

impl Iterator for QualIterator {
    type Item = (u8, u8, u8);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.total_count {
            let gap = self.i % GAP_QUALS.len();
            let indel_group = self.i % DEL_QUALS.len();
            let del = indel_group % DEL_QUALS.len();
            let ins = indel_group % DEL_QUALS.len();
            self.i += 1;
            Some((INS_QUALS[ins], DEL_QUALS[del], GAP_QUALS[gap]))
        } else {
            None
        }
    }
}
