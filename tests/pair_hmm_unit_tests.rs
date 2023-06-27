#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::allele_likelihoods::AlleleLikelihoods;


use lorikeet_genome::pair_hmm::pair_hmm::PairHMM;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::{
    PairHMMInputScoreImputator,
};



use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::base_utils::{BASES};
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::quality_utils::QualityUtils;
use lorikeet_genome::utils::simple_interval::{SimpleInterval};


use rayon::prelude::*;

use rust_htslib::bam::record::{Cigar, CigarString};
use std::cmp::{max, min};
use std::collections::{HashMap};




static ALLOW_READS_LONGER_THAN_HAPLOTYPE: bool = true;
static DEBUG: bool = false;
static EXTENSIVE_TESTING: bool = true;
static CONTEXT: &'static str = "ACGTAATGACGATTGCA";
static LEFT_FLANK: &'static str = "GATTTATCATCGAGTCTGC";
static RIGHT_FLANK: &'static str = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTTA";
static MASSIVE_QUAL: u8 = 100;

#[derive(Debug)]
struct BasicLikelihoodTestProvider {
    reference: String,
    next_reference: Option<String>,
    read: String,
    ref_bases_with_context: String,
    next_ref_bases_with_context: Option<String>,
    read_bases_with_context: String,
    base_qual: u8,
    ins_qual: u8,
    del_qual: u8,
    gcp: u8,
    expected_qual: usize,
    left: bool,
    right: bool,
}

impl BasicLikelihoodTestProvider {
    pub fn new(
        reference: String,
        next_reference: Option<String>,
        read: String,
        base_qual: u8,
        ins_qual: u8,
        del_qual: u8,
        expected_qual: usize,
        gcp: u8,
        left: bool,
        right: bool,
    ) -> Self {
        let ref_bases_with_context = as_bytes(&reference, left, right);
        let next_ref_bases_with_context = if next_reference.is_some() {
            Some(as_bytes(next_reference.as_ref().unwrap(), left, right))
        } else {
            None
        };
        let read_bases_with_context = as_bytes(&read, false, false);

        Self {
            reference,
            next_reference,
            read,
            ref_bases_with_context,
            next_ref_bases_with_context,
            read_bases_with_context,
            base_qual,
            ins_qual,
            del_qual,
            gcp,
            expected_qual,
            left,
            right,
        }
    }

    fn expected_log_likelihood(&self) -> f64 {
        println!(
            "expcted qual {} ref length {}",
            self.expected_qual,
            self.read_bases_with_context.len()
        );
        return (self.expected_qual as f64 / -10.0)
            + 0.03
            + (1.0 / self.ref_bases_with_context.len() as f64).log10();
    }

    fn get_tolerance(&self) -> f64 {
        0.2
    }

    fn tolerance_from_reference(&self) -> f64 {
        1e-3 // has to be very tolerant -- this approximation is quite approximate
    }

    fn tolerance_from_exact(&self) -> f64 {
        1e-9
    }

    fn calc_log10_likelihood(&self, anchor_indel: bool) -> f64 {
        let mut pair_hmm = PairHMM::quick_initialize(
            self.read_bases_with_context.len(),
            self.ref_bases_with_context.len(),
        );
        pair_hmm.do_not_use_tristate_correction();
        return pair_hmm.compute_read_likelihood_given_haplotype_log10(
            self.ref_bases_with_context.as_bytes(),
            self.read_bases_with_context.as_bytes(),
            &self.qual_as_bytes(self.base_qual, false, anchor_indel),
            &self.qual_as_bytes(self.ins_qual, true, anchor_indel),
            &self.qual_as_bytes(self.del_qual, true, anchor_indel),
            &self.qual_as_bytes(self.gcp, false, anchor_indel),
            true,
            match &self.next_ref_bases_with_context {
                Some(bases) => Some(bases.as_bytes()),
                None => None,
            },
            // Some(self.next_ref_bases_with_context.as_bytes()),
        );
    }

    fn qual_as_bytes(&self, phred_qual: u8, do_gop: bool, anchor_indel: bool) -> Vec<u8> {
        let mut phred_quals = vec![0; self.read_bases_with_context.len()];

        if anchor_indel {
            // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
            phred_quals = vec![MASSIVE_QUAL; self.read_bases_with_context.len()];

            // update just the bases corresponding to the provided micro read with the quality scores
            if do_gop {
                phred_quals[CONTEXT.len()] = phred_qual;
            } else {
                for i in 0..self.read.len() {
                    phred_quals[i + CONTEXT.len()] = phred_qual;
                }
            };
        } else {
            phred_quals = vec![phred_qual; self.read_bases_with_context.len()];
        }

        return phred_quals;
    }
}

fn as_bytes<'a>(bases: &'a str, left: bool, right: bool) -> String {
    format!(
        "{}{}{}{}{}",
        if left { LEFT_FLANK } else { "" },
        CONTEXT,
        bases,
        CONTEXT,
        if right { RIGHT_FLANK } else { "" }
    )
}

// Requires exact log pairhmm to be implemented. Will skip for now
#[test]
fn make_basic_likelihood_tests() {
    // context on either side is ACGTTGCA REF ACGTTGCA
    // test all combinations
    let base_quals = if EXTENSIVE_TESTING {
        vec![10, 20, 30, 40, 50]
    } else {
        vec![30]
    };
    let indel_quals = if EXTENSIVE_TESTING {
        vec![20, 30, 40, 50]
    } else {
        vec![40]
    };
    let gcps = if EXTENSIVE_TESTING {
        vec![8, 10, 20]
    } else {
        vec![10]
    };
    let sizes = if EXTENSIVE_TESTING {
        vec![2, 3, 4, 5, 7, 8, 9, 10, 20, 30, 35]
    } else {
        vec![2]
    };

    for base_qual in base_quals.iter() {
        for indel_qual in indel_quals.iter() {
            for gcp in gcps.iter() {
                for ref_base in BASES.iter() {
                    for read_base in BASES.iter() {
                        let reference = String::from_utf8(vec![*ref_base]).unwrap();
                        let read = String::from_utf8(vec![*read_base]).unwrap();
                        let expected = if ref_base == read_base {
                            0
                        } else {
                            *base_qual as usize
                        };
                        // runBasicLikelihoodTests uses calcLogLikelihood(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                        test_basic_likelihoods(BasicLikelihoodTestProvider::new(
                            reference,
                            None,
                            read,
                            *base_qual,
                            *indel_qual,
                            *indel_qual,
                            expected,
                            *gcp,
                            false,
                            false,
                        ));
                    }
                }

                // test insertions and deletions
                for size in sizes.iter() {
                    for base in BASES.iter() {
                        let expected =
                            (*indel_qual as i32 + (*size as i32 - 2) * *gcp as i32) as usize;

                        println!("expected {}, indel qual {} size {} gcp {} size - 2 {} expected recalc {}", expected, indel_qual, size, gcp, size-2, ((*indel_qual as i32 + (*size as i32 - 2) * *gcp as i32)));
                        for insertion_p in vec![true, false] {
                            let small = String::from_utf8(vec![*base]).unwrap();
                            let big = String::from_utf8(vec![*base; *size as usize]).unwrap();

                            let reference = if insertion_p {
                                small.clone()
                            } else {
                                big.clone()
                            };

                            let read = if insertion_p {
                                big.clone()
                            } else {
                                small.clone()
                            };

                            test_basic_likelihoods(BasicLikelihoodTestProvider::new(
                                reference.clone(),
                                None,
                                read.clone(),
                                *base_qual,
                                *indel_qual,
                                *indel_qual,
                                expected,
                                *gcp,
                                false,
                                false,
                            ));
                            test_basic_likelihoods(BasicLikelihoodTestProvider::new(
                                reference.clone(),
                                None,
                                read.clone(),
                                *base_qual,
                                *indel_qual,
                                *indel_qual,
                                expected,
                                *gcp,
                                true,
                                false,
                            ));
                            test_basic_likelihoods(BasicLikelihoodTestProvider::new(
                                reference.clone(),
                                None,
                                read.clone(),
                                *base_qual,
                                *indel_qual,
                                *indel_qual,
                                expected,
                                *gcp,
                                false,
                                true,
                            ));
                            test_basic_likelihoods(BasicLikelihoodTestProvider::new(
                                reference.clone(),
                                None,
                                read.clone(),
                                *base_qual,
                                *indel_qual,
                                *indel_qual,
                                expected,
                                *gcp,
                                true,
                                true,
                            ));
                        }
                    }
                }
            }
        }
    }
}

fn test_basic_likelihoods(cfg: BasicLikelihoodTestProvider) {
    println!(
        "CFG {:?}, read len {} ref len {}",
        &cfg,
        cfg.read.len(),
        cfg.reference.len()
    );
    if ALLOW_READS_LONGER_THAN_HAPLOTYPE || cfg.read.len() <= cfg.reference.len() {
        // let exact_log_l = cfg.calc_log10_likelihood(true);
        let actual_log_l = cfg.calc_log10_likelihood(true);
        let expected_log_l = cfg.expected_log_likelihood();

        // compare to our theoretical expectation
        println!("Actual {} Expected {}", actual_log_l, expected_log_l);
        assert!(
            relative_eq!(actual_log_l, expected_log_l, epsilon = cfg.get_tolerance()),
            "Failed with hmm calc {} -> {}",
            actual_log_l,
            expected_log_l
        );
        assert!(
            MathUtils::is_valid_log10_probability(actual_log_l),
            "Bad log likelihood {}",
            actual_log_l
        );
    }
}

#[test]
fn test_mismatch_in_every_position_in_the_read_with_centred_haplotype() {
    let haplotpe_1 =
        "TTCTCTTCTGTTGTGGCTGGTTTTCTCTTCTGTTGTGGCTGGTTTTCTCTTCTGTTGTGGCTGGTT".as_bytes();
    let match_qual = 90;
    let mismatch_qual = 20;
    let indel_qual = 80;
    let offset = 2;
    let gop = vec![indel_qual; haplotpe_1.len() - 2 * offset];
    let gcp = vec![indel_qual; haplotpe_1.len() - 2 * offset];

    let mut logless_hmm = PairHMM::quick_initialize(gop.len(), haplotpe_1.len());
    logless_hmm.do_not_use_tristate_correction();

    for k in 0..haplotpe_1.len() - 2 * offset {
        let mut quals = vec![match_qual; haplotpe_1.len() - 2 * offset];

        // one base mismatches the haplotype
        quals[k] = mismatch_qual;
        let mut m_read = haplotpe_1[offset..haplotpe_1.len() - offset].to_vec();
        // change single base at position k to C. If it's a C, change to T
        m_read[k] = if m_read[k] == b'C' { b'T' } else { b'C' };

        let res_1 = logless_hmm.compute_read_likelihood_given_haplotype_log10(
            haplotpe_1, &m_read, &quals, &gop, &gop, &gcp, true, None,
        );
        let expected = ((1.0 / haplotpe_1.len() as f64)
            * QualityUtils::qual_to_prob(match_qual).powf((m_read.len() - 1) as f64)
            * QualityUtils::qual_to_error_prob(mismatch_qual))
        .log10();

        assert!(
            relative_eq!(res_1, expected, epsilon = 1e-2),
            "Result {}, Expected {}",
            res_1,
            expected
        );
    }
}

#[test]
fn test_mismatch_in_every_position_in_the_read() {
    let haplotpe_1 = "TTCTCTTCTGTTGTGGCTGGTT".as_bytes();
    let match_qual = 90;
    let mismatch_qual = 20;
    let indel_qual = 80;
    let offset = 2;
    let gop = vec![indel_qual; haplotpe_1.len() - offset];
    let gcp = vec![indel_qual; haplotpe_1.len() - offset];

    let mut logless_hmm = PairHMM::quick_initialize(gop.len(), haplotpe_1.len());
    logless_hmm.do_not_use_tristate_correction();

    for k in 0..haplotpe_1.len() - offset {
        let mut quals = vec![match_qual; haplotpe_1.len() - offset];

        // one base mismatches the haplotype
        quals[k] = mismatch_qual;
        let mut m_read = haplotpe_1[offset..haplotpe_1.len()].to_vec();
        // change single base at position k to C. If it's a C, change to T
        m_read[k] = if m_read[k] == b'C' { b'T' } else { b'C' };

        let res_1 = logless_hmm.compute_read_likelihood_given_haplotype_log10(
            haplotpe_1, &m_read, &quals, &gop, &gop, &gcp, true, None,
        );
        let expected = ((1.0 / haplotpe_1.len() as f64)
            * QualityUtils::qual_to_prob(match_qual).powf((m_read.len() - 1) as f64)
            * QualityUtils::qual_to_error_prob(mismatch_qual))
        .log10();

        assert!(
            relative_eq!(res_1, expected, epsilon = 1e-2),
            "Result {}, Expected {}",
            res_1,
            expected
        );
    }
}

fn test_multiple_read_matches_in_haplotype(read_size: usize, ref_size: usize) {
    let read_bases = vec![b'A'; read_size];
    let ref_bases = format!("CC{}GGA", String::from_utf8(vec![b'A'; ref_size]).unwrap());
    let base_qual = 20;
    let base_quals = vec![base_qual; read_bases.len()];
    let ins_qual = 37;
    let ins_quals = vec![ins_qual; read_bases.len()];
    let del_qual = 37;
    let del_quals = vec![del_qual; read_bases.len()];
    let gcp = 10;
    let gcps = vec![gcp; read_bases.len()];

    let mut hmm = PairHMM::quick_initialize(read_bases.len(), ref_bases.len());
    hmm.do_not_use_tristate_correction();
    // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases

    let d = hmm.compute_read_likelihood_given_haplotype_log10(
        ref_bases.as_bytes(),
        &read_bases,
        &base_quals,
        &ins_quals,
        &del_quals,
        &gcps,
        true,
        None,
    );
    assert!(d <= 0.0);
}

fn test_all_matching_read(read_size: usize, ref_size: usize) {
    let read_bases = vec![b'A'; read_size];
    let ref_bases = vec![b'A'; ref_size];
    let base_qual = 20;
    let base_quals = vec![base_qual; read_bases.len()];
    let ins_qual = 100;
    let ins_quals = vec![ins_qual; read_bases.len()];
    let del_qual = 100;
    let del_quals = vec![del_qual; read_bases.len()];
    let gcp = 100;
    let gcps = vec![gcp; read_bases.len()];

    let mut hmm = PairHMM::quick_initialize(read_bases.len(), ref_bases.len());
    hmm.do_not_use_tristate_correction();
    // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
    let d = hmm.compute_read_likelihood_given_haplotype_log10(
        &ref_bases,
        &read_bases,
        &base_quals,
        &ins_quals,
        &del_quals,
        &gcps,
        true,
        None,
    );
    let expected =
        get_expected_matching_log_likelihood(&read_bases, &ref_bases, base_qual, ins_qual);
    assert!(
        relative_eq!(d, expected, epsilon = 1e-3),
        "Likelihoods should sum to just the error prob of the read but got Result {}, Expected {}",
        d,
        expected
    );
}

fn get_expected_matching_log_likelihood(
    read_bases: &[u8],
    ref_bases: &[u8],
    base_qual: u8,
    ins_qual: u8,
) -> f64 {
    let mut expected = 0.0;
    let initial_condition =
        (ref_bases.len() as f64 - read_bases.len() as f64 + 1.0).abs() / ref_bases.len() as f64;
    if read_bases.len() < ref_bases.len() {
        expected = (initial_condition
            * (QualityUtils::qual_to_prob(base_qual).powf(read_bases.len() as f64)))
        .log10();
    } else if read_bases.len() > ref_bases.len() {
        expected = (initial_condition
            * QualityUtils::qual_to_prob(base_qual).powf(ref_bases.len() as f64)
            * QualityUtils::qual_to_error_prob(ins_qual)
                .powf(read_bases.len() as f64 - ref_bases.len() as f64))
        .log10();
    }
    return expected;
}

fn test_read_same_as_haplotype(read_size: usize) {
    let read_bases = vec![b'A'; read_size];
    let base_qual = 20;
    let base_quals = vec![base_qual; read_bases.len()];
    let ins_qual = 37;
    let ins_quals = vec![ins_qual; read_bases.len()];
    let del_qual = 37;
    let del_quals = vec![del_qual; read_bases.len()];
    let gcp = 10;
    let gcps = vec![gcp; read_bases.len()];

    let mut hmm = PairHMM::quick_initialize(read_bases.len(), read_bases.len());
    hmm.do_not_use_tristate_correction();
    // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases

    let d = hmm.compute_read_likelihood_given_haplotype_log10(
        &read_bases,
        &read_bases,
        &base_quals,
        &ins_quals,
        &del_quals,
        &gcps,
        true,
        None,
    );
    assert!(d <= 0.0);
}

#[test]
fn hmm_provider_simple() {
    for read_size in vec![1, 2, 5, 10] {
        test_read_same_as_haplotype(read_size)
    }
}

#[test]
fn make_hmm_provider() {
    for read_size in vec![1, 2, 5, 10] {
        for ref_size in vec![1, 2, 5, 10] {
            if ref_size > read_size {
                test_multiple_read_matches_in_haplotype(read_size, ref_size);
                test_all_matching_read(read_size, ref_size);
            }
        }
    }
}

#[test]
fn make_big_read_hmm_provider() {
    let read_1 = "ACCAAGTAGTCACCGT";
    let ref_1 = "ACCAAGTAGTCACCGTAACG";

    for n_read_copies in vec![1, 2, 10, 20, 50] {
        for n_ref_copies in vec![1, 2, 10, 20, 100] {
            if n_ref_copies > n_read_copies {
                let read = std::iter::repeat(read_1)
                    .take(n_read_copies)
                    .collect::<String>();
                let reference = std::iter::repeat(ref_1)
                    .take(n_ref_copies)
                    .collect::<String>();
                test_really_big_reads(read, reference);
            }
        }
    }

    // vec![1, 2, 10, 20, 50, 100].into_par_iter().for_each(|n_read_copies| {
    //     vec![1, 2, 10, 20, 100, 200].into_par_iter().for_each(|n_ref_copies|{
    //         if n_ref_copies > n_read_copies {
    //             let read = std::iter::repeat(read_1).take(n_read_copies).collect::<String>();
    //             let reference = std::iter::repeat(ref_1).take(n_ref_copies).collect::<String>();
    //             test_really_big_reads(read, reference);
    //         }
    //     })
    // })
}

fn test_really_big_reads(read: String, reference: String) {
    let read_bases = read.as_bytes();
    let ref_bases = reference.as_bytes();

    let base_qual = 30;
    let ins_qual = 40;
    let del_qual = 40;
    let gcp = 10;
    let base_quals = vec![base_qual; read_bases.len()];
    let ins_quals = vec![ins_qual; read_bases.len()];
    let del_quals = vec![del_qual; read_bases.len()];
    let gcps = vec![gcp; read_bases.len()];

    let mut hmm = PairHMM::quick_initialize(read_bases.len(), ref_bases.len());
    hmm.do_not_use_tristate_correction();
    // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
    let _d = hmm.compute_read_likelihood_given_haplotype_log10(
        &ref_bases,
        &read_bases,
        &base_quals,
        &ins_quals,
        &del_quals,
        &gcps,
        true,
        None,
    );
}

#[test]
fn test_max_lengths_bigger_than_provided_read() {
    let read_bases = "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAACA";
    let ref_bases = "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAAAACTTCTGAGAAAAAAAAAAAAAATTAAATCAAACCCTGATTCCTTAAAGGTAGTAAAAAAACATCATTCTTTCTTAGTGGAATAGAAACTAGGTCAAAAGAACAGTGATTC";

    let quals = vec![
        35, 34, 31, 32, 35, 34, 32, 31, 36, 30, 31, 32, 36, 34, 33, 32, 32, 32, 33, 32, 30, 35, 33,
        35, 36, 36, 33, 33, 33, 32, 32, 32, 37, 33, 36, 35, 33, 32, 34, 31, 36, 35, 35, 35, 35, 33,
        34, 31, 31, 30, 28, 27, 26, 29, 26, 25, 29, 29,
    ];
    let ins_qual = vec![
        46, 46, 46, 46, 46, 47, 45, 46, 45, 48, 47, 44, 45, 48, 46, 43, 43, 42, 48, 48, 45, 47, 47,
        48, 48, 47, 48, 45, 38, 47, 45, 39, 47, 48, 47, 47, 48, 46, 49, 48, 49, 48, 46, 47, 48, 44,
        44, 43, 39, 32, 34, 36, 46, 48, 46, 44, 45, 45,
    ];
    let del_qual = vec![
        44, 44, 44, 43, 45, 44, 43, 42, 45, 46, 45, 43, 44, 47, 45, 40, 40, 40, 45, 46, 43, 45, 45,
        44, 46, 46, 46, 43, 35, 44, 43, 36, 44, 45, 46, 46, 44, 44, 47, 43, 47, 45, 45, 45, 46, 45,
        45, 46, 44, 35, 35, 35, 45, 47, 45, 44, 44, 43,
    ];
    let gcps = vec![10; del_qual.len()];
    let mut hmm = PairHMM::quick_initialize(read_bases.len(), ref_bases.len());
    hmm.do_not_use_tristate_correction();

    for _n_extra_max_size in 0..100 {
        hmm.compute_read_likelihood_given_haplotype_log10(
            ref_bases.as_bytes(),
            read_bases.as_bytes(),
            &quals,
            &ins_qual,
            &del_qual,
            &gcps,
            true,
            None,
        );
    }
}

#[test]
fn test_likelihoods_from_haplotypes() {
    let read_size = 10;
    let ref_size = 20;
    let read_bases = vec![b'A'; read_size];
    let ref_bases = vec![b'A'; ref_size];
    let base_qual = 20;
    let ins_qual = MASSIVE_QUAL;

    let ref_h: Haplotype<SimpleInterval> = Haplotype::new(ref_bases.as_slice(), true);
    let read_quals = vec![base_qual; read_bases.len()];
    let reads = vec![ArtificialReadUtils::create_artificial_read(
        read_bases.as_slice(),
        read_quals.as_slice(),
        CigarString::from(vec![Cigar::Match(read_bases.len() as u32)]),
    )];
    let mut sample_evidence_map = HashMap::new();
    sample_evidence_map.insert(0, reads.clone());

    let input_score_imputator = PairHMMInputScoreImputator::new(MASSIVE_QUAL);
    let mut hmm = PairHMM::quick_initialize(read_bases.len(), ref_bases.len());
    hmm.do_not_use_tristate_correction();
    let mut haplotype_mat = AlleleLikelihoods::new(
        vec![ref_h],
        vec!["sample_1".to_string()],
        sample_evidence_map,
    );
    hmm.compute_log10_likelihoods(0, &mut haplotype_mat, Vec::new(), &input_score_imputator);
    assert!(hmm.get_log_likelihood_array().len() == 0);

    hmm.compute_log10_likelihoods(0, &mut haplotype_mat, reads, &input_score_imputator);
    let expected = get_expected_matching_log_likelihood(
        read_bases.as_slice(),
        ref_bases.as_slice(),
        base_qual,
        ins_qual,
    );
    let la = hmm.get_log_likelihood_array();

    assert_eq!(la.len(), 1);
    assert!(
        relative_eq!(la[0], expected, epsilon = 1e-3),
        "Likelihoods should sum to just the error prob of the read but got Result {}, Expected {}",
        la[0],
        expected
    );
}

#[test]
fn test_find_first_position_where_haplotypes_differ() {
    for haplotype_size_1 in 10..30 {
        for haplotype_size_2 in 10..50 {
            let max_length = max(haplotype_size_1, haplotype_size_2);
            let min_length = min(haplotype_size_1, haplotype_size_2);
            for differing_site in 0..max_length + 1 {
                for one_is_diff in vec![true, false] {
                    let mut hap_1 = vec![b'A'; haplotype_size_1];
                    let mut hap_2 = vec![b'A'; haplotype_size_2];
                    let expected = if one_is_diff {
                        make_diff(&mut hap_1, differing_site, min_length)
                    } else {
                        make_diff(&mut hap_2, differing_site, min_length)
                    };

                    let actual =
                        PairHMM::find_first_position_where_haplotypes_differ(&hap_1, &hap_2);
                    assert_eq!(
                        actual,
                        expected,
                        "Bad differing site for {} vs. {}",
                        String::from_utf8(hap_1).unwrap(),
                        String::from_utf8(hap_2).unwrap()
                    );
                }
            }
        }
    }
}

fn make_diff(bytes: &mut [u8], site: usize, min_size: usize) -> usize {
    if site < bytes.len() {
        bytes[site] = b'C';
        return min(site, min_size);
    } else {
        return min_size;
    }
}

fn test_haplotype_indexing(root_1: &str, root_2: &str, root_3: &str, read: &str) {
    let TOLERANCE = 1e-9;

    let prefix = "AACCGGTTTTTGGGCCCAAACGTACGTACAGTTGGTCAACATCGATCAGGTTCCGGAGTAC";
    let _max_read_length = read.len();
    let _max_haplotype_length = prefix.len() + root_1.len();

    let max_read_length = read.len();
    let max_haplotype_length = prefix.len() + root_1.len();

    // the initialization occurs once, at the start of the evalution of reads
    let mut hmm = PairHMM::quick_initialize(max_read_length, max_haplotype_length);
    hmm.do_not_use_tristate_correction();

    for prefix_start in (0..=prefix.len()).into_iter().rev() {
        let my_prefix = &prefix[prefix_start..prefix.len()];
        let hap_1 = my_prefix.to_string() + root_1;
        let hap_2 = my_prefix.to_string() + root_2;
        let hap_3 = my_prefix.to_string() + root_3;

        let hap_start = PairHMM::find_first_position_where_haplotypes_differ(
            hap_1.as_bytes(),
            hap_2.as_bytes(),
        );

        // Run the HMM on the first haplotype, peaking ahead the second, to set up caching
        // Then run on the second haplotype in both cached and uncached mode, and verify that results are the same
        // When evaluating actual2, it is important that we both apply old caching from hap1 and
        // set up new caching for hap3, to ensure read/write operations do not cause conflicts
        let _actual_1 = test_haplotype_indexing_calc(&mut hmm, &hap_1, Some(&hap_2), read, 0, true);
        let actual_2 =
            test_haplotype_indexing_calc(&mut hmm, &hap_2, Some(&hap_3), read, hap_start, false);
        let expected_2 = test_haplotype_indexing_calc(&mut hmm, &hap_2, None, read, 0, true);
        assert!(relative_eq!(actual_2, expected_2, epsilon=TOLERANCE), "HMM caching calculation failed for read {} against hap with prefix {} expected {} got {} with hap_start {}", &read, &prefix, expected_2, actual_2, hap_start);
    }
}

fn test_haplotype_indexing_calc(
    hmm: &mut PairHMM,
    hap: &str,
    next_hap: Option<&str>,
    read: &str,
    _hap_start: usize,
    recache: bool,
) -> f64 {
    let read_bases = read.as_bytes();
    // if not peaking ahead to capture info for a future cache run, the next haplotype will be None, and this should be passed to HMM
    let next_hap_bases = match next_hap {
        None => None,
        Some(next_hap) => Some(next_hap.as_bytes()),
    };

    let base_quals = vec![30; read_bases.len()];
    let ins_quals = vec![45; read_bases.len()];
    let del_quals = vec![40; read_bases.len()];
    let gcp = vec![10; read_bases.len()];

    let d = hmm.compute_read_likelihood_given_haplotype_log10(
        hap.as_bytes(),
        read_bases,
        &base_quals,
        &ins_quals,
        &del_quals,
        &gcp,
        recache,
        next_hap_bases,
    );

    return d;
}

#[test]
fn make_haplotype_indexing_provider() {
    // First difference (root2, root3) is the base position immediately following first difference (root1, root2)
    let root_1 = "ACGTGTCAAACCGGGTT";
    let root_2 = "ACGTGTCACACTGGGTT"; // differs in two locations from root1
    let root_3 = "ACGTGTCACTCCGCGTT"; // differs in two locations from root2

    let read_1 = "ACGTGTCACACTGGATT"; // 1 diff from 2, 2 diff from root1, 2 diff from root3
    let read_2 = root_1;
    let read_3 = root_2;
    let read_4 = "ACGTGTCACACTGGATTCGAT";
    let read_5 = "CCAGTAACGTGTCACACTGGATTCGAT";
    for read in vec![read_1, read_2, read_3, read_4, read_5] {
        for read_length in 10..read.len() {
            let my_read = &read[0..read_length];
            test_haplotype_indexing(root_1, root_2, root_3, my_read);
        }
    }
}
