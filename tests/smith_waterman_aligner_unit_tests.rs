#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
extern crate rayon;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;
extern crate itertools;
extern crate rand;
extern crate term;

use itertools::Itertools;
use lorikeet_genome::genotype::genotype_builder::Genotype;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::genotype::genotype_likelihoods::GenotypeLikelihoods;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::allele_frequency_calculator::AlleleFrequencyCalculator;
use lorikeet_genome::model::allele_likelihoods::AlleleLikelihoods;
use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::{allele_list::AlleleList, variants::SPAN_DEL_ALLELE};
use lorikeet_genome::pair_hmm::pair_hmm::PairHMM;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::PairHMMInputScoreImputator;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::{
    CigarUtils, ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS,
};
use lorikeet_genome::smith_waterman::bindings::{SWOverhangStrategy, SWParameters};
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    SmithWatermanAligner, SmithWatermanAlignmentResult, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::base_utils::BaseUtils;
use lorikeet_genome::utils::math_utils::{MathUtils, LOG10_ONE_HALF};
use lorikeet_genome::utils::quality_utils::QualityUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use rand::rngs::ThreadRng;
use rand::seq::index::sample;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView, Seq};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::ops::Deref;
use std::sync::Mutex;

fn get_aligner() -> SmithWatermanAligner {
    SmithWatermanAligner::new()
}

fn print_alignment(
    reference: &[u8],
    read: &[u8],
    alignment: SmithWatermanAlignmentResult,
    overhang_strategy: SWOverhangStrategy,
) {
    let mut b_read = String::new();
    let mut b_ref = String::new();
    let mut matches = String::new();

    let mut i = 0;
    let mut j = 0;

    let offset = alignment.get_alignment_offset();
    let mut cigar = alignment.get_cigar();

    if overhang_strategy != SWOverhangStrategy::SoftClip {
        // we need to go through all the hassle below only if we do not do soft-clipping;
        // otherwise offset is never negative
        if offset < 0 {
            while j < (-offset) as usize {
                b_read.push(read[j] as char);
                b_ref.push(b' ' as char);
                matches.push(b' ' as char);
                j += 1;
            }

            // at negative offsets, our cigar's first element carries overhanging bases
            // that we have just printed above. Tweak the first element to
            // exclude those bases. Here we create a new list of cigar elements, so the original
            // list/original cigar are unchanged (they are unmodifiable anyway!)
            let mut tweaked = Vec::new();
            tweaked.extend(cigar.0.clone());
            tweaked.insert(
                0,
                CigarUtils::cigar_from_element_and_length(
                    &cigar.0[0],
                    (cigar.0[0].len() as i32 + offset) as u32,
                ),
            );
            cigar = CigarString::from(tweaked);
        }
    }

    if offset > 0 {
        // note: the way this implementation works, cigar will ever start from S *only* if read starts before the ref, i.e. offset = 0
        while i < offset as usize {
            b_ref.push(reference[i] as char);
            b_read.push(b' ' as char);
            matches.push(b' ' as char);
            i += 1;
        }
    };

    for e in cigar.0 {
        match e {
            Cigar::Match(len) => {
                for _ in 0..len {
                    b_ref.push(if i < reference.len() {
                        reference[i] as char
                    } else {
                        b' ' as char
                    });
                    b_read.push(if j < read.len() {
                        read[j] as char
                    } else {
                        b' ' as char
                    });
                    matches.push(if i < reference.len() && j < read.len() {
                        if reference[i] == read[j] {
                            b'.' as char
                        } else {
                            b'*' as char
                        }
                    } else {
                        b' ' as char
                    })
                }
            }
            Cigar::Ins(len) => {
                for _ in 0..len {
                    b_ref.push(b'-' as char);
                    b_read.push(read[j] as char);
                    matches.push(b'I' as char);
                }
            }
            Cigar::SoftClip(len) => {
                for _ in 0..len {
                    b_ref.push(b' ' as char);
                    b_read.push(read[j] as char);
                    matches.push(b'S' as char);
                }
            }
            Cigar::Del(len) => {
                for _ in 0..len {
                    b_ref.push(reference[i] as char);
                    b_read.push(b'-' as char);
                    matches.push(b'D' as char);
                }
            }
            _ => {
                panic!("Unexpected Cigar Element");
            }
        }
    }

    while i < reference.len() {
        b_ref.push(reference[i] as char);
        i += 1;
    }
    while j < read.len() {
        b_read.push(read[j] as char);
        j += 1;
    }

    let mut pos = 0;
    let max_length = vec![matches.len(), b_read.len(), b_ref.len()]
        .into_iter()
        .max()
        .unwrap();
    while pos < max_length {
        print_cautiously(&matches, pos, 100);
        print_cautiously(&b_read, pos, 100);
        print_cautiously(&b_ref, pos, 100);
        println!("");
        pos += 100;
    }
}

fn print_cautiously(s: &String, start: usize, width: usize) {
    if start >= s.len() {
        println!("");
    } else {
        let end = min(start + width, s.len());
        println!("{}", &s[start..end]);
    };
}

fn assert_alignment_matches_expected(
    reference: &str,
    read: &str,
    expected_start: i32,
    expected_cigar: &str,
    weights: SWParameters,
    strategy: SWOverhangStrategy,
) {
    let alignment =
        SmithWatermanAligner::align(reference.as_bytes(), read.as_bytes(), &weights, strategy);
    print_alignment(
        reference.as_bytes(),
        read.as_bytes(),
        alignment.clone(),
        strategy,
    );
    assert_eq!(alignment.get_alignment_offset(), expected_start);
    assert_eq!(
        alignment.get_cigar().to_string(),
        expected_cigar.to_string()
    );
}

#[test]
fn make_test_read_alignment_to_ref_complex_alignment() {
    test_read_alignment_to_ref_complex_alignment("AAAGGACTGACTG", "ACTGACTGACTG", 1, "12M")
}

fn test_read_alignment_to_ref_complex_alignment(
    reference: &str,
    read: &str,
    expected_start: i32,
    expected_cigar: &str,
) {
    assert_alignment_matches_expected(
        reference,
        read,
        expected_start,
        expected_cigar,
        *ORIGINAL_DEFAULT,
        SWOverhangStrategy::SoftClip,
    )
}

#[test]
fn make_test_odd_no_alignment() {
    let ref_1 = "AAAGACTACTG";
    let read_1 = "AACGGACACTG";

    test_odd_no_alignment(
        ref_1,
        read_1,
        SWParameters::new(50, -100, -220, -12),
        1,
        "2M2I3M1D4M",
    );
    test_odd_no_alignment(
        ref_1,
        read_1,
        SWParameters::new(200, -50, -300, -22),
        0,
        "11M",
    );
}

fn test_odd_no_alignment(
    reference: &str,
    read: &str,
    weights: SWParameters,
    expected_start: i32,
    expected_cigar: &str,
) {
    assert_alignment_matches_expected(
        reference,
        read,
        expected_start,
        expected_cigar,
        weights,
        SWOverhangStrategy::SoftClip,
    )
}

#[test]
fn test_indels_at_start_and_end() {
    let matc = "CCCCC";
    let reference = "AAA".to_string() + matc;
    let read = matc.to_string() + "GGG";
    let expected_start = 3;
    let expected_cigar = "5M3S";
    assert_alignment_matches_expected(
        reference.as_str(),
        read.as_str(),
        expected_start,
        expected_cigar,
        *ORIGINAL_DEFAULT,
        SWOverhangStrategy::SoftClip,
    );
}

#[test]
fn test_degenerate_alignment_with_indels_at_both_ends() {
    let reference = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
    let alt = "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
    let expected_start = 14;
    let expected_cigar = "31M20S";
    assert_alignment_matches_expected(
        reference,
        alt,
        expected_start,
        expected_cigar,
        *STANDARD_NGS,
        SWOverhangStrategy::SoftClip,
    );
}

#[test]
fn test_for_identical_alignments_with_differing_flank_lengths() {
    //This test is designed to ensure that the indels are correctly placed
    //if the region flanking these indels is extended by a varying amount.
    //It checks for problems caused by floating point rounding leading to different
    //paths being selected.

    //Create two versions of the same sequence with different flanking regions.
    let padded_ref = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT";
    let padded_hap = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA--GGGCC---------------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT".replace("-", "");
    let not_padded_ref = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA";
    let not_padded_hap = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA---------GGGCC--------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA".replace("-", "");

    //a simplified version of the getCigar routine in the haplotype caller to align these
    let SW_PAD = "NNNNNNNNNN";
    let paddeds_ref = format!("{}{}{}", SW_PAD, padded_ref, SW_PAD);
    let paddeds_hap = format!("{}{}{}", SW_PAD, padded_hap, SW_PAD);
    let not_paddeds_ref = format!("{}{}{}", SW_PAD, not_padded_ref, SW_PAD);
    let not_paddeds_hap = format!("{}{}{}", SW_PAD, not_padded_hap, SW_PAD);

    let mut padded_alignment = SmithWatermanAligner::align(
        paddeds_ref.as_bytes(),
        paddeds_hap.as_bytes(),
        &*NEW_SW_PARAMETERS,
        SWOverhangStrategy::SoftClip,
    );
    let mut not_padded_alignment = SmithWatermanAligner::align(
        not_paddeds_ref.as_bytes(),
        not_paddeds_hap.as_bytes(),
        &*NEW_SW_PARAMETERS,
        SWOverhangStrategy::SoftClip,
    );

    //Now verify that the two sequences have the same alignment and not match positions.
    let raw_padded = padded_alignment.get_cigar();
    let not_padded = padded_alignment.get_cigar();
    let padded_c = raw_padded.0.clone();
    let not_padded_c = not_padded.0.clone();

    assert_eq!(padded_c.len(), not_padded_c.len());
    for i in 0..not_padded_c.len() {
        let pc = &padded_c[i];
        let npc = &not_padded_c[i];
        if CigarUtils::cigar_elements_are_same_type(pc, &Some(Cigar::Match(0)))
            && CigarUtils::cigar_elements_are_same_type(npc, &Some(Cigar::Match(0)))
        {
            continue;
        };
        let l1 = pc.len();
        let l2 = npc.len();
        assert_eq!(l1, l2);
        assert!(CigarUtils::cigar_elements_are_same_type(
            pc,
            &Some(npc.clone())
        ));
    }
}

#[test]
fn get_substrings_match_tests() {
    test_substring_match(3, "5M", SWOverhangStrategy::SoftClip);
    test_substring_match(0, "3D5M", SWOverhangStrategy::Indel);
    test_substring_match(0, "3D5M", SWOverhangStrategy::LeadingIndel);
    test_substring_match(3, "5M", SWOverhangStrategy::Ignore);
}

fn test_substring_match(expected_start: i32, expected_cigar: &str, strategy: SWOverhangStrategy) {
    let matching_section = "CCCCC";
    let reference = "AAA".to_string() + matching_section;
    let read = matching_section;
    assert_alignment_matches_expected(
        reference.as_str(),
        read,
        expected_start,
        expected_cigar,
        *ORIGINAL_DEFAULT,
        strategy,
    )
}
