#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
extern crate rust_htslib;

use lorikeet_genome::haplotype::event_map::EventMap;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::variant_context_utils::VariantContextUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rust_htslib::bam::record::CigarString;
use std::convert::TryFrom;

// #[cfg(not(test))]
// use log::{info, debug, warn};
//
// #[cfg(test)]
// use std::{println as info, println as debug, println as warn};

static NAME: &str = "foo";

fn test_mnps(
    ref_bases: &str,
    haplotype_bases: &str,
    cigar: &str,
    max_mnp_distance: Vec<usize>,
    expected_alleles: Vec<Vec<&str>>,
) {
    let mut hap = Haplotype::new(haplotype_bases.as_bytes(), false);
    hap.set_cigar(CigarString::try_from(cigar).unwrap().0);
    let mut loc = SimpleInterval::new(0, 1, ref_bases.len());
    hap.set_genome_location(loc.clone());
    for max_dist in max_mnp_distance {
        let events = EventMap::new(
            &hap,
            ref_bases.as_bytes(),
            loc.clone(),
            NAME.to_string(),
            max_dist,
        );
        assert_eq!(events.get_number_of_events(), expected_alleles.len());
        let found_alleles = events.get_variant_contexts();
        for i in 0..events.get_number_of_events() {
            let actual = found_alleles[i];
            assert_eq!(
                actual.get_reference().get_bases(),
                expected_alleles[i][0].as_bytes()
            );
            assert_eq!(
                actual.get_alternate_alleles()[0].get_bases(),
                expected_alleles[i][1].as_bytes()
            );
        }
    }
}

fn test_get_overlapping_events(
    haplotype_bases: &str,
    cigar: &str,
    query_loc: usize,
    expected_ref: Option<&ByteArrayAllele>,
    expected_alt: Option<&ByteArrayAllele>,
) {
    // Parameters that are shared across test cases
    let ref_bases = "AAAAAAAAAACGGTCA";
    let hap_start_wrt_ref = 7; // zero-based index into the ref array where the 0th base of the hap lines up
    let ref_loc = SimpleInterval::new(0, 1, ref_bases.len());

    let mut hap: Haplotype<SimpleInterval> = Haplotype::new(haplotype_bases.as_bytes(), false);
    hap.set_alignment_start_hap_wrt_ref(hap_start_wrt_ref);
    hap.set_cigar(CigarString::try_from(cigar).unwrap().0);
    let event_map = EventMap::new(&hap, ref_bases.as_bytes(), ref_loc, NAME.to_string(), 1);

    let overlapping_events = event_map.get_overlapping_events(query_loc);

    let events_expected = expected_alt.is_some() || expected_ref.is_some();
    assert_eq!(
        overlapping_events.len(),
        if events_expected { 1 } else { 0 }
    );

    if events_expected {
        assert_eq!(overlapping_events[0].get_reference(), expected_ref.unwrap());
        assert_eq!(
            &overlapping_events[0].get_alternate_alleles()[0],
            expected_alt.unwrap()
        );
    }
}

fn test_make_blocks(
    first_alleles: Vec<&str>,
    second_alleles: Vec<&str>,
    expected_alleles: Vec<&str>,
) {
    let vc1 = VariantContextUtils::make_from_alleles("x".to_string(), 20, 10, first_alleles);
    let vc2 = VariantContextUtils::make_from_alleles("x".to_string(), 20, 10, second_alleles);
    let expected =
        VariantContextUtils::make_from_alleles("x".to_string(), 20, 10, expected_alleles);

    // println!("vc1 {:?} vc2 {:?} expected {:?}", &vc1, &vc2, &expected);
    let mut event_map = EventMap::state_for_testing(Vec::new());
    let block = EventMap::make_block(vc1, vc2);

    assert_eq!(block.loc.get_start(), expected.loc.get_start());
    assert_eq!(block.alleles, expected.alleles);
}

#[test]
fn run_mnp_tests() {
    test_mnps(
        "TTTGGGAAA",
        "TTTCCCAAA",
        "3M3X3M",
        vec![1, 2, 3, 5, 10],
        vec![vec!["GGG", "CCC"]],
    );
    test_mnps(
        "TTTGGGAAA",
        "TTTCCCAAA",
        "3M3X3M",
        vec![0],
        vec![vec!["G", "C"], vec!["G", "C"], vec!["G", "C"]],
    );
    test_mnps(
        "TTTGGGAAA",
        "TTTCCCAAA",
        "9M",
        vec![1, 2, 3, 5, 10],
        vec![vec!["GGG", "CCC"]],
    );
    test_mnps(
        "TTTGGGAAA",
        "TTTCCCAAA",
        "9M",
        vec![0],
        vec![vec!["G", "C"], vec!["G", "C"], vec!["G", "C"]],
    );
    test_mnps(
        "TTTTTTTTT",
        "ATATATATA",
        "9M",
        vec![2],
        vec![vec!["TTTTTTTTT", "ATATATATA"]],
    );
    test_mnps(
        "ACGT",
        "CGTA",
        "4M",
        vec![1, 2, 3, 5, 10],
        vec![vec!["ACGT", "CGTA"]],
    );
    test_mnps(
        "ACGT",
        "CGTA",
        "4M",
        vec![0],
        vec![
            vec!["A", "C"],
            vec!["C", "G"],
            vec!["G", "T"],
            vec!["T", "A"],
        ],
    );
    test_mnps(
        "ACTTGC",
        "CATTCG",
        "6M",
        vec![1, 2],
        vec![vec!["AC", "CA"], vec!["GC", "CG"]],
    );
    test_mnps(
        "ACTTGC",
        "CATTCG",
        "6M",
        vec![3, 5, 10],
        vec![vec!["ACTTGC", "CATTCG"]],
    );
    test_mnps(
        "ACTTGC",
        "CATTCG",
        "6M",
        vec![0],
        vec![
            vec!["A", "C"],
            vec!["C", "A"],
            vec!["G", "C"],
            vec!["C", "G"],
        ],
    );
}

#[test]
fn run_overlapping_events_tests() {
    let deletion_ref_allele = ByteArrayAllele::new("ACGG".as_bytes(), true);
    let deletion_alt_allele = ByteArrayAllele::new("A".as_bytes(), false);
    let insertion_ref_allele = ByteArrayAllele::new("G".as_bytes(), true);
    let insertion_alt_allele = ByteArrayAllele::new("GTT".as_bytes(), false);
    let snp_ref_allele = ByteArrayAllele::new("G".as_bytes(), true);
    let snp_alt_allele = ByteArrayAllele::new("A".as_bytes(), false);

    // hap1
    test_get_overlapping_events(
        "AAATTTCA",
        "3M3D2I3M",
        10,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATTTCA",
        "3M3D2I3M",
        11,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATTTCA",
        "3M3D2I3M",
        12,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATTTCA",
        "3M3D2I3M",
        13,
        Some(&insertion_ref_allele),
        Some(&insertion_alt_allele),
    );

    // hap 2
    test_get_overlapping_events(
        "AAATCA",
        "3M3D3M",
        10,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATCA",
        "3M3D3M",
        11,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATCA",
        "3M3D3M",
        12,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );
    test_get_overlapping_events(
        "AAATCA",
        "3M3D3M",
        13,
        Some(&deletion_ref_allele),
        Some(&deletion_alt_allele),
    );

    // hap 3
    test_get_overlapping_events("AAACGATCA", "9M", 10, None, None);
    test_get_overlapping_events("AAACGATCA", "9M", 11, None, None);
    test_get_overlapping_events("AAACGATCA", "9M", 12, None, None);
    test_get_overlapping_events(
        "AAACGATCA",
        "9M",
        13,
        Some(&snp_ref_allele),
        Some(&snp_alt_allele),
    );

    // hap 4
    test_get_overlapping_events("AAACGGTTTCA", "6M2I3M", 10, None, None);
    test_get_overlapping_events("AAACGGTTTCA", "6M2I3M", 11, None, None);
    test_get_overlapping_events("AAACGGTTTCA", "6M2I3M", 12, None, None);
    test_get_overlapping_events(
        "AAACGGTTTCA",
        "6M2I3M",
        13,
        Some(&insertion_ref_allele),
        Some(&insertion_alt_allele),
    );
}

#[test]
fn run_test_blocks() {
    test_make_blocks(vec!["A", "G"], vec!["AGT", "A"], vec!["AGT", "G"]);
    test_make_blocks(vec!["A", "G"], vec!["A", "AGT"], vec!["A", "GGT"]);

    test_make_blocks(vec!["AC", "A"], vec!["A", "AGT"], vec!["AC", "AGT"]);
    test_make_blocks(vec!["ACGTA", "A"], vec!["A", "AG"], vec!["ACGTA", "AG"]);
    test_make_blocks(vec!["AC", "A"], vec!["A", "AGCGT"], vec!["AC", "AGCGT"]);
    test_make_blocks(vec!["A", "ACGTA"], vec!["AG", "A"], vec!["AG", "ACGTA"]);
    test_make_blocks(vec!["A", "AC"], vec!["AGCGT", "A"], vec!["AGCGT", "AC"]);
}
