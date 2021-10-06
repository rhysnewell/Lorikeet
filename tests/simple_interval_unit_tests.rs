#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
extern crate rust_htslib;

use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reads::read_clipper::ReadClipper;
use lorikeet_genome::test_utils::read_clipper_test_utils::ReadClipperTestUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::ops::Deref;

#[test]
fn test_equality() {
    let i1 = SimpleInterval::new(1, 0, 100);
    let i2 = SimpleInterval::new(1, 0, 100);
    let i3 = SimpleInterval::new(1, 0, 100);
    let i4 = SimpleInterval::new(1, 1, 100);
    let i5 = SimpleInterval::new(1, 0, 200);

    assert!(i1 == i1);
    assert!(i1 == i2);
    assert!(i1 == i3);

    assert_ne!(i1, i4);
    assert_ne!(i1, i5);
    assert_ne!(i4, i5);
}

fn test_get_size(interval: SimpleInterval, expected_size: usize) {
    assert_eq!(
        interval.size(),
        expected_size,
        "Size incorrect for interval {:?} -> expected {} got {}",
        &interval,
        expected_size,
        interval.size()
    );
}

#[test]
fn test_interval_size() {
    test_get_size(SimpleInterval::new(1, 1, 1), 1);
    test_get_size(SimpleInterval::new(1, 1, 2), 2);
    test_get_size(SimpleInterval::new(1, 1, 10), 10);
    test_get_size(SimpleInterval::new(1, 2, 10), 9);
}

fn test_overlap(
    first_interval: &SimpleInterval,
    second_interval: SimpleInterval,
    expected_overlap_result: bool,
) {
    assert_eq!(
        first_interval.overlaps(&second_interval),
        expected_overlap_result,
        "Overlap returned incorrect result for intervals {:?} and {:?}",
        first_interval,
        second_interval
    )
}

#[test]
fn get_interval_overlap_data() {
    let standard_interval = SimpleInterval::new(1, 10, 20);
    let one_base_interval = SimpleInterval::new(1, 10, 10);

    test_overlap(&standard_interval, SimpleInterval::new(2, 10, 20), false);
    test_overlap(&standard_interval, SimpleInterval::new(1, 1, 5), false);
    test_overlap(&standard_interval, SimpleInterval::new(1, 1, 9), false);
    test_overlap(&standard_interval, SimpleInterval::new(1, 1, 10), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 1, 15), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 10, 10), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 10, 15), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 10, 20), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 15, 20), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 15, 25), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 20, 20), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 20, 25), true);
    test_overlap(&standard_interval, SimpleInterval::new(1, 21, 25), false);
    test_overlap(&standard_interval, SimpleInterval::new(1, 25, 30), false);
    test_overlap(&one_base_interval, SimpleInterval::new(2, 10, 10), false);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 1, 5), false);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 1, 9), false);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 1, 10), true);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 10, 10), true);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 10, 15), true);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 11, 15), false);
    test_overlap(&one_base_interval, SimpleInterval::new(1, 15, 20), false);
    test_overlap(&standard_interval, standard_interval.clone(), true)
}

#[test]
fn overlaps_with_margin() {
    let standard_interval = SimpleInterval::new(1, 10, 20);
    let middle_interval = SimpleInterval::new(1, 100, 200);

    test_overlaps_with_margin(
        &standard_interval,
        SimpleInterval::new(2, 10, 20),
        100,
        false,
    );
    test_overlaps_with_margin(&standard_interval, SimpleInterval::new(1, 1, 15), 0, true);
    test_overlaps_with_margin(&standard_interval, SimpleInterval::new(1, 30, 50), 9, false);
    test_overlaps_with_margin(&standard_interval, SimpleInterval::new(1, 30, 50), 10, true);
    test_overlaps_with_margin(&middle_interval, SimpleInterval::new(1, 50, 99), 0, false);
    test_overlaps_with_margin(&middle_interval, SimpleInterval::new(1, 50, 90), 9, false);
    test_overlaps_with_margin(&middle_interval, SimpleInterval::new(1, 50, 90), 10, true);
}

fn test_overlaps_with_margin(
    first_interval: &SimpleInterval,
    second_interval: SimpleInterval,
    margin: usize,
    expected_overlap_result: bool,
) {
    assert_eq!(first_interval.overlaps_with_margin(&second_interval, margin), expected_overlap_result, "Overlap with margin returned incorrect result for first {:?} and {:?} with margin {}, expected {}", &first_interval, &second_interval, margin, expected_overlap_result)
}

#[test]
fn get_interval_contains_data() {
    let containing_interval = SimpleInterval::new(1, 10, 20);
    test_contains(&containing_interval, SimpleInterval::new(2, 10, 20), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 1, 5), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 1, 10), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 5, 15), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 9, 10), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 9, 20), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 10, 10), true);
    test_contains(&containing_interval, SimpleInterval::new(1, 10, 15), true);
    test_contains(&containing_interval, SimpleInterval::new(1, 10, 20), true);
    test_contains(&containing_interval, SimpleInterval::new(1, 10, 21), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 15, 25), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 20, 20), true);
    test_contains(&containing_interval, SimpleInterval::new(1, 20, 21), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 20, 25), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 21, 25), false);
    test_contains(&containing_interval, SimpleInterval::new(1, 25, 30), false);
    test_contains(&containing_interval, containing_interval.clone(), true);
}

fn test_contains(
    first_interval: &SimpleInterval,
    second_interval: SimpleInterval,
    expected_contains_result: bool,
) {
    assert_eq!(
        first_interval.contains(&second_interval),
        expected_contains_result,
        "contains returned incorrect result for intervals {:?} and {:?}",
        first_interval,
        second_interval
    );
}

#[test]
fn test_not_contiguous_loc() {
    let loc1 = SimpleInterval::new(1, 10, 20);
    let loc2 = SimpleInterval::new(1, 22, 30);
    let loc3 = SimpleInterval::new(1, 31, 40);
    let loc4 = SimpleInterval::new(2, 20, 30);
    assert!(loc1.merge_with_contiguous(&loc2).is_err());
    assert!(loc1.merge_with_contiguous(&loc3).is_err());
    assert!(loc1.merge_with_contiguous(&loc4).is_err());
}

#[test]
fn test_merge_contiguous() {
    let loc1 = SimpleInterval::new(1, 10, 20);
    let loc2 = SimpleInterval::new(1, 20, 30);
    let loc3 = SimpleInterval::new(1, 21, 30);
    assert!(loc1.merge_with_contiguous(&loc2).is_ok());
    assert!(loc1.merge_with_contiguous(&loc3).is_ok());
}

#[test]
fn expand_within_contig() {
    let CONTIG_LENGTH = 10000;
    test_expand_within_contig(
        SimpleInterval::new(1, 5, 10),
        0,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 5, 10),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, 5, 10),
        1,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 4, 11),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, 1, 10),
        10,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 0, 20),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, 10, 20),
        10,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 0, 30),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, 10, 20),
        9,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 1, 29),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, 30, 40),
        5,
        CONTIG_LENGTH,
        SimpleInterval::new(1, 25, 45),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, CONTIG_LENGTH - 10, CONTIG_LENGTH),
        10,
        CONTIG_LENGTH,
        SimpleInterval::new(1, CONTIG_LENGTH - 20, CONTIG_LENGTH),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, CONTIG_LENGTH - 20, CONTIG_LENGTH - 10),
        11,
        CONTIG_LENGTH,
        SimpleInterval::new(1, CONTIG_LENGTH - 31, CONTIG_LENGTH),
    );
    test_expand_within_contig(
        SimpleInterval::new(1, CONTIG_LENGTH - 20, CONTIG_LENGTH - 10),
        10,
        CONTIG_LENGTH,
        SimpleInterval::new(1, CONTIG_LENGTH - 30, CONTIG_LENGTH),
    );
}

fn test_expand_within_contig(
    starting_interval: SimpleInterval,
    padding: usize,
    contig_len: usize,
    expected_interval: SimpleInterval,
) {
    assert_eq!(
        starting_interval.expand_within_contig(padding, contig_len),
        expected_interval
    );
}
