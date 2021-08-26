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
use lorikeet_genome::utils::simple_interval::Locatable;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::ops::Deref;

struct ReadClipperUnitTest {
    cigar_list: Vec<CigarString>,
    // 9 is the minimum necessary number to try all combinations of cigar types with guarantee of clipping an element with length = 2
    // and there are already 3 cigar elements on the basic cigar
    maximum_cigar_elements: usize,
}

impl ReadClipperUnitTest {
    pub fn new() -> ReadClipperUnitTest {
        let mut cigar_list = ReadClipperTestUtils::generate_cigar_list(6, false);
        let additional_cigar = CigarString::try_from("2M3I5M").unwrap();
        cigar_list.push(additional_cigar);

        Self {
            cigar_list,
            maximum_cigar_elements: 6,
        }
    }
}

#[test]
fn test_hard_clip_both_ends_by_reference() {
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let aln_start = read.get_start() as i64;
        let aln_end = read.get_end() as i64;
        let read_length = aln_start - aln_end;
        for i in (read_length / 2) + 1..=0 {
            let clipped_read = ReadClipper::new(&read)
                .hard_clip_both_ends_by_reference_coordinates(
                    (aln_start + i) as usize,
                    (aln_end - i) as usize,
                );
            assert!(
                clipped_read.get_start() >= (aln_start + i) as usize,
                "Clipped alignment is less than original read (minus {}): {:?} -> {:?} {} > {}",
                i,
                read.read.cigar().to_string(),
                clipped_read.read.cigar().to_string(),
                aln_start + 1,
                clipped_read.get_start()
            );
            assert!(
                clipped_read.get_end() <= (aln_end - i) as usize,
                "Clipped alignment is greater than original read (minus {}): {:?} -> {:?} {} < {}, start {} ref length {} original ref length {}",
                i,
                read.read.cigar().to_string(),
                clipped_read.read.cigar().to_string(),
                aln_end + i,
                clipped_read.get_end(),
                clipped_read.get_start(),
                CigarUtils::get_reference_length(clipped_read.read.cigar().deref()),
                CigarUtils::get_reference_length(read.read.cigar().deref())
            );
            // assert_unclipped_limits(&read, &clipped_read);
        }
    }
}

#[test]
fn test_hard_clip_by_reference_coordinates() {
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let start = read.get_soft_start().unwrap_or(0);
        let stop = read.get_soft_end();
        for i in start..=stop {
            let clip_left =
                ReadClipper::new(&read).hard_clip_by_reference_coordinates(None, Some(i));

            if !clip_left.is_empty() {
                assert!(
                    clip_left.get_start() >= min(read.get_end(), i),
                    "Clipped alignment start ({}) is less the expected ({}): {} -> {}",
                    clip_left.get_start(),
                    min(read.get_end(), i + 1),
                    read.read.cigar().to_string(),
                    clip_left.read.cigar().to_string()
                );
                assert_ref_alignment_consistent(&clip_left);
                assert_read_length_consistent(&clip_left);
            }

            let clip_right =
                ReadClipper::new(&read).hard_clip_by_reference_coordinates(Some(i), None);
            if !clip_right.is_empty() && clip_right.get_start() <= clip_right.get_end() {
                // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                assert!(
                    clip_right.get_end() <= max(read.get_start(), i),
                    "Clipped alignment end ({}) is greater than expected ({}): {} -> {}",
                    clip_right.get_end(),
                    max(read.get_start(), i - 1),
                    read.read.cigar().to_string(),
                    clip_right.read.cigar().to_string()
                );
                assert_ref_alignment_consistent(&clip_right);
                assert_read_length_consistent(&clip_right);
            }
        }
    }
}

#[test]
fn test_hard_clip_by_reference_coordinates_left_tail() {
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let aln_start = read.get_start();
        let aln_end = read.get_end();
        if read.get_soft_start().unwrap() == aln_start {
            for i in aln_start..=aln_end {
                let clip_left =
                    ReadClipper::new(&read).hard_clip_by_reference_coordinates(None, Some(i));
                if !clip_left.is_empty() {
                    assert!(
                        clip_left.get_start() >= i,
                        "Clipped alignment start ({}) is less the expected ({}): {} -> {}",
                        clip_left.get_start(),
                        i,
                        read.read.cigar().to_string(),
                        clip_left.read.cigar().to_string()
                    );
                    assert_ref_alignment_consistent(&clip_left);
                    assert_read_length_consistent(&clip_left);
                }
            }
        }
    }
}

#[test]
fn test_hard_clip_by_reference_coordinates_right_tail() {
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let aln_start = read.get_start();
        let aln_end = read.get_end();
        if read.get_soft_end() == aln_start {
            for i in aln_start..=aln_end {
                let clip_right =
                    ReadClipper::new(&read).hard_clip_by_reference_coordinates(Some(i), None);
                if !clip_right.is_empty() && clip_right.get_start() <= clip_right.get_end() {
                    assert!(
                        clip_right.get_end() <= i,
                        "Clipped alignment end ({}) is greater than the expected ({}): {} -> {}",
                        clip_right.get_start(),
                        i - 1,
                        read.read.cigar().to_string(),
                        clip_right.read.cigar().to_string()
                    );
                    assert_ref_alignment_consistent(&clip_right);
                    assert_read_length_consistent(&clip_right);
                }
            }
        }
    }
}

#[test]
fn test_hard_clip_low_qual_ends() {
    let LOW_QUAL: u8 = 2;
    let HIGH_QUAL: u8 = 30;

    // create a read for every cigar permutation
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let mut read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let read_length = read.len();
        let mut quals = vec![0; read_length];

        let qname = read.read.qname().to_vec();
        let seq = read.read.seq().as_bytes();
        let cigar = read.read.cigar().take();

        for n_low_qual_bases in 0..read_length {
            // create a read with n_low_qual_bases in the tail
            quals = vec![HIGH_QUAL; read_length];
            for add_left in 0..n_low_qual_bases {
                quals[add_left] = LOW_QUAL;
            }

            read.read.set(&qname, Some(&cigar), &seq, &quals);

            let clip_left = ReadClipper::new(&read).hard_clip_low_qual_ends(LOW_QUAL);
            assert_no_low_qual_bases(&clip_left, LOW_QUAL);

            // create a read with n_low_qual_bases in the right tail
            quals = vec![HIGH_QUAL; read_length];
            for add_right in 0..n_low_qual_bases {
                quals[read_length - add_right - 1] = LOW_QUAL;
            }

            read.read.set(&qname, Some(&cigar), &seq, &quals);

            let clip_right = ReadClipper::new(&read).hard_clip_low_qual_ends(LOW_QUAL);
            assert_no_low_qual_bases(&clip_right, LOW_QUAL);

            // create a read with n_low_qual_bases in both tails
            if n_low_qual_bases <= read_length / 2 {
                // create a read with n_low_qual_bases in the right tail
                quals = vec![HIGH_QUAL; read_length];
                for add_both in 0..n_low_qual_bases {
                    quals[add_both] = LOW_QUAL;
                    quals[read_length - add_both - 1] = LOW_QUAL;
                }

                read.read.set(&qname, Some(&cigar), &seq, &quals);

                let clip_right = ReadClipper::new(&read).hard_clip_low_qual_ends(LOW_QUAL);
                assert_no_low_qual_bases(&clip_right, LOW_QUAL);
            }
        }
    }
}

#[test]
fn test_hard_clip_soft_clipped_bases() {
    // create a read for every cigar permutation
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let clipped_read = ReadClipper::new(&read).hard_clip_soft_clipped_bases();

        if !clipped_read.is_empty() {
            let original = CigarCounter::new(&read);
            let clipped = CigarCounter::new(&clipped_read);

            assert_unclipped_limits(&read, &clipped_read);
            original.assert_hard_clipping_soft_clips(&clipped);
        }
    }
}

#[test]
fn test_revert_soft_clipped_bases() {
    // create a read for every cigar permutation
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let leading_soft_clips = leading_cigar_element_length(&cigar, Some(Cigar::SoftClip(0)));
        let tail_soft_clips = leading_cigar_element_length(
            &CigarString(cigar.0.iter().cloned().rev().collect::<Vec<Cigar>>()),
            Some(Cigar::SoftClip(0)),
        );
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let unclipped = ReadClipper::new(&read).revert_soft_clipped_bases();

        assert_unclipped_limits(&read, &unclipped);

        if leading_soft_clips > 0 || tail_soft_clips > 0 {
            let expected_start = read.get_start() as i64 - leading_soft_clips as i64;
            let expected_end = read.get_end() + tail_soft_clips as usize;
            assert_eq!(unclipped.get_start() as i64, expected_start);
            assert_eq!(
                unclipped.get_end(),
                expected_end,
                "{} -> {}",
                read.read.cigar().to_string(),
                unclipped.read.cigar().to_string()
            );
        } else {
            assert_eq!(
                read.read.cigar().to_string(),
                unclipped.read.cigar().to_string()
            );
        }
    }
}

#[test]
fn test_revert_entirely_soft_clipped_reads() {
    let read = ReadClipperTestUtils::make_read_from_str("2H1S3H", 0);
    let clipped_read = ReadClipper::new(&read).revert_soft_clipped_bases();

    assert_eq!(clipped_read.get_start(), read.get_soft_start().unwrap());
}

#[test]
fn test_soft_clip_both_ends_by_reference_coordinates() {
    // create a read for every cigar permutation
    let mut read_clipper_unit_test = ReadClipperUnitTest::new();

    for cigar in read_clipper_unit_test.cigar_list {
        let read = ReadClipperTestUtils::make_read_from_cigar(cigar, 0);
        let aln_start = read.get_start() as i64;
        let aln_end = read.get_end() as i64;
        let read_length = aln_start - aln_end;
        for i in (read_length / 2) + 1..=0 {
            let clipped_read = ReadClipper::new(&read)
                .soft_clip_both_ends_by_reference_coordinates(
                    (aln_start + i) as usize,
                    (aln_end - i) as usize,
                );
            assert!(
                clipped_read.get_start() >= (aln_start + i) as usize,
                "Clipped alignment is less than original read (minus {}): {:?} -> {:?}",
                i,
                read.read.cigar().to_string(),
                clipped_read.read.cigar().to_string()
            );
            assert!(
                clipped_read.get_end() <= (aln_end - i) as usize,
                "Clipped alignment is greater than original read (minus {}): {:?} -> {:?}",
                i,
                read.read.cigar().to_string(),
                clipped_read.read.cigar().to_string()
            );
            // assert_unclipped_limits(&read, &clipped_read);
        }
    }
}

#[test]
fn make_revert_soft_clips_before_contig() {
    for soft_start in vec![-10, -1, 0] {
        for alignment_start in vec![1, 10] {
            test_revert_soft_clipped_bases_before_start_of_contig(soft_start, alignment_start)
        }
    }
}

fn test_revert_soft_clipped_bases_before_start_of_contig(soft_start: i64, alignment_start: i64) {
    let n_matches = 10;
    let n_soft = -1 * (soft_start - alignment_start);
    let cigar = format!("{}S{}M", n_soft, n_matches);
    let mut read = ReadClipperTestUtils::make_read_from_str(cigar.as_str(), 0);
    read.read.set_pos(alignment_start);
    assert_eq!(read.get_soft_start_i64(), soft_start);
    assert_eq!(read.get_start() as i64, alignment_start);
    assert_eq!(read.read.cigar().to_string(), cigar);

    let reverted = ReadClipper::new(&read).revert_soft_clipped_bases();
    let expected_alignment_start = 0;
    let expected_cigar = format!("{}H{}M", 1 - soft_start, read.get_end());
    assert_eq!(reverted.get_soft_start_i64(), expected_alignment_start);
    assert_eq!(reverted.get_start() as i64, expected_alignment_start);
    assert_eq!(reverted.read.cigar().to_string(), expected_cigar);
}

fn leading_cigar_element_length(cigar: &CigarString, operator: Option<Cigar>) -> u32 {
    for cigar_element in cigar.iter() {
        if CigarUtils::cigar_elements_are_same_type(cigar_element, &operator) {
            return cigar_element.len();
        }

        match cigar_element {
            Cigar::HardClip(_) => continue,
            _ => return 0,
        }
    }
    return 0;
}

fn assert_no_low_qual_bases(read: &BirdToolRead, low_qual: u8) {
    if !read.is_empty() {
        let quals = read.read.qual();
        for i in 0..quals.len() {
            assert!(
                quals[i] > low_qual,
                "Found low qual ({}) base after hard clipping. Position: {} -- {}",
                low_qual,
                i,
                read.read.cigar().to_string()
            );
        }
    }
}

/**
 * Asserts that clipping doesn't change the getUnclippedStart / getUnclippedEnd
 *
 * @param original original read
 * @param clipped clipped read
 */
fn assert_unclipped_limits(original: &BirdToolRead, clipped: &BirdToolRead) {
    if original
        .read
        .cigar()
        .iter()
        .any(|el| !CigarUtils::is_clipping(el))
    {
        assert_eq!(
            original.get_unclipped_start(),
            clipped.get_unclipped_start(),
            "Unclipped start locations did not match {} -> {}",
            original.read.cigar().to_string(),
            clipped.read.cigar().to_string()
        );
        assert_eq!(
            original.get_unclipped_end(),
            clipped.get_unclipped_end(),
            "Unclipped end locations did not match {} -> {}",
            original.read.cigar().to_string(),
            clipped.read.cigar().to_string()
        );
    }
}

/**
* Asserts that the length of alignment on the reference is consistent with the CIGAR
* after clipping
*
* @param clippedRead input read
* */
fn assert_ref_alignment_consistent(clipped_read: &BirdToolRead) {
    let cigar_ref_length = CigarUtils::get_reference_length(clipped_read.read.cigar().deref());
    let read_ref_length = if clipped_read.read.is_unmapped() {
        0
    } else {
        clipped_read.get_length_on_reference()
    };

    assert_eq!(
        cigar_ref_length,
        read_ref_length as u32,
        "cigar {}, start {} end {} unmapped? {}",
        clipped_read.read.cigar().to_string(),
        clipped_read.get_start(),
        clipped_read.get_end(),
        clipped_read.read.is_unmapped()
    );
}

/**
 * Asserts that the length of the read is consistent with the CIGAR after clipping
 *
 * @param clippedRead input read
 * */
fn assert_read_length_consistent(clipped_read: &BirdToolRead) {
    let cigar_read_length = CigarUtils::get_read_length(clipped_read.read.cigar().deref());
    let read_read_length = clipped_read.len();
    assert_eq!(
        cigar_read_length,
        read_read_length as u32,
        "Read {:?} start {} stop {} seq {:?}",
        clipped_read.read.cigar(),
        clipped_read.get_start(),
        clipped_read.get_end(),
        clipped_read.read.seq().as_bytes()
    );
}

struct CigarCounter {
    counter: HashMap<Cigar, u32>,
}

impl CigarCounter {
    fn new(read: &BirdToolRead) -> Self {
        let mut counter = HashMap::new();

        for cigar_element in read.read.cigar().iter() {
            let count = counter
                .entry(CigarUtils::cigar_from_element_and_length(cigar_element, 0))
                .or_insert(0);
            *count += cigar_element.len();
        }

        Self { counter }
    }

    fn assert_hard_clipping_soft_clips(&self, clipped: &Self) {
        let zero = 0;
        for cigar_operator in self.counter.keys() {
            if CigarUtils::is_clipping(cigar_operator) {
                let counter_total = (self.counter.get(&Cigar::HardClip(0)).unwrap_or(&zero)
                    + self.counter.get(&Cigar::SoftClip(0)).unwrap_or(&zero));
                let clipped_hard = *clipped.counter.get(&Cigar::HardClip(0)).unwrap_or(&zero);
                let clipped_soft = *clipped.counter.get(&Cigar::SoftClip(0)).unwrap_or(&zero);

                assert_eq!(counter_total, clipped_hard);
                assert_eq!(clipped_soft, 0);
            } else {
                assert_eq!(
                    self.counter
                        .get(&CigarUtils::cigar_from_element_and_length(
                            cigar_operator,
                            0
                        ))
                        .unwrap_or(&zero),
                    clipped
                        .counter
                        .get(&CigarUtils::cigar_from_element_and_length(
                            cigar_operator,
                            0
                        ))
                        .unwrap_or(&zero),
                )
            }
        }
    }
}
