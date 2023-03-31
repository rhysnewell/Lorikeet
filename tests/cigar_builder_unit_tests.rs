#![allow(
    non_upper_case_globals,
    non_snake_case
)]







use lorikeet_genome::reads::cigar_builder::CigarBuilder;


use rust_htslib::bam::record::{CigarString};
use std::convert::TryFrom;

fn test_simple_concatenation(cigar_element_strings: Vec<&str>) {
    let mut builder = CigarBuilder::new(true);
    for element_str in cigar_element_strings.iter() {
        builder
            .add(CigarString::try_from(*element_str).unwrap().0[0])
            .unwrap();
    }

    let expected = cigar_element_strings.join("");
    assert_eq!(
        builder
            .make(false)
            .unwrap_or_else(|err| panic!("Unhandled error in cigar builder: {:?}", err))
            .to_string(),
        expected.to_string()
    );
}

#[test]
fn cigar_algebra() {
    let leading_clips = vec![Vec::new(), vec!["10H"], vec!["10S"], vec!["10H", "10S"]];

    let middle_operators = vec![
        vec!["10M"],
        vec!["10M", "10I", "10M"],
        vec!["10M", "10D", "10M"],
    ];

    let trailing_clips = vec![Vec::new(), vec!["10H"], vec!["10S"], vec!["10S", "10H"]];

    for leading_string in leading_clips.iter() {
        for middle_string in middle_operators.iter() {
            for trailing_string in trailing_clips.iter() {
                let mut all_elements = Vec::new();
                all_elements.extend(leading_string);
                all_elements.extend(middle_string);
                all_elements.extend(trailing_string);
                test_simple_concatenation(all_elements)
            }
        }
    }
}

fn test_initial_and_final_deletion(cigar_elements_strings: Vec<&str>, expected: &str) {
    let mut builder = CigarBuilder::new(true);
    for element_string in cigar_elements_strings {
        builder
            .add(CigarString::try_from(element_string).unwrap().0[0])
            .unwrap();
    }

    assert_eq!(
        builder
            .make(false)
            .unwrap_or_else(|err| panic!("{:?}", err))
            .to_string(),
        expected.to_string()
    );
}

#[test]
fn initial_and_final_deletion() {
    test_initial_and_final_deletion(vec!["10M", "10D"], "10M");
    test_initial_and_final_deletion(vec!["10D", "10M"], "10M");
    test_initial_and_final_deletion(vec!["10H", "10D", "10M"], "10H10M");
    test_initial_and_final_deletion(vec!["10S", "10D", "10M"], "10S10M");
    test_initial_and_final_deletion(vec!["10S", "10D", "10M", "10S"], "10S10M10S");
    test_initial_and_final_deletion(vec!["10M", "10D", "10S"], "10M10S");
    test_initial_and_final_deletion(vec!["10M", "10D", "10H"], "10M10H");
    test_initial_and_final_deletion(vec!["10S", "10M", "10D", "10H"], "10S10M10H");
}

fn test_retain_deletions(cigar_elements_strings: Vec<&str>, expected: &str) {
    let mut builder = CigarBuilder::new(false);
    for element_string in cigar_elements_strings {
        builder
            .add(CigarString::try_from(element_string).unwrap().0[0])
            .unwrap();
    }

    assert_eq!(
        builder
            .make(false)
            .unwrap_or_else(|err| panic!("{:?}", err))
            .to_string(),
        expected.to_string()
    );
}

#[test]
fn retain_deletions() {
    test_retain_deletions(vec!["10M", "10D"], "10M10D");
    test_retain_deletions(vec!["10D", "10M"], "10D10M");
    test_retain_deletions(vec!["10H", "10D", "10M"], "10H10D10M");
    test_retain_deletions(vec!["10S", "10D", "10M"], "10S10D10M");
    test_retain_deletions(vec!["10S", "10D", "10M", "10S"], "10S10D10M10S");
    test_retain_deletions(vec!["10M", "10D", "10S"], "10M10D10S");
    test_retain_deletions(vec!["10M", "10D", "10H"], "10M10D10H");
    test_retain_deletions(vec!["10S", "10M", "10D", "10H"], "10S10M10D10H");
}

fn test_merge_consecutive(cigar_elements_strings: Vec<&str>, expected: &str) {
    let mut builder = CigarBuilder::new(true);
    for element_string in cigar_elements_strings {
        builder.add(CigarString::try_from(element_string).unwrap().0[0]);
    }

    assert_eq!(
        builder
            .make(false)
            .unwrap_or_else(|err| panic!("{:?}", err))
            .to_string(),
        expected.to_string()
    );
}

#[test]
fn merge_consecutive() {
    test_merge_consecutive(vec!["10H", "10H", "10M"], "20H10M");
    test_merge_consecutive(vec!["10S", "10M", "10M"], "10S20M");
    test_merge_consecutive(vec!["10S", "10M", "10S", "10S"], "10S10M20S");
    test_merge_consecutive(
        vec!["10S", "10M", "10I", "10I", "10I", "10S", "10H"],
        "10S10M30I10S10H",
    );
    test_merge_consecutive(
        vec!["10S", "10S", "10M", "10M", "10I", "10I", "10S", "10H"],
        "20S20M20I10S10H",
    );
}

#[test]
fn tricky() {
    test_merge_consecutive(vec!["10H", "10H", "10D", "10D", "10M"], "20H10M");
}

#[test]
fn indel_sandwich() {
    test_merge_consecutive(vec!["10M", "10I", "10D", "10M"], "10M10D10I10M");
    test_merge_consecutive(vec!["10M", "10D", "10I", "10M"], "10M10D10I10M");
    test_merge_consecutive(vec!["10M", "10I", "10D", "10I", "10M"], "10M10D20I10M");
    test_merge_consecutive(
        vec!["10M", "10I", "10D", "10I", "10D", "10I", "10M"],
        "10M20D30I10M",
    );
    test_merge_consecutive(
        vec!["10M", "10I", "10D", "10I", "10M", "10D", "10I", "10M"],
        "10M10D20I10M10D10I10M",
    );

    //does the indel sandwich logic interfere with removing leading/trailing deletions
    test_merge_consecutive(vec!["10D", "10I", "10M"], "10I10M");
    test_merge_consecutive(vec!["10M", "10I", "10D"], "10M10I");
    test_merge_consecutive(vec!["10M", "10D", "10I"], "10M10I");
    test_merge_consecutive(vec!["10M", "10D", "10I", "10S"], "10M10I10S");
    test_merge_consecutive(vec!["10S", "10D", "10I", "10M"], "10S10I10M");
    test_merge_consecutive(vec!["10S", "10I", "10D", "10I", "10M"], "10S20I10M");
}

fn test_invalid(cigar_elements_strings: Vec<&str>) {
    let mut builder = CigarBuilder::new(true);
    for element_string in cigar_elements_strings {
        builder.add(CigarString::try_from(element_string).unwrap().0[0]);
    }

    assert!(builder.make(false).is_err());
}

#[test]
fn invalid() {
    // completely soft-clipped
    test_invalid(vec!["10S"]);
    test_invalid(vec!["10S", "10S"]);

    // also completely clipped
    test_invalid(vec!["10S", "10D"]);
    test_invalid(vec!["10S", "10D", "10S"]);
    test_invalid(vec!["10S", "10D", "10D", "10S"]);

    // wrong order of hard and soft clips
    test_invalid(vec!["10S", "10H", "10M"]);
    test_invalid(vec!["10M", "10H", "10S"]);

    // clipping in middle of read
    test_invalid(vec!["10M", "10H", "10M"]);
    test_invalid(vec!["10M", "10S", "10M"]);
}

fn test_removed_deletions(
    cigar_elements_strings: Vec<&str>,
    removed_leading: u32,
    removed_trailing: u32,
) {
    let mut builder = CigarBuilder::new(true);
    for element_string in cigar_elements_strings {
        builder
            .add(CigarString::try_from(element_string).unwrap().0[0])
            .unwrap();
    }

    builder.make(false);
    assert_eq!(
        builder.get_leading_deletion_bases_removed(),
        removed_leading
    );
    assert_eq!(
        builder.get_trailing_deletion_bases_removed(),
        removed_trailing
    )
}

#[test]
fn removed_deletions() {
    test_removed_deletions(vec!["10M"], 0, 0);
    test_removed_deletions(vec!["10S", "10M"], 0, 0);
    test_removed_deletions(vec!["10M", "10S"], 0, 0);
    test_removed_deletions(vec!["10M", "10I", "10D", "10M"], 0, 0);
    test_removed_deletions(vec!["10M", "10D", "10I", "10M"], 0, 0);

    test_removed_deletions(vec!["10D", "10I", "10M"], 10, 0);
    test_removed_deletions(vec!["10D", "10D", "10I", "10M"], 20, 0);
    test_removed_deletions(vec!["10D", "10D", "10I", "10D", "10M"], 30, 0);
    test_removed_deletions(vec!["10S", "10D", "10D", "10I", "10D", "10M"], 30, 0);

    test_removed_deletions(vec!["10M", "10I", "10D"], 0, 10);
    test_removed_deletions(vec!["10M", "10D", "10I"], 0, 10);
    test_removed_deletions(vec!["10M", "10D", "10I", "10D"], 0, 20);
    test_removed_deletions(vec!["10M", "10D", "10I", "10D", "10S", "10H"], 0, 20);

    test_removed_deletions(
        vec![
            "10H", "10S", "10D", "10M", "10D", "10I", "10D", "10S", "10H",
        ],
        10,
        20,
    );
}

fn test_removed_deletions_two_makes(
    cigar_elements_strings1: Vec<&str>,
    cigar_elements_strings2: Vec<&str>,
    removed_leading: u32,
    removed_trailing: u32,
) {
    let mut builder = CigarBuilder::new(true);
    for element_string in cigar_elements_strings1 {
        builder
            .add(CigarString::try_from(element_string).unwrap().0[0])
            .unwrap();
    }

    builder.make(false);

    for element_string in cigar_elements_strings2 {
        builder
            .add(CigarString::try_from(element_string).unwrap().0[0])
            .unwrap();
    }

    builder.make(false);
    assert_eq!(
        builder.get_leading_deletion_bases_removed(),
        removed_leading
    );
    assert_eq!(
        builder.get_trailing_deletion_bases_removed(),
        removed_trailing
    )
}

#[test]
fn removed_deletions_two_makes() {
    test_removed_deletions_two_makes(vec!["10M"], vec!["10M"], 0, 0);
    test_removed_deletions_two_makes(vec!["10M", "10I"], vec!["10D", "10M"], 0, 0);
    test_removed_deletions_two_makes(vec!["10M", "10D"], vec!["10I", "10M"], 0, 0);

    test_removed_deletions_two_makes(vec!["10D", "10I"], vec!["10M"], 10, 0);
    test_removed_deletions_two_makes(vec!["10D", "10D", "10I"], vec!["10D", "10M"], 30, 0);

    test_removed_deletions_two_makes(
        vec!["10H", "10S", "10D", "10M"],
        vec!["10D", "10I", "10D", "10S", "10H"],
        10,
        20,
    );
}
