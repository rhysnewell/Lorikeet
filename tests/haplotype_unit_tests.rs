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
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rust_htslib::bam::record::{Cigar, CigarString};
use std::convert::TryFrom;

// fn basic_insert_test(reference: &str, alt: &str, loc: usize, cigar: Vec<Cigar>, hap: &str, new_hap: &str) {
//     let h = Haplotype::new(hap.as_bytes(), false);
//     let h1_ref_allele = ByteArrayAllele::new(reference.as_bytes(), true);
//     let h1_alt_allele = ByteArrayAllele::new(alt.as_bytes(), false);
//
//     let mut alleles = vec![h1_ref_allele, h1_alt_allele];
//     let mut vc = VariantContext::build(1, loc, loc + alleles[0].get_bases().len() - 1, alleles);
//
//     h.set_alignment_start_hap_wrt_ref(0);
//     h.set_cigar(cigar);
//     let h1 = h.inser
// }

fn test_trim(
    full: Haplotype<SimpleInterval>,
    trim_to: SimpleInterval,
    expected: Option<Haplotype<SimpleInterval>>,
) {
    let actual = full.trim(trim_to.clone());
    match expected {
        None => {
            assert!(actual
                .unwrap_or_else(|_| panic!("Unhandled error in cigar builder"))
                .is_none())
        }
        Some(expected) => {
            match actual.unwrap_or_else(|_| panic!("Unhandled error in cigar builder")) {
                None => assert!(false),
                Some(actual) => {
                    assert_eq!(actual.get_bases(), expected.get_bases());
                    assert_eq!(actual.get_start_position(), trim_to.get_start());
                    assert_eq!(actual.get_stop_position(), trim_to.get_end());
                    assert_eq!(actual.get_cigar(), expected.get_cigar());
                    assert_eq!(
                        actual.get_alignment_start_hap_wrt_ref(),
                        expected.get_alignment_start_hap_wrt_ref()
                    );
                }
            }
        }
    }
}

fn make_hcf_for_cigar(bases: &str, cigar: &str) -> Haplotype<SimpleInterval> {
    let mut h = Haplotype::new(bases.as_bytes(), false);
    h.set_cigar(CigarString::try_from(cigar).unwrap().0);
    return h;
}

fn test_trim_leading_and_trailing_insertions(
    bases: &str,
    cigar: &str,
    start: usize,
    stop: usize,
    expected_trimmed_cigar: &str,
    expected_trimmed_bases: &str,
) {
    let offset = 10;
    let untrimmed_cigar = CigarString::try_from(cigar).unwrap();
    let mut haplotype = Haplotype::new(bases.as_bytes(), false);
    haplotype.set_genome_location(SimpleInterval::new(
        0,
        offset,
        offset + CigarUtils::get_reference_length(&untrimmed_cigar) as usize - 1,
    ));
    haplotype.set_alignment_start_hap_wrt_ref(offset);
    haplotype.set_cigar(untrimmed_cigar.0);

    let trimmed = haplotype
        .trim(SimpleInterval::new(0, offset + start, offset + stop))
        .unwrap_or_else(|_| panic!("Unhandled error in cigar builder"))
        .unwrap();
    assert_eq!(
        trimmed.get_cigar().to_string(),
        expected_trimmed_cigar.to_string()
    );
    assert_eq!(trimmed.get_bases(), expected_trimmed_bases.as_bytes());
}

#[test]
fn test_consolidate_cigar() {
    assert_eq!(
        make_hcf_for_cigar("AGCT", "4M")
            .get_consolidated_padded_cigar(0)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("4M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "4M")
            .get_consolidated_padded_cigar(1)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("5M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1M")
            .get_consolidated_padded_cigar(0)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M2I1M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1M")
            .get_consolidated_padded_cigar(1)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M2I2M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1M")
            .get_consolidated_padded_cigar(2)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M2I3M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1I")
            .get_consolidated_padded_cigar(0)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M3I").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1I")
            .get_consolidated_padded_cigar(1)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M3I1M").unwrap()
    );
    assert_eq!(
        make_hcf_for_cigar("AGCT", "1M1I1I1I")
            .get_consolidated_padded_cigar(2)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder")),
        CigarString::try_from("1M3I2M").unwrap()
    );
}

#[test]
fn test_trimming() {
    let loc = SimpleInterval::new(0, 10, 20);
    let full_bases = "ACGTAACCGGT";
    for trim_start in loc.get_start()..loc.get_end() {
        for trim_stop in trim_start..=loc.get_end() {
            let start = trim_start - loc.get_start();
            let stop = start + (trim_stop - trim_start) + 1;
            let trimmed_loc =
                SimpleInterval::new(0, start + loc.get_start(), stop + loc.get_start() - 1);
            let expected_bases = &full_bases.as_bytes()[start..stop];
            let mut full = Haplotype::new(full_bases.as_bytes(), false);
            full.set_genome_location(loc.clone());
            let mut trimmed = Haplotype::new(expected_bases, false);
            trimmed.set_genome_location(trimmed_loc.clone());

            let hap_start = 10;
            full.set_alignment_start_hap_wrt_ref(hap_start);
            full.set_cigar(
                CigarString::try_from(format!("{}M", full.len()).as_bytes())
                    .unwrap()
                    .0,
            );

            trimmed.set_alignment_start_hap_wrt_ref(hap_start + start);
            trimmed.set_cigar(
                CigarString::try_from(format!("{}M", trimmed.len()).as_bytes())
                    .unwrap()
                    .0,
            );
            test_trim(full, trimmed_loc, Some(trimmed));
        }
    }

    let mut full = Haplotype::new("ACT".as_bytes(), false);
    full.set_genome_location(SimpleInterval::new(0, 10, 14));
    full.set_alignment_start_hap_wrt_ref(10);
    full.set_cigar(CigarString::try_from("1M2D2M").unwrap().0);
    test_trim(full.clone(), SimpleInterval::new(0, 11, 12), None);
    test_trim(full.clone(), SimpleInterval::new(0, 10, 12), None);
    test_trim(full, SimpleInterval::new(0, 11, 13), None);
}

#[test]
fn test_trimming_data_with_leading_and_trailing_insertions() {
    // order is: String bases, String cigar, int start, int stop, String expectedTrimmedCigar, String expectedTrimmedBases
    // start (inclusive ie the trimmed haplotype includes it) and stop (inclusive ie the trimmed haplotype includes it) are
    // zero-based ref coordinates relative to the alignment start of the untrimmed haplotype
    // the first GT is an insertion in the following examples
    test_trim_leading_and_trailing_insertions("ACGTACGT", "2M2I4M", 1, 5, "1M2I4M", "CGTACGT"); // no leading insertion
    test_trim_leading_and_trailing_insertions("ACGTACGT", "2M2I4M", 2, 5, "4M", "ACGT"); // leading insertion removed
    test_trim_leading_and_trailing_insertions("ACGTACGT", "2M2I4M", 3, 5, "3M", "CGT"); // no leading insertion
    test_trim_leading_and_trailing_insertions("ACGTACGT", "2M2I4M", 0, 2, "2M2I1M", "ACGTA"); // no trailing insertion
    test_trim_leading_and_trailing_insertions("ACGTACGT", "2M2I4M", 0, 1, "2M", "AC");
    // trailing insertion removed
}

#[test]
fn test_bad_trim_loc() {
    let loc = SimpleInterval::new(0, 10, 20);
    let mut hap = Haplotype::new("ACGTAACCGGT".as_bytes(), false);
    hap.set_genome_location(loc);
    assert!(hap.trim(SimpleInterval::new(0, 1, 20)).is_err());
}

#[test]
fn test_bad_trim_no_loc() {
    let mut hap: Haplotype<SimpleInterval> = Haplotype::new("ACGTAACCGGT".as_bytes(), false);
    assert!(hap.trim(SimpleInterval::new(0, 1, 20)).is_err());
}
