#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use gkl::smithwaterman::{OverhangStrategy, Parameters};
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::AVXMode;
use lorikeet_genome::reads::cigar_utils::CigarUtils;

use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    NEW_SW_PARAMETERS,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use std::convert::TryFrom;

lazy_static! {
    static ref HAPLOTYPE_TO_REFERENCE_SW_PARAMETERS: Parameters = *NEW_SW_PARAMETERS;
}

fn test_compute_cigar(s1: &str, s2: &str, expected_cigar: &str) {
    let actual_cigar = CigarUtils::calculate_cigar(
        s1.as_bytes(),
        s2.as_bytes(),
        OverhangStrategy::InDel,
        &*HAPLOTYPE_TO_REFERENCE_SW_PARAMETERS,
        AVXMode::detect_mode(),
    );
    let expected = CigarString::try_from(expected_cigar).unwrap();

    assert_eq!(actual_cigar.unwrap(), expected);
}

#[test]
fn make_test_compute_cigar_data() {
    test_compute_cigar("ATGGAGGGGC", "ATGGTGGGGC", "10M");
    test_compute_cigar("ATGGAGGGGC", "ATGGAAAATGGGGC", "5M4I5M");
    test_compute_cigar("ATGGAGGGGC", "ATGGAAAAAAAAATGGGGC", "5M9I5M");
    test_compute_cigar("ATGGAAAAAGGGGC", "ATGGTGGGGC", "4M4D6M");
    test_compute_cigar("ATGGAAAAAGGGGC", "ATGGAAAATGGGGC", "14M");
    test_compute_cigar("ATGGAAAAAGGGGC", "ATGGAAAAAAAAATGGGGC", "9M5I5M");
    test_compute_cigar("ATGGAAAAAAAAAAGGGGC", "ATGGTGGGGC", "4M9D6M");
    test_compute_cigar("ATGGAAAAAAAAAAGGGGC", "ATGGAAAATGGGGC", "4M5D10M");
    test_compute_cigar("ATGGAAAAAAAAAAGGGGC", "ATGGAAAAAAAAATGGGGC", "19M");
    test_compute_cigar(
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "51M",
    );
    test_compute_cigar(
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "NNNACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "3M6I48M",
    );
    test_compute_cigar(
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "51M",
    );
    test_compute_cigar(
        "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "NNNACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN",
        "3M6I48M",
    );
    test_compute_cigar("TCCCCCGGGT", "TCCCCCGGGT", "10M");
    test_compute_cigar("TCCCCCGGGT", "TAAACCCCCT", "1M3I5M3D1M");
    test_compute_cigar("G", "C", "1M");
    test_compute_cigar("G", "", "1D");
    test_compute_cigar("", "C", "1I");
    test_compute_cigar("AAA", "CGT", "3M");
    test_compute_cigar("TAT", "CAC", "3M");
    test_compute_cigar("GCTG", "GTCG", "4M");
    test_compute_cigar("AAAAA", "", "5D");
    test_compute_cigar("", "AAAAA", "5I");
    test_compute_cigar("AAAAACC", "CCGGGGGG", "5D2M6I");
    test_compute_cigar("GX", "CX", "2M");
    test_compute_cigar("GX", "X", "1D1M");
    test_compute_cigar("X", "CX", "1I1M");
    test_compute_cigar("AAAX", "CGTX", "4M");
    test_compute_cigar("TATX", "CACX", "4M");
    test_compute_cigar("GCTGX", "GTCGX", "5M");
    test_compute_cigar("AAAAAX", "X", "5D1M");
    test_compute_cigar("X", "AAAAAX", "5I1M");
    test_compute_cigar("AAAAACCX", "CCGGGGGGX", "5D2M6I1M");
    test_compute_cigar("GXXXXXXXXXXXXX", "CXXXXXXXXXXXXX", "14M");
    test_compute_cigar("GXXXXXXXXXXXXX", "XXXXXXXXXXXXX", "1D13M");
    test_compute_cigar("XXXXXXXXXXXXX", "CXXXXXXXXXXXXX", "1I13M");
    test_compute_cigar("AAAXXXXXXXXXXXXX", "CGTXXXXXXXXXXXXX", "16M");
    test_compute_cigar("TATXXXXXXXXXXXXX", "CACXXXXXXXXXXXXX", "16M");
    test_compute_cigar("GCTGXXXXXXXXXXXXX", "GTCGXXXXXXXXXXXXX", "17M");
    test_compute_cigar("AAAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXX", "5D13M");
    test_compute_cigar("XXXXXXXXXXXXX", "AAAAAXXXXXXXXXXXXX", "5I13M");
    test_compute_cigar("AAAAACCXXXXXXXXXXXXX", "CCGGGGGGXXXXXXXXXXXXX", "5D2M6I13M");
    test_compute_cigar("XG", "XC", "2M");
    test_compute_cigar("XG", "X", "1M1D");
    test_compute_cigar("X", "XC", "1M1I");
    test_compute_cigar("XAAA", "XCGT", "4M");
    test_compute_cigar("XTAT", "XCAC", "4M");
    test_compute_cigar("XGCTG", "XGTCG", "5M");
    test_compute_cigar("XAAAAA", "X", "1M5D");
    test_compute_cigar("X", "XAAAAA", "1M5I");
    test_compute_cigar("XAAAAACC", "XCCGGGGGG", "1M5D2M6I");
    test_compute_cigar("XGX", "XCX", "3M");
    test_compute_cigar("XGX", "XX", "1M1D1M");
    test_compute_cigar("XX", "XCX", "1M1I1M");
    test_compute_cigar("XAAAX", "XCGTX", "5M");
    test_compute_cigar("XTATX", "XCACX", "5M");
    test_compute_cigar("XGCTGX", "XGTCGX", "6M");
    test_compute_cigar("XAAAAAX", "XX", "1M5D1M");
    test_compute_cigar("XX", "XAAAAAX", "1M5I1M");
    test_compute_cigar("XAAAAACCX", "XCCGGGGGGX", "1M5D2M6I1M");
    test_compute_cigar("XGXXXXXXXXXXXXX", "XCXXXXXXXXXXXXX", "15M");
    test_compute_cigar("XGXXXXXXXXXXXXX", "XXXXXXXXXXXXXX", "1M1D13M");
    test_compute_cigar("XXXXXXXXXXXXXX", "XCXXXXXXXXXXXXX", "1M1I13M");
    test_compute_cigar("XAAAXXXXXXXXXXXXX", "XCGTXXXXXXXXXXXXX", "17M");
    test_compute_cigar("XTATXXXXXXXXXXXXX", "XCACXXXXXXXXXXXXX", "17M");
    test_compute_cigar("XGCTGXXXXXXXXXXXXX", "XGTCGXXXXXXXXXXXXX", "18M");
    test_compute_cigar("XAAAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXXX", "1M5D13M");
    test_compute_cigar("XXXXXXXXXXXXXX", "XAAAAAXXXXXXXXXXXXX", "1M5I13M");
    test_compute_cigar(
        "XAAAAACCXXXXXXXXXXXXX",
        "XCCGGGGGGXXXXXXXXXXXXX",
        "1M5D2M6I13M",
    );
    test_compute_cigar("XXXXXXXXXXXXXG", "XXXXXXXXXXXXXC", "14M");
    test_compute_cigar("XXXXXXXXXXXXXG", "XXXXXXXXXXXXX", "13M1D");
    test_compute_cigar("XXXXXXXXXXXXX", "XXXXXXXXXXXXXC", "13M1I");
    test_compute_cigar("XXXXXXXXXXXXXAAA", "XXXXXXXXXXXXXCGT", "16M");
    test_compute_cigar("XXXXXXXXXXXXXTAT", "XXXXXXXXXXXXXCAC", "16M");
    test_compute_cigar("XXXXXXXXXXXXXGCTG", "XXXXXXXXXXXXXGTCG", "17M");
    test_compute_cigar("XXXXXXXXXXXXXAAAAA", "XXXXXXXXXXXXX", "13M5D");
    test_compute_cigar("XXXXXXXXXXXXX", "XXXXXXXXXXXXXAAAAA", "13M5I");
    test_compute_cigar("XXXXXXXXXXXXXAAAAACC", "XXXXXXXXXXXXXCCGGGGGG", "13M5D2M6I");
    test_compute_cigar("XXXXXXXXXXXXXGX", "XXXXXXXXXXXXXCX", "15M");
    test_compute_cigar("XXXXXXXXXXXXXGX", "XXXXXXXXXXXXXX", "13M1D1M");
    test_compute_cigar("XXXXXXXXXXXXXX", "XXXXXXXXXXXXXCX", "13M1I1M");
    test_compute_cigar("XXXXXXXXXXXXXAAAX", "XXXXXXXXXXXXXCGTX", "17M");
    test_compute_cigar("XXXXXXXXXXXXXTATX", "XXXXXXXXXXXXXCACX", "17M");
    test_compute_cigar("XXXXXXXXXXXXXGCTGX", "XXXXXXXXXXXXXGTCGX", "18M");
    test_compute_cigar("XXXXXXXXXXXXXAAAAAX", "XXXXXXXXXXXXXX", "13M5D1M");
    test_compute_cigar("XXXXXXXXXXXXXX", "XXXXXXXXXXXXXAAAAAX", "13M5I1M");
    test_compute_cigar(
        "XXXXXXXXXXXXXAAAAACCX",
        "XXXXXXXXXXXXXCCGGGGGGX",
        "13M5D2M6I1M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXGXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXCXXXXXXXXXXXXX",
        "27M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXGXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXX",
        "13M1D13M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXCXXXXXXXXXXXXX",
        "13M1I13M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXAAAXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXCGTXXXXXXXXXXXXX",
        "29M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXTATXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXCACXXXXXXXXXXXXX",
        "29M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXGCTGXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXGTCGXXXXXXXXXXXXX",
        "30M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXAAAAAXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXX",
        "13M5D13M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXAAAAAXXXXXXXXXXXXX",
        "13M5I13M",
    );
    test_compute_cigar(
        "XXXXXXXXXXXXXAAAAACCXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXCCGGGGGGXXXXXXXXXXXXX",
        "13M5D2M6I13M",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGTTTCCTCCCCAAAAAAAAAAAATGGCCGCCCC",
        "32M4I",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGTTTCCTCCCCAAAAAAAAAAAATGGCCG",
        "32M",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGAAAATTTCCTCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC",
        "5M4I22M4I5M4I",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGAAAATTTCCTCCCCAAAAAAAAAAAAGGGGTGGCCG",
        "5M4I22M4I5M",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC",
        "5M9I22M9I5M4I",
    );
    test_compute_cigar(
        "ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG",
        "5M9I22M9I5M",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGTTTCCTCCCCCCCCAAAAAAAAAAAATGGCCGCCCC",
        "4M4D26M4D6M4I",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGTTTCCTCCCCCCCCAAAAAAAAAAAATGGCCG",
        "4M4D26M4D6M",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC",
        "44M4I",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGTGGCCG",
        "44M",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC",
        "9M5I30M5I5M4I",
    );
    test_compute_cigar(
        "ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG",
        "9M5I30M5I5M",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGTTTCCTCCCCCCCCCCCCCAAAAAAAAAAAATGGCCGCCCC",
        "4M9D31M9D6M4I",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGTTTCCTCCCCCCCCCCCCCAAAAAAAAAAAATGGCCG",
        "4M9D31M9D6M",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC",
        "4M5D35M5D10M4I",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGTGGCCG",
        "4M5D35M5D10M",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC",
        "59M4I",
    );
    test_compute_cigar(
        "ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG",
        "ATGGAAAAAAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG",
        "59M",
    );
}

fn test_clip_cigar(
    original: &str,
    start: u32,
    stop: u32,
    expected_soft_clip: &str,
    expected_hard_clip: &str,
) {
    let original_cigar = CigarString::try_from(original).unwrap().into_view(0);
    assert_eq!(
        CigarUtils::clip_cigar(&original_cigar, start, stop, Cigar::SoftClip(0)).to_string(),
        expected_soft_clip.to_string()
    );
    assert_eq!(
        CigarUtils::clip_cigar(&original_cigar, start, stop, Cigar::HardClip(0)).to_string(),
        expected_hard_clip.to_string()
    );
}

#[test]
fn clip_cigar() {
    //simple cases
    test_clip_cigar("10M", 0, 5, "5S5M", "5H5M");
    test_clip_cigar("10M", 5, 10, "5M5S", "5M5H");
    test_clip_cigar("10H10M", 0, 5, "10H5S5M", "15H5M");
    test_clip_cigar("10H10M", 5, 10, "10H5M5S", "10H5M5H");
    test_clip_cigar("10M10H", 0, 5, "5S5M10H", "5H5M10H");

    // clipping into insertion
    test_clip_cigar("10M10I10M", 0, 5, "5S5M10I10M", "5H5M10I10M");
    test_clip_cigar("10M10I10M", 0, 15, "15S5I10M", "15H5I10M");
    test_clip_cigar("10M10I10M", 15, 30, "10M5I15S", "10M5I15H");

    // clipping into a soft clip
    test_clip_cigar("10S10M10S", 0, 5, "10S10M10S", "5H5S10M10S");
    test_clip_cigar("10S10M10S", 25, 30, "10S10M10S", "10S10M5S5H");
    test_clip_cigar("10S10M10S", 0, 15, "15S5M10S", "15H5M10S");

    // clipping over a deletion
    test_clip_cigar("10M10D10M", 0, 10, "10S10M", "10H10M");
    test_clip_cigar("10M10D10M", 0, 15, "15S5M", "15H5M");
    test_clip_cigar("10M10D10M", 5, 20, "5M15S", "5M15H");

    // removing leading deletions
    test_clip_cigar("10D10M", 0, 5, "5S5M", "5H5M");
}

fn test_alignment_start_shift(cigar_string: &str, num_clips: i64, expected_result: i64) {
    let cigar = CigarString::try_from(cigar_string).unwrap().into_view(0);
    let actual_result = CigarUtils::alignment_start_shift(&cigar, num_clips);
    assert_eq!(actual_result, expected_result);
}

#[test]
fn alignment_start_shift() {
    test_alignment_start_shift("70M", 10, 10);
    test_alignment_start_shift("70M", 0, 0);

    test_alignment_start_shift("30M10D30M", 29, 29);
    test_alignment_start_shift("30M10D30M", 30, 40);
    test_alignment_start_shift("30M10D30M", 31, 41);

    test_alignment_start_shift("30M10D30M", 29, 29);
    test_alignment_start_shift("30M10I30M", 30, 30);
    test_alignment_start_shift("30M10I30M", 31, 30);
    test_alignment_start_shift("30M10I30M", 40, 30);
    test_alignment_start_shift("30M10I30M", 41, 31);

    test_alignment_start_shift("10H10M", 5, 5);
    test_alignment_start_shift("10S10M", 5, 0);
    test_alignment_start_shift("10S10M", 5, 0);
}
