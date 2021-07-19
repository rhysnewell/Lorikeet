use smith_waterman::bindings::*;
use bio_types::alignment::Alignment;
use rust_htslib::bam::record::{CigarString, Cigar};
use reads::alignment_utils::AlignmentUtils;

lazy_static! {
    static ref ORIGINAL_DEFAULT: SWParameters = SWParameters::new(3, -1, -4, -3);
    static ref STANDARD_NGS: SWParameters = SWParameters::new(25, -50, -110, -6);
}


pub struct SmithWatermanAligner {}

impl SmithWatermanAligner {

    pub fn new() -> SmithWatermanAligner {
        // let score = |a: u8, b: u8| if a == b { 3i32 } else { -1i32 };

        SmithWatermanAligner {}
    }

    /**
     *  perform a Smith-Waterman alignment of alt against ref
     *
     * @param ref bases to align to, values must be the byte equivalent of uppercase chars
     * @param alt bases to align against ref, values must be the byte equivalent of uppercase chars
     * @param parameters a set of weights to use when performing the alignment
     * @param overhangStrategy how to treat overhangs during alignment
     */
    pub fn align(
        reference: &[u8], alternate: &[u8], parameters: SWParameters, overhang_strategy: SWOverhangStrategy
    ) -> SmithWatermanAlignmentResult {
        assert!(reference.len() > 0 && alternate.len() > 0, "non-empty sequences are required for the Smith-Waterman calculation");

        let mut match_index = None;
        if overhang_strategy == SWOverhangStrategy::SoftClip
            || overhang_strategy == SWOverhangStrategy::Ignore {
            // Use a substring search to find an exact match of the alternate in the reference
            // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
            match_index = AlignmentUtils::last_index_of(reference, alternate);
        }

        let mut alignment_result;

        match match_index {
            None => {
                alignment_result = SmithWatermanAlignmentResult::new(CigarString::from(vec![Cigar::Match(alternate.len() as u32)]), match_index);
            },
            Some(match_index) => {
                // run full smith-waterman
                
            }
        }
    }
}

pub struct SmithWatermanAlignmentResult {
    pub(crate) cigar: CigarString,
    pub(crate) alignment_offset: Option<usize>
}

impl SmithWatermanAlignmentResult {
    pub fn new(cigar: CigarString, alignment_offset: Option<usize>) -> Self {
        Self {
            cigar,
            alignment_offset
        }
    }
}

enum State {
    Match,
    Insertion,
    Deletion,
    Clip
}