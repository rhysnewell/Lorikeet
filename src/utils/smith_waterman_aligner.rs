use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE, MatchParams};
use bio_types::alignment::Alignment;
use rust_htslib::bam::record::CigarString;

lazy_static! {
    static ref ORIGINAL_DEFAULT: Scoring<MatchParams> = Scoring::from_scores(-4, -3, 3, -1).xclip(MIN_SCORE).yclip(MIN_SCORE);
    static ref STANDARD_NGS: Scoring<MatchParams> = Scoring::from_scores(-110, -6, 25, -50).xclip(MIN_SCORE).yclip(MIN_SCORE);
}


pub struct SmithWatermanAligner {
    aligner: Aligner<MatchParams>
}

impl SmithWatermanAligner {

    pub fn new() -> SmithWatermanAligner {
        // let score = |a: u8, b: u8| if a == b { 3i32 } else { -1i32 };

        SmithWatermanAligner {
            aligner: Aligner::with_scoring(*ORIGINAL_DEFAULT)
        }
    }

    /**
     *  perform a Smith-Waterman alignment of alt against ref
     *
     * @param ref bases to align to, values must be the byte equivalent of uppercase chars
     * @param alt bases to align against ref, values must be the byte equivalent of uppercase chars
     * @param parameters a set of weights to use when performing the alignment
     * @param overhangStrategy how to treat overhangs during alignment
     */
    pub fn align(ref_seq: &[u8], alt_seq: &[u8], parameters: Scoring<MatchParams>) -> Alignment {
        // avoid running full Smith-Waterman if there is an exact match of alternate in reference
        let mut match_index = -1;

        let mut aligner = Aligner::with_scoring(parameters);
        let alignment = aligner.custom(ref_seq, alt_seq);

        return alignment
    }
}

pub struct SWNativeResultWrapper {
    cigar: CigarString,
    alignment_offset: i64
}

impl SWNativeResultWrapper {
    // pub fn new(alignment: Alignment)
}