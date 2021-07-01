use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE, MatchFunc};
use test::TestName::AlignedTestName;
use bio_types::alignment::Alignment;
use rust_htslib::bam::record::CigarString;

lazy_static! {
    static ORIGINAL_DEFAULT: Scoring = Scoring::new(-4, -3, 3, -1).xclip(MIN_SCORE).yclip(MIN_SCORE);
    static STANDARD_NGS: Scoring = Scoring::new(-110, -6, 25, -50).xclip(MIN_SCORE).yclip(MIN_SCORE);
}

pub struct SmithWatermanAligner {
    aligner: Aligner<MatchFunc>
}

impl SmithWatermanAligner {

    pub fn new() -> SmithWatermanAligner {
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
    pub fn align<F: MatchFunc>(ref_seq: &[u8], alt_seq: &[u8], parameters: Scoring<F>) -> Alignment {
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
    pub fn new(alignment: Alignment)
}