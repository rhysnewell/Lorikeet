/**
 * How overhangs should be treated during Smith-Waterman alignment
 */
pub enum SWOverhangStrategy {
    /*
     * Add softclips for the overhangs
     */
    SoftClip,

    /*
     * Treat the overhangs as proper insertions/deletions
     */
    InDel,

    /*
     * Treat the overhangs as proper insertions/deletions for leading (but not trailing) overhangs.
     * This is useful e.g. when we want to merge dangling tails in an assembly graph: because we don't
     * expect the dangling tail to reach the end of the reference path we are okay ignoring trailing
     * deletions - but leading indels are still very much relevant.
     */
    LeadingInDel,

    /*
     * Just ignore the overhangs
     */
    Ignore
}

/**
 * a set of parameters to configure Smith-Waterman assembly
 */
pub struct SWParameters {
    pub(crate) match_value: i32,
    pub(crate) mismatch_penalty: i32,
    pub(crate) gap_open_penalty: i32,
    pub(crate) gap_extend_penalty: i32,
}

impl SWParameters {
    /**
     * Create a new set of parameters for Smith-Waterman alignment
     * @param matchValue how much to reward a match during alignment >= 0
     * @param mismatchPenalty how much to penalize a mismatch during alignment <= 0
     * @param gapOpenPenalty how much penalize the creation of a new gap in the alignment <= 0
     * @param gapExtendPenalty how much to penalize extending an already open gap in the alignment <= 0
     */
    pub fn new(
        match_value: i32, mismatch_penalty: i32,
        gap_open_penalty: i32, gap_extend_penalty: i32,
    ) -> Self {
        assert!(match_value >= 0, "matchValue must be >= 0");
        assert!(mismatch_penalty <= 0, "mismatchPenalty must be <= 0");
        assert!(gap_open_penalty <= 0, "gapOpenPenalty must be <= 0");
        assert!(gap_extend_penalty <= 0, "gapExtendPenalty must be <= 0");
        Self {
            match_value,
            mismatch_penalty,
            gap_open_penalty,
            gap_extend_penalty,
        }
    }

}