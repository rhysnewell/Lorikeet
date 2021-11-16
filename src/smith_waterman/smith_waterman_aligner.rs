use ndarray::{Array2, Axis};
use reads::alignment_utils::AlignmentUtils;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
// use smith_waterman::bindings::*;
use std::cmp::max;
use gkl::smithwaterman::{OverhangStrategy, Parameters};

lazy_static! {
    pub static ref ORIGINAL_DEFAULT: Parameters = Parameters::new(3, -1, -4, -3);
    pub static ref STANDARD_NGS: Parameters = Parameters::new(25, -50, -110, -6);
        // FROM GATK COMMENTS:
    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    pub static ref NEW_SW_PARAMETERS: Parameters = Parameters::new(200, -150, -260, -11);
    // FROM GATK COMMENTS:
    // In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
    // The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
    // Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
    // would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
    // Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
    // we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
    pub static ref ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS: Parameters = Parameters::new(10, -15, -30, -5);
}

pub struct SmithWatermanAligner {}

impl SmithWatermanAligner {
    const MATRIX_MIN_CUTOFF: i32 = -100000000;

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
        reference: &[u8],
        alternate: &[u8],
        parameters: &Parameters,
        overhang_strategy: OverhangStrategy,
    ) -> SmithWatermanAlignmentResult {
        assert!(
            reference.len() > 0 && alternate.len() > 0,
            "non-empty sequences are required for the Smith-Waterman calculation"
        );

        let mut match_index = None;
        match &overhang_strategy {
            OverhangStrategy::SoftClip
            | OverhangStrategy::Ignore => {
                // Use a substring search to find an exact match of the alternate in the reference
                // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
                match_index = AlignmentUtils::last_index_of(reference, alternate);
            },
            _ => {
                // pass
            }
        }

        let mut alignment_result;

        match match_index {
            None => {
                // run full smith-waterman
                let n = reference.len() + 1;
                let m = alternate.len() + 1;
                let mut sw = Array2::<i32>::zeros((n, m));
                let mut btrack = Array2::<i32>::zeros((n, m));

                Self::calculate_matrix(
                    reference,
                    alternate,
                    &mut sw,
                    &mut btrack,
                    &overhang_strategy,
                    &parameters,
                );
                alignment_result = Self::calculate_cigar(&sw, &btrack, &overhang_strategy);
            }
            Some(match_index) => {
                alignment_result = SmithWatermanAlignmentResult::new(
                    CigarString::from(vec![Cigar::Match(alternate.len() as u32)]),
                    match_index as i32,
                )
            }
        }

        return alignment_result;
    }

    /**
     * Calculates the SW matrices for the given sequences
     * @param reference  ref sequence
     * @param alternate  alt sequence
     * @param sw         the Smith-Waterman matrix to populate
     * @param btrack     the back track matrix to populate
     * @param overhangStrategy    the strategy to use for dealing with overhangs
     * @param parameters the set of weights to use to configure the alignment
     */
    fn calculate_matrix(
        reference: &[u8],
        alternate: &[u8],
        sw: &mut Array2<i32>,
        btrack: &mut Array2<i32>,
        overhang_strategy: &OverhangStrategy,
        parameters: &Parameters,
    ) {
        if reference.len() == 0 || alternate.len() == 0 {
            panic!("Non-null, non-empty sequences are required for the Smith-Waterman calculation");
        };

        let ncol = sw.ncols();
        let nrow = sw.nrows();

        let low_init_value = std::i32::MIN / 2;
        let mut best_gap_v = vec![low_init_value; ncol + 1];
        let mut gap_size_v = vec![0; ncol + 1];
        let mut best_gap_h = vec![low_init_value; nrow + 1];
        let mut gap_size_h = vec![0; nrow + 1];

        // we need to initialize the SW matrix with gap penalties if we want to keep track of indels at the edges of alignments
        match overhang_strategy {
            OverhangStrategy::InDel
            | OverhangStrategy::LeadingInDel => {
                let current_value = parameters.gap_open_penalty;

                // initialize the first row
                {
                    let mut top_row = sw.row_mut(0);
                    // .accumulate_axis_inplace(Axis(1), |&prev, curr| *curr = prev + current_value);
                    let mut current_value = parameters.gap_open_penalty;
                    top_row[1] = current_value;
                    for i in 2..top_row.len() {
                        current_value += parameters.gap_extend_penalty;
                        top_row[i] = current_value;
                    }
                }
                debug!("Top row values {:?}", &sw.row(0));

                // initialize the first column
                {
                    sw[[1, 0]] = parameters.gap_open_penalty;
                    let mut current_value = parameters.gap_open_penalty;
                    for i in 2..sw.nrows() {
                        current_value += parameters.gap_extend_penalty;
                        sw[[i, 0]] = current_value;
                    }
                }
                // sw.column_mut(0)
                //     .accumulate_axis_inplace(Axis(0), |&prev, curr| *curr = prev + current_value);
                debug!("first column values {:?}", &sw.column(0));
            },
            _ => {
                // pass
            }
        }

        // build smith-waterman matrix and keep backtrack info:

        //access is pricey if done enough times so we extract those out
        let w_open = parameters.gap_open_penalty;
        let w_extend = parameters.gap_extend_penalty;
        let w_match = parameters.match_value;
        let w_mismatch = parameters.mismatch_penalty;

        for i in 1..nrow {
            let a_base = reference[i - 1]; // letter in a at the current pos
                                           // let mut cur_row = sw.row_mut(i);
            let mut cur_back_track_row = btrack.row_mut(i);

            for j in 1..ncol {
                let b_base = alternate[j - 1]; // letter in b at the current pos
                                               // in other words, step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                let step_diag = sw[[i - 1, j - 1]]
                    + if a_base == b_base {
                        w_match
                    } else {
                        w_mismatch
                    };

                // optimized "traversal" of all the matrix cells above the current one (i.e. traversing
                // all 'step down' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                // if a gap (length 1) was just opened above, this is the cost of arriving to the current cell:
                let mut prev_gap = sw[[i - 1, j]] + w_open;
                best_gap_v[j] += w_extend; // for the gaps that were already opened earlier, extending them by 1 costs w_extend
                if prev_gap > best_gap_v[j] {
                    // opening a gap just before the current cell results in better score than extending by one
                    // the best previously opened gap. This will hold for ALL cells below: since any gap
                    // once opened always costs w_extend to extend by another base, we will always get a better score
                    // by arriving to any cell below from the gap we just opened (prev_gap) rather than from the previous best gap
                    best_gap_v[j] = prev_gap;
                    gap_size_v[j] = 1; // remember that the best step-down gap from above has length 1 (we just opened it)
                } else {
                    // previous best gap is still the best, even after extension by another base, so we just record that extension:
                    gap_size_v[j] += 1;
                }

                let step_down = best_gap_v[j];
                let kd = gap_size_v[j];

                // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
                // all 'step right' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                prev_gap = sw[[i, j - 1]] + w_open;
                best_gap_h[i] += w_extend;

                if prev_gap > best_gap_h[i] {
                    // newly opened gap is better (score-wise) than any previous gap with the same row index i; since
                    // gap penalty is linear with k, this new gap location is going to remain better than any previous ones
                    best_gap_h[i] = prev_gap;
                    gap_size_h[i] = 1;
                } else {
                    gap_size_h[i] += 1;
                }

                let step_right = best_gap_h[i];
                let ki = gap_size_h[i];

                //priority here will be step diagonal, step right, step down
                let diag_highest_or_equal = step_diag >= step_down && step_diag >= step_right;

                if diag_highest_or_equal {
                    sw[[i, j]] = max(Self::MATRIX_MIN_CUTOFF, step_diag);
                    cur_back_track_row[j] = 0;
                } else if step_right >= step_down {
                    //moving right is the highest
                    sw[[i, j]] = max(Self::MATRIX_MIN_CUTOFF, step_right);
                    cur_back_track_row[j] = -ki; // negative = horizontal
                } else {
                    sw[[i, j]] = max(Self::MATRIX_MIN_CUTOFF, step_down);
                    cur_back_track_row[j] = kd; // positive = vertical
                }
            }
        }
    }

    /**
     * Calculates the CIGAR for the alignment from the back track matrix
     *
     * @param sw                   the Smith-Waterman matrix to use
     * @param btrack               the back track matrix to use
     * @param overhangStrategy    the strategy to use for dealing with overhangs
     * @return non-null SWPairwiseAlignmentResult object
     */
    fn calculate_cigar(
        sw: &Array2<i32>,
        btrack: &Array2<i32>,
        overhang_strategy: &OverhangStrategy,
    ) -> SmithWatermanAlignmentResult {
        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order
        let mut p1 = 0;
        let mut p2 = 0;

        let ref_length = sw.nrows() - 1;
        let alt_length = sw.ncols() - 1;

        let mut max_score = std::i32::MIN; // sw scores are allowed to be negative
        let mut segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        // if we want to consider overhangs as legitimate operators, then just start from the corner of the matrix
        match overhang_strategy {
            OverhangStrategy::InDel => {
                p1 = ref_length;
                p2 = alt_length;
            },
            _ => {
                // look for the largest score on the rightmost column. we use >= combined with the traversal direction
                // to ensure that if two scores are equal, the one closer to diagonal gets picked
                //Note: this is not technically smith-waterman, as by only looking for max values on the right we are
                //excluding high scoring local alignments
                p2 = alt_length;

                for i in 1..sw.nrows() {
                    let cur_score = sw[[i, alt_length]];
                    if cur_score >= max_score {
                        p1 = i;
                        max_score = cur_score;
                    };
                }

                // now look for a larger score on the bottom-most row
                match overhang_strategy {
                    OverhangStrategy::LeadingInDel => {
                        // pass
                    },
                    _ => {
                        let bottom_row = sw.row(ref_length);
                        for j in 1..bottom_row.len() {
                            let cur_score = bottom_row[j];
                            // data_offset is the offset of [n][j]
                            if cur_score > max_score
                                || (cur_score == max_score
                                && (ref_length as i32 - j as i32).abs() < (p1 as i32 - p2 as i32).abs())
                            {
                                p1 = ref_length;
                                p2 = j;
                                max_score = cur_score;
                                segment_length = alt_length as i32 - j as i32; // end of sequence 2 is overhanging; we will just record it as 'M' segment
                            }
                        }
                    }
                }
            }
        }

        let mut lce = Vec::new();
        if segment_length > 0 {
            match overhang_strategy {
                OverhangStrategy::SoftClip => {
                    lce.push(Cigar::SoftClip(segment_length as u32));
                    segment_length = 0;
                },
                _ => {
                    // pass
                }
            }
        };

        // we will be placing all insertions and deletions into sequence b, so the states are named w/regard
        // to that sequence
        let mut p1 = p1 as i32; // type conversion since p1 and p2 will have to go negative
        let mut p2 = p2 as i32;
        let mut state = State::Match;
        loop {
            let btr = btrack[[p1 as usize, p2 as usize]];
            let mut new_state;
            let mut step_length = 1;
            if btr > 0 {
                new_state = State::Deletion;
                new_state = State::Deletion;
                step_length = btr;
            } else if btr < 0 {
                new_state = State::Insertion;
                step_length = -btr;
            } else {
                new_state = State::Match; // and step_length =1, already set above
            }

            // move to next best location in the sw matrix:
            match new_state {
                State::Match => {
                    // move back along the diag in the sw matrix
                    p1 -= 1;
                    p2 -= 1;
                }
                State::Insertion => {
                    p2 -= step_length;
                }
                State::Deletion => {
                    p1 -= step_length;
                }
                _ => {
                    // do nothing
                }
            };

            // now let's see if the state actually changed:
            if new_state == state {
                segment_length += step_length;
            } else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                if segment_length > 0 {
                    lce.push(Self::make_element(state, segment_length as u32))
                }

                segment_length = step_length;
                state = new_state;
            }
            // next condition is equivalent to  while ( sw[p1][p2] != 0 ) (with modified p1 and/or p2:
            if p1 <= 0 || p2 <= 0 {
                break;
            }
        }

        // post-process the last segment we are still keeping;
        // NOTE: if reads "overhangs" the ref on the left (i.e. if p2>0) we are counting
        // those extra bases sticking out of the ref into the first cigar element if DO_SOFTCLIP is false;
        // otherwise they will be soft-clipped. For instance,
        // if read length is 5 and alignment starts at offset -2 (i.e. read starts before the ref, and only
        // last 3 bases of the read overlap with/align to the ref), the cigar will be still 5M if
        // DO_SOFTCLIP is false or 2S3M if DO_SOFTCLIP is true.
        // The consumers need to check for the alignment offset and deal with it properly.
        let mut alignment_offset = 0;
        match overhang_strategy {
            OverhangStrategy::SoftClip => {
                lce.push(Self::make_element(state, segment_length as u32));
                if p2 > 0 {
                    lce.push(Self::make_element(State::Clip, p2 as u32));
                }
                alignment_offset = p1;
            },
            OverhangStrategy::Ignore => {
                lce.push(Self::make_element(state, (segment_length + p2) as u32));
                alignment_offset = p1 - p2;
            },
            _ => {
                // overhangStrategy == OverhangStrategy.INDEL || overhangStrategy == OverhangStrategy.LEADING_INDEL
                // take care of the actual alignment
                lce.push(Self::make_element(state, segment_length as u32));
                // take care of overhangs at the beginning of the alignment
                if p1 > 0 {
                    lce.push(Self::make_element(State::Deletion, p1 as u32));
                } else if p2 > 0 {
                    lce.push(Self::make_element(State::Insertion, p2 as u32));
                }

                alignment_offset = 0;
            }
        }

        // reverse the cigar string
        lce.reverse();
        return SmithWatermanAlignmentResult::new(CigarString::from(lce), alignment_offset);
    }

    fn make_element(state: State, length: u32) -> Cigar {
        return match state {
            State::Match => Cigar::Match(length),
            State::Insertion => Cigar::Ins(length),
            State::Deletion => Cigar::Del(length),
            State::Clip => Cigar::SoftClip(length),
        };
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmithWatermanAlignmentResult {
    pub(crate) cigar: CigarString,
    pub(crate) alignment_offset: i32,
}

impl SmithWatermanAlignmentResult {
    pub fn new(cigar: CigarString, alignment_offset: i32) -> Self {
        Self {
            cigar,
            alignment_offset,
        }
    }

    pub fn get_alignment_offset(&self) -> i32 {
        self.alignment_offset
    }

    pub fn get_cigar(&self) -> CigarString {
        self.cigar.clone()
    }
}

#[derive(Debug, Eq, PartialEq)]
enum State {
    Match,
    Insertion,
    Deletion,
    Clip,
}
