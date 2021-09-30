use ndarray::parallel::prelude::*;
use ndarray::{Array2, ArrayViewMut1, Axis};
use num::traits::real::Real;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use utils::math_utils::MathUtils;
use utils::quality_utils::QualityUtils;

lazy_static! {
    /**
      * Convenient ln10 constant.
      */
    static ref LN10: f64 = 10.0_f64.ln();

    /**
      * Convenient (ln10)^-1 constant.
      */
    static ref INV_LN10: f64 = 1.0 / *LN10;
}

/**
 * Helper class that implement calculations required to implement the PairHMM Finite State Automation (FSA) model.
 */
pub struct PairHMMModel {
    /**
     * Holds pre-calculated the matchToMath probability values in linear scale.
     *
     * <p/>
     * This is a triangular matrix stored in a unidimentional array like so:
     * <p/>
     * (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), (0,3) ... ({@link QualityUtils#MAX_QUAL},{@link QualityUtils#MAX_QUAL})
     */
    match_to_match_prob: Vec<f64>,
    /**
     * Holds pre-calculated the matchToMatch log10-probabilities.
     *
     * <p/>
     * This is a triangular matrix stored in a unidimentional array like so:
     * <p/>
     * (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), (0,3) ... ({@link QualityUtils#MAX_QUAL},{@link QualityUtils#MAX_QUAL})
     */
    match_to_match_log10: Vec<f64>,
}

impl PairHMMModel {
    pub fn new() -> PairHMMModel {
        let mut match_to_match_prob = vec![
            0.0;
            ((QualityUtils::MAX_QUAL as usize + 1)
                * (QualityUtils::MAX_QUAL as usize + 2))
                >> 1
        ];
        let mut match_to_match_log10 = vec![
            0.0;
            ((QualityUtils::MAX_QUAL as usize + 1)
                * (QualityUtils::MAX_QUAL as usize + 2))
                >> 1
        ];
        let mut offset = 0;
        for i in 0..=(QualityUtils::MAX_QUAL as usize) {
            for j in 0..=i {
                let log10_sum =
                    MathUtils::approximate_log10_sum_log10(-0.1 * (i as f64), -0.1 * (j as f64));
                let log10_sum_pow: f64 = 10.0.powf(log10_sum);
                match_to_match_log10[offset + j] =
                    (-std::cmp::min(OrderedFloat(1.0), OrderedFloat(log10_sum_pow))).ln_1p()
                        * *INV_LN10;
                match_to_match_prob[offset + j] = 10.0.powf(match_to_match_log10[offset + j]);
            }
            offset += i + 1;
        }

        PairHMMModel {
            match_to_match_log10,
            match_to_match_prob,
        }
    }

    /**
     * Length of the standard transition probability array.
     */
    pub const TRANS_PROB_ARRAY_LENGTH: usize = 6;

    /**
     * Position in the transition probability array for the Match-to-Match transition.
     */
    pub const match_to_match: usize = 0;

    /**
     * Position in the transition probability array for the Indel-to-Match transition.
     */
    pub const indel_to_match: usize = 1;

    /**
     * Position in the transition probability array for the Match-to-Insertion transition.
     */
    pub const match_to_insertion: usize = 2;

    /**
     * Position in the transition probability array for the Insertion-to-Insertion transition.
     */
    pub const insertion_to_insertion: usize = 3;

    /**
     * Position in the transition probability array for the Match-to-Deletion transition.
     */
    pub const match_to_deletion: usize = 4;

    /**
     * Position in the transition probability array for the Deletion-to-Deletion transition.
     */
    pub const deletion_to_deletion: usize = 5;

    /**
     * Fills a transition probability array given the different quality scores affecting a read site
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     */
    pub fn qual_to_trans_probs_with_vec(
        &self,
        dest: &mut Vec<f64>,
        ins_qual: u8,
        del_qual: u8,
        gcp: u8,
    ) {
        dest[Self::match_to_match] = self.match_to_match_prob(ins_qual as usize, del_qual as usize);
        dest[Self::match_to_insertion] = QualityUtils::qual_to_error_prob(ins_qual);
        dest[Self::match_to_deletion] = QualityUtils::qual_to_error_prob(del_qual);
        dest[Self::indel_to_match] = QualityUtils::qual_to_prob(gcp);
        let tmp = QualityUtils::qual_to_error_prob(gcp);
        dest[Self::insertion_to_insertion] = tmp;
        dest[Self::deletion_to_deletion] = tmp;
    }

    pub fn qual_to_trans_probs_with_array1(
        &self,
        dest: &mut ArrayViewMut1<f64>,
        ins_qual: u8,
        del_qual: u8,
        gcp: u8,
    ) {
        dest[Self::match_to_match] = self.match_to_match_prob(ins_qual as usize, del_qual as usize);
        dest[Self::match_to_insertion] = QualityUtils::qual_to_error_prob(ins_qual);
        dest[Self::match_to_deletion] = QualityUtils::qual_to_error_prob(del_qual);
        dest[Self::indel_to_match] = QualityUtils::qual_to_prob(gcp);
        let tmp = QualityUtils::qual_to_error_prob(gcp);
        dest[Self::insertion_to_insertion] = tmp;
        dest[Self::deletion_to_deletion] = tmp;
    }

    /**
     * Returns a transition probability array given the different quality scores affecting a read site.
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     *
     * @return never {@code null}. An array of length {@link #TRANS_PROB_ARRAY_LENGTH}.
     */
    pub fn qual_to_trans_probs_return_vec(&self, ins_qual: u8, del_qual: u8, gcp: u8) -> Vec<f64> {
        let mut dest = vec![0.0; Self::TRANS_PROB_ARRAY_LENGTH];
        self.qual_to_trans_probs_with_vec(&mut dest, ins_qual, del_qual, gcp);
        dest
    }

    /**
     * Fills ax matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param dest the matrix to update
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays, matrices is {@code null} or any entry in {@code dest} is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} or any of its elements is not large enough to contain the
     *  transition  matrix.
     */
    pub fn qual_to_trans_probs_with_array(
        &self,
        dest: &mut Array2<f64>,
        ins_quals: &[u8],
        del_quals: &[u8],
        gcps: &[u8],
    ) {
        let read_length = ins_quals.len();
        if del_quals.len() != read_length {
            panic!("deletion quality array length does not match insertion quality array length");
        };

        if gcps.len() != read_length {
            panic!("deletion quality array length does not match insertion quality array length");
        };

        if dest.shape()[0] < read_length + 1 {
            panic!("destination is not enough for the read length");
        };

        let model = &self;
        dest.axis_iter_mut(Axis(0))
            .into_par_iter()
            .skip(1)
            .enumerate()
            .for_each(|(i, mut row)| {
                if i < read_length {
                    model.qual_to_trans_probs_with_array1(
                        &mut row,
                        ins_quals[i],
                        del_quals[i],
                        gcps[i],
                    );
                }
            });
    }

    /**
     * Returns a matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     *
     * @return never {@code null}, an matrix of the dimensions explained above.
     */
    pub fn qual_to_trans_probs_return_array(
        &self,
        ins_quals: &[u8],
        del_quals: &[u8],
        gcps: &[u8],
    ) -> Array2<f64> {
        let mut dest = Self::create_transition_matrix(ins_quals.len());
        self.qual_to_trans_probs_with_array(&mut dest, ins_quals, del_quals, gcps);
        dest
    }

    /**
     * Creates a transition probability matrix large enough to work with sequences of a particular length.
     *
     * @param maxReadLength the maximum read length for the transition matrix.
     *
     * @return never {@code null}. A matrix of {@code maxReadLength + 1} by {@link #TRANS_PROB_ARRAY_LENGTH} positions.
     */
    pub fn create_transition_matrix(max_read_length: usize) -> Array2<f64> {
        Array2::zeros((max_read_length + 1, Self::TRANS_PROB_ARRAY_LENGTH))
    }

    /**
     * Fills a transition probability array given the different quality scores affecting a read site
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     */
    pub fn qual_to_trans_probs_log10_with_vec(
        &self,
        dest: &mut [f64],
        ins_qual: u8,
        del_qual: u8,
        gcp: u8,
    ) {
        dest[Self::match_to_match] =
            self.match_to_match_prob_log10(ins_qual as usize, del_qual as usize);
        dest[Self::match_to_insertion] = QualityUtils::qual_to_error_prob_log10(ins_qual);
        dest[Self::match_to_deletion] = QualityUtils::qual_to_error_prob_log10(del_qual);
        dest[Self::indel_to_match] = QualityUtils::qual_to_prob_log10(gcp);
        let tmp = QualityUtils::qual_to_error_prob_log10(gcp);
        dest[Self::insertion_to_insertion] = tmp;
        dest[Self::deletion_to_deletion] = tmp;
    }

    pub fn qual_to_trans_probs_log10_with_array1(
        &self,
        dest: &mut ArrayViewMut1<f64>,
        ins_qual: u8,
        del_qual: u8,
        gcp: u8,
    ) {
        dest[Self::match_to_match] =
            self.match_to_match_prob_log10(ins_qual as usize, del_qual as usize);
        dest[Self::match_to_insertion] = QualityUtils::qual_to_error_prob_log10(ins_qual);
        dest[Self::match_to_deletion] = QualityUtils::qual_to_error_prob_log10(del_qual);
        dest[Self::indel_to_match] = QualityUtils::qual_to_prob_log10(gcp);
        let tmp = QualityUtils::qual_to_error_prob_log10(gcp);
        dest[Self::insertion_to_insertion] = tmp;
        dest[Self::deletion_to_deletion] = tmp;
    }

    /**
     * Returns a transition probability array given the different quality scores affecting a read site.
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     *
     * @return never {@code null}. An array of length {@link #TRANS_PROB_ARRAY_LENGTH}.
     */
    pub fn qual_to_trans_probs_log10_return_vec(
        &self,
        ins_qual: u8,
        del_qual: u8,
        gcp: u8,
    ) -> Vec<f64> {
        let mut dest = vec![0.0; Self::TRANS_PROB_ARRAY_LENGTH];
        self.qual_to_trans_probs_log10_with_vec(&mut dest, ins_qual, del_qual, gcp);
        dest
    }

    /**
     * Fills ax matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param dest the matrix to update
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays, matrices is {@code null} or any entry in {@code dest} is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} or any of its elements is not large enough to contain the
     *  transition  matrix.
     */
    pub fn qual_to_trans_probs_log10_with_array(
        &self,
        dest: &mut Array2<f64>,
        ins_quals: Vec<u8>,
        del_quals: Vec<u8>,
        gcps: Vec<u8>,
    ) {
        let read_length = ins_quals.len();
        if del_quals.len() != read_length {
            panic!("deletion quality array length does not match insertion quality array length");
        };

        if gcps.len() != read_length {
            panic!("deletion quality array length does not match insertion quality array length");
        };

        if dest.shape()[0] < read_length + 1 {
            panic!("destination is not enough for the read length");
        };

        let model = &self;
        dest.axis_iter_mut(Axis(0))
            .into_par_iter()
            .enumerate()
            .for_each(|(i, mut row)| {
                model.qual_to_trans_probs_log10_with_array1(
                    &mut row,
                    ins_quals[i],
                    del_quals[i],
                    gcps[i],
                );
            });
    }

    /**
     * Returns a matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     *
     * @return never {@code null}, an matrix of the dimensions explained above.
     */
    pub fn qual_to_trans_probs_log10_return_array(
        &self,
        ins_quals: Vec<u8>,
        del_quals: Vec<u8>,
        gcps: Vec<u8>,
    ) -> Array2<f64> {
        let mut dest = Self::create_transition_matrix(ins_quals.len());
        self.qual_to_trans_probs_log10_with_array(&mut dest, ins_quals, del_quals, gcps);
        dest
    }

    /**
     * Returns the probability that neither of two events, insertion and deletion, takes place.
     * <p/>
     *
     * We assume that both event never occur together and that delQual is the conditional probability
     * (qual. encoded) of the second event, given the first event didn't took place. So that the
     * probability of no event is: <br/>
     *
     * <code>1 - ProbErr(insQual) - ProbErr(delQual)</code> <br/>
     *
     * @param insQual PhRED scaled quality/probability of an insertion.
     * @param delQual PhRED scaled quality/probability of a deletion.
     * @return a value between 0 and 1.
     */
    pub fn match_to_match_prob(&self, ins_qual: usize, del_qual: usize) -> f64 {
        let mut min_qual;
        let mut max_qual;
        if ins_qual <= del_qual {
            min_qual = ins_qual;
            max_qual = del_qual;
        } else {
            min_qual = del_qual;
            max_qual = ins_qual;
        }

        if (QualityUtils::MAX_QUAL as usize) < max_qual {
            1.0 - 10.0.powf(MathUtils::approximate_log10_sum_log10(
                -0.1 * min_qual as f64,
                -0.1 * max_qual as f64,
            ))
        } else {
            self.match_to_match_prob[((max_qual * (max_qual + 1)) >> 1) + min_qual]
        }
    }

    pub fn match_to_match_prob_static(ins_qual: u8, del_qual: u8) -> f64 {
        let mut min_qual;
        let mut max_qual;
        if ins_qual <= del_qual {
            min_qual = ins_qual;
            max_qual = del_qual;
        } else {
            min_qual = del_qual;
            max_qual = ins_qual;
        }

        1.0 - 10.0.powf(MathUtils::approximate_log10_sum_log10(
            -0.1 * min_qual as f64,
            -0.1 * max_qual as f64,
        ))
    }

    /**
     * Returns the log-probability that neither of two events takes place.
     *
     *
     * We assume that both events never occur together and that delQual is the conditional probability (qual. encoded)
     * of the second event, given the first event didn't took place. So that the probability of no event is:
     *
     * 1 - ProbErr(insQual) - ProbErr(delQual)
     *
     * @param insQual PhRED scaled quality/probability of an insertion.
     * @param delQual PhRED scaled quality/probability of a deletion.
     *
     * @return a value between 0 and -Inf.
     */
    pub fn match_to_match_prob_log10(&self, ins_qual: usize, del_qual: usize) -> f64 {
        let mut min_qual;
        let mut max_qual;
        if ins_qual <= del_qual {
            min_qual = ins_qual;
            max_qual = del_qual;
        } else {
            min_qual = del_qual;
            max_qual = ins_qual;
        }

        if (QualityUtils::MAX_QUAL as usize) < max_qual {
            (-std::cmp::min(
                OrderedFloat(1.0),
                OrderedFloat(10.0.powf(MathUtils::approximate_log10_sum_log10(
                    -0.1 * min_qual as f64,
                    -0.1 * max_qual as f64,
                ))),
            )
            .into_inner())
            .ln_1p()
                * *INV_LN10
        } else {
            self.match_to_match_log10[((max_qual * (max_qual + 1)) >> 1) + min_qual]
        }
    }
}
