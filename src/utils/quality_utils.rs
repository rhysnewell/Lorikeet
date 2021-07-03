use utils::math_utils::MathUtils;
use num::traits::Float;
use rayon::prelude::*;
use std::cmp::{min, max};

lazy_static! {
    static ref MIN_LOG10_SCALED_QUAL: f64 = (std::f64::MIN).log10();
    static ref MIN_PHRED_SCALED_QUAL: f64 = -10.0 * *MIN_LOG10_SCALED_QUAL;
}

pub struct QualityUtils {}

impl QualityUtils {


    const MAX_REASONABLE_Q_SCORE: u8 = 60;

    pub const MAX_QUAL: u8 = 254;

    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * The calculation is extremely efficient
     *
     * @param qual a phred-scaled quality score encoded as a double
     * @return log of probability (0.0-1.0)
     */
    pub fn qual_to_error_prob_log10(qual: f64) -> f64 {
        qual * -0.1
    }

    /**
     * Calculate the sum of phred scores.
     * @param phreds the phred score values.
     * @return the sum.
     */
    pub fn phred_sum(phreds: &[f64]) -> f64 {
        match phreds.len() {
            0 => std::f64::MAX,
            1 => phreds[0],
            2 => -10.0 * MathUtils::log10_sum_log10_two_values(phreds[0] * -0.1, phreds[1] * -0.1),
            3 => -10.0 * MathUtils::log10_sum_log10_three_values(phreds[0] * -0.1, phreds[1] * -0.1, phreds[2] * -0.1),
            _ => {
                let log10_vals = phreds.par_iter().map(|p| {
                    *p * -0.1
                }).collect::<Vec<f64>>();
                -10.0 * MathUtils::log10_sum_log10(&log10_vals, 0, log10_vals.len())
            }
        }
    }

    /**
     * Convert a phred-scaled quality score to its probability of being true (Q30 => 0.999)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a double value, this function must call Math.pow so can be quite expensive
     *
     * @param qual a phred-scaled quality score encoded as a double.  Can be non-integer values (30.5)
     * @return a probability (0.0-1.0)
     */
    pub fn qual_to_prob(qual: f64) -> f64 {
        if qual < 0.0 {
            panic!("Qual must be >= 0.0 but got {}", qual)
        }

        return 1.0 - QualityUtils::qual_to_error_prob(qual)
    }

    /**
     * Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a double value, this function must call Math.pow so can be quite expensive
     *
     * @param qual a phred-scaled quality score encoded as a double.  Can be non-integer values (30.5)
     * @return a probability (0.0-1.0)
     */
    pub fn qual_to_error_prob(qual: f64) -> f64 {
        if qual < 0.0 {
            panic!("Qual must be >= 0.0 but got {}", qual)
        }

        10.0.powf(qual / -10.0)
    }

    // ----------------------------------------------------------------------
    //
    // Functions to convert a probability to a phred-scaled quality score
    //
    // ----------------------------------------------------------------------

    /**
     * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
     *
     * Note, this function caps the resulting quality score by the public static value MAX_SAM_QUAL_SCORE
     * and by 1 at the low-end.
     *
     * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
     * @return a quality score (0-MAX_SAM_QUAL_SCORE)
     */
    pub fn error_prob_to_qual(error_rate: f64) -> u8 {
        return Self::error_prob_to_qual_with_max_qual(error_rate, Self::MAX_QUAL)
    }

    /**
     * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
     *
     * Note, this function caps the resulting quality score by the public static value MIN_REASONABLE_ERROR
     * and by 1 at the low-end.
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
     * @return a quality score (0-maxQual)
     */
    pub fn error_prob_to_qual_with_max_qual(error_rate: f64, max_qual: u8) -> u8 {
        assert!(MathUtils::is_valid_probability(error_rate), "Error rate is not a valid probability");
        let d = (-10.0 * error_rate.log10()).round();
        return Self::bound_qual(d as u8, max_qual)
    }

    /**
     * Return a quality score that bounds qual by maxQual and 1
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual the uncapped quality score as an integer.  Can be < 0 (which may indicate an error in the
     *             client code), which will be brought back to 1, but this isn't an error, as some
     *             routines may use this functionality (BaseRecalibrator, for example)
     * @param maxQual the maximum quality score, must be less < 255
     * @return the bounded quality score
     */
    pub fn bound_qual(qual: u8, max_qual: u8) -> u8 {
        return max(min(qual, max_qual), 1)
    }
}