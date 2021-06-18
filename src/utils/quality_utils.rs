use utils::math_utils::MathUtils;
use num::traits::Float;
use rayon::prelude::*;

lazy_static! {
    static ref MIN_LOG10_SCALED_QUAL: f64 = (std::f64::MIN).log10();
    static ref MIN_PHRED_SCALED_QUAL: f64 = -10.0 * *MIN_LOG10_SCALED_QUAL;
}

pub struct QualityUtils {}

impl QualityUtils {


    const MAX_REASONABLE_Q_SCORE: u64 = 60;

    pub const MAX_QUAL: u64 = 254;

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
}