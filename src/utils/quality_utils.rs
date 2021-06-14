pub struct QualityUtils {}

impl QualityUtils {
    const MIN_LOG10_SCALED_QUAL: f64 = (std::f64::MIN).log10();
    const MIN_PHRED_SCALED_QUAL: f64 = -10.0 * QualityUtils::MIN_LOG10_SCALED_QUAL;

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
    pub fn qual_to_error_prob_log10<T: Float + Copy>(qual: T) -> T {
        qual * (-0.1 as T)
    }
}