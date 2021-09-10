use libm;
use num::traits::Float;
use rayon::prelude::*;
use utils::quality_utils::QualityUtils;

lazy_static! {
    pub static ref LOG_ONE_HALF: f64 = 0.5.ln();
    pub static ref LOG_ONE_THIRD: f64 = (1.0 / 3.0).ln();
    static ref LOG1MEXP_THRESHOLD: f64 = 0.5.ln();
    static ref PHRED_TO_ERROR_PROB_FACTOR: f64 = -(10.0.ln()) / 10.0;
    static ref qual_to_log_prob_cache: Vec<f64> = (0..QualityUtils::MAX_QUAL + 1)
        .into_iter()
        .map(|n| { NaturalLogUtils::log1mexp(NaturalLogUtils::qual_to_log_error_prob(n)) })
        .collect::<Vec<f64>>();
}

pub struct NaturalLogUtils {}

impl NaturalLogUtils {
    /**
     * Calculates {@code log(1-exp(a))} without losing precision.
     *
     * <p>
     *     This is based on the approach described in:
     *
     * </p>
     * <p>
     *     Maechler M, Accurately Computing log(1-exp(-|a|)) Assessed by the Rmpfr package, 2012 <br/>
     *     <a ref="http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf">Online document</a>.
     *
     * </p>
     *
     * @param a the input exponent.
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    pub fn log1mexp(a: f64) -> f64 {
        if a > 0.0 {
            return std::f64::NAN;
        }
        if a == 0.0 {
            return std::f64::NEG_INFINITY;
        }

        if a < *LOG1MEXP_THRESHOLD {
            NaturalLogUtils::log1p(-(a).exp())
        } else {
            (-(a.exp_m1())).ln()
        }
    }

    #[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
    pub fn log1p(x: f64) -> f64 {
        libm::log1p(x)
    }

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
    pub fn qual_to_log_error_prob(qual: u8) -> f64 {
        return (qual as f64) * *PHRED_TO_ERROR_PROB_FACTOR;
    }

    pub fn qual_to_log_prob(qual: u8) -> f64 {
        return qual_to_log_prob_cache[qual as usize]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }
}
