use utils::quality_utils::QualityUtils;
use num::traits::Float;
use rayon::prelude::*;
use libm;

pub struct NaturalLogUtils {}

impl NaturalLogUtils {
    pub const LOG_ONE_HALF: f64 = 0.5.ln();
    pub const LOG_ONE_THIRD: f64 = (1.0/3.0).ln();
    const LOG1MEXP_THRESHOLD: f64 = 0.5.ln();
    const PHRED_TO_ERROR_PROB_FACTOR: f64 = -(10.0.ln())/10.0;

    const qual_to_log_prob_cache: Vec<f64> = (0..QualityUtils::MAX_QUAL as usize).into_par_iter().map(|n| {
        NaturalLogUtils::log1mexp(NaturalLogUtils::qual_to_log_error_prob(n as f64))
    }).collect::<Vec<f64>>();

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
            return std::f64::NAN
        }
        if a == 0.0 {
            return std::f64::NEG_INFINITY
        }

        if a < NaturalLogUtils::LOG1MEXP_THRESHOLD {
            NaturalLogUtils::log1p(-(a).exp())
        } else {
            (-(a.exp_m1())).ln()
        }
    }

    #[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
    pub fn log1p(x: f64) -> f64 {
        libm::log1p(x)
    }

    pub fn qual_to_log_error_prob(qual: f64) -> f64 {
        if qual < 0.0 {
            panic!("Qual must be >= 0 but got {}", qual)
        }

        return qual * NaturalLogUtils::PHRED_TO_ERROR_PROB_FACTOR
    }
}