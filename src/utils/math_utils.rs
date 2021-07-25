use ordered_float::{NotNan, OrderedFloat};
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use std::clone::Clone;
use std::ops::{Add, AddAssign, Sub};
use utils::natural_log_utils::NaturalLogUtils;

lazy_static! {
    static ref cache: Vec<f64> =
        (0..((JacobianLogTable::MAX_TOLERANCE / JacobianLogTable::TABLE_STEP) + 1.0) as usize)
            .into_par_iter()
            .map(|k| { 1.0 + (10.0_f64).powf(-(k as f64) * JacobianLogTable::INV_STEP) })
            .collect::<Vec<f64>>();
    pub static ref LOG10_ONE_HALF: f64 = (0.5 as f64).log10();
    pub static ref LOG10_ONE_THIRD: f64 = -((3.0 as f64).log10());
    pub static ref LOG_ONE_THIRD: f64 = -((3.0 as f64).ln());
    pub static ref INV_LOG_2: f64 = (1.0 as f64) / (2.0 as f64).ln();
    static ref LOG_10: f64 = (10. as f64).ln();
    static ref INV_LOG_10: f64 = (1.0) / *LOG_10;
    pub static ref LOG10_E: f64 = std::f64::consts::E.log10();
    static ref ROOT_TWO_PI: f64 = (2.0 * std::f64::consts::PI).sqrt();
}

pub struct MathUtils {}

impl MathUtils {
    pub const LOG10_P_OF_ZERO: f64 = -1000000.0;

    // const LOG_10_CACHE: Log10Cache
    // const LOG_10_FACTORIAL_CACHE: Log10FactorialCache
    // const DIGAMMA_CACHE: DiGammaCache

    pub fn normalize_pls(pls: &[f64]) -> Vec<f64> {
        let mut new_pls = vec![0.0; pls.len()];
        let smallest = *pls
            .par_iter()
            .min_by_key(|x| OrderedFloat(**x))
            .unwrap_or(&std::f64::NAN);
        new_pls.par_iter_mut().enumerate().for_each(|(i, pl)| {
            *pl = pls[i] - smallest;
        });

        return new_pls;
    }

    /**
     * Element by elemnt addition of two vectors in place
     */
    pub fn ebe_add_in_place<T: Add + Copy + AddAssign>(a: &mut [T], b: &[T]) {
        a.iter_mut().enumerate().for_each(|(i, val)| *val += b[i]);
    }

    /**
     * Element by elemnt addition of two vectors
     */
    pub fn ebe_add<T: Add + Copy + Add<Output = T>>(a: &[T], b: &[T]) -> Vec<T> {
        let mut z = Vec::with_capacity(a.len());
        for (i, (aval, bval)) in a.iter().zip(b).enumerate() {
            z[i] = *aval + *bval;
        }
        z
    }

    /**
     * Element by elemnt subtraction of two vectors
     */
    pub fn ebe_subtract<T: Sub + Copy + Sub<Output = T>>(a: &[T], b: &[T]) -> Vec<T> {
        let mut z = Vec::with_capacity(a.len());
        for (i, (aval, bval)) in a.iter().zip(b).enumerate() {
            z[i] = *aval - *bval;
        }
        z
    }

    /**
     * Converts LN to LOG10
     * @param ln log(x)
     * @return log10(x)
     */
    pub fn log_to_log10(ln: f64) -> f64 {
        ln * *LOG10_E
    }

    /**
     * @see #binomialCoefficient(int, int) with log10 applied to result
     */
    pub fn log10_binomial_coeffecient(n: f64, k: f64) -> f64 {
        return MathUtils::log10_factorial(n)
            - MathUtils::log10_factorial(k)
            - MathUtils::log10_factorial(n - k);
    }

    pub fn log10_factorial(n: f64) -> f64 {
        ln_gamma(n) * *LOG10_E
    }

    /**
     * Gets the maximum element's index of an array of f64 values
     * Rather convoluted due to Rust not allowing proper comparisons between floats
     */
    pub fn max_element_index(array: &[f64], start: usize, finish: usize) -> usize {
        let non_nan_floats: Vec<_> = array[start..finish]
            .par_iter()
            .cloned()
            .map(NotNan::new)
            .filter_map(Result::ok)
            .collect::<Vec<_>>();

        let max = non_nan_floats.par_iter().max().unwrap();
        let index = non_nan_floats
            .par_iter()
            .position_first(|f| f == max)
            .unwrap();

        return index;
    }

    pub fn normalize_log10(mut array: Vec<f64>, take_log10_of_output: bool) -> Vec<f64> {
        let log10_sum = MathUtils::log10_sum_log10(&array, 0, array.len());
        array.par_iter_mut().for_each(|x| *x = *x - log10_sum);
        if take_log10_of_output {
            array.par_iter_mut().for_each(|x| *x = 10.0_f64.powf(*x))
        }
        return array;
    }

    pub fn log10_sum_log10(log10_values: &[f64], start: usize, finish: usize) -> f64 {
        if start >= finish {
            return std::f64::NEG_INFINITY;
        }

        let max_element_index = MathUtils::max_element_index(log10_values, start, finish);

        let max_value = log10_values[max_element_index];

        if max_value == std::f64::NEG_INFINITY {
            return max_value;
        }

        let sum_tot = log10_values[start..finish]
            .par_iter()
            .enumerate()
            .filter(|(index, value)| {
                *index != max_element_index || **value != std::f64::NEG_INFINITY
            })
            .map(|(_, value)| value)
            .sum::<f64>();

        if sum_tot == std::f64::NAN || sum_tot == std::f64::INFINITY {
            panic!("log10 p: Values must be non-infinite and non-NAN")
        }

        return max_value + (if sum_tot != 1.0 { sum_tot.log10() } else { 0.0 });
    }

    pub fn log10_sum_log10_two_values(a: f64, b: f64) -> f64 {
        if a > b {
            a + (1. + 10.0_f64.powf(b - a))
        } else {
            b + (1. + 10.0_f64.powf(a - b))
        }
    }

    /**
     * Do the log-sum trick for three double values.
     * @param a
     * @param b
     * @param c
     * @return the sum... perhaps NaN or infinity if it applies.
     */
    pub fn log10_sum_log10_three_values(a: f64, b: f64, c: f64) -> f64 {
        if a >= b && a >= c {
            return a + (1.0 + 10.0_f64.powf(b - a)).log10() + 10.0_f64.powf(c - a);
        } else if b >= c {
            return b + (1.0 + 10.0_f64.powf(a - b)).log10() + 10.0_f64.powf(c - b);
        } else {
            return c + (1.0 + 10.0_f64.powf(a - c)).log10() + 10.0_f64.powf(b - c);
        }
    }

    /**
     * Given an array of log space (log or log10) values, subtract all values by the array maximum so that the max element in log space
     * is zero.  This is equivalent to dividing by the maximum element in real space and is useful for avoiding underflow/overflow
     * when the array's values matter only up to an arbitrary normalizing factor, for example, an array of likelihoods.
     *
     * @param array
     * @return the scaled-in-place array
     */
    pub fn scale_log_space_array_for_numeric_stability(array: &[f64]) -> Vec<f64> {
        let max_value: f64 = *array
            .par_iter()
            .max_by_key(|x| OrderedFloat(**x))
            .unwrap_or(&std::f64::NAN);
        let result = array
            .par_iter()
            .map(|x| *x - max_value)
            .collect::<Vec<f64>>();
        return result;
    }

    /**
     * See #normalizeFromLog10 but with the additional option to use an approximation that keeps the calculation always in log-space
     *
     * @param array
     * @param takeLog10OfOutput
     * @param keepInLogSpace
     *
     * @return array
     */
    //TODO: Check that this works
    pub fn normalize_from_log10(
        array: &[f64],
        take_log10_of_output: bool,
        keep_in_log_space: bool,
    ) -> Vec<f64> {
        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        let max_value: f64 = *array
            .par_iter()
            .max_by_key(|x| OrderedFloat(**x))
            .unwrap_or(&std::f64::NAN);

        // we may decide to just normalize in log space without converting to linear space
        if keep_in_log_space {
            let array: Vec<f64> = array
                .par_iter()
                .map(|x| *x - max_value)
                .collect::<Vec<f64>>();
            return array;
        }
        // default case: go to linear space
        let mut normalized: Vec<f64> = (0..array.len())
            .into_par_iter()
            .map(|i| 10.0_f64.powf(array[i] - max_value))
            .collect::<Vec<f64>>();

        let sum: f64 = normalized.par_iter().sum::<f64>();

        normalized.par_iter_mut().enumerate().for_each(|(i, x)| {
            *x = *x / sum;
            if take_log10_of_output {
                *x = x.log10();
                if *x < MathUtils::LOG10_P_OF_ZERO || x.is_infinite() {
                    *x = array[i] - max_value
                }
            }
        });

        return normalized;
    }

    pub fn is_valid_log10_probability(result: f64) -> bool {
        result <= 0.0
    }

    pub fn log10_to_log(log10: f64) -> f64 {
        return log10 * (*LOG_10);
    }
    /**
     * Calculates {@code log10(1-10^a)} without losing precision.
     *
     * @param a the input exponent.
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    pub fn log10_one_minus_pow10(a: f64) -> f64 {
        if a > 0.0 {
            return std::f64::NAN;
        }
        if a == 0.0 {
            return std::f64::NEG_INFINITY;
        }

        let b = a * *LOG_10;
        return NaturalLogUtils::log1mexp(b) * *INV_LOG_10;
    }

    pub fn approximate_log10_sum_log10(a: f64, b: f64) -> f64 {
        // this code works only when a <= b so we flip them if the order is opposite
        if a > b {
            return MathUtils::approximate_log10_sum_log10(b, a);
        } else if a == std::f64::NEG_INFINITY {
            return b;
        } else {
            // if |b-a| < tol we need to compute log(e^a + e^b) = log(e^b(1 + e^(a-b))) = b + log(1 + e^(-(b-a)))
            // we compute the second term as a table lookup with integer quantization
            // we have pre-stored correction for 0,0.1,0.2,... 10.0
            let diff = b - a;

            return b + if diff < JacobianLogTable::MAX_TOLERANCE {
                JacobianLogTable::get(diff)
            } else {
                0.0
            };
        }
    }

    pub fn well_formed_f64(val: f64) -> bool {
        return !val.is_nan() && !val.is_infinite();
    }

    /**
     * Calculate f(x) = Normal(x | mu = mean, sigma = sd)
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    pub fn normal_distribution(mean: f64, sd: f64, x: f64) -> f64 {
        assert!(sd >= 0.0, "Standard deviation must be >= 0.0");
        assert!(
            Self::well_formed_f64(mean) && Self::well_formed_f64(sd) && Self::well_formed_f64(x),
            "mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)"
        );

        return ((-((x - mean).powf(2.0))) / (2.0 * sd.powf(2.0))).exp() / (sd * *ROOT_TWO_PI);
    }

    /**
     * normalizes the real-space probability array.
     *
     * Does not assume anything about the values in the array, beyond that no elements are below 0.  It's ok
     * to have values in the array of > 1, or have the sum go above 0.
     *
     * @param array the array to be normalized
     * @return a newly allocated array corresponding the normalized values in array
     */
    pub fn normalize_sum_to_one(mut array: Vec<f64>) -> Vec<f64> {
        if array.len() == 0 {
            return array;
        }

        let sum = array.par_iter().sum::<f64>();
        assert!(
            sum >= 0.0,
            "Values in probability array sum to a negative number"
        );
        array.par_iter_mut().for_each(|x| *x = *x / sum);

        return array;
    }

    /**
     * Computes the entropy -p*ln(p) - (1-p)*ln(1-p) of a Bernoulli distribution with success probability p
     * using an extremely fast Pade approximation that is very accurate for all values of 0 <= p <= 1.
     *
     * See http://www.nezumi.demon.co.uk/consult/logx.htm
     */
    pub fn fast_bernoulli_entropy(p: f64) -> f64 {
        let product = p * (1.0 - p);
        return product * (11.0 + 33.0 * product) / (2.0 + 20.0 * product);
    }

    pub fn is_valid_probability(result: f64) -> bool {
        return result >= 0.0 && result <= 1.0;
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RunningAverage {
    mean: f64,
    s: f64,
    obs_count: usize,
}

impl RunningAverage {
    pub fn new() -> RunningAverage {
        RunningAverage {
            mean: 0.0,
            s: 0.0,
            obs_count: 0,
        }
    }

    pub fn add(&mut self, obs: f64) {
        self.obs_count += 1;
        let old_mean = self.mean;
        self.mean += (obs - self.mean) / self.obs_count as f64;
        self.s += (obs - old_mean) * (obs - old_mean)
    }

    pub fn add_all(&mut self, col: &[f64]) {
        for obs in col {
            self.add(*obs)
        }
    }

    pub fn mean(&self) -> f64 {
        self.mean
    }

    pub fn stddev(&self) -> f64 {
        (self.s / (self.obs_count - 1) as f64).sqrt()
    }

    pub fn var(&self) -> f64 {
        self.s / (self.obs_count - 1) as f64
    }

    pub fn obs_count(&self) -> usize {
        self.obs_count
    }
}

/**
 * Encapsulates the second term of Jacobian log identity for differences up to MAX_TOLERANCE
 */
struct JacobianLogTable {}

impl JacobianLogTable {
    // if log(a) - log(b) > MAX_TOLERANCE, b is effectively treated as zero in approximateLogSumLog
    // MAX_TOLERANCE = 8.0 introduces an error of at most one part in 10^8 in sums
    pub const MAX_TOLERANCE: f64 = 8.0;

    //  Phred scores Q and Q+1 differ by 0.1 in their corresponding log-10 probabilities, and by
    // 0.1 * log(10) in natural log probabilities.  Setting TABLE_STEP to an exact divisor of this
    // quantity ensures that approximateSumLog in fact caches exact values for integer phred scores
    pub const TABLE_STEP: f64 = 0.0001;
    pub const INV_STEP: f64 = (1.0) / JacobianLogTable::TABLE_STEP;

    pub fn get(difference: f64) -> f64 {
        let index = (difference * JacobianLogTable::INV_STEP).abs() as usize;
        return cache[index];
    }
}
