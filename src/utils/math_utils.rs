use rayon::prelude::*;
use std::ops::{Add, Sub, Mul};
use mathru::special;

pub struct MathUtils {

}

impl MathUtils {
    pub const LOG10_P_OF_ZERO: f64 = -1000000.0;
    pub const LOG10_ONE_HALF: f64 = (0.5 as f64).log10();
    pub const LOG10_ONE_THIRD: f64 = -((3.0 as f64).log10());
    pub const LOG_ONE_THIRD: f64 = -((3.0 as f64).ln());
    pub const INV_LOG_2: f64 = (1.0 as f64) / (2.0 as f64).ln();
    const LOG_10: f64 = (10. as f64).ln();
    const INV_LOG_10: f64 = (1.0) / MathUtils::LOG_10;
    pub const LOG10_E: f64 = std::f64::consts::E.log10();
    c

    const ROOT_TWO_PI: f64 = (2.0 * std::f64::consts::PI).sqrt();

    // const LOG_10_CACHE: Log10Cache
    // const LOG_10_FACTORIAL_CACHE: Log10FactorialCache
    // const DIGAMMA_CACHE: DiGammaCache

    pub fn normalize_pls<T: Sub>(pls: &[T]) -> Vec<T> {
        let mut new_pls = vec![0 as T; pls.len()];
        let smallest = pls.min();
        new_pls.par_iter_mut().enumerate()
            .for_each(|(i, pl)| {
                pl = pls[i] - smallest;
            });

        return new_pls
    }

    /**
    * Element by elemnt addition of two vectors
    */
    pub fn ebe_add<T: Add>(a: &[T], b: &[T]) -> Vec<T> {
        let mut z = Vec::with_capacity(a.len());
        for (i, (aval, bval)) in a.iter().zip(&b).enumerate() {
            z[i] = aval + bval;
        }
        z
    }

    /**
    * Element by elemnt subtraction of two vectors
    */
    pub fn ebe_subtract<T: Sub>(a: &[T], b: &[T]) -> Vec<T> {
        let mut z = Vec::with_capacity(a.len());
        for (i, (aval, bval)) in a.iter().zip(&b).enumerate() {
            z[i] = aval - bval;
        }
        z
    }

    /**
     * Converts LN to LOG10
     * @param ln log(x)
     * @return log10(x)
     */
    pub fn log_to_log10<T: Float>(ln: T) -> T {
        ln * MathUtils::LOG10_E
    }

    /**
     * @see #binomialCoefficient(int, int) with log10 applied to result
     */
    pub fn log10_binomial_coeffecient<T: Float + Copy>(n: T, k: T) -> T {
        return MathUtils::log10_factorial(n)
            - MathUtils::log10_factorial(k)
            - MathUtils::log10_factorial(n - k)
    }

    pub fn log10_factorial<T: Float>(n: T) -> T {
        n.log_gamma() * MathUtils::LOG10_E
    }

    pub fn max_element_index<T: PartialOrd + PartialEq>(
        array: &[T],
        start: usize,
        finish: usize
    ) -> usize {

        let result = array[start..finish]
            .par_iter()
            .enumerate()
            .max_by(|&(_, item)| item);

        return result
    }

    pub fn normalize_log10(
        mut array: Vec<f64>,
        take_log10_of_output: bool,
    ) -> Vec<f64> {
        let log10_sum = MathUtils::log10_sum_log10(&array, 0, array.len());
        array
            .par_iter_mut()
            .for_each(|x| x = x - log10_sum);
        if take_log10_of_output {
            array
                .par_iter_mut()
                .for_each(|x| x = (10.0).powf(x))
        }
        return array
    }

    pub fn log10_sum_log10(log10_values: &Vec<f64>, start: usize, finish: usize) -> f64 {
        if start >= finish {
            return std::f64::NEG_INFINITY
        }

        let max_element_index = MathUtils::max_element_index(
            log10_values,
            start,
            finish
        );

        let max_value = log10_values[max_element_index];

        if max_value == std::f64::NEG_INFINITY {
            return max_value
        }

        let sum_tot = log10_values[start..finish]
            .par_iter()
            .enumerate()
            .filter(|&(index, value)| {
                i != max_element_index || value != std::f64::NEG_INFINITY
            })
            .map(|&(_, value)| {
                value
            }).sum();

        if sum_tot == std::f64::NAN || sum == std::f64::INFINITY {
            panic!("log10 p: Values must be non-infinite and non-NAN")
        }

        return max_value + (if sum_tot != 1.0 { sum.log10() } else { 0.0 })

    }
}

#[Derive(Debug, Clone, Copy)]
pub struct RunningAverage<T: Float + Copy> {
    mean: T,
    s: T,
    obs_count: usize,
}

impl<T: Float + Copy> RunningAverage<T> {
    pub fn new() -> RunningAverage<T> {
        RunningAverage {
            mean: 0.0,
            s: 0.0,
            obs_count: 0
        }
    }

    pub fn add<T: Float + Copy>(&mut self, obs: T) {
        self.obs_count += 1;
        let old_mean = self.mean;
        self.mean += (obs - self.mean) / self.obs_count;
        self.s += (obs - old_mean) * (obs - old_mean)
    }

    pub fn add_all<T: Float + Copy>(&mut self, col: &[T]) {
        for obs in col {
            self.add(obs)
        }
    }

    pub fn mean<T: Float + Copy>(&self) -> T {
        self.mean
    }

    pub fn stddev<T: Float + Copy>(&self) -> T {
        (self.s / (self.obs_count - 1)).sqrt()
    }

    pub fn var<T: Float + Copy>(&self) -> T {
        self.s / (self.obs_count - 1)
    }

    pub fn obs_count(&self) -> usize {
        self.obs_count
    }
}