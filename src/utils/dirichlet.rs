use statrs::function::gamma;
use utils::math_utils::MathUtils;
use num::traits::Float;
use rayon::prelude::*;

pub struct Dirichlet {
    alpha: Vec<f64>
}

impl Dirichlet{
    pub fn new(alpha: &[f64]) -> Dirichlet {
        Dirichlet {
            alpha: alpha.to_vec()
        }
    }

    /**
     * Create a symmetric distribution Dir(a/K, a/K, a/K . . .) where K is the number of states and
     * a is the concentration.
     */
    pub fn symmetric_dirichlet(num_states: usize, concentration: f64) -> Dirichlet {
        if num_states <= 0 {
            panic!("Must have at least one state")
        }
        if concentration <= 0.0 {
            panic!("Concentration must be positive")
        }

        Dirichlet::new(&vec![concentration / (num_states as f64); num_states])
    }

    // in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
    // Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exp of the expected log weights.
    // note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
    // amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
    // normalized weights.  In that case the calling code must normalize the weights.
    pub fn effective_multinomial_weights(&self) -> Vec<f64> {
        let digamma_of_sum = gamma::digamma(self.alpha.par_iter().sum::<f64>());
        let result = self.alpha.par_iter().map(|a| {
            (gamma::digamma(*a) - digamma_of_sum).exp()
        }).collect::<Vec<f64>>();

        return result
    }

    pub fn effective_log10_multinomial_weights(&self) -> Vec<f64> {
        let digamma_of_sum = gamma::digamma(self.alpha.par_iter().sum::<f64>());
        let result = self.alpha.par_iter().map(|a| {
            (gamma::digamma(*a) - digamma_of_sum) * MathUtils::LOG10_E
        }).collect::<Vec<f64>>();

        return result
    }

    pub fn effective_log_multinomial_weights(&self) -> Vec<f64> {
        let digamma_of_sum = gamma::digamma(self.alpha.par_iter().sum::<f64>());
        let result = self.alpha.par_iter().map(|a| {
            gamma::digamma(*a) - digamma_of_sum
        }).collect::<Vec<f64>>();

        return result
    }

    pub fn mean_weights(&self) -> Vec<f64> {
        let sum = self.alpha.par_iter().sum::<f64>();
        let result = self.alpha.par_iter().map(|x| {
            *x / sum
        }).collect::<Vec<f64>>();

        return result
    }

    pub fn log10_mean_weights(&self) -> Vec<f64> {
        let sum = self.alpha.par_iter().sum::<f64>();
        let result = self.alpha.par_iter().map(|x| {
            (*x / sum).log10()
        }).collect::<Vec<f64>>();

        return result
    }

    pub fn size(&self) -> usize {
        self.alpha.len()
    }
}