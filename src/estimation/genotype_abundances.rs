use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use model::variants;

#[derive(Debug)]
pub struct Genotype {
    // The variants in a genotype, not order is random but consistent with variant_weights
    // pub variants: Vec<variants::Variant>,
    // The genotype index
    pub index: usize,
    // The expected weights of each variant that decide the strain abundance
    pub variant_weights: Vec<f64>,
    // The geometric weighting of this strains abundance
    // Strain abundance would be calculated as abundance_weight * total_abundance
    // denoted as `theta`
    pub abundance_weight: f64,
}

impl Genotype {
    pub fn new(index: usize) -> Genotype {
        Genotype {
            // variants: vec!(),
            index,
            variant_weights: vec!(),
            abundance_weight: 1.,
        }
    }
}

pub fn calculate_abundances(
    sample_genotypes: &mut Vec<Genotype>,
) {
    // Placeholder stopping criterion vector
    // Going to check for small changes between theta_prev and theta_curr
    let mut theta_prev_mean: f64 = 1.;
    let mut theta_curr_mean: f64 = 2.;
    let mut theta_curr = vec![1.; sample_genotypes.len()];
    let eps = 0.001;

    while (theta_curr_mean - theta_prev_mean).abs() > eps {
        // Update theta values
        theta_prev_mean = theta_curr_mean;

        for (index, genotype) in sample_genotypes.iter_mut().enumerate() {
            // Step 1: update variant weights
            genotype.variant_weights =
                genotype.variant_weights
                    .par_iter()
                    .map(|w| w*genotype.abundance_weight)
                    .collect();

            // Step 2: update abundance weight based on mean of variant weights
            genotype.abundance_weight =
                genotype.variant_weights.par_iter().sum::<f64>() / genotype.variant_weights.len() as f64;

            // Update list of theta_curr
            theta_curr[index] = genotype.abundance_weight;
        }

        // Update theta_curr_mean
        theta_curr_mean = theta_curr.par_iter().sum::<f64>() / theta_curr.len() as f64;
    }
}