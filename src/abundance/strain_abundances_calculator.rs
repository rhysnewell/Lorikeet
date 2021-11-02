use rayon::prelude::*;

/// Implementation of the EM algorithm used by centrifuge
///
/// @author Rhys Newell <rhys.newell@.hdr.qut.edu.au>
#[derive(Debug, Clone)]
pub struct StrainAbundanceCalculator {
    // The genotype index
    pub index: usize,
    // The expected weights of each variant that decide the strain abundance
    pub variant_weights: Vec<f64>,
    // The indices of all the genotypes that share this variant
    pub variant_genotype_ids: Vec<Vec<usize>>,
    // The geometric weighting of this strains abundance
    // Strain abundance would be calculated as abundance_weight * total_abundance
    // denoted as `theta`
    pub abundance_weight: f64,
}

impl StrainAbundanceCalculator {
    pub fn new(index: usize) -> Self {
        Self {
            index,
            variant_weights: vec![],
            variant_genotype_ids: vec![],
            abundance_weight: 1.,
        }
    }

    pub fn calculate_abundances(sample_genotypes: &mut Vec<Self>, eps: f64) {
        // the difference between theta curr and theta prev
        let mut omega = 1.;

        let mut theta_prev;
        let mut theta_curr = vec![1.; sample_genotypes.len()];
        let mut n = 0;

        while omega > eps {
            // Update theta values
            theta_prev = theta_curr.clone();

            for (index, genotype) in sample_genotypes.iter_mut().enumerate() {
                if (genotype.abundance_weight - eps).abs() <= f64::EPSILON || genotype.abundance_weight.is_infinite() {
                    continue;
                }

                // Step 1: update variant weights
                genotype.variant_weights = genotype
                    .variant_weights
                    .par_iter()
                    .enumerate()
                    .map(|(variant_index, w)| {
                        let mut pooled_weights = genotype.variant_genotype_ids[variant_index]
                            .iter()
                            .map(|genotype_index| theta_curr[*genotype_index])
                            .sum::<f64>();

                        if pooled_weights <= f64::EPSILON {
                            pooled_weights = 1.0;
                        }
                        debug!(
                            "Variant index {} weight {} pooled weight {}",
                            variant_index, w, pooled_weights
                        );
                        // if variant weights are between 0 and 1
                        let w = w * (genotype.abundance_weight / pooled_weights);
                        w
                    })
                    .collect();

                // Step 2: update abundance weight based on mean of variant weights
                genotype.abundance_weight = genotype.variant_weights.iter().sum::<f64>()
                    / genotype.variant_weights.len() as f64;

                debug!(
                    "Index {} abundance weight {} variant weights",
                    index,
                    genotype.abundance_weight, //&genotype.variant_weights
                );

                if genotype.abundance_weight.is_nan()
                    // || genotype.variant_weights.contains(&0.0)
                    || genotype.abundance_weight.is_infinite()
                {
                    genotype.abundance_weight = 0.;
                    genotype.variant_weights = vec![];
                }

                // Update list of theta_curr
                theta_curr[index] = genotype.abundance_weight;
            }

            // Update omega
            omega = theta_curr.iter().zip(theta_prev.iter()).map(|(curr, prev)| (curr - prev).abs()).sum::<f64>();

            debug!(
                "Theta Current {:?} Prev {:?} Omega {}",
                &theta_curr, &theta_prev, &omega,
            );
            n += 1;
        }
        debug!("EM Algorithm Finished in {} iterations", n);
    }
}
