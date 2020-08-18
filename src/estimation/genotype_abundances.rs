use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct Genotype {
    // The variants in a genotype, not order is random but consistent with variant_weights
    // pub variants: Vec<variants::Variant>,
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

impl Genotype {
    pub fn new(index: usize) -> Genotype {
        Genotype {
            index,
            variant_weights: vec!(),
            variant_genotype_ids: vec!(),
            abundance_weight: 1.,
        }
    }
}

pub fn calculate_abundances(
    sample_genotypes: &mut Vec<Genotype>,
) {
    // Placeholder stopping criterion vector
    // Going to check for small changes between theta_prev and theta_curr
    let mut theta_prev_mean: f64;
    let mut theta_curr_mean: f64 = 1.;
    // the difference between theta curr and theta prev
    let mut omega = 1.;
    // The minimum value of omega before stopping
    let eps = 0.001;


    let mut theta_curr = vec![1.; sample_genotypes.len()];
    let mut n = 0;

    while omega > eps {
        // Update theta values
        theta_prev_mean = theta_curr_mean;

        for (index, genotype) in sample_genotypes.iter_mut().enumerate() {
            if genotype.variant_weights.len() == 0 {
                break
            }

            // Step 1: update variant weights
            genotype.variant_weights =
                genotype.variant_weights
                    .par_iter()
                    .enumerate()
                    .map(|(variant_index, w)| {
                        let pooled_weights = genotype.variant_genotype_ids[variant_index]
                            .par_iter()
                            .fold_with(0., |acc, genotype_index| acc + theta_curr[*genotype_index])
                            .sum::<f64>();
                        w*(genotype.abundance_weight / pooled_weights)
                    })
                    .collect();

            // Step 2: update abundance weight based on mean of variant weights
            genotype.abundance_weight =
                genotype.variant_weights.par_iter().sum::<f64>() / genotype.variant_weights.len() as f64;

            // Update list of theta_curr
            theta_curr[index] = genotype.abundance_weight;
        }

        // Update theta_curr_mean and omega
        theta_curr_mean = theta_curr.par_iter().sum::<f64>() / theta_curr.len() as f64;
        omega = (theta_curr_mean - theta_prev_mean).abs();

        debug!(
            "Theta Current {:?} \
            theta curr mean {} \
            Omega {}",
            &theta_curr,
            &theta_prev_mean,
            &omega,
        );
        n += 1;
    }
    debug!("EM Algorithm Finished in {} iterations", n);
}