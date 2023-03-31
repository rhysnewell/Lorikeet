use ordered_float::OrderedFloat;

use crate::utils::simple_interval::SimpleInterval;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 */
#[derive(Debug)]
pub struct RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    pub genotype_likelihoods: Vec<f64>,
    pub final_phred_scaled_genotype_likelihoods: Vec<i32>,
    pub ref_depth: i32,
    pub non_ref_depth: i32,
    pub read_counts: i32,
    pub loc: SimpleInterval,
}

impl RefVsAnyResult {
    pub fn new(likelihood_capacity: usize, pos: usize, tid: usize) -> RefVsAnyResult {
        RefVsAnyResult {
            genotype_likelihoods: vec![0.0; likelihood_capacity],
            final_phred_scaled_genotype_likelihoods: vec![0; likelihood_capacity],
            ref_depth: 0,
            non_ref_depth: 0,
            read_counts: 0,
            loc: SimpleInterval::new(tid, pos, pos),
        }
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    pub fn get_dp(&self) -> i32 {
        return self.ref_depth + self.non_ref_depth;
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    pub fn get_ad(&self) -> Vec<i32> {
        return vec![self.ref_depth, self.non_ref_depth];
    }

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihood_capacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihood_capacity} is negative.
     */
    pub fn ref_vs_any_result(&mut self, likelihood_capacity: usize) {
        self.genotype_likelihoods = vec![0.0; likelihood_capacity];
        self.final_phred_scaled_genotype_likelihoods = vec![0; likelihood_capacity];
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     */
    pub fn get_genotype_likelihoods_capped_by_hom_ref_likelihood(&self) -> Vec<f64> {
        let mut output = vec![0.0; self.genotype_likelihoods.len()];

        for i in 0..self.genotype_likelihoods.len() {
            output[i] = std::cmp::min(
                OrderedFloat(self.genotype_likelihoods[i]),
                OrderedFloat(self.genotype_likelihoods[0]),
            )
            .into()
        }

        return output;
    }
}
