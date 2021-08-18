use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use haplotype::homogenous_ploidy_model::PloidyModel;
use model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use model::allele_likelihoods::AlleleLikelihoods;
use model::allele_list::AlleleList;
use model::byte_array_allele::ByteArrayAllele;

pub struct IndependentSamplesGenotypesModel {
    cache_allele_count_capacity: usize,
    cache_ploidy_capacity: usize,
    likelihood_calculators: Vec<Vec<Option<GenotypeLikelihoodCalculator>>>,
    calculators: GenotypeLikelihoodCalculators,
}

impl IndependentSamplesGenotypesModel {
    const DEFAULT_CACHE_PLOIDY_CAPACITY: usize = 10;
    const DEFAULT_CACHE_ALLELE_CAPACITY: usize = 50;

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    pub fn default() -> Self {
        Self::new(
            Self::DEFAULT_CACHE_PLOIDY_CAPACITY,
            Self::DEFAULT_CACHE_ALLELE_CAPACITY,
        )
    }

    pub fn new(
        calculator_cache_ploidy_capacity: usize,
        calculator_cache_allele_capacity: usize,
    ) -> Self {
        Self {
            cache_ploidy_capacity: calculator_cache_ploidy_capacity,
            cache_allele_count_capacity: calculator_cache_allele_capacity,
            likelihood_calculators: vec![
                vec![None; calculator_cache_allele_capacity];
                calculator_cache_ploidy_capacity
            ],
            calculators: GenotypeLikelihoodCalculators::build_empty(),
        }
    }

    pub fn calculate_likelihoods<P: PloidyModel>(
        &mut self,
        genotyping_alleles: &AlleleList<ByteArrayAllele>,
        read_likelihoods: &AlleleLikelihoods,
        ploidy_model: &P,
        padded_reference: &[u8],
        offset_for_into_event: usize,
    ) -> Vec<GenotypeLikelihoods> {
        let permutation = read_likelihoods
            .get_allele_list()
            .permutation(genotyping_alleles.clone());
        let mut allele_likelihood_matrix_mapper = AlleleLikelihoodMatrixMapper::new(permutation);

        let sample_count = read_likelihoods.samples.len();
        assert!(sample_count > 0, "Sample count must be greater than 0");
        let mut genotype_likelihoods = Vec::with_capacity(sample_count);
        let allele_count = genotyping_alleles.len();

        let mut likelihoods_calculator =
            self.get_likelihood_calculator(ploidy_model.sample_ploidy(0), allele_count);
        for i in 0..sample_count {
            let sample_ploidy = ploidy_model.sample_ploidy(i);
            let sample_likelihoods = &read_likelihoods.values_by_sample_index[i];

            match likelihoods_calculator {
                None => {
                    let mut likelihoods_calculator =
                        Self::get_uncached_likelihood_calculator(sample_ploidy, allele_count);
                    genotype_likelihoods.push(likelihoods_calculator.genotype_likelihoods(
                        sample_likelihoods,
                        &allele_likelihood_matrix_mapper,
                    ));
                }
                Some(ref mut likelihoods_calculator) => {
                    genotype_likelihoods.push(likelihoods_calculator.genotype_likelihoods(
                        sample_likelihoods,
                        &allele_likelihood_matrix_mapper,
                    ));
                }
            };
        }

        return genotype_likelihoods;
    }

    fn get_likelihood_calculator(
        &mut self,
        sample_ploidy: usize,
        allele_count: usize,
    ) -> Option<&mut GenotypeLikelihoodCalculator> {
        if sample_ploidy >= self.cache_ploidy_capacity
            || allele_count >= self.cache_allele_count_capacity
        {
            return None;
        }

        if self.likelihood_calculators[sample_ploidy][allele_count].is_some() {
            return self.likelihood_calculators[sample_ploidy][allele_count].as_mut();
        } else {
            let new_one = GenotypeLikelihoodCalculators::get_instance(sample_ploidy, allele_count);
            self.likelihood_calculators[sample_ploidy][allele_count] = Some(new_one);
            return self.likelihood_calculators[sample_ploidy][allele_count].as_mut();
        }
        // match &self.likelihood_calculators[sample_ploidy][allele_count] {
        //     Some(result) => {
        //         return Some(result)
        //     },
        //     None => {
        //         // Put this here to satisfy the borrow checker
        //         let new_one = GenotypeLikelihoodCalculators::get_instance(sample_ploidy, allele_count);
        //         self.likelihood_calculators[sample_ploidy][allele_count] = Some(new_one);
        //         return self.likelihood_calculators[sample_ploidy][allele_count].as_ref();
        //     }
        // };
    }

    fn get_uncached_likelihood_calculator(
        sample_ploidy: usize,
        allele_count: usize,
    ) -> GenotypeLikelihoodCalculator {
        return GenotypeLikelihoodCalculators::get_instance(sample_ploidy, allele_count);
    }
}
