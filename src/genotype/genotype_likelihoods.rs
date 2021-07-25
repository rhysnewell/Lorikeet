use statrs::function::factorial::binomial;
use std::collections::HashMap;
use utils::math_utils::MathUtils;

#[derive(Debug, Clone)]
pub struct GenotypeLikelihoods {
    num_likelihood_cache: GenotypeNumLikelihoodsCache,
    //
    // There are two objects here because we are lazy in creating both representations
    // for this object: a vector of log10 Probs and the PL phred-scaled string.  Supports
    // having one set during initializating, and dynamic creation of the other, if needed
    //
    log10_likelihoods: Vec<f64>,
    // likelihoods_as_string_pls: String,
}

impl GenotypeLikelihoods {
    pub const MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED: usize = 50;

    pub fn calc_num_likelihoods(num_alleles: usize, ploidy: usize) -> i64 {
        binomial((num_alleles + ploidy - 1) as u64, ploidy as u64) as i64
    }

    pub fn new() -> GenotypeLikelihoods {
        GenotypeLikelihoods {
            num_likelihood_cache: GenotypeNumLikelihoodsCache::new_empty(),
            log10_likelihoods: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.log10_likelihoods.len() == 0
    }

    pub fn from_log10_likelihoods(log10_likelihoods: Vec<f64>) -> GenotypeLikelihoods {
        GenotypeLikelihoods {
            num_likelihood_cache: GenotypeNumLikelihoodsCache::new_empty(),
            log10_likelihoods,
        }
    }

    pub fn get_gq_log10_from_likelihoods(&self, i_of_chosen_genotype: usize) -> f64 {
        let mut qual = std::f64::NEG_INFINITY;

        for i in (0..self.log10_likelihoods.len()).into_iter() {
            if i == i_of_chosen_genotype {
                continue;
            } else if self.log10_likelihoods[i] >= qual {
                qual = self.log10_likelihoods[i]
            }
        }

        // qual contains now max(likelihoods[k]) for all k != bestGTguess
        qual = self.log10_likelihoods[i_of_chosen_genotype] - qual;

        if qual < 0.0 {
            let normalized =
                MathUtils::normalize_from_log10(&self.log10_likelihoods[..], false, false);
            let chosen_genotype = normalized[i_of_chosen_genotype];

            (1.0 - chosen_genotype).log10()
        } else {
            -1. * qual
        }
    }

    pub fn get_gq_log10_from_likelihoods_on_the_fly(
        i_of_chosen_genotype: usize,
        log10_likelihoods: &[f64],
    ) -> f64 {
        let mut qual = std::f64::NEG_INFINITY;

        for i in (0..log10_likelihoods.len()).into_iter() {
            if i == i_of_chosen_genotype {
                continue;
            } else if log10_likelihoods[i] >= qual {
                qual = log10_likelihoods[i]
            }
        }

        // qual contains now max(likelihoods[k]) for all k != bestGTguess
        qual = log10_likelihoods[i_of_chosen_genotype] - qual;

        if qual < 0. {
            let normalized = MathUtils::normalize_from_log10(&log10_likelihoods[..], false, false);
            let chosen_genotype = normalized[i_of_chosen_genotype];

            (1.0 - chosen_genotype).log10()
        } else {
            -1.0 * qual
        }
    }

    /**
     * Compute how many likelihood elements are associated with the given number of alleles
     * Equivalent to asking in how many ways N non-negative integers can add up to P is S(N,P)
     * where P = ploidy (number of chromosomes) and N = total # of alleles.
     * Each chromosome can be in one single state (0,...,N-1) and there are P of them.
     * Naive solution would be to store N*P likelihoods, but this is not necessary because we can't distinguish chromosome states, but rather
     * only total number of alt allele counts in all chromosomes.
     *
     * For example, S(3,2) = 6: For alleles A,B,C, on a diploid organism we have six possible genotypes:
     * AA,AB,BB,AC,BC,CC.
     * Another way of expressing is with vector (#of A alleles, # of B alleles, # of C alleles)
     * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
     * In general, for P=2 (regular biallelic), then S(N,2) = N*(N+1)/2
     *
     * Note this method caches the value for most common num Allele / ploidy combinations for efficiency
     *
     * For non-cached values, the result is calculated via a call to calcNumLikelihoods,
     * which uses the Apache Commons CombinatoricsUtils class
     * using the formula (numAlleles + ploidy - 1) choose ploidy
     *
     *   @param  numAlleles      Number of alleles (including ref)
     *   @param  ploidy          Ploidy, or number of chromosomes in set
     *   @return    Number of likelihood elements we need to hold.
     */
    pub fn num_likelihoods(&mut self, num_alleles: i64, ploidy: i64) -> i64 {
        self.num_likelihood_cache.get(num_alleles, ploidy)
    }

    pub fn get_as_vector(&mut self) -> &mut Vec<f64> {
        &mut self.log10_likelihoods
    }

    pub fn get_likelihoods(&self) -> &Vec<f64> {
        &self.log10_likelihoods
    }

    pub fn len(&self) -> usize {
        self.log10_likelihoods.len()
    }
}

#[derive(Debug, Clone)]
pub struct GenotypeNumLikelihoodsCache {
    static_cache: Vec<Vec<i64>>,
    dynamic_cache: HashMap<CacheKey, i64>,
}

impl GenotypeNumLikelihoodsCache {
    const DEFAULT_N_ALLELES: i64 = 5;
    const DEFAULT_PLOIDY: i64 = 5;

    pub fn new_empty() -> GenotypeNumLikelihoodsCache {
        GenotypeNumLikelihoodsCache {
            static_cache: Vec::new(),
            dynamic_cache: HashMap::new(),
        }
    }

    pub fn default_values(&mut self) {
        self.add(
            GenotypeNumLikelihoodsCache::DEFAULT_N_ALLELES,
            GenotypeNumLikelihoodsCache::DEFAULT_PLOIDY,
        )
    }

    pub fn add(&mut self, num_alleles: i64, ploidy: i64) {
        self.static_cache = vec![vec![0; ploidy as usize]; num_alleles as usize];

        self.fill_cache();
    }

    fn fill_cache(&mut self) {
        for num_alleles in (0..self.static_cache.len()).into_iter() {
            for ploidy in (0..self.static_cache[num_alleles].len()).into_iter() {
                self.static_cache[num_alleles][ploidy] =
                    GenotypeLikelihoods::calc_num_likelihoods(num_alleles + 1, ploidy + 1);
            }
        }
    }

    fn put(&mut self, num_alleles: i64, ploidy: i64, num_likelihoods: i64) {
        self.dynamic_cache
            .insert(CacheKey::build(num_alleles, ploidy), num_likelihoods);
    }

    /**
     * Returns the number of likelihoods for the specified allele count and ploidy
     * Values not present yet in the cache will be calculated and cached when get is called
     * @param numAlleles
     * @param ploidy
     * @return number of likelihoods
     */
    fn get(&mut self, num_alleles: i64, ploidy: i64) -> i64 {
        if num_alleles <= 0 || ploidy <= 0 {
            panic!(
                "num_alleles and ploidy must both exceed 0, but they are numAlleles {}, ploidy {}",
                num_alleles, ploidy
            )
        }
        if (num_alleles as usize) < self.static_cache.len()
            && (ploidy as usize) < self.static_cache[num_alleles as usize].len()
        {
            self.static_cache[(num_alleles - 1) as usize][(ploidy - 1) as usize]
        } else {
            let cached_value = self
                .dynamic_cache
                .get(&CacheKey::build(num_alleles, ploidy));
            match cached_value {
                Some(value) => return *value,
                None => {
                    let new_value = GenotypeLikelihoods::calc_num_likelihoods(
                        num_alleles as usize,
                        ploidy as usize,
                    );
                    self.put(num_alleles, ploidy, new_value);

                    return new_value;
                }
            }
        }
    }
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct CacheKey {
    num_alleles: i64,
    ploidy: i64,
}

impl CacheKey {
    pub fn build(num_alleles: i64, ploidy: i64) -> CacheKey {
        CacheKey {
            num_alleles,
            ploidy,
        }
    }

    pub fn equals(&self, other_cache: &CacheKey) -> bool {
        self == other_cache
    }

    pub fn hash_code(&self) -> i64 {
        self.num_alleles * 31 + self.ploidy
    }
}
