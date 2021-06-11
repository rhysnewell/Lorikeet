use model::allele_list::AlleleList;
use model::variants::Allele;
use genotype::genotype_builder::Genotype;
use model::variant_context::VariantContext;
use clap::ArgMatches;
use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use utils::math_utils::MathUtils;
use ordered_float::OrderedFloat;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;

pub struct AlleleFrequencyCalculator {
    pub ref_pseudo_count: f64,
    pub snp_pseudo_count: f64,
    pub indel_pseudo_count: f64,
    pub default_ploidy: usize,
}

impl AlleleFrequencyCalculator {
    const GL_CALCS: GenotypeLikeliHoodCalculators = GenotypeLikelihoodCalculators::build_empty();
    const THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE: f64 = 0.1;
    const HOM_REF_GENOTYPE_INDEX: usize = 0;

    pub fn new(
        ref_pseudo_count: f64,
        snp_pseudo_count: f64,
        indel_pseudo_count: f64,
        default_ploidy: usize
    ) -> AlleleFrequencyCalculator {
        AlleleFrequencyCalculator {
            ref_pseudo_count,
            snp_pseudo_count,
            indel_pseudo_count,
            default_ploidy
        }
    }

    pub fn make_calculator(args: &ArgMatches) -> AlleleFrequencyCalculator {
        let snp_het = args.value_of("snp-heterozygosity").unwrap().parse::<f64>().unwrap();
        let ind_het = args.value_of("indel-heterozygosity").unwrap().parse::<f64>().unwrap();
        let het_std = args.value_of("heterozygosity-stdev").unwrap().parse::<f64>().unwrap();
        let ploidy: usize = m.value_of("ploidy").unwrap().parse().unwrap();

        let ref_pseudo_count = snp_het / (het_std.powf(2.));
        let snp_pseudo_count = snp_het * ref_pseudo_count;
        let indel_pseudo_count = ind_het * ref_pseudo_count;

        AlleleFrequencyCalculator::new(
            ref_pseudo_count,
            snp_pseudo_count,
            indel_pseudo_count,
            ploidy
        )
    }

    fn log10_normalized_genotype_posteriors(
        &mut self,
        g: &mut Genotype,
        gl_calc: &mut GenotypeLikelihoodCalculator,
        log10_allele_frequencies: Vec<f64>,
    ) -> Vec<f64> {
        let mut log10_likelihoods = g.get_likelihoods();
        let log10_posteriors =
            (0..gl_calc.genotype_count as usize).iter().map(|genotype_index| {
                let mut gac = gl_calc.genotype_allele_counts_at(genotype_index);
                let result = gac.log10_combination_count()
                    + log10_likelihoods[genotype_index]
                    + gac.sum_over_allele_indices_and_counts(
                    |index: usize, count: usize| {
                        (count as f64) * &log10_allele_frequencies[index]
                    }
                );
                result
            }).collect::<Vec<f64>>();

        return MathUtils::normalize_log10(log10_posteriors, true);
    }

    /**
     * Calculate the posterior probability that a single biallelic genotype is non-ref
     *
     * The nth genotype (n runs from 0 to the sample ploidy, inclusive) contains n copies of the alt allele
     * @param log10GenotypeLikelihoods
     */
    pub fn calculate_single_sample_biallelic_non_ref_posterior(
        &self,
        log10_genotype_likelihoods: &Vec<OrderedFloat<f64>>,
        return_zero_if_ref_is_max: bool,
    ) -> f64 {
        if return_zero_if_ref_is_max
            && log10_genotype_likelihoods.iter()
            .position(|&item| item == *(log10_genotype_likelihoods.iter().max().unwrap())).unwrap() == 0
            && (log10_genotype_likelihoods[0] != OrderedFloat(0.5) && log10_genotype_likelihoods.len() == 2){
            return 0.
        }

        let ploidy = log10_genotype_likelihoods.len() - 1;

        let log10_unnormalized_posteriors = (0..ploidy + 1)
            .iter()
            .map(|n| {
                log10_genotype_likelihoods[n]
                    + MathUtils::log10_binomial_coeffecient(ploidy as f64, n as f64)
                    + MathUtils::log_to_log10(
                    (n as f64 + self.snp_pseudo_count).log_gamma()
                        + ((ploidy - n) as f64 + self.ref_pseudo_count).log_gamma()
                )
            });

        return if return_zero_if_ref_is_max
            && MathUtils::max_element_index(log10_unnormalized_posteriors, 0, log10_unnormalized_posteriors.len()) == 0 {
            0.0
        } else {
            1 - MathUtils::normalize_log10(log10_unnormalized_posteriors, false)[0]
        }
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @return result (for programming convenience)
     */
    pub fn calculate(vc: VariantContext, )
}