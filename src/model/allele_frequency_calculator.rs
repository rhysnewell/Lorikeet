use clap::ArgMatches;
use genotype::genotype_builder::{Genotype, GenotypeType};
use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use model::allele_frequency_calculator_result::AFCalculationResult;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::VariantContext;
use model::variants::SPAN_DEL_ALLELE;
use num::traits::Float;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, Mutex};
use utils::dirichlet::Dirichlet;
use utils::math_utils::MathUtils;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use std::rc::Rc;
use std::cell::RefCell;

lazy_static!{
    //from the genotype likelihoods equations assuming the SNP ref conf model with no mismatches
    //PL[2] = GQ; scaleFactor = PL[3]/GQ ~ -10 * DP * log10(P_error) / (-10 * DP * log10(1/ploidy)) where BASE_QUALITY = -10 * log10(P_error)
    static ref PLOIDY_2_HOM_VAR_SCALE_FACTOR: i64 = (AlleleFrequencyCalculator::TYPICAL_BASE_QUALITY as f64/-10.0/(0.5.log10())).round() as i64;
}

#[derive(Debug, Clone)]
pub struct AlleleFrequencyCalculator {
    pub ref_pseudo_count: f64,
    pub snp_pseudo_count: f64,
    pub indel_pseudo_count: f64,
    pub default_ploidy: usize,
}

impl AlleleFrequencyCalculator {
    const THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE: f64 = 0.1;
    const HOM_REF_GENOTYPE_INDEX: usize = 0;
    const TYPICAL_BASE_QUALITY: i64 = 30;

    pub fn new(
        ref_pseudo_count: f64,
        snp_pseudo_count: f64,
        indel_pseudo_count: f64,
        default_ploidy: usize,
    ) -> AlleleFrequencyCalculator {
        AlleleFrequencyCalculator {
            ref_pseudo_count,
            snp_pseudo_count,
            indel_pseudo_count,
            default_ploidy,
        }
    }

    pub fn make_calculator(args: &ArgMatches) -> AlleleFrequencyCalculator {
        let snp_het = args
            .value_of("snp-heterozygosity")
            .unwrap()
            .parse::<f64>()
            .unwrap();
        let ind_het = args
            .value_of("indel-heterozygosity")
            .unwrap()
            .parse::<f64>()
            .unwrap();
        let het_std = args
            .value_of("heterozygosity-stdev")
            .unwrap()
            .parse::<f64>()
            .unwrap();
        let ploidy: usize = args.value_of("ploidy").unwrap().parse().unwrap();

        let ref_pseudo_count = snp_het / (het_std.powf(2.));
        let snp_pseudo_count = snp_het * ref_pseudo_count;
        let indel_pseudo_count = ind_het * ref_pseudo_count;

        AlleleFrequencyCalculator::new(
            ref_pseudo_count,
            snp_pseudo_count,
            indel_pseudo_count,
            ploidy,
        )
    }

    fn log10_normalized_genotype_posteriors(
        &mut self,
        g: &Genotype,
        gl_calc: &mut GenotypeLikelihoodCalculator,
        log10_allele_frequencies: &mut [f64],
    ) -> Vec<f64> {
        let mut log10_likelihoods;

        if g.has_likelihoods() {
            log10_likelihoods = g.get_likelihoods().get_as_vector();
        } else if g.genotype_type == Some(GenotypeType::NoCall) || g.genotype_type == Some(GenotypeType::HomRef) {
            if g.ploidy != 2 {
                panic!("Likelihoods are required to calculate posteriors for hom-refs with ploidy != 2");
            };

            if g.has_gq() {
                //for a hom-ref, as long as we have GQ we can make a very accurate QUAL calculation
                // since the hom-var likelihood should make a minuscule contribution

                let mut per_sample_indexes_of_relevant_alleles = vec![1; log10_allele_frequencies.len()];
                per_sample_indexes_of_relevant_alleles[0] = 0; // ref still maps to ref

                let gq = g.gq;
                let ploidy = g.ploidy;

                //use these values for diploid ref/ref, ref/alt, alt/alt likelihoods
                let mut approx_likelihoods = vec![0, gq, *PLOIDY_2_HOM_VAR_SCALE_FACTOR * gq];

                //map likelihoods for any other alts to biallelic ref/alt likelihoods above
                let genotype_index_map_by_ploidy = GenotypeLikelihoodCalculators::get_instance(ploidy, log10_allele_frequencies.len()).genotype_index_map(&per_sample_indexes_of_relevant_alleles);
                let mut PLs = vec![0; genotype_index_map_by_ploidy.len()];
                for i in 0..PLs.len() {
                    PLs[i] = approx_likelihoods[genotype_index_map_by_ploidy[i]]
                }

                log10_likelihoods = GenotypeLikelihoods::from_pls(PLs).get_as_vector();
            } else {
                panic!("Genotype does not contain the GQ necessary to calculate posteriors")
            }
        } else {
            panic!("Genotype does not contain likelihoods necessary to calculate posteriors")
        }
        let log10_posteriors = (0..gl_calc.genotype_count as usize)
            .into_iter()
            .map(|genotype_index| {
                let mut gac = gl_calc.genotype_allele_counts_at(genotype_index);

                let result = gac.log10_combination_count()
                    + log10_likelihoods[genotype_index]
                    + gac.sum_over_allele_indices_and_counts(|index: usize, count: usize| {
                        (count as f64) * &log10_allele_frequencies[index]
                    });

                result
            })
            .collect::<Vec<f64>>();

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
        log10_genotype_likelihoods: &Vec<f64>,
        return_zero_if_ref_is_max: bool,
    ) -> f64 {
        if return_zero_if_ref_is_max
            && MathUtils::max_element_index(
                log10_genotype_likelihoods,
                0,
                log10_genotype_likelihoods.len(),
            ) == 0
        {
            return 0.;
        }

        let ploidy = log10_genotype_likelihoods.len() - 1;

        let log10_unnormalized_posteriors = (0..ploidy + 1)
            .into_iter()
            .map(|n| {
                log10_genotype_likelihoods[n]
                    + MathUtils::log10_binomial_coeffecient(ploidy as f64, n as f64)
                    + MathUtils::log_to_log10(
                        ln_gamma(n as f64 + self.snp_pseudo_count)
                            + ln_gamma(ploidy as f64 - n as f64 + self.ref_pseudo_count),
                    )
            })
            .collect::<Vec<f64>>();

        return if return_zero_if_ref_is_max
            && MathUtils::max_element_index(
                &log10_unnormalized_posteriors,
                0,
                log10_unnormalized_posteriors.len(),
            ) == 0
        {
            0.0
        } else {
            1.0 - MathUtils::normalize_log10(log10_unnormalized_posteriors, false)[0]
        };
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @return result (for programming convenience)
     */
    pub fn calculate(&mut self, vc: VariantContext, default_ploidy: usize) -> AFCalculationResult {
        let num_alleles = vc.get_n_alleles();

        let alleles = vc.get_alleles();
        if num_alleles <= 1 {
            panic!("Variant context has only a single reference allele, but get_log10_p_non_ref requires at least one at all {:?}", vc);
        }
        let prior_pseudo_counts = alleles
            .iter()
            .map(|a| {
                if a.is_reference() {
                    self.ref_pseudo_count
                } else if a.length() == vc.get_reference().length() {
                    self.snp_pseudo_count
                } else {
                    self.indel_pseudo_count
                }
            })
            .collect::<Vec<f64>>();

        let mut allele_counts = vec![0.0; num_alleles];
        let flat_log10_allele_frequency = -((num_alleles as f64).log10());
        let mut log10_allele_frequencies = vec![flat_log10_allele_frequency; num_alleles];

        let mut allele_counts_maximum_difference = std::f64::INFINITY;
        while allele_counts_maximum_difference
            > AlleleFrequencyCalculator::THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE
        {
            let new_allele_counts =
                self.effective_allele_counts(&vc, &mut log10_allele_frequencies);
            allele_counts_maximum_difference =
                MathUtils::ebe_subtract(&allele_counts, &new_allele_counts)
                    .into_iter()
                    .map(|x| x.abs())
                    .max_by_key(|x| OrderedFloat(*x))
                    .unwrap_or(std::f64::NAN);
            allele_counts = new_allele_counts;

            let posterior_pseudo_counts = MathUtils::ebe_add(&prior_pseudo_counts, &allele_counts);
            // first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts gives such a low
            // effective allele frequency that it overwhelms the genotype likelihood of a real variant
            // basically, we want a chance to get non-zero pseudocounts before using a prior that's biased against a variant
            log10_allele_frequencies =
                Dirichlet::new(&posterior_pseudo_counts).log10_mean_weights();
        }

        let mut log10_p_of_zero_counts_by_allele = vec![0.0; num_alleles];
        let mut log10_p_no_variant = 0.0;

        let spanning_deletion_present = alleles.iter().any(|allele| allele == &*SPAN_DEL_ALLELE);

        let mut non_variant_indices_by_ploidy = BTreeMap::new();

        // re-usable buffers of the log10 genotype posteriors of genotypes missing each allele
        let mut log10_absent_posteriors = Rc::new(RefCell::new(vec![Vec::new(); num_alleles]));

        for (i, g) in vc.get_genotypes().genotypes().iter().enumerate() {
            if !g.genotype_usable_for_af_calculation()
            {
                continue;
            }

            let ploidy = if g.get_ploidy() == 0 {
                default_ploidy
            } else {
                g.get_ploidy()
            };

            let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, num_alleles);

            let log10_genotype_posteriors = self.log10_normalized_genotype_posteriors(
                &g,
                &mut gl_calc,
                &mut log10_allele_frequencies,
            );

            if !spanning_deletion_present {
                log10_p_no_variant +=
                    log10_genotype_posteriors[AlleleFrequencyCalculator::HOM_REF_GENOTYPE_INDEX];
            } else {
                let non_variant_indices = non_variant_indices_by_ploidy
                    .entry(ploidy)
                    .or_insert_with(|| {
                        AlleleFrequencyCalculator::genotype_indices_with_only_ref_and_span_del(
                            ploidy, alleles,
                        )
                    });

                let non_variant_log10_posteriors = non_variant_indices
                    .iter()
                    .map(|n| log10_genotype_posteriors[*n])
                    .collect::<Vec<f64>>();
                // when the only alt allele is the spanning deletion the probability that the site is non-variant
                // may be so close to 1 that finite precision error in log10SumLog10 yields a positive value,
                // which is bogus.  Thus we cap it at 0.
                log10_p_no_variant += std::cmp::min(
                    OrderedFloat(0.0),
                    OrderedFloat(MathUtils::log10_sum_log10(
                        &non_variant_log10_posteriors,
                        0,
                        non_variant_log10_posteriors.len(),
                    )),
                )
                .into_inner();
            }

            // if the VC is biallelic the allele-specific qual equals the variant qual
            if num_alleles == 2 && !spanning_deletion_present {
                continue;
            }

            // for each allele, we collect the log10 probabilities of genotypes in which the allele is absent, then add (in log space)
            // to get the log10 probability that the allele is absent in this sample
            {
                let mut log10_absent_posteriors = log10_absent_posteriors.borrow_mut();
                log10_absent_posteriors
                    .iter_mut()
                    .for_each(|arr| arr.clear());
            }
            for genotype in (0..gl_calc.genotype_count as usize).into_iter() {
                let log10_genotype_posterior = log10_genotype_posteriors[genotype];
                let log10_absent_posteriors = &log10_absent_posteriors;
                gl_calc
                    .genotype_allele_counts_at(genotype)
                    .for_each_absent_allele_index(
                        |a| {
                            let mut log10_absent_posteriors =
                                log10_absent_posteriors.borrow_mut();
                            log10_absent_posteriors[a].push(log10_genotype_posterior)
                        },
                        num_alleles,
                    );
            }

            let log10_absent_posteriors = log10_absent_posteriors.borrow_mut();
            let log10_p_no_allele = log10_absent_posteriors
                .iter()
                .map(|buffer| {
                    let mut result = MathUtils::log10_sum_log10(&buffer, 0, buffer.len());
                    result = std::cmp::min(OrderedFloat(0.0), OrderedFloat(result)).into_inner();
                    result
                })
                .collect::<Vec<f64>>();

            // multiply the cumulative probabilities of alleles being absent, which is addition of logs
            MathUtils::ebe_add_in_place(&mut log10_p_of_zero_counts_by_allele, &log10_p_no_allele);
        }

        // for biallelic the allele-specific qual equals the variant qual, and we short-circuited the calculation above
        if num_alleles == 2 && !spanning_deletion_present {
            log10_p_of_zero_counts_by_allele[1] = log10_p_no_variant
        }

        let int_allele_counts: Vec<i64> = allele_counts
            .iter()
            .map(|n| n.round() as i64)
            .collect::<Vec<i64>>();

        let current_ref_index = alleles.iter().position(|a| a.is_reference()).unwrap_or(0);

        let int_alt_allele_counts: Vec<i64> = int_allele_counts
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != current_ref_index)
            .map(|(_, v)| *v)
            .collect::<Vec<i64>>();

        let log10_p_ref_by_allele: HashMap<ByteArrayAllele, f64> = alleles
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != current_ref_index)
            .map(|(i, a)| (a.clone(), log10_p_of_zero_counts_by_allele[i]))
            .collect();

        return AFCalculationResult::new(
            int_alt_allele_counts,
            alleles.to_vec(),
            log10_p_no_variant,
            log10_p_ref_by_allele,
        );
    }

    fn genotype_indices_with_only_ref_and_span_del(
        ploidy: usize,
        alleles: &Vec<ByteArrayAllele>,
    ) -> Vec<usize> {
        let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles.len());

        let spanning_deletion_present = alleles.contains(&*SPAN_DEL_ALLELE);

        if !spanning_deletion_present {
            vec![AlleleFrequencyCalculator::HOM_REF_GENOTYPE_INDEX]
        } else {
            let span_del_index = alleles
                .iter()
                .position(|allele| allele == &*SPAN_DEL_ALLELE)
                .unwrap();
            let result = (0..ploidy + 1)
                .into_iter()
                .map(|n| gl_calc.allele_counts_to_index(&[0, ploidy - n, span_del_index, n]))
                .collect::<Vec<usize>>();

            return result;
        }
    }

    /**
     * effectiveAlleleCounts[allele a] = SUM_{genotypes g} (posterior_probability(g) * num_copies of a in g), which we denote as SUM [n_g p_g]
     * for numerical stability we will do this in log space:
     * count = SUM 10^(log (n_g p_g)) = SUM 10^(log n_g + log p_g)
     * thanks to the log-sum-exp trick this lets us work with log posteriors alone
     */
    fn effective_allele_counts(
        &mut self,
        vc: &VariantContext,
        log10_allele_frequencies: &mut [f64],
    ) -> Vec<f64> {
        let num_alleles = vc.get_n_alleles();
        let mut log10_result = Rc::new(RefCell::new(vec![std::f64::NEG_INFINITY; num_alleles]));
        for (i, g) in vc.get_genotypes().genotypes().iter().enumerate() {
            if !g.genotype_usable_for_af_calculation() {
                continue;
            }
            let mut gl_calc =
                GenotypeLikelihoodCalculators::get_instance(g.get_ploidy(), num_alleles);

            let log10_genotype_posteriors = self.log10_normalized_genotype_posteriors(
                g,
                &mut gl_calc,
                log10_allele_frequencies,
            );

            for genotype_index in (0..gl_calc.genotype_count as usize).into_iter() {
                let log10_result = &log10_result;
                gl_calc
                    .genotype_allele_counts_at(genotype_index)
                    .for_each_allele_index_and_count(|allele_index: usize, count: usize| {
                        let mut log10_result = log10_result.borrow_mut();
                        log10_result[allele_index] = MathUtils::log10_sum_log10_two_values(
                            log10_result[allele_index],
                            log10_genotype_posteriors[genotype_index] + (count as f64).log10(),
                        );
                    });
            }
        }

        let mut log10_result = log10_result.borrow_mut();
        log10_result.iter_mut().for_each(|x| *x = (10.0).powf(*x));

        return log10_result.to_vec();
    }
}
