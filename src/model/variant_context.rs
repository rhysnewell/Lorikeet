use model::variants::Allele;
use itertools::Itertools;
use genotype::genotype_builder::{GenotypesContext, GenotypeAssignmentMethod, Genotype};
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use utils::math_utils::MathUtils;
use rayon::prelude::*;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use utils::vcf_constants::VCFConstants;
use utils::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use utils::simple_interval::SimpleInterval;
use std::collections::{HashMap, HashSet};
use num::traits::Float;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct VariantContext {
    pub loc: SimpleInterval,
    // variant alleles
    pub alleles: Vec<Allele>,
    // per sample likelihoods
    pub genotypes: GenotypesContext,
    pub source: String,
    pub log10_p_error: f64,
    pub filters: HashSet<String>,
    pub attributes: HashMap<String, Vec<f64>>
}

impl VariantContext {
    pub const MAX_ALTERNATE_ALLELES: usize = 180;
    pub const SUM_GL_THRESH_NOCALL: f64 = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    pub fn build(tid: usize, start: usize, end: usize, alleles: Vec<Allele>) -> VariantContext {
        VariantContext {
            loc: SimpleInterval::new(tid, start, end),
            alleles,
            genotypes: GenotypesContext::empty(),
            source: "".to_string(),
            log10_p_error: 0.0,
            filters: HashSet::new(),
            attributes: HashMap::new(),
        }
    }

    pub fn build_from_vc(vc: &VariantContext) -> VariantContext {
        VariantContext {
            loc: vc.loc,
            alleles: Vec::new(),
            genotypes: GenotypesContext::empty(),
            source: "".to_string(),
            log10_p_error: 0.0,
            filters: HashSet::new(),
            attributes: HashMap::new(),
        }
    }

    pub fn attributes(&mut self, attributes: HashMap<String, Vec<f64>>) {
        self.attributes.merge(attributes)
    }

    pub fn filter(&mut self, filter: String) {
        self.filters.insert(filter)
    }

    pub fn get_log10_p_error(&self) -> f64 {
        self.log10_p_error
    }

    pub fn get_phred_scaled_qual(&self) -> f64 {
        (self.log10_p_error * (-10.0)) + 0.0
    }

    pub fn log10_p_error(&mut self, error: f64) {
        self.log10_p_error = error
    }

    pub fn add_alleles(&mut self, alleles: Vec<Allele>) {
        self.alleles = alleles
    }

    pub fn add_source(&mut self, source: String) {
        self.source = source
    }

    pub fn get_n_samples(&self) -> usize {
        self.genotypes.len()
    }

    pub fn get_n_alleles(&self) -> usize {
        self.alleles.len()
    }

    pub fn add_genotypes(
        &mut self,
        genotypes: Vec<Genotype>
    ) {
        for genotype in genotypes.iter() {
            if genotype.pl.len() != self.alleles.len() {
                panic!(
                    "Number of likelihoods does not match number of alleles at position {} on tid {}",
                    self.start, self.tid
                )
            }
        }
        self.genotypes = GenotypesContext::new(genotypes);
    }

    pub fn has_non_ref_allele(&self) -> bool {
        self.alleles.iter().any(|a| !a.is_reference())
    }

    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get into trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    pub fn has_too_many_alternative_alleles(&self) -> bool {
        if self.get_n_alleles() <= GenotypeLikelihoods::MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED {
            false
        } else {
            true
        }
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested {@link GenotypeAssignmentMethod}
     *
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param genotypeLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles with respect to which the likelihoods are defined
     */
    pub fn make_genotype_call<T: Float + Copy>(
        ploidy: usize,
        gb: &mut Genotype,
        assignment_method: &GenotypeAssignmentMethod,
        genotype_likelihoods: Option<&[T]>,
        alleles_to_use: &Vec<Allele>,
        original_gt: &Vec<Allele>,
        gpc: &GenotypePriorCalculator
    ) {
        match genotype_likelihoods {
            Some(genotype_likelihoods) => {
                match assignment_method {
                    &GenotypeAssignmentMethod::SetToNoCall => {
                        gb.no_call_alleles(ploidy);
                        gb.no_qg();
                    },
                    &GenotypeAssignmentMethod::UsePLsToAssign => {
                        if !VariantContext::is_informative(genotype_likelihoods) {
                            gb.no_call_alleles(ploidy);
                            gb.no_qg();
                        } else {
                            let max_likelihood_index = MathUtils::max_element_index(genotype_likelihoods, 0, genotype_likelihoods.len());
                            let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles_to_use.len());
                            let allele_counts = gl_calc.genotype_allele_counts_at(max_likelihood_index);
                            let final_alleles = allele_counts.as_allele_list(alleles_to_use);
                            if final_alleles.contains(&Allele::NON_REF_ALLELE) {
                                gb.no_call_alleles(ploidy);
                                gb.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![0.0 as T; genotype_likelihoods.len()]);
                            } else {
                                gb.alleles = final_alleles;
                            }
                            let num_alt_alleles = alleles_to_use.len() - 1;
                            if num_alt_alleles > 0 {
                                gb.log10_p_error(
                                    GenotypeLikelihoods::get_gq_log10_from_likelihoods_on_the_fly(max_likelihood_index, genotype_likelihoods)
                                );
                            }
                        }
                    },
                    &GenotypeAssignmentMethod::SetToNoCallNoAnnotations => {
                        gb.no_call_alleles(ploidy);
                        gb.no_annotations();
                    },
                    &GenotypeAssignmentMethod::BestMatchToOriginal => {
                        let mut best = Vec::new();
                        let refr = alleles_to_use[0];
                        for original_allele in original_gt.iter() {
                            if alleles_to_use.contains(original_allele) || original_allele.is_no_call() {
                                best.push(original_allele.clone())
                            } else {
                                best.push(refr)
                            }
                        }
                        gb.alleles = best;
                    },
                    &GenotypeAssignmentMethod::UsePosteriorProbabilities => {
                        // Calculate posteriors.
                        let gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles_to_use.len());
                        let log10_priors = gpc.get_log10_priors(gl_calc, alleles_to_use);
                        let log10_posteriors = MathUtils::ebe_add(&log10_priors, &genotype_likelihoods);
                        let normalized_log10_posteriors = MathUtils::scale_log_space_array_for_numeric_stability(
                            &log10_posteriors
                        );
                        // Update GP and PG annotations:
                        gb.attribute(
                            VCFConstants::GENOTYPE_POSTERIOR_KEY,
                            normalized_log10_posteriors.par_iter().map(|v| {
                                if v == 0.0 { // the reason for the == 0.0 is to avoid a signed 0 output "-0.0"
                                    0.0
                                } else {
                                    v * -10.
                                }
                            }).collect_vec()
                        );

                        gb.attribute(
                            VCFConstants::GENOTYPE_PRIOR_KEY,
                            log10_priors.par_iter().map(|v| {
                                if v == 0.0 {
                                    0.0
                                } else {
                                    v * -10.
                                }
                            }).collect_vec()
                        );
                        // Set the GQ accordingly
                        let max_posterior_index = MathUtils::max_element_index(&log10_posteriors, 0, log10_posteriors.len());
                        if alleles_to_use.len() > 0 {
                            gb.log10_p_error(VariantContext::get_gq_log10_from_posteriors(
                                max_posterior_index, &normalized_log10_posteriors,
                            ));
                        }
                        // Finally we update the genotype alleles.
                        gb.alleles(gl_calc.genotype_allele_counts_at(max_posterior_index).as_allele_list(alleles_to_use));
                    },
                    _ => panic!("Unknown GenotypeAssignmentMethod")
                }
            },
            _ => {
                gb.no_call_alleles(ploidy);
                gb.no_qg();
            }
        }
    }

    fn get_gq_log10_from_posteriors<T: Float + Copy>(best_genotype_index: usize, log10_posteriors: &[T]) -> T {
        match log10_posteriors.len() {
            0 | 1 => {
                1.0 as T
            },
            2 => {
                if best_genotype_index == 0 {
                    log10_posteriors[1]
                } else {
                    log10_posteriors[0]
                }
            },
            3 => {
                std::cmp::min(0. as T, MathUtils::log10_sum_log10_two_values(
                    log10_posteriors[
                            if best_genotype_index == 0 {
                                2
                            } else {
                                best_genotype_index - 1
                            }
                        ],
                    log10_posteriors[
                            if best_genotype_index == 2 {
                                0
                            } else {
                                best_genotype_index + 1
                            }
                        ]
                ))
            },
            _ => {
                if best_genotype_index == 0 {
                    MathUtils::log10_sum_log10(log10_posteriors, 1, log10_posteriors.len())
                } else if best_genotype_index == (log10_posteriors.len() - 1) {
                    MathUtils::log10_sum_log10(log10_posteriors, 0, best_genotype_index);
                } else {
                    std::cmp::min(
                        0. as T,
                        MathUtils::log10_sum_log10_two_values(
                            MathUtils::log10_sum_log10(
                                log10_posteriors, 0, best_genotype_index
                            ),
                            MathUtils::log10_sum_log10(
                                log10_posteriors, best_genotype_index + 1, log10_posteriors.len()
                            )
                        )
                    )
                }
            }
        }
    }

    pub fn is_informative<T: Float + Copy>(gls: &[T]) -> bool {
        gls.par_iter().sum() < (VariantContext::SUM_GL_THRESH_NOCALL as T)
    }

    /**
     * Subset the samples in VC to reference only information with ref call alleles
     *
     * Preserves DP if present
     *
     * @param vc the variant context to subset down to
     * @param defaultPloidy defaultPloidy to use if a genotype doesn't have any alleles
     * @return a GenotypesContext
     */
    pub fn subset_to_ref_only(vc: &mut VariantContext, default_ploidy: usize) -> GenotypesContext {
        if default_ploidy < 1 {
            panic!("default_ploidy must be >= 1, got {}", default_ploidy)
        } else {
            let old_gts = vc.get_genotypes();
            if old_gts.is_empty() { return old_gts }

            let mut new_gts = GenotypesContext::create(old_gts.size());

            let ref_allele = vc.get_reference();
            let diploid_ref_alleles = vec![ref_allele; 2];

            // create the new genotypes
            for g in old_gts.genotypes() {
                let g_ploidy = if g.get_ploidy() == 0 { default_ploidy } else { g.get_ploidy() };
                let ref_alleles = if g_ploidy == 2 { diploid_ref_alleles.clone() } else { vec![ref_allele; g_ploidy] };

                let genotype = Genotype::build_from_alleles(ref_alleles);
                new_gts.add(genotype);
            }

            return new_gts
        }
    }

    pub fn get_genotypes(&self) -> GenotypesContext {
        self.genotypes.clone()
    }

    pub fn get_alleles(&self) -> &Vec<Allele> {
        self.alleles
    }

    pub fn get_reference(&self) -> &Allele {
        self.alleles[0]
    }

    pub fn get_dp(&self) -> i64 {
        self.genotypes.get_dp()
    }

    pub fn get_alternate_alleles(&self) -> &Vec<Allele> {
        self.variants[1..]
    }
}