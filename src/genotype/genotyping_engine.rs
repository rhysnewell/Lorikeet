use hashlink::LinkedHashMap;
use ordered_float::OrderedFloat;
use std::collections::{BinaryHeap, HashSet};

use crate::assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use crate::genotype::genotype_builder::{
    AttributeObject, Genotype, GenotypeAssignmentMethod, GenotypesContext,
};
use crate::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use crate::genotype::genotype_likelihoods::GenotypeLikelihoods;
use crate::genotype::genotype_prior_calculator::GenotypePriorCalculator;
use crate::model::allele_frequency_calculator::AlleleFrequencyCalculator;
use crate::model::allele_frequency_calculator_result::AFCalculationResult;
use crate::model::allele_subsetting_utils::AlleleSubsettingUtils;
use crate::model::byte_array_allele::{Allele, ByteArrayAllele};
use crate::model::variant_context::VariantContext;
use crate::model::variants::{Filter, NON_REF_ALLELE};
use crate::utils::math_utils::MathUtils;
use crate::utils::quality_utils::QualityUtils;
use crate::utils::simple_interval::{Locatable, SimpleInterval};
use crate::utils::vcf_constants::*;

#[derive(Debug, Clone)]
pub struct GenotypingEngine {
    pub(crate) allele_frequency_calculator: AlleleFrequencyCalculator,
    // pub(crate) number_of_genomes: usize,
    pub(crate) samples: Vec<String>,
    // the top of the queue is the upstream deletion that ends first
    // note that we can get away with ordering just by the end and not the contig as long as we preserve the invariant
    // that everything in this queue belongs to the same contig
    upstream_deletions_loc: BinaryHeap<SimpleInterval>,
    do_allele_specific_calcs: bool,
    genotype_assignment_method: GenotypeAssignmentMethod,
    use_posterior_probabilities_to_calculate_qual: bool,
    annotate_number_of_alleles_discovered: bool,
}

impl GenotypingEngine {
    pub const TOO_LONG_PL: usize = 100000;

    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param doAlleleSpecificCalcs Whether the AS_QUAL key should be calculated and added to newly genotyped variants.
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} is {@code null}.
     */
    pub fn make(
        args: &clap::ArgMatches,
        samples: Vec<String>,
        do_allele_specific_calcs: bool,
        _sample_ploidy: usize,
    ) -> GenotypingEngine {
        GenotypingEngine {
            allele_frequency_calculator: AlleleFrequencyCalculator::make_calculator(args),
            // number_of_genomes: samples.len() * sample_ploidy,
            samples,
            do_allele_specific_calcs,
            upstream_deletions_loc: BinaryHeap::new(),
            genotype_assignment_method: GenotypeAssignmentMethod::from_args(args),
            use_posterior_probabilities_to_calculate_qual: args
                .get_flag("use-posteriors-to-calculate-qual"),
            annotate_number_of_alleles_discovered: args
                .get_flag("annotate-with-num-discovered-alleles"),
        }
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GL's that is shared across genotypers (namely UG and HC).
     *
     * Completes a variant context with genotype calls and associated annotations given the genotype likelihoods and
     * the model that need to be applied.
     *
     * @param self                               Input variant context to complete.
     * @return                                   VC with assigned genotypes
     */
    pub fn calculate_genotypes<'b>(
        &'b mut self,
        mut vc: VariantContext,
        ploidy: usize,
        gpc: &'b GenotypePriorCalculator,
        given_alleles: &'b Vec<VariantContext>,
        stand_min_conf: f64,
    ) -> Option<VariantContext> {
        if vc.has_too_many_alternative_alleles() || vc.get_n_samples() == 0 {
            return None;
        }

        let mut reduced_vc: VariantContext;
        if VariantContext::MAX_ALTERNATE_ALLELES < (vc.get_alternate_alleles().len()) {
            let alleles_to_keep = AlleleSubsettingUtils::calculate_most_likely_alleles(
                &vc,
                ploidy,
                VariantContext::MAX_ALTERNATE_ALLELES,
            );

            let reduced_genotypes = if alleles_to_keep.len() == 1 {
                VariantContext::subset_to_ref_only(&vc, ploidy)
            } else {
                AlleleSubsettingUtils::subset_alleles(
                    vc.get_genotypes().clone(),
                    ploidy,
                    vc.get_alleles(),
                    &alleles_to_keep,
                    gpc,
                    &GenotypeAssignmentMethod::SetToNoCall,
                    vc.get_dp(),
                    true,
                )
            };
            reduced_vc = vc.clone();
            reduced_vc.alleles = alleles_to_keep;
            reduced_vc.genotypes = reduced_genotypes;
        } else {
            reduced_vc = vc.clone();
        }

        //Calculate the expected total length of the PL arrays for this VC to warn the user in the case that they will be exceptionally large
        let max_pl_length =
            GenotypeLikelihoods::calc_num_likelihoods(reduced_vc.get_n_alleles(), ploidy);
        // println!("Max pl length {}", max_pl_length);

        if max_pl_length >= GenotypingEngine::TOO_LONG_PL as i64 {
            warn!("Length of PL arrays for this Variant Context \
            (position: {}, allles: {}, ploidy: {}) is likely to reach {} so processing may take a long time",
                  vc.loc.get_start(), vc.get_n_alleles(), vc.genotypes.get_max_ploidy(ploidy), max_pl_length)
        }

        let af_result = self
            .allele_frequency_calculator
            .calculate(reduced_vc, ploidy);
        debug!("AFresult {:?}", &af_result);
        debug!("VC {:?}", &vc);
        let given_alleles_empty = given_alleles.is_empty();
        let output_alternative_alleles =
            self.calculate_output_allele_subset(&af_result, &vc, given_alleles, stand_min_conf);
        debug!("Ouput alt alleles {:?}", &output_alternative_alleles);
        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        let log10_confidence = if !output_alternative_alleles.site_is_monomorphic {
            af_result.log10_prob_only_ref_allele_exists() + 0.0
        } else {
            af_result.log10_prob_variant_present() + 0.0
        };

        // Add 0.0 removes -0.0 occurrences.
        // if log10_confidence.is_finite() {
        debug!(
            "pos {} log10 {} phred {}",
            vc.loc.start,
            log10_confidence,
            (-10.0 * log10_confidence) + 0.0
        );
        // }
        let phred_scaled_confidence = (-10.0 * log10_confidence) + 0.0;

        // return a None if we don't pass the confidence cutoff or the most likely allele frequency is zero
        // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
        if !GenotypingEngine::passes_emit_threshold(
            phred_scaled_confidence,
            stand_min_conf,
            output_alternative_alleles.site_is_monomorphic,
        ) && GenotypingEngine::no_alleles_or_first_allele_is_not_non_ref(
            &output_alternative_alleles.alleles,
        ) && given_alleles_empty
        {
            debug!(
                "Did not pass emit threshold {} {}",
                phred_scaled_confidence, stand_min_conf
            );
            return None;
        }

        // start constructing the resulting VC
        let output_alleles = output_alternative_alleles.output_alleles(vc.get_reference());

        self.record_deletions(&vc, &output_alleles);

        let mut builder = VariantContext::build(
            vc.loc.get_contig(),
            vc.loc.get_start(),
            vc.loc.get_end(),
            output_alleles.clone(),
        );

        builder.log10_p_error(log10_confidence);
        if !GenotypingEngine::passes_call_threshold(phred_scaled_confidence, stand_min_conf) {
            builder.filter(Filter::from(&(*LOW_QUAL_FILTER_NAME)))
        }

        // create the genotypes
        //TODO: omit subsetting if output alleles is not a proper subset of vc.getAlleles
        let mut genotypes = if builder.alleles.len() == 1 {
            VariantContext::subset_to_ref_only(&vc, ploidy)
        } else {
            AlleleSubsettingUtils::subset_alleles(
                vc.get_genotypes().clone(),
                ploidy,
                vc.get_alleles(),
                &output_alleles,
                gpc,
                &self.genotype_assignment_method,
                vc.get_dp(),
                true,
            )
        };

        if self.use_posterior_probabilities_to_calculate_qual
            && GenotypingEngine::has_posteriors(&genotypes)
        {
            let log10_no_variant_posterior =
                GenotypingEngine::phred_no_variant_posterior_probability(
                    &output_alleles,
                    &mut genotypes,
                ) * (-0.1);

            let qual_update = if !output_alternative_alleles.site_is_monomorphic {
                // TODO: add || self.annotate_all_sites_with_pls
                log10_no_variant_posterior + 0.0
            } else {
                MathUtils::log10_one_minus_pow10(log10_no_variant_posterior) + 0.0
            };

            if !qual_update.is_nan() {
                builder.log10_p_error(qual_update)
            }
        }

        // calculating strand bias involves overwriting data structures, so we do it last
        let attributes = self.compose_call_attributes(
            &vc,
            &output_alternative_alleles.alternative_allele_mle_counts(),
            &af_result,
            &output_alternative_alleles.output_alleles(vc.get_reference()),
            &genotypes,
        );

        builder.attributes(attributes);
        builder.genotypes = genotypes;

        return Some(builder);
    }

    fn phred_no_variant_posterior_probability(
        alleles: &Vec<ByteArrayAllele>,
        gc: &mut GenotypesContext,
    ) -> f64 {
        gc.genotypes_mut()
            .iter_mut()
            .map(|mut gt| GenotypingEngine::extract_p_no_alt(alleles, &mut gt))
            .filter(|d| !d.is_nan())
            .fold(std::f64::NAN, |a, b| {
                if a.is_nan() {
                    b
                } else if b.is_nan() {
                    a
                } else {
                    a + b
                }
            })
    }

    fn extract_p_no_alt(alleles: &Vec<ByteArrayAllele>, gt: &mut Genotype) -> f64 {
        let gp_array = gt.get_attribute(&*GENOTYPE_POSTERIORS_KEY);

        match gp_array {
            Some(array) => {
                GenotypingEngine::extract_p_no_alt_with_posteriors(alleles, &gt, array.to_ref_vec())
            }
            None => std::f64::NAN,
        }
    }

    fn extract_p_no_alt_with_posteriors(
        alleles: &Vec<ByteArrayAllele>,
        gt: &Genotype,
        posteriors: &[f64],
    ) -> f64 {
        if !alleles
            .iter()
            .any(|allele| VCFConstants::is_spanning_deletion(allele))
        {
            let reducer: f64 = std::cmp::max(
                OrderedFloat(0.0),
                OrderedFloat(QualityUtils::phred_sum(posteriors)),
            )
            .into();
            return posteriors[0] - reducer;
        } else {
            // here we need to get indices of genotypes composed of REF and * alleles
            let ploidy = gt.ploidy;
            let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles.len());
            let span_del_index = alleles
                .iter()
                .position(|allele| VCFConstants::is_spanning_deletion(allele))
                .unwrap();
            // allele counts are in the GenotypeLikelihoodCalculator format of {ref index, ref count, span del index, span del count}
            let mut non_variant_log10_posteriors = (0..ploidy)
                .into_iter()
                .map(|n| {
                    gl_calc.allele_counts_to_index(&[0, ploidy - n, span_del_index, n]);
                    posteriors[n]
                })
                .collect::<Vec<f64>>();

            // when the only alt allele is the spanning deletion the probability that the site is non-variant
            // may be so close to 1 that finite precision error in log10SumLog10 (called by phredSum) yields a positive value,
            // which is bogus.  Thus we cap it at 0. See AlleleFrequencyCalculator.
            return (std::cmp::max(
                OrderedFloat(0.0),
                OrderedFloat(QualityUtils::phred_sum(&mut non_variant_log10_posteriors)),
            ) - std::cmp::max(
                OrderedFloat(0.0),
                OrderedFloat(QualityUtils::phred_sum(posteriors)),
            ))
            .into();
        }
    }

    /**
     *  Record emitted deletions in order to remove downstream spanning deletion alleles that are not covered by any emitted deletion.
     *  In addition to recording new deletions, this method culls previously-recorded deletions that end before the current variant
     *  context.  This assumes that variants are traversed in order.
     *
     * @param vc                VariantContext, potentially multiallelic and potentially containing one or more deletion alleles
     * @param emittedAlleles    The subset of the variant's alt alleles that are actually emitted
     */
    fn record_deletions(&mut self, vc: &VariantContext, emitted_alleles: &Vec<ByteArrayAllele>) {
        while !self.upstream_deletions_loc.is_empty() {
            let condition_met;
            match self.upstream_deletions_loc.peek() {
                Some(upstream) => {
                    if !upstream.contigs_match(&vc.loc) || upstream.get_end() < vc.loc.get_start() {
                        condition_met = true
                    } else {
                        break;
                    }
                }
                None => break,
            }
            if condition_met {
                self.upstream_deletions_loc.pop();
            }
        }

        for allele in emitted_alleles.iter() {
            let deletion_size = if !allele.is_symbolic {
                vc.get_reference().length() - allele.length()
            } else {
                0
            };

            if deletion_size > 0 {
                let genome_loc = SimpleInterval::new(
                    vc.loc.get_contig(),
                    vc.loc.get_start(),
                    vc.loc.get_start() + deletion_size,
                );
                self.upstream_deletions_loc.push(genome_loc);
            }
        }
    }

    fn no_alleles_or_first_allele_is_not_non_ref(alt_alleles: &Vec<ByteArrayAllele>) -> bool {
        alt_alleles.is_empty() || alt_alleles[0] != *NON_REF_ALLELE
    }

    pub fn passes_emit_threshold(conf: f64, min_conf: f64, best_guess_is_ref: bool) -> bool {
        !best_guess_is_ref && GenotypingEngine::passes_call_threshold(conf, min_conf)
    }

    pub fn passes_call_threshold(conf: f64, min_conf: f64) -> bool {
        conf >= min_conf
    }

    /**
     * Provided the exact mode computations it returns the appropriate subset of alleles that progress to genotyping.
     * @param afCalculationResult the allele fraction calculation result.
     * @param vc the variant context
     * @return information about the alternative allele subsetting {@code null}.
     */
    fn calculate_output_allele_subset<'b>(
        &'b self,
        af_calculation_result: &'b AFCalculationResult,
        vc: &'b VariantContext,
        given_alleles: &'b Vec<VariantContext>,
        stand_min_conf: f64,
    ) -> OutputAlleleSubset {
        let mut output_alleles = Vec::new();
        let mut mle_counts = Vec::new();

        let mut site_is_monomorphic = true;
        let alleles = af_calculation_result.get_alleles_used_in_genotyping();
        let alternative_allele_count = alleles.len() - 1;
        let mut _reference_size;

        let forced_alleles: HashSet<&ByteArrayAllele> =
            AssemblyBasedCallerUtils::get_alleles_consistent_with_given_alleles(given_alleles, vc);
        for allele in alleles.iter() {
            if allele.is_reference() {
                _reference_size = allele.length();
            } else {
                // we want to keep the NON_REF symbolic allele but only in the absence of a non-symbolic allele, e.g.
                // if we combined a ref / NON_REF gVCF with a ref / alt gVCF
                let is_non_ref_which_is_lone_alt_allele =
                    alternative_allele_count == 1 && allele.eq(&*NON_REF_ALLELE);

                debug!(
                    "is non ref which is lone alt_allele {}",
                    is_non_ref_which_is_lone_alt_allele
                );
                let is_plausible = af_calculation_result.passes_threshold(allele, stand_min_conf);
                debug!(
                    "plausible {} {} {}",
                    is_plausible,
                    af_calculation_result.get_log10_posterior_of_allele_absent(allele),
                    QualityUtils::qual_to_error_prob_log10(stand_min_conf as u8)
                );

                //it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion allele
                let is_spurious_spanning_deletion = VCFConstants::is_spanning_deletion(allele)
                    || self.is_vc_covered_by_deletion(vc);
                debug!("is spurious spanning del {}", is_spurious_spanning_deletion);

                let to_output = (is_plausible
                    || is_non_ref_which_is_lone_alt_allele
                    || forced_alleles.contains(&allele))
                    && !is_spurious_spanning_deletion;
                debug!("to output {}", to_output);
                site_is_monomorphic =
                    site_is_monomorphic & !(is_plausible && !is_spurious_spanning_deletion);

                if to_output {
                    output_alleles.push(allele.clone());
                    mle_counts.push(af_calculation_result.get_allele_count_at_mle(allele));
                }
            }
        }

        return OutputAlleleSubset::new(output_alleles, mle_counts, site_is_monomorphic);
    }

    fn is_vc_covered_by_deletion(&self, vc: &VariantContext) -> bool {
        // note: the code below seems like it's duplicating Locatable.overlaps, but here if the upstream deletion
        // has the same start as the vc we don't want to count it
        return !self.upstream_deletions_loc.is_empty()
            && self.upstream_deletions_loc.iter().any(|loc| {
                loc.get_contig() == vc.loc.tid
                    && loc.get_start() < vc.loc.start
                    && vc.loc.start <= loc.get_end()
            });
    }

    fn has_posteriors(gc: &GenotypesContext) -> bool {
        gc.genotypes()
            .iter()
            .any(|genotype| genotype.has_attribute(&(*GENOTYPE_POSTERIORS_KEY).to_string()))
    }

    fn compose_call_attributes<'b>(
        &self,
        vc: &'b VariantContext,
        allele_counts_of_mle: &'b Vec<i64>,
        af_result: &'b AFCalculationResult,
        all_alleles_to_use: &'b Vec<ByteArrayAllele>,
        genotypes: &'b GenotypesContext,
    ) -> LinkedHashMap<String, AttributeObject> {
        let mut attributes = LinkedHashMap::new();

        // add the MLE AC and AF annotations
        if !allele_counts_of_mle.is_empty() {
            attributes.insert(
                MLE_ALLELE_COUNT_KEY.to_string(),
                AttributeObject::Vecf64(
                    allele_counts_of_mle
                        .iter()
                        .map(|count| *count as f64)
                        .collect::<Vec<f64>>(),
                ),
            );
            let mle_frequencies = GenotypingEngine::calculate_mle_allele_frequencies(
                &allele_counts_of_mle,
                genotypes,
            );
            attributes.insert(
                MLE_ALLELE_FREQUENCY_KEY.to_string(),
                AttributeObject::Vecf64(mle_frequencies),
            );
        }

        if self.do_allele_specific_calcs {
            let mut per_allele_quals = Vec::new();
            //Per-allele quals are not calculated for biallelic sites
            if af_result.get_alleles_used_in_genotyping().len() > 2 {
                for a in all_alleles_to_use.iter() {
                    if !a.is_reference() {
                        //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                        per_allele_quals
                            .push(af_result.get_log10_posterior_of_allele_absent(a) * -10.0);
                    }
                }
            } else {
                //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                per_allele_quals.push(af_result.log10_prob_only_ref_allele_exists() * -10.0);
            }
            attributes.insert(
                AS_QUAL_KEY.to_string(),
                AttributeObject::Vecf64(per_allele_quals),
            );
        }

        if self.annotate_number_of_alleles_discovered {
            attributes.insert(
                NUMBER_OF_DISCOVERED_ALLELES_KEY.to_string(),
                AttributeObject::Vecf64(vec![vc.get_alternate_alleles().len() as f64]),
            );
        }

        return attributes;
    }

    fn calculate_mle_allele_frequencies(
        allele_counts_of_mle: &[i64],
        genotypes: &GenotypesContext,
    ) -> Vec<f64> {
        let an = genotypes
            .genotypes()
            .iter()
            .flat_map(|g| g.alleles.iter())
            .filter(|a| a.is_called())
            .count();

        return allele_counts_of_mle
            .iter()
            .map(|ac| {
                std::cmp::min(OrderedFloat(0.0), OrderedFloat((*ac as f64) / (an as f64))).into()
            })
            .collect::<Vec<f64>>();
    }
}

#[derive(Debug)]
struct OutputAlleleSubset {
    alleles: Vec<ByteArrayAllele>,
    site_is_monomorphic: bool,
    mle_counts: Vec<i64>,
}

impl OutputAlleleSubset {
    pub fn new(
        alleles: Vec<ByteArrayAllele>,
        mle_counts: Vec<i64>,
        site_is_monomorphic: bool,
    ) -> OutputAlleleSubset {
        OutputAlleleSubset {
            alleles,
            mle_counts,
            site_is_monomorphic,
        }
    }

    pub fn output_alleles(&self, reference_allele: &ByteArrayAllele) -> Vec<ByteArrayAllele> {
        let mut result = self.alleles.clone();
        result.insert(0, reference_allele.clone());
        return result;
    }

    pub fn alternative_allele_mle_counts(&self) -> &Vec<i64> {
        &self.mle_counts
    }
}
