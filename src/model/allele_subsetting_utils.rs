use genotype::genotype_builder::{Genotype, GenotypeAssignmentMethod, GenotypesContext};
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use itertools::Itertools;
use model::allele_list::{AlleleList, AlleleListPermutation};
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::VariantContext;
use model::variants::NON_REF_ALLELE;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use utils::math_utils::MathUtils;
use utils::vcf_constants::*;

pub struct AlleleSubsettingUtils {}

impl AlleleSubsettingUtils {
    /**
     * Returns the new set of alleles to use based on a likelihood score: alleles' scores are the sum of their counts in
     * sample genotypes, weighted by the confidence in the genotype calls.
     *
     * In the case of ties, the alleles will be chosen from lowest index to highest index.
     *
     * @param vc target variant context.
     * @param numAltAllelesToKeep number of alt alleles to keep.
     * @return the list of alleles to keep, including the reference and {@link Allele#NonRefAllele} if present
     *
     */
    pub fn calculate_most_likely_alleles(
        vc: &VariantContext,
        ploidy: usize,
        num_alt_alleles_to_keep: usize,
    ) -> Vec<ByteArrayAllele> {
        let has_symbolic_non_ref = vc.has_non_ref_allele();
        let number_of_alleles_that_arent_proper_alt = if has_symbolic_non_ref { 2 } else { 1 };
        let number_of_proper_alts = vc.get_n_alleles() - number_of_alleles_that_arent_proper_alt;
        let all_hom_ref_data = vc.get_genotypes().genotypes().iter().all(|g| {
            if g.has_likelihoods() {
                g.pl[0] == 0
            } else {
                false
            }
        });

        if num_alt_alleles_to_keep >= number_of_proper_alts {
            return vc.get_alleles().to_vec();
        }

        let likelihood_sums =
            AlleleSubsettingUtils::calculate_likelihood_sums(vc, ploidy, all_hom_ref_data);

        return AlleleSubsettingUtils::filter_to_max_number_of_alt_alleles_based_on_scores(
            num_alt_alleles_to_keep,
            vc.get_alleles().to_vec(),
            likelihood_sums,
        );
    }

    /**
     * @param alleles a list of alleles including the reference and possible the NON_REF
     * @return a list of the best proper alt alleles based on the likelihood sums, keeping the reference allele and {@link Allele#NonRefAllele}
     * the list will include no more than {@code numAltAllelesToKeep + 2} alleles and will maintain the order of the original alleles in {@code vc}
     *
     */
    pub fn filter_to_max_number_of_alt_alleles_based_on_scores(
        num_alt_alleles_to_keep: usize,
        alleles: Vec<ByteArrayAllele>,
        likelihood_sums: Vec<f64>,
    ) -> Vec<ByteArrayAllele> {
        let non_ref_alt_allele_index = alleles.iter().position(|a| a == &*NON_REF_ALLELE);
        let num_alleles = alleles.len();
        let mut indices = (1..num_alleles).collect_vec();
        indices.sort_by_key(|n| OrderedFloat(likelihood_sums[*n]));

        let mut proper_alt_indexes_to_keep = indices
            .into_iter()
            .filter(|n| Some(*n) != non_ref_alt_allele_index)
            .collect::<Vec<usize>>();
        proper_alt_indexes_to_keep.sort_by_key(|n| -OrderedFloat(likelihood_sums[*n]));
        let proper_alt_indexes_to_keep = proper_alt_indexes_to_keep
            .into_iter()
            .take(num_alt_alleles_to_keep)
            .collect::<HashSet<usize>>();

        let result = alleles
            .into_iter()
            .take(num_alleles)
            .enumerate()
            .filter(|(i, a)| {
                *i == 0
                    || Some(*i) == non_ref_alt_allele_index
                    || proper_alt_indexes_to_keep.contains(i)
            })
            .map(|(_, a)| a)
            .collect::<Vec<ByteArrayAllele>>();

        return result;
    }

    /** the likelihood sum for an alt allele is the sum over all samples whose likeliest genotype contains that allele of
     * the GL difference between the most likely genotype and the hom ref genotype
     *
     * Since GLs are log likelihoods, this quantity is thus
     * SUM_{samples whose likeliest genotype contains this alt allele} log(likelihood alt / likelihood hom ref)
     */
    pub fn calculate_likelihood_sums(
        vc: &VariantContext,
        ploidy: usize,
        all_hom_ref_data: bool,
    ) -> Vec<f64> {
        let mut likelihood_sums = vec![0.; vc.get_n_alleles()];

        for genotype in vc.get_genotypes().genotypes().iter() {
            if genotype.pl.is_empty() {
                continue;
            }
            let gls = genotype.get_likelihoods();
            let index_of_most_likely_genotype = MathUtils::max_element_index(
                gls.get_likelihoods(),
                if all_hom_ref_data { 1 } else { 0 },
                genotype.get_likelihoods().len(),
            );

            let difference_between_best_and_ref =
                (gls.get_likelihoods()[index_of_most_likely_genotype] - gls.get_likelihoods()[0])
                    .abs();
            let ploidy = if genotype.get_ploidy() > 0 {
                genotype.get_ploidy()
            } else {
                ploidy
            };

            let allele_counts =
                GenotypeLikelihoodCalculators::get_instance(ploidy, vc.get_n_alleles())
                    .genotype_allele_counts_at(index_of_most_likely_genotype)
                    .allele_counts_by_index(vc.get_n_alleles() - 1);

            for allele in 1..allele_counts.len() {
                if allele_counts[allele] > 0 {
                    likelihood_sums[allele] += difference_between_best_and_ref;
                }
            }
        }

        return likelihood_sums;
    }

    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     *
     * Will reorder subsetted alleles according to the ordering provided by the list allelesToKeep
     *
     * @param originalGs               the original GenotypesContext
     * @param originalAlleles          the original alleles
     * @param allelesToKeep            the subset of alleles to use with the new Genotypes
     * @param assignmentMethod         assignment strategy for the (subsetted) PLs
     * @param depth                    the original variant DP or 0 if there was no DP
     * @return                         a new non-null GenotypesContext
     */
    pub fn subset_alleles<'b>(
        mut original_gs: GenotypesContext,
        default_ploidy: usize,
        original_alleles: &'b Vec<ByteArrayAllele>,
        alleles_to_keep: &'b Vec<ByteArrayAllele>,
        gpc: &'b GenotypePriorCalculator,
        assignment_method: &'b GenotypeAssignmentMethod,
        depth: i64,
        emit_empty_pls: bool,
    ) -> GenotypesContext {
        if alleles_to_keep.len() == 0 {
            panic!("alleles_to_keep is empty")
        }
        if !alleles_to_keep[0].is_reference() {
            panic!("First allele must be reference allele")
        }

        let mut new_gts = GenotypesContext::create(original_gs.size());
        let allele_permutation =
            AlleleList::new(&original_alleles).permutation(AlleleList::new(&alleles_to_keep));
        let mut subsetted_likelihood_indices_by_ploidy = BTreeMap::new();

        for g in original_gs.genotypes_mut().iter_mut() {
            let ploidy = if g.get_ploidy() > 0 {
                g.get_ploidy()
            } else {
                default_ploidy
            };

            let subsetted_likelihoods_indices = subsetted_likelihood_indices_by_ploidy
                .entry(ploidy)
                .or_insert(AlleleSubsettingUtils::subsetted_pl_indices(
                    ploidy,
                    original_alleles,
                    alleles_to_keep,
                    g,
                ));

            let expected_num_likelihoods =
                g.num_likelihoods(original_alleles.len() as i64, ploidy as i64);

            let mut new_likelihoods: Option<Vec<f64>> = None;
            let mut new_log10_gq = std::f64::NEG_INFINITY;

            if g.has_likelihoods() {
                let original_likelihoods = g.get_likelihoods().get_as_vector();
                new_likelihoods = if original_likelihoods.len() == expected_num_likelihoods {
                    let mut subsetted_likeihoods = subsetted_likelihoods_indices
                        .iter()
                        .map(|i| original_likelihoods[*i])
                        .collect::<Vec<f64>>();
                    MathUtils::scale_log_space_array_for_numeric_stability(
                        &mut subsetted_likeihoods,
                    );
                    Some(subsetted_likeihoods)
                } else {
                    None
                };
            } else if g.has_gq() {
                new_log10_gq = -0.1 * g.gq as f64
            }

            match &new_likelihoods {
                Some(new_likelihoods) => {
                    let pl_index =
                        MathUtils::max_element_index(new_likelihoods, 0, new_likelihoods.len());
                    new_log10_gq = GenotypeLikelihoods::get_gq_log10_from_likelihoods(
                        pl_index,
                        new_likelihoods,
                    );
                }
                _ => {}
            }

            let use_new_likelihoods = match &new_likelihoods {
                Some(new_likelihoods) => {
                    if depth != 0 || VariantContext::is_informative(new_likelihoods) {
                        true
                    } else {
                        false
                    }
                }
                _ => false,
            } || emit_empty_pls;

            let mut gb = Genotype::from(g.clone());

            if new_log10_gq != std::f64::NEG_INFINITY {
                gb.log10_p_error(new_log10_gq);
            };
            if use_new_likelihoods {
                gb.pl(GenotypeLikelihoods::from_log10_likelihoods(
                    new_likelihoods.clone().unwrap(),
                )
                .as_pls());
            }

            gb.attributes.remove(&*PHRED_SCALED_POSTERIORS_KEY);
            gb.attributes.remove(&*GENOTYPE_POSTERIORS_KEY);
            gb.attributes.remove(&*GENOTYPE_POSTERIORS_KEY);

            VariantContext::make_genotype_call(
                g.ploidy,
                &mut gb,
                &assignment_method,
                new_likelihoods,
                alleles_to_keep,
                &g.alleles,
                &gpc,
            );

            if g.has_ad() {
                let old_ad = g.get_ad();
                let mut new_ad = (0..alleles_to_keep.len())
                    .into_iter()
                    .map(|n| old_ad[allele_permutation.from_index(n)])
                    .collect::<Vec<i64>>();
                let non_ref_index = alleles_to_keep.iter().position(|p| p == &*NON_REF_ALLELE);

                match non_ref_index {
                    Some(index) => {
                        new_ad[index] = 0;
                    }
                    _ => {
                        // pass
                    }
                }
                gb.ad = new_ad;
            }

            new_gts.add(gb)
        }
        return new_gts;
    }

    /**
     * Given a list of original alleles and a subset of new alleles to retain, find the array of old PL indices that correspond
     * to new PL indices i.e. result[7] = old PL index of genotype containing same alleles as the new genotype with PL index 7.
     *
     * This method is written in terms f indices rather than subsetting PLs directly in order to produce output that can be
     * recycled from sample to sample, provided that the ploidy is the same.
     *
     * @param ploidy                Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param newAlleles            New alleles -- must be a subset of {@code originalAlleles}
     * @return                      old PL indices of new genotypes
     */
    pub fn subsetted_pl_indices(
        ploidy: usize,
        original_alleles: &Vec<ByteArrayAllele>,
        new_alleles: &Vec<ByteArrayAllele>,
        g: &mut Genotype,
    ) -> Vec<usize> {
        let mut result =
            vec![0; g.num_likelihoods(new_alleles.len() as i64, ploidy as i64) as usize];

        let allele_permutation =
            AlleleList::new(original_alleles).permutation(AlleleList::new(new_alleles));

        let mut gl_calc =
            GenotypeLikelihoodCalculators::get_instance(ploidy, original_alleles.len());

        for old_pl_index in (0..gl_calc.genotype_count as usize).into_iter() {
            let old_allele_counts = gl_calc.genotype_allele_counts_at(old_pl_index);
            let contains_only_new_alleles = (0..old_allele_counts.distinct_allele_count())
                .map(|i| old_allele_counts.allele_index_at(i))
                .all(|b| allele_permutation.is_kept(b));

            if contains_only_new_alleles {
                // make an array in the format described in {@link GenotypeAlleleCounts}:
                // [(new) index of first allele, count of first allele, (new) index of second allele, count of second allele. . .]
                let new_allele_counts = (0..new_alleles.len())
                    .into_iter()
                    .flat_map(|new_allele_index| {
                        vec![
                            new_allele_index,
                            old_allele_counts
                                .allele_count_for(allele_permutation.from_index(new_allele_index)),
                        ]
                    })
                    .collect_vec();

                let new_pl_index = gl_calc.allele_counts_to_index(&new_allele_counts);

                result[new_pl_index] = old_pl_index
            }
        }
        return result;
    }
}
