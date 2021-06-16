use model::variants::Allele;
use itertools::Itertools;
use model::variant_context::VariantContext;
use rayon::prelude::*;
use std::collections::{HashSet, BTreeMap};
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use genotype::genotype_builder::{GenotypesContext, GenotypeAssignmentMethod, Genotype};
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use model::allele_list::{AlleleList, AlleleListPermutation};
use utils::math_utils::MathUtils;
use ordered_float::OrderedFloat;
use genotype::genotype_likelihoods::GenotypeLikelihoods;

pub struct AlleleSubsettingUtils {}

impl AlleleSubsettingUtils {
    const GL_CALCS: GenotypeLikelihoodCalculators = GenotypeLikelihoodCalculators::build_empty();

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
        vc: &mut VariantContext,
        ploidy: usize,
        num_alt_alleles_to_keep: usize,
    ) -> Vec<Allele> {
        let has_symbolic_non_ref = vc.has_non_ref_allele();
        let number_of_alleles_that_arent_proper_alt = if has_symbolic_non_ref { 2 } else { 1 };
        let number_of_proper_alts = vc.get_n_alleles() - number_of_alleles_that_arent_proper_alt;

        if num_alt_alleles_to_keep >= number_of_proper_alts {
            return vc.get_alleles()
        }

        let likelihood_sums = AlleleSubsettingUtils::calculate_likelihood_sums(vc, ploidy);

        return AlleleSubsettingUtils::filter_to_max_number_of_alt_alleles_based_on_scores(
            num_alt_alleles_to_keep,
            vc.get_alleles(),
            likelihood_sums,
        )
    }

    /**
     * @param alleles a list of alleles including the reference and possible the NON_REF
     * @return a list of the best proper alt alleles based on the likelihood sums, keeping the reference allele and {@link Allele#NonRefAllele}
     * the list will include no more than {@code numAltAllelesToKeep + 2} alleles and will maintain the order of the original alleles in {@code vc}
     *
     */
    pub fn filter_to_max_number_of_alt_alleles_based_on_scores(
        num_alt_alleles_to_keep: usize,
        alleles: Vec<Allele>,
        likelihood_sums: Vec<f64>
    ) -> Vec<Allele> {
        let non_ref_alt_allele_index = alleles.par_iter().position(|&a| a == Allele::NON_REF_ALLELE).unwrap();
        let num_alleles = alleles.len();

        let proper_alt_indexes_to_keep =
            (1..num_alleles).into_par_iter()
            .filter(|n| n != non_ref_alt_allele_index)
            .par_sort_by_key(|n| likelihood_sums[n])
            .reverse()
            .take(num_alt_alleles_to_keep)
            .collect::<HashSet<usize>>();

        let result = (0..num_alleles).into_par_iter()
            .filter(|i| {i == 0 || i == non_ref_alt_allele_index || proper_alt_indexes_to_keep.contains(i)})
            .map(|i| {alleles[i].clone()}).collect::<Vec<Allele>>();

        return result
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
    ) -> Vec<f64> {
        let mut likelihood_sums = vec![0.; vc.get_n_alleles()];

        for genotype in vc.genotypes.iter() {
            let index_of_most_likely_genotype = genotype.likelihoods.iter()
                .position(|&item| item == *(genotype.likelihoods.iter().max().unwrap())).unwrap();

            let difference_between_best_and_ref =
                genotype.likelihoods[index_of_most_likely_genotype] - genotype.likelihoods[0];
            let ploidy = if genotype.get_ploidy() > 0 { genotype.get_ploidy() } else { ploidy };

            let allele_counts = GenotypeLikelihoodCalculators::get_instance(
                ploidy,
                vc.get_n_alleles()
            ).genotype_allele_counts_at(
                index_of_most_likely_genotype
            ).allele_counts_by_index(
                vc.get_n_alleles() - 1
            );

            for allele in 1..allele_counts.len() {
                if allele_counts[allele] > 0 {
                    likelihood_sums[allele] += difference_between_best_and_ref;
                }
            }
        }

        return likelihood_sums
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
    pub fn subset_alleles(
        original_gs: &GenotypesContext,
        default_ploidy: usize,
        original_alleles: &Vec<Allele>,
        alleles_to_keep: &Vec<Allele>,
        gpc: &GenotypePriorCalculator,
        assignment_method: &GenotypeAssignmentMethod,
        depth: i64,
    ) -> GenotypesContext {
        if alleles_to_keep.len() == 0 {
            panic!("alleles_to_keep is empty")
        }
        if !alleles_to_keep[0].is_reference() {
            panic!("First allele must be reference allele")
        }

        let mut new_gts = GenotypesContext::create(original_gs.size());
        let allele_permutation = AlleleList::new(&original_alleles).permuation(AlleleList::new(&alleles_to_keep));
        let mut subsetted_likelihood_indices_by_ploidy = BTreeMap::new();

        for g in original_gs.genotypes().iter_mut() {
            let ploidy = if g.get_ploidy() > 0 { g.get_ploidy() } else { default_ploidy };

            let mut subsetted_likelihoods_indices = subsetted_likelihood_indices_by_ploidy.entry(ploidy)
                .or_insert(
                    AlleleSubsettingUtils::subsetted_pl_indices(
                        ploidy,
                        &original_alleles,
                        &alleles_to_keep,
                        g,
                    )
                );

            let expected_num_likelihoods = g.num_likelihoods(original_alleles.len() as i64, ploidy as i64);
            let mut new_likelihoods: Option<Vec<OrderedFloat<f64>>>;
            let mut new_log10_gq = -1;

            let mut original_likelihoods = g.get_likelihoods();
            new_likelihoods =
                if original_likelihoods.len() == expected_num_likelihoods {
                    let mut subsetted_likeihoods =
                        subsetted_likelihoods_indices.iter().map(|i| original_likelihoods[i]).collect_vec();
                    MathUtils::scale_log_space_array_for_numeric_stability(&mut subsetted_likeihoods);
                    Some(subsetted_likeihoods)
                } else {
                    None
                };
            match new_likelihoods {
                Some(&new_likelihoods) => {
                    let pl_index = MathUtils::max_element_index(new_likelihoods, 0, new_likelihoods.len());
                    new_log10_gq = g.get_likelihoods().get_gq_log10_from_likelihoods(pl_index);
                },
                _ => {}
            }

            let use_new_likelihoods = match new_likelihoods {
                Some(&new_likelihoods) => {
                    if depth != 0 || VariantContext::is_informative(new_likelihoods) {
                        true
                    } else {
                        false
                    }
                },
                _ => false,
            };

            let mut gb = Genotype::from(g.clone());
            if use_new_likelihoods {
                gb.pl(GenotypeLikelihoods::from_log10_likelihoods(new_likelihoods.unwrap()));
                // Don't worry about extended attributes here, we won't use them
            }

            VariantContext::make_genotype_call(
                g.ploidy,
                &mut gb,
                &assignment_method,
                new_likelihoods,
                &alleles_to_keep,
                &g.alleles,
                &gpc
            );

            if g.has_ad() {
                let old_ad = g.get_ad();
                let new_ad = (0..alleles_to_keep.len()).into_par_iter().map(|n| {
                    old_ad[allele_permutation.from_index(n)]
                }).collect_vec();

                let non_ref_index = alleles_to_keep.par_iter().position(|&p| p == Allele::NON_REF_ALLELE);

                match non_ref_index {
                    Some(index) => {
                        new_ad[index] = 0;
                    },
                    _ => {
                        // pass
                    },
                }
                gb.ad = new_ad.clone();
            }

            new_gts.add(gb)
        }
        return new_gts
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
        original_alleles: &Vec<Allele>,
        new_alleles: &Vec<Allele>,
        g: &mut Genotype
    ) -> Vec<usize> {
        let mut result = vec![0; g.num_likelihoods(new_alleles.len() as i64, ploidy as i64) as usize];

        let allele_permutation = AlleleList::new(original_alleles).permuation(AlleleList::new(new_alleles));

        let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, original_alleles.len());

        for old_pl_index in (0..gl_calc).into_iter() {
            let old_allele_counts = gl_calc.genotype_allele_counts_at(old_pl_index);
            let contains_only_new_alleles = (0..old_allele_counts.distinct_allele_count())
                .map(|i| {
                    old_allele_counts.allele_index_at(i)
                }).all(|b| allele_permutation.is_kept(b));

            if contains_only_new_alleles {
                // make an array in the format described in {@link GenotypeAlleleCounts}:
                // [(new) index of first allele, count of first allele, (new) index of second allele, count of second allele. . .]
                let new_allele_counts = (0..new_alleles.len()).into_iter()
                    .flat_map(|new_allele_index| {
                        new_allele_index..old_allele_counts.allele_count_for(
                            allele_permutation.from_index(new_allele_index)
                        )
                    }).collect_vec();

                let new_pl_index = gl_calc.allele_counts_to_index(&new_allele_counts);

                result[new_pl_index] = old_pl_index
            }
        }
        return result
    }

}