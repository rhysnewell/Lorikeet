use bio::stats::LogProb;
use ordered_float::{NotNan, OrderedFloat};
use model::genotype_builder::{Genotype, GenotypeLikelihoodCalculator};
use model::variants::Allele;
use itertools::Itertools;


pub struct VariantContext {
    // contig id
    pub tid: usize,
    // context start
    pub start: usize,
    // context end
    pub end: usize,
    // reference allele
    pub refr: Allele,
    // variant alleles
    pub variants: Vec<Allele>,
    // per sample likelihoods
    pub genotypes: Vec<Genotype>,
}

impl VariantContext {
    pub fn build(tid: usize, start: usize, end: usize, alleles: Vec<Alleles>) -> VariantContext {
        VariantContext {
            tid,
            start,
            end,
            refr: alleles[0].clone(),
            variants: alleles[1..].iter().cloned().collect_vec(),
            genotypes: Vec::new(),
        }
    }

    pub fn get_n_samples(&self) -> usize {
        self.genotypes.len()
    }

    pub fn get_n_alleles(&self) -> usize {
        self.variants.len() + 1
    }

    pub fn add_genotypes(
        &mut self,
        genotypes: Vec<Genotype>
    ) {
        for genotype in genotypes.iter() {
            if genotype.likelihoods.len() != self.variants.len() + 1 {
                panic!(format!(
                    "Number of likelihoods does not match number of alleles at position {} on tid {}",
                    self.start, self.tid
                )
                )
            }
        }
        self.genot = likelihoods
    }

    pub fn has_non_ref_allele(&self) -> bool {
        for alt_allele in self.variants.iter() {
            if alt_allele != self.refr {
                return true
            }
        }

        false
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
        if self.get_n_alleles() <= 180 {
            false
        } else {
            true
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
    pub fn calculate_genotypes(
        &mut self,
        ploidy: f32,
        gpc: Option<f32>,
        given_alleles: Vec<VariantContext>,

    ) -> Option<VariantContext> {
        if self.has_too_many_alternative_alleles() || self.get_n_samples() == 0 {
            return None
        }

        if 180 < self.variants.len() {
            let alleles_to_keep = self.calculate_most_likely_alleles(ploidy, 180);
        }

        return Some(self)
    }

    pub fn get_alleles(&self) -> Vec<Allele> {
        let mut all_alleles = vec![self.refr.clone()];
        all_alleles.extend(self.variants.clone());

        return all_alleles
    }

    /**
     * Returns the new set of alleles to use based on a likelihood score: alleles' scores are the sum of their counts in
     * sample genotypes, weighted by the confidence in the genotype calls.
     *
     * In the case of ties, the alleles will be chosen from lowest index to highest index.
     *
     * @param vc target variant context.
     * @param numAltAllelesToKeep number of alt alleles to keep.
     * @return the list of alleles to keep, including the reference and {@link Allele#NON_REF_ALLELE} if present
     *
     */
    pub fn calculate_most_likely_alleles(
        &mut self,
        ploidy: f32,
        num_alt_alleles_to_keep: usize,
    ) -> Vec<Allele> {
        let has_symbolic_non_ref = self.has_non_ref_allele();
        let number_of_alleles_that_arent_proper_alt = if has_symbolic_non_ref { 2 } else { 1 };
        let number_of_proper_alts = self.get_n_alleles() - number_of_alleles_that_arent_proper_alt;

        if num_alt_alleles_to_keep >= number_of_proper_alts {
            return self.get_alleles()
        }

        let likelihood_sums = self.calculate_likelihood_sums(ploidy);

        return self.f
    }

    /** the likelihood sum for an alt allele is the sum over all samples whose likeliest genotype contains that allele of
     * the GL difference between the most likely genotype and the hom ref genotype
     *
     * Since GLs are log likelihoods, this quantity is thus
     * SUM_{samples whose likeliest genotype contains this alt allele} log(likelihood alt / likelihood hom ref)
     */
    pub fn calculate_likelihood_sums(
        &self,
        ploidy: f32,
    ) -> Vec<f32> {
        let mut likelihood_sums = vec![0.; self.get_n_alleles()];

        for genotype in self.genotypes.iter() {
            let index_of_most_likely_genotype = genotype.likelihoods.iter()
                .position(|&item| item == *(genotype.likelihoods.iter().max().unwrap())).unwrap();

            let difference_between_best_and_ref =
                genotype.likelihoods[index_of_most_likely_genotype] - genotype.likelihoods[0];
            let ploidy = if genotype.get_ploidy() > 0 { genotype.get_ploidy() } else { ploidy };

            let allele_counts = GenotypeLikelihoodCalculator::get_instance(
                ploidy,
                self.get_n_alleles() as i32
            )
        }
    }
}