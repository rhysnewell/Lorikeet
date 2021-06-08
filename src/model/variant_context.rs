use bio::stats::LogProb;
use ordered_float::{NotNan, OrderedFloat};
use model::genotype_builder::{Genotype, GenotypesContext};
use model::variants::Allele;
use itertools::Itertools;
use model::allele_subsetting_utils::AlleleSubsettingUtils;


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
    pub genotypes: GenotypesContext,
}

impl VariantContext {
    pub const MAX_ALTERNATE_ALLELES: usize = 180;

    pub fn build(tid: usize, start: usize, end: usize, alleles: Vec<Alleles>) -> VariantContext {
        VariantContext {
            tid,
            start,
            end,
            refr: alleles[0].clone(),
            variants: alleles[1..].iter().cloned().collect_vec(),
            genotypes: GenotypesContext::empty(),
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
                panic!(
                    "Number of likelihoods does not match number of alleles at position {} on tid {}",
                    self.start, self.tid
                )
            }
        }
        self.genotypes = GenotypesContext::new(genotypes);
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
        mut vc: VariantContext,
        ploidy: usize,
        gpc: Option<f64>,
        given_alleles: Vec<VariantContext>,

    ) -> Option<VariantContext> {
        if vc.has_too_many_alternative_alleles() || vc.get_n_samples() == 0 {
            return None
        }

        if VariantContext::MAX_ALTERNATE_ALLELES < vc.variants.len() {
            let alleles_to_keep = AlleleSubsettingUtils::calculate_most_likely_alleles(
                &mut vc,
                ploidy,
                VariantContext::MAX_ALTERNATE_ALLELES
            );

            let reduced_genotypes = if alleles_to_keep.len() == 1 {
                VariantContext::subset_to_ref_only(vc, ploidy)
            } else {

            }
        }

        return Some(vc)
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
    pub fn subset_to_ref_only(vc: VariantContext, default_ploidy: usize) -> GentypesContext {
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

    pub fn get_alleles(&self) -> Vec<Allele> {
        let mut all_alleles = vec![self.refr.clone()];
        all_alleles.extend(self.variants.clone());

        return all_alleles
    }

    pub fn get_reference(&self) -> Allele {
        self.refr.clone()
    }
}