use model::byte_array_allele::{Allele, ByteArrayAllele};
use rayon::prelude::*;
use std::collections::HashMap;
use utils::math_utils::MathUtils;
use utils::quality_utils::QualityUtils;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
#[derive(Debug)]
pub struct AFCalculationResult {
    log10_posterior_of_no_variant: f64,
    log10_p_ref_by_allele: HashMap<ByteArrayAllele, f64>,
    allele_counts_of_mle: Vec<i64>,
    alleles_used_in_genotyping: Vec<ByteArrayAllele>,
}

impl AFCalculationResult {
    const EPSILON: f64 = 1.0e-10;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    pub fn new(
        allele_counts_of_mle: Vec<i64>,
        alleles_used_in_genotyping: Vec<ByteArrayAllele>,
        log10_posterior_of_no_variant: f64,
        log10_p_ref_by_allele: HashMap<ByteArrayAllele, f64>,
    ) -> AFCalculationResult {
        if !MathUtils::is_valid_log10_probability(log10_posterior_of_no_variant) {
            panic!("log10 posterior must be a valid log probability")
        }
        if alleles_used_in_genotyping.len() == 0 {
            panic!("allelesUsedInGenotyping must be non-null list of at least 1 value")
        }
        if allele_counts_of_mle.len() != alleles_used_in_genotyping.len() - 1 {
            panic!("Allele_counts_of_mle.len() does not equal alleles_used_genotyping.len() - 1")
        }
        if log10_p_ref_by_allele.len() != alleles_used_in_genotyping.len() - 1 {
            panic!("log10_p_ref_by_allele wrong size")
        }
        // if !log10_p_ref_by_allele.keys

        AFCalculationResult {
            log10_p_ref_by_allele,
            log10_posterior_of_no_variant,
            alleles_used_in_genotyping,
            allele_counts_of_mle,
        }
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * This method returns a copy of the internally-stored array.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    pub fn get_allele_counts_of_mle(&self) -> Vec<i64> {
        self.allele_counts_of_mle.clone()
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    pub fn get_allele_count_at_mle(&self, allele: &ByteArrayAllele) -> i64 {
        if allele.is_reference() {
            panic!(
                "Cannot get the alt allele index for reference allele {:?}",
                allele
            )
        }
        let index_in_all_alleles_including_ref = self
            .alleles_used_in_genotyping
            .par_iter()
            .position_first(|a| a == allele)
            .unwrap();
        let index_in_alt_alleles = index_in_all_alleles_including_ref - 1;
        return self.allele_counts_of_mle[index_in_alt_alleles];
    }

    /**
     * Get the list of alleles actually used in genotyping, which may be smaller than the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    pub fn get_alleles_used_in_genotyping(&self) -> &Vec<ByteArrayAllele> {
        &self.alleles_used_in_genotyping
    }

    pub fn log10_prob_only_ref_allele_exists(&self) -> f64 {
        self.log10_posterior_of_no_variant
    }

    pub fn log10_prob_variant_present(&self) -> f64 {
        MathUtils::log10_one_minus_pow10(self.log10_posterior_of_no_variant)
    }

    pub fn passes_threshold(
        &self,
        allele: &ByteArrayAllele,
        phred_scale_qual_threshold: f64,
    ) -> bool {
        debug!(
            "log 10 posterior {} qual thresh {}",
            (self.get_log10_posterior_of_allele_absent(allele) + AFCalculationResult::EPSILON),
            QualityUtils::qual_to_error_prob_log10(phred_scale_qual_threshold as u8)
        );
        (self.get_log10_posterior_of_allele_absent(allele) + AFCalculationResult::EPSILON)
            < QualityUtils::qual_to_error_prob_log10(phred_scale_qual_threshold as u8)
    }

    /**
     * Returns the log10 probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log10 probability that allele is not segregating at this site
     */
    pub fn get_log10_posterior_of_allele_absent(&self, allele: &ByteArrayAllele) -> f64 {
        debug!(
            "allele {:?} available {:?}",
            allele, &self.log10_p_ref_by_allele
        );
        let log10_p_non_ref = self.log10_p_ref_by_allele.get(allele).unwrap();
        return *log10_p_non_ref;
    }
}
