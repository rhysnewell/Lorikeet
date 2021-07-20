use utils::math_utils::MathUtils;
use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use model::variants::{Allele, Variant};
use rayon::prelude::*;
use ordered_float::OrderedFloat;

lazy_static! {
    static ref NUMBER_OF_ALLELE_TYPES: usize = 4;

    // A snp can go to 3 different bases (standard-nucs - 1), so we normalize SNP lks accordingly. Here is the
    // log10 constant used for that:
    static ref LOG10_SNP_NORMALIZATION_CONSTANT: f64 = (3. as f64).log10();
}

#[derive(Debug, PartialEq, Eq, Ordinalize)]
enum AlleleType {
    REF,
    SNP,
    INDEL,
    OTHER,
}

/**
 * Class to compose genotype prior probability calculators.
 *
 * <p>
 *     Contains a collection of static method to create calculators based on different
 *     assumptions and source of knowledge a priori (e.g. {@link #assumingHW(double, double) assumingHW}
 *     or {@link #givenDragstrParams(DragstrParams, int, int, double, double) givenDragstrParams}).
 * </p>
 *
 * <p>
 *     Such priors are obtained by invoking {@link #getLog10Priors(GenotypeLikelihoodCalculator, List).
 *     This method takes on the list of alleles for that variant, an a reference to the genotype likelihood calculator witch determines the ploidy.
 * </p>
 * assumptions
 */
pub struct GenotypePriorCalculator {
    het_values: Vec<f64>,
    hom_values: Vec<f64>,
    diff_values: Vec<f64>
}

impl GenotypePriorCalculator {

    fn genotype_prior_calculator(
        snp_het: f64,
        snp_hom: f64,
        indel_het: f64,
        indel_hom: f64,
        other_het: f64,
        other_hom: f64
    ) -> GenotypePriorCalculator {
        let mut het_values = vec![0.; *NUMBER_OF_ALLELE_TYPES];
        let mut hom_values = vec![0.; *NUMBER_OF_ALLELE_TYPES];

        // REF
        // by convention ref log10 likelihoods are set to 0.
        // so they are already set.

        // SNPs: normalized for all possible mutations (number of nucs (4) - 1)
        het_values[AlleleType::SNP.ordinal() as usize] = snp_het - *LOG10_SNP_NORMALIZATION_CONSTANT;
        hom_values[AlleleType::SNP.ordinal() as usize] = snp_hom - *LOG10_SNP_NORMALIZATION_CONSTANT;
        // INDELs:
        het_values[AlleleType::INDEL.ordinal() as usize] = indel_het;
        hom_values[AlleleType::INDEL.ordinal() as usize] = indel_hom;
        // Others:
        het_values[AlleleType::OTHER.ordinal() as usize] = other_het;
        hom_values[AlleleType::OTHER.ordinal() as usize] = other_hom;

        let diff_values = MathUtils::ebe_subtract(&hom_values, &het_values);

        GenotypePriorCalculator {
            het_values,
            hom_values,
            diff_values
        }
    }

    /**
     * Calculate priors based on fix heterozygosities (per event type) and het to hom-var prior ratio.
     *
     * @param log10SnpHet snp heterozygosity in log10 scale.
     * @param log10IndelHet indel heterozygosity in log10 scale.
     * @param log10OtherHet heterozygosity for other type of variants in log10 scale.
     * @param hetHomRatio ratio between the het-var and hom-var genotype priors for the same allele in linear scale.
     * @return never {@code null}.
     */
    pub fn given_het_to_hom_ratio(
        log10_snp_het: f64,
        log10_indel_het: f64,
        log10_other_het: f64,
        het_hom_ratio: f64
    ) -> GenotypePriorCalculator {
        let log10_ratio = het_hom_ratio.log10();

        GenotypePriorCalculator::genotype_prior_calculator(
            log10_snp_het, log10_snp_het - log10_ratio,
            log10_indel_het, log10_indel_het - log10_ratio,
            log10_other_het, log10_other_het - log10_ratio
        )
    }

    /**
     * Composes a calculator based on Hardy-Weinberg equilibrium so that only the het-priors
     * are need to calculate the rest.
     * @param snpHet the prior for an SNP alternative allele in log10 scale.
     * @param indelHet the prior for an INDEL alternative allele in log10 scale.
     * @return never {@code null}.
     */
    pub fn assuming_hw(
        snp_het: f64,
        indel_het: f64,
        other_het: Option<f64>
    ) -> GenotypePriorCalculator {
        match other_het {
            Some(other) => {
                GenotypePriorCalculator::genotype_prior_calculator(
                    snp_het, snp_het * 2.,
                    indel_het, indel_het * 2.,
                    other, other * 2.,
                )
            },
            _ => {
                GenotypePriorCalculator::genotype_prior_calculator(
                    snp_het, snp_het * 2.,
                    indel_het, indel_het * 2.,
                    std::cmp::max(OrderedFloat(snp_het), OrderedFloat(indel_het)).into_inner(), std::cmp::max(OrderedFloat(snp_het), OrderedFloat(indel_het)).into_inner() * 2.0
                )
            }
        }

    }

    pub fn make(args: &clap::ArgMatches) -> GenotypePriorCalculator {
        let snp_het = args.value_of("snp-heterozygosity").unwrap().parse::<f64>().unwrap();
        let ind_het = args.value_of("indel-heterozygosity").unwrap().parse::<f64>().unwrap();

        GenotypePriorCalculator::assuming_hw(snp_het, ind_het, None)
    }

    /**
     * Calculates the priors given the alleles to genetype and a likelihood calculator that determines the ploidy
     * of the sample at that site.
     * @param likelihood_calculator the input calculator
     * @param alleles the input alleles.
     * @throws IllegalArgumentException if either input is {@code null} or the calculator maximum number of supported alleles is less that the input allele size.
     * @return never {@code null}, the array will have as many positions as necessary to hold the priors of all possible
     * unphased genotypes as per the number of input alleles and the input calculator's ploidy.
     */
    pub fn get_log10_priors(&self, likelihood_calculator: &mut GenotypeLikelihoodCalculator, alleles: &Vec<Allele>) -> Vec<f64> {
        if likelihood_calculator.allele_count < alleles.len() {
            panic!("the number of alleles in the input calculator must be at least as large as the number of alleles in the input list")
        }

        let allele_types = GenotypePriorCalculator::calculate_allele_types(alleles);
        let mut number_of_genotypes = likelihood_calculator.genotype_count;
        if number_of_genotypes == -1 {
            number_of_genotypes = 0
        }
        let number_of_genotypes = number_of_genotypes as usize;
        let mut result = vec![0.; number_of_genotypes];

        for g in (1..number_of_genotypes).into_iter() {
            let gac = likelihood_calculator.genotype_allele_counts_at(g);
            result[g] = gac.sum_over_allele_indices_and_counts(
                |idx: usize, cnt: usize| {
                    if cnt == 2 {
                        self.hom_values[allele_types[idx] as usize]
                    } else {
                        self.het_values[allele_types[idx] as usize] + self.diff_values[allele_types[idx] as usize] * (cnt - 1) as f64
                    }
                }
            );
        }

        return result
    }

    fn calculate_allele_types(alleles: &Vec<Allele>) -> Vec<i64> {
        let ref_allele = &alleles[0];
        if !ref_allele.is_reference() {
            panic!("The first allele must be a valid reference")
        }
        // let ref_length = ref_allele.length();

        let result = alleles.par_iter().map(|allele| {
            if allele.is_reference() {
                return AlleleType::REF.ordinal() as i64
            } else if allele.is_called() && !allele.is_symbolic() {
                match allele.variant() {
                    Variant::Insertion(_) | Variant::Deletion(_) => {
                        AlleleType::INDEL.ordinal() as i64
                    },
                    Variant::SNV(_) => {
                        AlleleType::SNP.ordinal() as i64
                    },
                    _ => AlleleType::OTHER.ordinal() as i64
                }
            } else if allele.is_called() && allele.is_symbolic() {
                panic!("Cannot handle symbolic structural variants at the moment")
            } else {
                AlleleType::OTHER.ordinal() as i64
            }
        }).collect::<Vec<i64>>();

        result
    }

}

