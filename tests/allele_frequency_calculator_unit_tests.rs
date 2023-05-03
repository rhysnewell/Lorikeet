#![allow(
    non_upper_case_globals,
    non_snake_case
)]

use lorikeet_genome::genotype::genotype_builder::Genotype;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::genotype::genotype_likelihoods::GenotypeLikelihoods;
use lorikeet_genome::model::allele_frequency_calculator::AlleleFrequencyCalculator;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::variants::SPAN_DEL_ALLELE;
use lorikeet_genome::utils::math_utils::{MathUtils, LOG10_ONE_HALF};
use rayon::prelude::*;

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

lazy_static! {
    static ref A: ByteArrayAllele = ByteArrayAllele::new("A".as_bytes(), true);
    static ref B: ByteArrayAllele = ByteArrayAllele::new("C".as_bytes(), false);
    static ref C: ByteArrayAllele = ByteArrayAllele::new("G".as_bytes(), false);
}

static EPS: f64 = 1e-3;
static DIPLOID: usize = 2;
static TRIPLOID: usize = 3;
static BIALLELIC: usize = 2;
static TRIALLELIC: usize = 3;

static EXTREMELY_CONFIDENT_PL: i64 = 1000;
static FAIRLY_CONFIDENT_PL: i64 = 20;

static DEFAULT_PLOIDY: usize = 2;

#[test]
fn test_symmetries() {
    let mut sample = 0;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 0.1, 0.1, DEFAULT_PLOIDY);
    let alleles = vec![A.clone(), B.clone(), C.clone()];
    let AA =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![0, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let BB =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![1, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let CC =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![2, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let AB = genotype_with_obvious_call(
        DIPLOID,
        TRIALLELIC,
        vec![0, 1, 1, 1],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;
    let AC = genotype_with_obvious_call(
        DIPLOID,
        TRIALLELIC,
        vec![0, 1, 2, 1],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;

    let BBB = genotype_with_obvious_call(
        TRIPLOID,
        TRIALLELIC,
        vec![1, 3],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;
    let CCC = genotype_with_obvious_call(
        TRIPLOID,
        TRIALLELIC,
        vec![2, 3],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    // sample += 1;

    let switch_b_with_c_pairs = vec![
        (
            make_vc(alleles.clone(), vec![AA.clone(), BB.clone()]),
            make_vc(alleles.clone(), vec![AA.clone(), CC.clone()]),
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AB.clone()]),
            make_vc(alleles.clone(), vec![AA.clone(), AC.clone()]),
        ),
        (
            make_vc(alleles.clone(), vec![AB.clone(), AB.clone()]),
            make_vc(alleles.clone(), vec![AC.clone(), AC.clone()]),
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AA.clone(), BB.clone()]),
            make_vc(alleles.clone(), vec![AA.clone(), AA.clone(), CC.clone()]),
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AB.clone(), AB.clone()]),
            make_vc(alleles.clone(), vec![AA.clone(), AC.clone(), AC.clone()]),
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), BBB.clone()]),
            make_vc(alleles.clone(), vec![AA.clone(), CCC.clone()]),
        ),
    ];

    for (_, pair) in switch_b_with_c_pairs.into_iter().enumerate() {
        let vc1 = pair.0;
        let vc2 = pair.1;
        let result1 = af_calc.calculate(vc1, DEFAULT_PLOIDY);
        let result2 = af_calc.calculate(vc2, DEFAULT_PLOIDY);

        assert!(
            relative_eq!(
                result1.log10_prob_only_ref_allele_exists(),
                result2.log10_prob_only_ref_allele_exists(),
                epsilon = EPS
            ),
            "Symmetry failed prob only ref {} -> {}",
            result1.log10_prob_only_ref_allele_exists(),
            result2.log10_prob_only_ref_allele_exists()
        );
        assert!(
            relative_eq!(
                result1.get_log10_posterior_of_allele_absent(&B),
                result2.get_log10_posterior_of_allele_absent(&C),
                epsilon = EPS
            ),
            "Symmetry failed {} -> {}",
            result1.get_log10_posterior_of_allele_absent(&B),
            result2.get_log10_posterior_of_allele_absent(&C)
        );
    }
}

#[test]
fn test_MLE_counts() {
    let mut sample = 0;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 1.0, 1.0, DEFAULT_PLOIDY);
    let alleles = vec![A.clone(), B.clone(), C.clone()];
    let AA =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![0, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let BB =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![1, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    // let CC =
    //     genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![2, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let AB = genotype_with_obvious_call(
        DIPLOID,
        TRIALLELIC,
        vec![0, 1, 1, 1],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;
    let AC = genotype_with_obvious_call(
        DIPLOID,
        TRIALLELIC,
        vec![0, 1, 2, 1],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;

    let BBB = genotype_with_obvious_call(
        TRIPLOID,
        TRIALLELIC,
        vec![1, 3],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    sample += 1;
    let CCC = genotype_with_obvious_call(
        TRIPLOID,
        TRIALLELIC,
        vec![2, 3],
        FAIRLY_CONFIDENT_PL,
        sample,
    );
    // sample += 1;

    let vc_with_expected_counts = vec![
        (
            make_vc(alleles.clone(), vec![AA.clone(), BB.clone()]),
            vec![2, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AB.clone()]),
            vec![1, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AB.clone(), AB.clone()]),
            vec![2, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AA.clone(), BB.clone()]),
            vec![2, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AB.clone(), AB.clone()]),
            vec![2, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), BBB.clone()]),
            vec![3, 0],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), BBB.clone(), CCC.clone()]),
            vec![3, 3],
        ),
        (
            make_vc(alleles.clone(), vec![AA.clone(), AB.clone(), AC.clone()]),
            vec![1, 1],
        ),
        (
            make_vc(
                alleles.clone(),
                vec![AA.clone(), AB.clone(), AC.clone(), BBB.clone(), CCC.clone()],
            ),
            vec![4, 4],
        ),
    ];

    for (_, pair) in vc_with_expected_counts.into_iter().enumerate() {
        let vc = pair.0;
        let expected = pair.1;
        let actual = af_calc
            .calculate(vc, DEFAULT_PLOIDY)
            .get_allele_counts_of_mle();
        assert_eq!(actual, expected);
    }
}

// many samples with low confidence should yield a non-zero MLE, in contrast to the old exact model
#[test]
fn test_many_samples_with_low_confidence() {
    // prior corresponding to 1000 observations of ref, 1 of a SNP
    // for this test, we want many pseudocounts in the prior because the new AF calculator learns the allele frequency
    // and we don't want the complication of the posterior being differetn from the prior
    let mut af_calc = AlleleFrequencyCalculator::new(1000.0, 1.0, 1.0, DEFAULT_PLOIDY);
    let alleles = vec![A.clone(), B.clone()];

    // for FAIRLY_CONFIDENT_PL = 20, this genotype has about 100 times greater likelihood to be het than hom ref
    // with our prior giving 1000 times as much weight to ref, this implies a 1 in 5 chance of each sample having a copy of the alt allele
    // (that is, 100/1000 times the combinatorial factor of 2).  Thus the MLE for up to 2 samples should be zero
    // for five samples we should have one
    // for ten samples we will have more than twice as many as for five since the counts fromt he samples start to influence
    // the estimated allele frequency
    let AB =
        genotype_with_obvious_call(DIPLOID, BIALLELIC, vec![0, 1, 1, 1], FAIRLY_CONFIDENT_PL, 0);

    let vcs_with_different_numbers_of_samples = (1..11)
        .into_par_iter()
        .map(|n| make_vc(alleles.clone(), vec![AB.clone(); n]))
        .collect::<Vec<VariantContext>>();

    let counts = vcs_with_different_numbers_of_samples
        .into_iter()
        .map(|vc| {
            let result = af_calc.calculate(vc, DEFAULT_PLOIDY);
            result.get_allele_count_at_mle(&B)
        })
        .collect::<Vec<i64>>();
    assert_eq!(counts[0], 0);
    assert_eq!(counts[1], 0);
    assert_eq!(counts[4], 2);
    assert!(counts[8] >= 3);
}

#[test]
fn test_many_very_confident_samples() {
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 1.0, 1.0, DEFAULT_PLOIDY);
    let alleles = vec![A.clone(), B.clone(), C.clone()];

    let AC = genotype_with_obvious_call(
        DIPLOID,
        TRIALLELIC,
        vec![0, 1, 2, 1],
        EXTREMELY_CONFIDENT_PL,
        0,
    );
    for num_samples in vec![100, 1000] {
        let vc = make_vc(alleles.clone(), vec![AC.clone(); num_samples]);
        let result = af_calc.calculate(vc, DEFAULT_PLOIDY);
        assert_eq!(result.get_allele_count_at_mle(&B), 0);
        assert_eq!(result.get_allele_count_at_mle(&C), num_samples as i64);

        assert!(relative_eq!(
            result.log10_prob_only_ref_allele_exists(),
            result.get_log10_posterior_of_allele_absent(&C),
            epsilon = num_samples as f64 * 0.01
        ));
        // with a large number of samples all with the AC genotype, the calculator will learn that the frequencies of the A and C alleles
        // are 1/2, while the frequency of the B allele is 0.  Thus the only genotypes with appreciable priors are AA, AC, and CC
        // with priors of 1/4, 1/2, and 1/4 and relative likelihoods of 1, 10^(PL/10), and 1
        // the posterior probability of each sample having the C allele is thus
        // (1 + 2*10^(PL/10))/(1 + 2*10^(PL/10) + 1) = (1 + x/2)/(1 + x), where x = 10^(-PL/10)

        // to first-order in x, which is an extremely good approximation, this is 1 - x/2
        // thus the probability that N identical samples don't have the C allele is (x/2)^N, and the log-10 probability of this is
        // N * [log_10(1/2) - PL/10]
        let expected_log10_probability_of_no_C_allele =
            num_samples as f64 * (*LOG10_ONE_HALF - EXTREMELY_CONFIDENT_PL as f64 / 10.0);
        assert!(relative_eq!(
            result.get_log10_posterior_of_allele_absent(&C),
            expected_log10_probability_of_no_C_allele,
            epsilon = num_samples as f64 * 0.01
        ));
    }
}

#[test]
fn test_approximate_multiplicative_confidence() {
    let mut sample = 0;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 1.0, 1.0, DEFAULT_PLOIDY); //flat prior -- we will choose genotypes such that the posterior remains flat
    let alleles = vec![A.clone(), B.clone()];
    let AA =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![0, 2], FAIRLY_CONFIDENT_PL, sample);
    sample += 1;
    let BB =
        genotype_with_obvious_call(DIPLOID, TRIALLELIC, vec![1, 2], FAIRLY_CONFIDENT_PL, sample);
    // sample += 1;

    let mut vcs_with_different_numbers_of_samples = Vec::new();
    let mut genotype_list = Vec::new();

    for _ in 0..10 {
        genotype_list.push(AA.clone());
        genotype_list.push(BB.clone()); //adding both keeps the flat prior.  Thus the posterior will equal the likelihood
        vcs_with_different_numbers_of_samples.push(make_vc(alleles.clone(), genotype_list.clone()));
    }

    // since we maintain a flat allele frequency distribution, the probability of being ref as each successive sample is added
    // is multiplied by the probability of any one.  Thus we get an arithmetic series in log space
    let log10_p_refs = vcs_with_different_numbers_of_samples
        .into_iter()
        .map(|vc| {
            af_calc
                .calculate(vc, DEFAULT_PLOIDY)
                .log10_prob_only_ref_allele_exists()
        })
        .collect::<Vec<f64>>();

    for n in 0..9 {
        assert!(relative_eq!(
            log10_p_refs[n + 1] - log10_p_refs[n],
            log10_p_refs[0],
            epsilon = 0.01
        ));
    }
}

#[test]
fn test_many_ref_samples_dont_kill_good_variant() {
    let mut sample = 0;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 0.1, 0.1, DEFAULT_PLOIDY);
    let alleles = vec![A.clone(), B.clone()];
    let AA = genotype_with_obvious_call(
        DIPLOID,
        BIALLELIC,
        vec![0, 2],
        EXTREMELY_CONFIDENT_PL,
        sample,
    );
    sample += 1;
    let AB = genotype_with_obvious_call(
        DIPLOID,
        BIALLELIC,
        vec![0, 1, 1, 1],
        EXTREMELY_CONFIDENT_PL,
        sample,
    );
    // sample += 1;

    for num_ref in vec![1, 10, 100, 1000, 10000, 100000] {
        let mut genotype_list = vec![AA.clone(); num_ref];
        genotype_list.push(AB.clone());
        let vc = make_vc(alleles.clone(), genotype_list);
        let log10_ref = af_calc
            .calculate(vc, DEFAULT_PLOIDY)
            .log10_prob_only_ref_allele_exists();
        assert!(
            log10_ref < (-(EXTREMELY_CONFIDENT_PL as f64) / 10.0) + (num_ref as f64).log10() + 1.0
        );
    }
}

#[test]
fn test_spanning_deletion_is_not_considered_variant() {
    let ploidy = 2;
    let mut sample = 0;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 0.1, 0.1, ploidy);
    let alleles = vec![A.clone(), B.clone(), SPAN_DEL_ALLELE.clone()];

    // some pls that have high likelihood for span del allele but not for the SNP (B)
    let span_del_pls = vec![50, 100, 100, 0, 100, 100];

    // some pls that weakly support the SNP
    let low_qual_snp_pls = vec![10, 0, 40, 100, 70, 300];

    let span_del = make_genotype(ploidy, sample, span_del_pls);
    sample += 1;
    let low_qual_snp = make_genotype(ploidy, sample, low_qual_snp_pls);
    // sample += 1;

    // first test the span del genotype alone.  Its best PL containing the SNP is 100, so we expect a variant probability
    // of about 10^(-100/10) -- a bit less due to the prior bias in favor of the reference
    let vc_span_del = make_vc(alleles.clone(), vec![span_del.clone()]);
    let log10_p_variant = af_calc
        .calculate(vc_span_del, ploidy)
        .log10_prob_variant_present();
    assert!(log10_p_variant < -10.0);

    // now test a realistic situation of two samples, one with a low-quality SNP and one with the spanning deletion
    // we want to find that the spanning deletion has little effect on the qual
    // In fact, we also want to check that it *decreases* the qual, because it's essentially one more hom ref sample
    // Furthermore, to be precise it should be really behave almost identically to a hom ref *haploid* sample,
    // so we check that, too
    let vc_low_qual_snp = make_vc(alleles.clone(), vec![low_qual_snp.clone()]);
    let low_qual_snp_score = af_calc
        .calculate(vc_low_qual_snp, ploidy)
        .log10_prob_variant_present();
    let vc_both = make_vc(
        alleles.clone(),
        vec![low_qual_snp.clone(), span_del.clone()],
    );

    let both_qual_score = af_calc
        .calculate(vc_both, ploidy)
        .log10_prob_variant_present();
    assert!(relative_eq!(
        low_qual_snp_score,
        both_qual_score,
        epsilon = 0.1
    ));
    assert!(both_qual_score < low_qual_snp_score);

    let haploid_ref_pls = vec![0, 100, 100];
    let haploid_ref = make_genotype(1, sample, haploid_ref_pls);
    sample += 1;

    let vc_low_qual_snp_and_haploid_ref = make_vc(
        alleles.clone(),
        vec![low_qual_snp.clone(), haploid_ref.clone()],
    );
    let low_qual_snp_and_haploid_ref_qual_score = af_calc
        .calculate(vc_low_qual_snp_and_haploid_ref, ploidy)
        .log10_prob_variant_present();
    assert!(relative_eq!(
        low_qual_snp_and_haploid_ref_qual_score,
        both_qual_score,
        epsilon = 1e-5
    ));

    // as a final test, we check that getting rid of the spanning deletion allele, in the sense that
    // REF / SPAN_DEL --> haploid REF; REF / SNP --> REF / SNP
    // does not affect the qual score
    let haploid_ref_pls_without_span_del = vec![0, 100];
    let snp_pls_without_span_del = vec![10, 0, 40];
    let vc_no_span_del = make_vc(
        vec![A.clone(), B.clone()],
        vec![
            make_genotype(ploidy, sample, snp_pls_without_span_del),
            make_genotype(1, sample + 1, haploid_ref_pls_without_span_del),
        ],
    );
    // sample += 2;
    let no_span_del_qual_score = af_calc
        .calculate(vc_no_span_del, ploidy)
        .log10_prob_variant_present();
    assert!(relative_eq!(
        no_span_del_qual_score,
        both_qual_score,
        epsilon = 1e-6
    ));
}

#[test]
fn test_presence_of_unlikely_spanning_deletion_doesnt_affect_results() {
    let ploidy = 2;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 0.1, 0.1, ploidy);
    let alleles_without_span_del = vec![A.clone(), B.clone()];
    let alleles_with_span_del = vec![A.clone(), B.clone(), SPAN_DEL_ALLELE.clone()];

    // make PLs that support an A/B genotype
    let pls_without_span_del = vec![50, 0, 50];
    let pls_with_span_del = vec![50, 0, 50, 100, 100, 100];

    let genotype_without_span_del = make_genotype(ploidy, 0, pls_without_span_del);
    let genotype_with_span_del = make_genotype(ploidy, 0, pls_with_span_del);
    let vc_without_span_del = make_vc(alleles_without_span_del, vec![genotype_without_span_del]);
    let vc_with_span_del = make_vc(alleles_with_span_del, vec![genotype_with_span_del]);

    let log10_p_variant_without_span_del = af_calc
        .calculate(vc_without_span_del, ploidy)
        .log10_prob_variant_present();
    let log10_p_variant_with_span_del = af_calc
        .calculate(vc_with_span_del, ploidy)
        .log10_prob_variant_present();
    assert!(relative_eq!(
        log10_p_variant_with_span_del,
        log10_p_variant_without_span_del,
        epsilon = 1e-4
    ));
}

// test that a finite precision bug for span del sites with a very unlikely alt allele doesn't occur
#[test]
fn test_spanning_deletion_with_very_unlikely_alt_allele() {
    let ploidy = 4;
    let mut af_calc = AlleleFrequencyCalculator::new(1.0, 0.1, 0.1, ploidy);
    let alleles = vec![A.clone(), SPAN_DEL_ALLELE.clone(), B.clone()];

    // make PLs that don't support the alt allele
    let pls = vec![
        0, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,
        10000, 10000,
    ];
    let vc = make_vc(alleles, vec![make_genotype(ploidy, 0, pls)]);
    let _log10_p_variant = af_calc.calculate(vc, ploidy).log10_prob_variant_present();
}

#[test]
fn test_single_sample_biallelic_shortcut() {
    // in the haploid case, if the AF calc has equal pseudocounts of ref and alt, the posterior is proportional to the likelihoods:
    for pseudo_count in vec![1.0, 5.0, 10.0] {
        let af_calc = AlleleFrequencyCalculator::new(
            pseudo_count,
            pseudo_count,
            pseudo_count,
            DEFAULT_PLOIDY,
        );
        for pl in vec![10, 100, 1000] {
            let log_10_likelihoods = vec![0.0, (pl as f64) / 10.0];
            let result = af_calc
                .calculate_single_sample_biallelic_non_ref_posterior(&log_10_likelihoods, false);
            let expected = MathUtils::normalize_log10(log_10_likelihoods, false)[1];
            assert!(relative_eq!(result, expected, epsilon = 1e-10));
        }
    }

    // in the diploid case, we roughly multiply the prior by the likelihoods -- it's not exact because the allele frequency is a random variable and not a single value
    for heterozygosity in vec![0.1, 0.01, 0.001] {
        let af_calc = AlleleFrequencyCalculator::new(
            100.0,
            100.0 * heterozygosity,
            100.0 * heterozygosity,
            DEFAULT_PLOIDY,
        );
        for pl in vec![10, 100, 1000] {
            let log10_likelihood = (pl as f64) / 10.0;
            let log10_likelihoods = vec![0.0, log10_likelihood, -100.0];
            let log10_priors = vec![
                (1.0 - heterozygosity).powf(2.0).log10(),
                (2.0 * heterozygosity * (1.0 - heterozygosity)).log10(),
                heterozygosity.powf(2.0),
            ];
            let result = af_calc
                .calculate_single_sample_biallelic_non_ref_posterior(&log10_likelihoods, false);
            let expected = 1.0
                - MathUtils::normalize_log10(
                    MathUtils::ebe_add(&log10_likelihoods, &log10_priors),
                    false,
                )[0];
            assert!(relative_eq!(result, expected, epsilon = 0.3));
        }
    }
}

// make PLs that correspond to an obvious call i.e. one PL is relatively big and the rest are zero
// alleleCounts is the GenotypeAlleleCounts format for the obvious genotype, with repeats but in no particular order
fn PLs_for_obvious_call(
    ploidy: usize,
    num_alleles: usize,
    allele_counts: Vec<usize>,
    PL: i64,
) -> Vec<i64> {
    let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, num_alleles);
    let mut result = vec![PL; gl_calc.genotype_count as usize];
    result[gl_calc.allele_counts_to_index(&allele_counts)] = 0;
    return result;
}

fn genotype_with_obvious_call(
    ploidy: usize,
    num_alleles: usize,
    alleles: Vec<usize>,
    PL: i64,
    sample: usize,
) -> Genotype {
    return make_genotype(
        ploidy,
        sample,
        PLs_for_obvious_call(ploidy, num_alleles, alleles, PL),
    );
}

fn make_genotype(ploidy: usize, sample: usize, pls: Vec<i64>) -> Genotype {
    let mut g = Genotype::build_from_alleles(
        vec![ByteArrayAllele::no_call(); ploidy],
        format!("sample_{}", sample),
    );
    g.pl(GenotypeLikelihoods::from_pls(pls).as_pls());
    return g;
}

fn make_vc(alleles: Vec<ByteArrayAllele>, genotypes: Vec<Genotype>) -> VariantContext {
    let mut vc = VariantContext::build(0, 1, 1, alleles);
    vc.add_genotypes(genotypes);
    return vc;
}
