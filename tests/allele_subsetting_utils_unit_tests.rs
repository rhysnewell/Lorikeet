#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::genotype::genotype_builder::{
    AttributeObject, Genotype, GenotypeAssignmentMethod, GenotypesContext,
};
use lorikeet_genome::genotype::genotype_likelihoods::GenotypeLikelihoods;
use lorikeet_genome::genotype::genotype_prior_calculator::GenotypePriorCalculator;
use lorikeet_genome::model::allele_subsetting_utils::AlleleSubsettingUtils;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::variant_context_utils::VariantContextUtils;
use lorikeet_genome::model::variants::NON_REF_ALLELE;
use lorikeet_genome::test_utils::variant_context_test_utils::VariantContextTestUtils;
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::vcf_constants::*;

lazy_static! {
    static ref A_ref: ByteArrayAllele = ByteArrayAllele::new(b"A", true);
    static ref C: ByteArrayAllele = ByteArrayAllele::new(b"C", false);
    static ref G: ByteArrayAllele = ByteArrayAllele::new(b"G", false);
}

fn test_update_pls_and_ad_data(
    original_vc: VariantContext,
    selected_vc: VariantContext,
    expected_genotypes: Vec<Genotype>,
) {
    let gpc = GenotypePriorCalculator::empty();
    let mut selected_vc_with_gts = VariantContext::build_from_vc(&selected_vc);
    selected_vc_with_gts.add_genotypes(original_vc.get_genotypes().genotypes().clone());

    let old_gs = selected_vc_with_gts.get_genotypes();
    let actual = if selected_vc_with_gts.get_n_alleles() == original_vc.get_n_alleles() {
        old_gs.clone()
    } else {
        AlleleSubsettingUtils::subset_alleles(
            old_gs.clone(),
            0,
            original_vc.get_alleles(),
            selected_vc_with_gts.get_alleles(),
            &gpc,
            &GenotypeAssignmentMethod::DoNotAssignGenotypes,
            original_vc.get_attribute_as_int(&DEPTH_KEY, 0) as i32,
            false,
        )
    };

    assert_eq!(actual.len(), expected_genotypes.len());
    let dummy = Vec::new();
    for expected in expected_genotypes {
        let actual_gt = actual.get(0); // this should be the sample index relating to expected.sample
        VariantContextTestUtils::assert_genotypes_are_equal(&actual_gt, &expected, &dummy)
    }
}

#[test]
fn make_update_pls_sacs_and_ad_data() {
    let AA = vec![A_ref.clone(), A_ref.clone()];
    let AC = vec![A_ref.clone(), C.clone()];
    let CC = vec![C.clone(), C.clone()];
    // let CG = vec![C.clone(), G.clone()];
    let AG = vec![A_ref.clone(), G.clone()];
    // let GG = vec![G.clone(), G.clone()];
    let ACG = vec![A_ref.clone(), C.clone(), G.clone()];
    let mut vc_base = VariantContext::build(20, 10, 10, AC.clone());

    let hom_ref_pl = MathUtils::normalize_sum_to_one(vec![0.9, 0.09, 0.01]);
    let het_pl = MathUtils::normalize_sum_to_one(vec![0.09, 0.9, 0.01]);
    let hom_var_pl = MathUtils::normalize_sum_to_one(vec![0.01, 0.09, 0.9]);
    let uninformative = vec![0.0, 0.0, 0.0];

    // let mut base = Genotype::build_from_alleles(Vec::new(), "NA12878".to_string());

    // the simple case where no selection occurs
    let mut aa_gt = Genotype::build_from_alleles(AA, 0);
    aa_gt.pl = GenotypeLikelihoods::from_log10_likelihoods(hom_ref_pl).as_pls();
    aa_gt.ad = vec![10, 10];
    aa_gt.attribute(
        STRAND_COUNT_BY_SAMPLE_KEY.to_string(),
        AttributeObject::VecUnsize(vec![5, 10, 15, 20]),
    );
    aa_gt.gq = 8;

    let mut ac_gt = Genotype::build_from_alleles(AC.clone(), 0);
    ac_gt.pl = GenotypeLikelihoods::from_log10_likelihoods(het_pl).as_pls();
    ac_gt.ad = vec![10, 2];
    ac_gt.attribute(
        STRAND_COUNT_BY_SAMPLE_KEY.to_string(),
        AttributeObject::VecUnsize(vec![5, 10, 15, 20]),
    );
    ac_gt.gq = 8;

    let mut cc_gt = Genotype::build_from_alleles(CC.clone(), 0);
    cc_gt.pl = GenotypeLikelihoods::from_log10_likelihoods(hom_var_pl).as_pls();
    cc_gt.ad = vec![10, 2];
    cc_gt.attribute(
        STRAND_COUNT_BY_SAMPLE_KEY.to_string(),
        AttributeObject::VecUnsize(vec![5, 10, 15, 20]),
    );
    cc_gt.gq = 8;

    let mut vc_aa = vc_base.clone();
    vc_aa.add_genotypes(vec![aa_gt.clone()]);
    vc_base.alleles = AC.clone();
    test_update_pls_and_ad_data(vc_aa, vc_base.clone(), vec![aa_gt.clone()]);

    let mut vc_ac = vc_base.clone();
    vc_ac.add_genotypes(vec![ac_gt.clone()]);
    vc_base.alleles = AC.clone();
    test_update_pls_and_ad_data(vc_ac, vc_base.clone(), vec![ac_gt.clone()]);

    let mut vc_cc = vc_base.clone();
    vc_cc.add_genotypes(vec![cc_gt.clone()]);
    vc_base.alleles = AC.clone();
    test_update_pls_and_ad_data(vc_cc, vc_base.clone(), vec![cc_gt.clone()]);

    // uninformative test cases
    let mut uninformative_gt = Genotype::build_from_alleles(CC, 0);
    uninformative_gt.pl =
        GenotypeLikelihoods::from_log10_likelihoods(uninformative).as_pls();
    uninformative_gt.gq = 0;
    let mut vc_uninformative = vc_base.clone();
    vc_uninformative.add_genotypes(vec![uninformative_gt.clone()]);
    test_update_pls_and_ad_data(
        vc_uninformative,
        vc_base.clone(),
        vec![uninformative_gt],
    );

    let empty_gt = Genotype::build_from_alleles(
        VariantContextUtils::no_call_alleles(2),
        0,
    );
    let mut vc_empty = vc_base.clone();
    vc_empty.add_genotypes(vec![empty_gt.clone()]);
    test_update_pls_and_ad_data(vc_empty, vc_base.clone(), vec![empty_gt]);

    // actually subsetting down from multiple alt values
    // let homRef3AllelesPL = vec![0.0, -30.0, -60.0, -30.0, -60.0, -60.0];
    // let hetRefC3AllelesPL = vec![-20.0, 0.0, -20.0, -30.0, -40.0, -60.0];
    // let homC3AllelesPL = vec![-50.0, -30.0, 0.0, -70.0, -30.0, -70.0];
    // let hetRefG3AllelesPL = vec![-50.0, -30.0, -70.0, 0.0, -30.0, -20.0];
    // let homG3AllelesPL = vec![-50.0, -70.0, -70.0, -30.0, -30.0, 0.0]; // AA, AC, CC, AG, CG, GG;

    let homRef3AllelesAD = vec![20, 0, 1];
    let hetRefC3AllelesAD = vec![14, 7, 1];
    let homC3AllelesAD = vec![0, 20, 1];
    let hetRefG3AllelesAD = vec![14, 0, 7];
    let homG3AllelesAD = vec![0, 1, 21]; // AA, AC, CC, AG, CG, GG;

    let haploidRef3AllelesPL = vec![0.0, -50.0, -50.0];
    let haploidAltC3AllelesPL = vec![-30.0, 0.0, -60.0];
    let haploidAltG3AllelesPL = vec![-40.0, -70.0, 0.0];

    // for P=3 and N=2, the ordering is 000, 001, 011, 111, 002, 012, 112, 022, 122, 222
    let triploidRef3AllelesPL = vec![
        0.0, -30.0, -60.0, -90.0, -30.0, -60.0, -90.0, -60.0, -90.0, -90.0,
    ];
    let triploidAltC3AllelesPL = vec![
        -20.0, 0.0, -20.0, -50.0, -40.0, -70.0, -90.0, -90.0, -100.0, -100.0,
    ];
    let triploidAltG3AllelesPL = vec![
        -20.0, -40.0, -90.0, -100.0, 0.0, -70.0, -100.0, -20.0, -90.0, -50.0,
    ];

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG.clone();
    let mut genotype_a_ref =
        Genotype::build_from_alleles(vec![A_ref.clone()], 0);
    genotype_a_ref.ad = homRef3AllelesAD.clone();
    genotype_a_ref.pl = GenotypeLikelihoods::from_log10_likelihoods(haploidRef3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_a_ref]);

    let mut vc_ac = vc_base.clone();
    vc_ac.alleles = AC.clone();

    let mut genotype_a_ref_expected =
        Genotype::build_from_alleles(vec![A_ref.clone()], 0);
    genotype_a_ref_expected.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![0.0, -50.0]).as_pls();
    genotype_a_ref_expected.ad = vec![20, 0];
    genotype_a_ref_expected.gq = 500;
    test_update_pls_and_ad_data(vc_acg, vc_ac, vec![genotype_a_ref_expected]);

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG.clone();
    let mut genotype_c = Genotype::build_from_alleles(vec![C.clone()], 0);
    genotype_c.ad = homC3AllelesAD;
    genotype_c.pl = GenotypeLikelihoods::from_log10_likelihoods(haploidAltC3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_c]);

    let mut vc_ac = vc_base.clone();
    vc_ac.alleles = AC.clone();

    let mut genotype_c_expected =
        Genotype::build_from_alleles(vec![C.clone()], 0);
    genotype_c_expected.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![-30.0, 0.0]).as_pls();
    genotype_c_expected.ad = vec![0, 20];
    genotype_c_expected.gq = 300;
    test_update_pls_and_ad_data(vc_acg, vc_ac, vec![genotype_c_expected]);

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG.clone();
    let mut genotype_g = Genotype::build_from_alleles(vec![G.clone()], 0);
    genotype_g.ad = homG3AllelesAD;
    genotype_g.pl = GenotypeLikelihoods::from_log10_likelihoods(haploidAltG3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_g]);

    let mut vc_ag = vc_base.clone();
    vc_ag.alleles = AG.clone();

    let mut genotype_g_expected =
        Genotype::build_from_alleles(vec![G.clone()], 0);
    genotype_g_expected.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![-40.0, 0.0]).as_pls();
    genotype_g_expected.ad = vec![0, 21];
    genotype_g_expected.gq = 400;
    test_update_pls_and_ad_data(vc_acg, vc_ag, vec![genotype_g_expected]);

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG.clone();
    let mut genotype_3a = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), A_ref.clone()],
        0,
    );
    genotype_3a.ad = homRef3AllelesAD;
    genotype_3a.pl = GenotypeLikelihoods::from_log10_likelihoods(triploidRef3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_3a]);

    let mut vc_ac = vc_base.clone();
    vc_ac.alleles = AC.clone();

    let mut genotype_3a_expected = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), A_ref.clone()],
        0,
    );
    genotype_3a_expected.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![0.0, -30.0, -60.0, -90.0]).as_pls();
    genotype_3a_expected.ad = vec![20, 0];
    genotype_3a_expected.gq = 300;
    test_update_pls_and_ad_data(vc_acg, vc_ac, vec![genotype_3a_expected]);

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG.clone();
    let mut genotype_3a = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), C.clone()],
        0,
    );
    genotype_3a.ad = hetRefC3AllelesAD;
    genotype_3a.pl = GenotypeLikelihoods::from_log10_likelihoods(triploidAltC3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_3a]);

    let mut vc_ac = vc_base.clone();
    vc_ac.alleles = AC;

    let mut genotype_3a_expected = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), C.clone()],
        0,
    );
    genotype_3a_expected.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![-20.0, 0.0, -20.0, -50.0]).as_pls();
    genotype_3a_expected.ad = vec![14, 7];
    genotype_3a_expected.gq = 200;
    test_update_pls_and_ad_data(vc_acg, vc_ac, vec![genotype_3a_expected]);

    let mut vc_acg = vc_base.clone();
    vc_acg.alleles = ACG;
    let mut genotype_3a = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), G.clone()],
        0,
    );
    genotype_3a.ad = hetRefG3AllelesAD;
    genotype_3a.pl = GenotypeLikelihoods::from_log10_likelihoods(triploidAltG3AllelesPL).as_pls();
    vc_acg.genotypes = GenotypesContext::new(vec![genotype_3a]);

    let mut vc_ac = vc_base;
    vc_ac.alleles = AG;

    let mut genotype_3a_expected = Genotype::build_from_alleles(
        vec![A_ref.clone(), A_ref.clone(), G.clone()],
        0,
    );
    genotype_3a_expected.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![-20.0, 0.0, -20.0, -50.0]).as_pls();
    genotype_3a_expected.ad = vec![14, 7];
    genotype_3a_expected.gq = 200;
    test_update_pls_and_ad_data(vc_acg, vc_ac, vec![genotype_3a_expected]);
}

fn test_that_filtering_works_correctly(
    num_to_keep: usize,
    alleles: Vec<ByteArrayAllele>,
    scores: Vec<f64>,
    expected: Vec<ByteArrayAllele>,
) {
    assert_eq!(
        AlleleSubsettingUtils::filter_to_max_number_of_alt_alleles_based_on_scores(
            num_to_keep,
            alleles,
            scores
        ),
        expected
    );
}

#[test]
fn get_alleles_with_scores() {
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), G.clone()],
        vec![0.0, 5.0, 2.0],
        vec![A_ref.clone(), C.clone()],
    ); //first is best
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), G.clone()],
        vec![0.0, 2.0, 5.0],
        vec![A_ref.clone(), G.clone()],
    ); //second is best
    println!("first tie");
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), G.clone()],
        vec![0.0, 1.0, 1.0],
        vec![A_ref.clone(), C.clone()],
    ); //tie chooses first
    println!("second tie");
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), G.clone()],
        vec![5.0, 1.0, 1.0],
        vec![A_ref.clone(), C.clone()],
    ); //ref score is ignored chooses first
    test_that_filtering_works_correctly(
        2,
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone(), G.clone()],
        vec![0.0, 5.0, 0.0, 2.0],
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone(), G.clone()],
    ); //keep NON_REF in order
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone(), G.clone()],
        vec![0.0, 5.0, 0.0, 2.0],
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone()],
    ); //keep NON_REF in order when trimming
    test_that_filtering_works_correctly(
        1,
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone(), G.clone()],
        vec![0.0, 5.0, 0.0, 2.0],
        vec![A_ref.clone(), C.clone(), NON_REF_ALLELE.clone()],
    ); //keep NON_REF in order when trimming
}

#[test]
fn test_calculate_likelihood_sums() {
    // diploid, biallelic, three samples
    let two_alleles = vec![A_ref.clone(), C.clone()];

    // the canonical ordering of genotypes is 0/0, 0/1, 1/1

    // sample 1 has GLs: {1.1, 0.1, 2.3}, so that its likeliest genotype has two copies of the alt allele
    // and its GL difference weight (see javadoc of the tested method) is 2.3 - 1.1 = 1.2

    // sample 2 will have GLs: {3.1, 0.1, 2.3}, so that its likeliest genotype is hom ref
    // and it contributes nothing to the likelihood sum

    // sample 3 will have GLs: {1.1, 4.1, 2.3}, so that its likeliest genotype has one copy of the alt allele
    // and its GL difference weight is 4.3 - 1.1 = 3.0

    // the total likelihood sum is thus 1.2 + 0.0 + 3.0 = 4.2
    let mut g1 = Genotype::build_from_alleles(two_alleles.clone(), 0);
    g1.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![1.1, 0.1, 2.3]).as_pls();
    let mut g2 = Genotype::build_from_alleles(two_alleles.clone(), 1);
    g2.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![3.1, 0.1, 2.3]).as_pls();
    let mut g3 = Genotype::build_from_alleles(two_alleles.clone(), 2);
    g3.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![1.1, 4.1, 2.3]).as_pls();
    let g4 = Genotype::build_from_alleles(two_alleles.clone(), 3);

    let mut vc1 = VariantContext::build(0, 1, 1, two_alleles);
    vc1.add_genotypes(vec![g1, g2, g3, g4]);
    assert!(
        relative_eq!(
            AlleleSubsettingUtils::calculate_likelihood_sums(&vc1, 2, false)[1],
            4.2,
            epsilon = 1.0e-8
        ),
        "value {} expected {}",
        AlleleSubsettingUtils::calculate_likelihood_sums(&vc1, 2, false)[1],
        4.2
    );

    // diploid, triallelic, three samples
    let three_alleles = vec![A_ref.clone(), C.clone(), G.clone()];
    // the canonical ordering of genotypes is 0/0, 0/1, 1/1, 0/2, 1/2, 2/2

    // sample 1 has GLs: {0.0, 0.0, 0.0, 0.0, 3.1, 0.0}, so that its likeliest genotype has one copy of each
    // alt allele and its GL difference weight is 3.1
    // thus it contributes 3.1 to each alt allele

    // sample 2 has GLs: {0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, so that its likeliest genotype has one copy of the first
    // alt allele and its GL difference weight is 1.0
    // thus it contributes 1.0 to the first alt allele

    // the likelihood sums are 4.1 and 3.1
    let mut g4 = Genotype::build_from_alleles(vec![C.clone(), G.clone()], 0);
    g4.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![0.0, 0.0, 0.0, 0.0, 3.1, 0.0]).as_pls();
    let mut g5 =
        Genotype::build_from_alleles(vec![A_ref.clone(), C.clone()], 1);
    g5.pl =
        GenotypeLikelihoods::from_log10_likelihoods(vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0]).as_pls();

    let mut vc2 = VariantContext::build(0, 1, 1, three_alleles);
    vc2.add_genotypes(vec![g4, g5]);
    let likelihood_sums = AlleleSubsettingUtils::calculate_likelihood_sums(&vc2, 2, false);
    assert!(
        relative_eq!(likelihood_sums[1], 4.1, epsilon = 1.0e-8),
        "value {} expected {}",
        likelihood_sums[1],
        4.1
    );
    assert!(
        relative_eq!(likelihood_sums[2], 3.1, epsilon = 1.0e-8),
        "value {} expected {}",
        likelihood_sums[2],
        3.1
    );
}
