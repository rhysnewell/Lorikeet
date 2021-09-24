use genotype::genotype_builder::Genotype;
use model::variant_context::VariantContext;
use utils::simple_interval::Locatable;

pub struct VariantContextTestUtils {}

impl VariantContextTestUtils {
    /**
     * VariantContext comparison function for testing
     * @param actual    vc derived from running test command
     * @param expected  vc we're hoping to get
     * @param attributesToIgnore    attributes (INFO or FORMAT) that may or may not exist in actual or expected
     * @param attributesWithJitter  attributes (INFO or FORMAT) that should existing in actual and expected, but may not match in value
     */
    pub fn assert_variant_contexts_are_equal(
        actual: &VariantContext,
        expected: &VariantContext,
        attributes_to_ignore: Vec<&str>,
        attributes_with_jitter: Vec<&str>,
    ) {
        assert_eq!(actual.loc.get_contig(), expected.loc.get_contig());
        assert_eq!(actual.loc.get_start(), expected.loc.get_start());
        assert_eq!(actual.loc.get_end(), expected.loc.get_end());
        assert_eq!(
            actual.get_alleles(),
            expected.get_alleles(),
            "Alleles not equal or wrong order?"
        );

        assert!(relative_eq!(
            actual.get_phred_scaled_qual(),
            expected.get_phred_scaled_qual(),
            epsilon = std::f64::EPSILON,
        ))
    }

    pub fn assert_genotypes_are_equal(
        actual: &Genotype,
        expected: &Genotype,
        extended_attributes_to_ignore: &Vec<String>,
    ) {
        assert_eq!(&actual.sample_name, &expected.sample_name, "genotype names");
        assert_eq!(&actual.alleles, &expected.alleles, "Genotype alleles");
        assert_eq!(
            &actual.genotype_type, &actual.genotype_type,
            "Genotype types"
        );

        // inline attributes
        println!("actual pl {:?} expected {:?}", &actual.pl, &expected.pl);
        assert_eq!(actual.has_dp(), expected.has_dp(), "Genotype hasDP");
        println!("actual dp {:?} expected {:?}", &actual.dp, &expected.dp);
        assert_eq!(&actual.dp, &expected.dp, "Genotype dp");
        assert_eq!(actual.has_ad(), expected.has_ad(), "Genotype hasAD");
        println!("actual ad {:?} expected {:?}", &actual.ad, &expected.ad);
        assert_eq!(&actual.ad, &expected.ad, "Genotype AD");
        assert_eq!(actual.has_gq(), expected.has_gq(), "Genotype hasGQ");
        assert_eq!(&actual.gq, &expected.gq, "Genotype gq");
        assert_eq!(&actual.pl, &expected.pl, "Genotype PL");

        assert_eq!(actual.has_likelihoods(), expected.has_likelihoods());
        assert_eq!(actual.get_likelihoods(), expected.get_likelihoods());

        assert_eq!(actual.is_phased, expected.is_phased);
        assert_eq!(actual.ploidy, expected.ploidy);
    }
}
