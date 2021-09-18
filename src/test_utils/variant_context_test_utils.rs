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
}
