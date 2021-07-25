use haplotype::haplotype::Haplotype;
use model::variant_context::VariantContext;
use std::collections::BTreeMap;
use utils::simple_interval::Locatable;

// lazy_static! {
//     pub static ref SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele
// }

/**
 * Extract simple VariantContext events from a single haplotype
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct EventMap<L: Locatable> {
    // haplotype: Haplotype<'a, L>,
    // reference: &'a [u8],
    reference_loc: L,
    source_name_to_add: String,
    map: BTreeMap<usize, VariantContext>,
}

impl<L: Locatable> EventMap<L> {
    const MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION: usize = 3;
    const MAX_EVENT_PER_HAPLOTYPE: usize = 3;
    const MAX_INDELS_PER_HAPLOTYPE: usize = 3;

    pub fn new(
        reference_loc: L,
        source_name_to_add: String,
        max_mnp_distance: usize,
    ) -> EventMap<L> {
        let mut result = EventMap {
            // haplotype,
            // reference,
            reference_loc,
            source_name_to_add,
            map: BTreeMap::new(),
        };

        return result;
    }
}
