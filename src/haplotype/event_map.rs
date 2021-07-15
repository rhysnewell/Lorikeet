use haplotype::haplotype::Haplotype;
use utils::simple_interval::Locatable;
use std::collections::BTreeMap;
use model::variant_context::VariantContext;
use model::variants::Allele;

// lazy_static! {
//     pub static ref SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele
// }

/**
 * Extract simple VariantContext events from a single haplotype
 */
pub struct EventMap<'a, L: Locatable> {
    haplotype: Haplotype<L>,
    reference: &'a [u8],
    reference_loc: L,
    source_name_to_add: String,
    map: BTreeMap<usize, VariantContext>,
}

impl<'a, L: Locatable> EventMap<'a, L> {
    const MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION: usize = 3;
    const MAX_EVENT_PER_HAPLOTYPE: usize = 3;
    const MAX_INDELS_PER_HAPLOTYPE: usize = 3;

    pub fn new(haplotype: Haplotype<L>, reference: &'a [u8], reference_loc: L, source_name_to_add: String, max_mnp_distance: usize) -> EventMap<'a, L> {
        let mut result = EventMap {
            haplotype,
            reference,
            reference_loc,
            source_name_to_add,
            map: BTreeMap::new(),
        };
    }
}