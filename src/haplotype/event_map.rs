use haplotype::haplotype::Haplotype;
use model::variant_context::VariantContext;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};
use utils::simple_interval::{Locatable, SimpleInterval};

// lazy_static! {
//     pub static ref SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele
// }

/**
 * Extract simple VariantContext events from a single haplotype
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct EventMap {
    // haplotype: Haplotype<'a, L>,
    // reference: &'a [u8],
    pub(crate) reference_loc: SimpleInterval,
    pub(crate) source_name_to_add: String,
    pub(crate) map: BTreeMap<usize, VariantContext>,
}

impl EventMap {
    const MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION: usize = 3;
    const MAX_EVENT_PER_HAPLOTYPE: usize = 3;
    const MAX_INDELS_PER_HAPLOTYPE: usize = 3;

    pub fn new(
        reference_loc: SimpleInterval,
        source_name_to_add: String,
        max_mnp_distance: usize,
    ) -> EventMap {
        let mut result = EventMap {
            // haplotype,
            // reference,
            reference_loc,
            source_name_to_add,
            map: BTreeMap::new(),
        };

        return result;
    }

    /**
     * Build event maps for each haplotype, returning the sorted set of all of the starting positions of all
     * events across all haplotypes
     *
     * @param haplotypes a list of haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     * @param debug if true, we'll emit debugging information during this operation
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occuring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @return a sorted set of start positions of all events among all haplotypes
     */
    pub fn build_event_maps_for_haplotypes<L: Locatable>(
        haplotypes: &mut Vec<Haplotype<L>>,
        ref_loc: SimpleInterval,
        max_mnp_distance: usize,
    ) {
        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        let mut hap_number = 0;
        debug!("=== Best Haplotypes ===");
        for h in haplotypes.iter_mut() {
            // Walk along the alignment and turn any difference from the reference into an event
            h.event_map = Some(EventMap::new(
                ref_loc.clone(),
                format!("HC{}", hap_number),
                max_mnp_distance,
            ));
            hap_number += 1;
            // Assert that all of the events discovered have 2 alleles
            h.event_map
                .as_ref()
                .unwrap()
                .map
                .values()
                .par_bridge()
                .for_each(|vc| {
                    if vc.get_alleles().len() == 2 {
                        // pass
                    } else {
                        panic!("Error Haplotype event map Variant Context has too many alleles")
                    }
                });

            debug!("{:?}", &h.genome_location.as_ref().unwrap());
            debug!("> Cigar {:?}", &h.cigar);
            debug!(">> Events = {:?}", &h.event_map.as_ref().unwrap().map)
        }
    }

    /**
     * Get the starting positions of events in this event map
     * @return
     */
    pub fn get_start_position(&self) -> Vec<&usize> {
        self.map.keys().collect::<Vec<&usize>>()
    }

    /**
     * Get the variant contexts in order of start position in this event map
     * @return
     */
    pub fn get_variant_contexts(&self) -> Vec<&VariantContext> {
        self.map.values().collect::<Vec<&VariantContext>>()
    }
}
