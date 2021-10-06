use haplotype::haplotype::Haplotype;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::{VariantContext, VariantType};
use rayon::prelude::*;
use rust_htslib::bam::record::Cigar;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use utils::base_utils::BaseUtils;
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

    pub fn empty() -> Self {
        Self {
            reference_loc: SimpleInterval::new(0, 0, 0),
            source_name_to_add: "".to_string(),
            map: BTreeMap::new(),
        }
    }

    pub fn new<L: Locatable>(
        haplotype: &Haplotype<L>,
        reference: &[u8],
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

        result.process_cigar_for_initial_events(haplotype, reference, max_mnp_distance);

        return result;
    }

    pub fn state_for_testing(vcs: Vec<VariantContext>) -> EventMap {
        let mut result = EventMap {
            reference_loc: SimpleInterval::new(0, 0, 0),
            source_name_to_add: " ".to_string(),
            map: BTreeMap::new(),
        };

        for vc in vcs {
            result.add_vc(vc, true);
        }

        return result;
    }

    pub fn get_number_of_events(&self) -> usize {
        self.map.len()
    }

    /**
     *
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     */
    fn process_cigar_for_initial_events<L: Locatable>(
        &mut self,
        haplotype: &Haplotype<L>,
        reference: &[u8],
        max_mnp_distance: usize,
    ) {
        let cigar = &haplotype.cigar;
        let alignment = haplotype.get_bases();

        let mut ref_pos = haplotype.alignment_start_hap_wrt_ref;

        let mut proposed_events = Vec::new();

        let mut alignment_pos = 0;
        for cigar_index in 0..cigar.0.len() {
            let ce = &cigar.0[cigar_index];
            debug!("Cigar {:?}", ce);

            match ce {
                Cigar::Ins(len) => {
                    if ref_pos > 0 {
                        // protect against trying to create insertions/deletions at the beginning of a contig
                        let mut insertion_alleles = Vec::new();
                        let insertion_start = self.reference_loc.start + ref_pos - 1;
                        let ref_byte = reference[ref_pos - 1];
                        if BaseUtils::is_regular_base(ref_byte) {
                            insertion_alleles
                                .push(ByteArrayAllele::new(vec![ref_byte].as_slice(), true));
                        };
                        if cigar_index == 0 || cigar_index == cigar.0.len() - 1 {
                            // if the insertion isn't completely resolved in the haplotype, skip it
                            // note this used to emit SYMBOLIC_UNASSEMBLED_EVENT_ALLELE but that seems dangerous
                        } else {
                            let mut insertion_bases = Vec::new();
                            insertion_bases.push(reference[ref_pos - 1]); // add the padding base
                            insertion_bases.extend_from_slice(
                                &alignment[alignment_pos..(alignment_pos + *len as usize)],
                            );
                            if BaseUtils::is_all_regular_base(insertion_bases.as_slice()) {
                                insertion_alleles
                                    .push(ByteArrayAllele::new(insertion_bases.as_slice(), false));
                            };
                        };
                        if insertion_alleles.len() == 2 {
                            // found a proper ref and alt allele
                            let mut proposed_vc = VariantContext::build(
                                self.reference_loc.tid,
                                insertion_start,
                                insertion_start,
                                insertion_alleles,
                            );
                            proposed_vc.source = self.source_name_to_add.clone();
                            proposed_vc.get_type();
                            proposed_events.push(proposed_vc);
                        }
                    }
                    alignment_pos += *len as usize;
                }
                Cigar::SoftClip(len) => {
                    alignment_pos += *len as usize;
                }
                Cigar::Del(len) => {
                    if ref_pos > 0 {
                        // protect against trying to create insertions/deletions at the beginning of a contig
                        let mut deletion_bases = &reference[ref_pos - 1..(ref_pos + *len as usize)];
                        let mut deletion_alleles = Vec::new();
                        let deletion_start = self.reference_loc.start + ref_pos - 1;
                        let ref_byte = reference[ref_pos - 1];
                        if BaseUtils::is_regular_base(ref_byte)
                            && BaseUtils::is_all_regular_base(deletion_bases)
                        {
                            deletion_alleles.push(ByteArrayAllele::new(deletion_bases, true));
                            deletion_alleles
                                .push(ByteArrayAllele::new(vec![ref_byte].as_slice(), false));
                            let mut proposed_vc = VariantContext::build(
                                self.reference_loc.tid,
                                deletion_start,
                                deletion_start + *len as usize,
                                deletion_alleles,
                            );
                            proposed_vc.source = self.source_name_to_add.clone();
                            proposed_vc.get_type();
                            proposed_events.push(proposed_vc);
                        }
                    }
                    ref_pos += *len as usize;
                }
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let mut mismatch_offsets = VecDeque::new();
                    for offset in 0..(*len as usize) {
                        let ref_byte = reference[ref_pos + offset];
                        let alt_byte = alignment[alignment_pos + offset];
                        let mismatch = ref_byte != alt_byte
                            && BaseUtils::is_regular_base(ref_byte)
                            && BaseUtils::is_regular_base(alt_byte);

                        // debug!(
                        //     "Cigar {:?} ref_byte {} alt_byte {} mismatch {}",
                        //     ce, ref_byte, alt_byte, mismatch
                        // );
                        if mismatch {
                            mismatch_offsets.push_back(offset);
                        };
                    }

                    while !mismatch_offsets.is_empty() {
                        let start = mismatch_offsets.pop_front().unwrap();
                        let mut end = start;
                        while !mismatch_offsets.is_empty() {
                            match mismatch_offsets.front() {
                                None => break,
                                Some(val) => {
                                    if val - end <= max_mnp_distance {
                                        end = mismatch_offsets.pop_front().unwrap();
                                    } else {
                                        break;
                                    };
                                }
                            }
                            // while *mismatch_offsets.front().unwrap() - end <= max_mnp_distance {
                            //     end = mismatch_offsets.pop_front().unwrap();
                            // };
                            // break;
                        }
                        let ref_allele = ByteArrayAllele::new(
                            &reference[(ref_pos + start)..(ref_pos + end + 1)],
                            true,
                        );
                        let alt_allele = ByteArrayAllele::new(
                            &alignment[(alignment_pos + start)..(alignment_pos + end + 1)],
                            false,
                        );
                        let mut proposed_vc = VariantContext::build(
                            self.reference_loc.tid,
                            self.reference_loc.start + ref_pos + start,
                            self.reference_loc.start + ref_pos + end,
                            vec![ref_allele, alt_allele],
                        );
                        proposed_vc.source = self.source_name_to_add.clone();
                        proposed_vc.get_type();
                        proposed_events.push(proposed_vc);
                    }

                    // move refPos and alignmentPos forward to the end of this cigar element
                    ref_pos += *len as usize;
                    alignment_pos += *len as usize;
                }
                Cigar::RefSkip(_) | Cigar::Pad(_) | Cigar::HardClip(_) => {
                    panic!(
                        "Unsupported cigar operator created during SW alignment: {:?}",
                        ce
                    );
                }
            };
        }

        debug!("Found {} events", proposed_events.len());
        for proposed_event in proposed_events {
            debug!("Adding event {:?}", &proposed_event);
            self.add_vc(proposed_event, true)
        }
    }

    /**
     * Add VariantContext vc to this map
     * @param vc the variant context to add
     * @param merge should we attempt to merge it with an already existing element, or should we throw an error in that case?
     */
    pub fn add_vc(&mut self, vc: VariantContext, merge: bool) {
        if self.map.contains_key(&vc.loc.start) {
            if merge {
                let prev = self.map.remove(&vc.loc.start).unwrap();
                self.map.insert(vc.loc.start, Self::make_block(prev, vc));
            }
        } else {
            self.map.insert(vc.loc.start, vc);
        };
    }

    /**
     * Create a block substitution out of two variant contexts that start at the same position
     *
     * vc1 can be SNP, and vc2 can then be either a insertion or deletion.
     * If vc1 is an indel, then vc2 must be the opposite type (vc1 deletion => vc2 must be an insertion)
     *
     * @param vc1 the first variant context we want to merge
     * @param vc2 the second
     * @return a block substitution that represents the composite substitution implied by vc1 and vc2
     */
    pub fn make_block(mut vc1: VariantContext, mut vc2: VariantContext) -> VariantContext {
        assert!(
            vc1.loc.start == vc2.loc.start,
            "vc1 and 2 must have the same start but got got {} and {}",
            vc1.loc.start,
            vc2.loc.start
        );
        assert!(vc1.is_biallelic(), "vc1 must be biallelic");
        if !vc1.is_snp() {
            assert!(
                (vc1.is_simple_deletion() && vc2.is_simple_insertion())
                    || (vc1.is_simple_insertion() && vc2.is_simple_deletion()),
                "can only merge single insertion with deletion"
            )
        } else {
            assert!(
                !vc2.is_snp(),
                "vc1 and vc2 are both SNPs, which implies a terrible bug in the cigar"
            );
        }

        let mut reference;
        let mut alt;

        let mut b = VariantContext::build_from_vc(&vc1);
        if vc1.is_snp() {
            // we have to repair the first base, so SNP case is special cased
            if vc1.get_reference() == vc2.get_reference() {
                // we've got an insertion, so we just update the alt to have the prev alt
                reference = vc1.alleles.remove(0);
                let mut bases = vc1.alleles.remove(0).bases;
                bases.extend_from_slice(&vc2.get_alternate_alleles()[0].bases[1..]);
                alt = ByteArrayAllele::new(bases.as_slice(), false);
            } else {
                // we're dealing with a deletion, so we patch the ref
                reference = vc2.alleles.remove(0);
                alt = vc1.alleles.remove(1);
                b.loc.end = vc2.loc.end;
            }
        } else {
            // let mut insertion = if vc1.is_simple_deletion() {
            //     vc1
            // } else {
            //     vc2
            // };
            // let mut deletion = if vc1.is_simple_deletion() {
            //     vc2
            // } else {
            //     vc1
            // };
            //
            if vc1.is_simple_insertion() {
                let mut insertion = vc1;
                let mut deletion = vc2;
                reference = deletion.alleles.remove(0);
                alt = insertion.alleles.remove(1);
                b.loc.end = deletion.loc.end;
            } else {
                let mut insertion = vc2;
                let mut deletion = vc1;
                reference = deletion.alleles.remove(0);
                alt = insertion.alleles.remove(1);
                b.loc.end = deletion.loc.end;
            }
        }

        b.get_type();
        b.alleles = vec![reference, alt];
        return b;
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
        reference: &[u8],
        ref_loc: &SimpleInterval,
        max_mnp_distance: usize,
    ) -> BTreeSet<usize> {
        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        let mut hap_number = 0;
        debug!("=== Best Haplotypes ===");

        let mut start_pos_key_set = BTreeSet::new();
        for h in haplotypes.iter_mut() {
            // Walk along the alignment and turn any difference from the reference into an event
            h.event_map = Some(EventMap::new(
                &h,
                reference,
                ref_loc.clone(),
                format!("HC{}", hap_number),
                max_mnp_distance,
            ));
            hap_number += 1;
            start_pos_key_set.par_extend(h.event_map.as_ref().unwrap().get_start_positions());
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

        return start_pos_key_set;
    }

    /**
     * Get the starting positions of events in this event map
     * @return
     */
    pub fn get_start_positions(&self) -> Vec<usize> {
        self.map.keys().cloned().collect::<Vec<usize>>()
    }

    /**
     * Get the variant contexts in order of start position in this event map
     * @return
     */
    pub fn get_variant_contexts<'b>(&'b self) -> Vec<&'b VariantContext> {
        self.map.values().collect::<Vec<&'b VariantContext>>()
    }

    /**
     * Returns any events in the map that overlap loc, including spanning deletions and events that start at loc.
     */
    pub fn get_overlapping_events<'b>(&'b self, loc: usize) -> Vec<&'b VariantContext> {
        let mut overlapping_events = self
            .map
            .range(..=loc)
            .filter(|(k, v)| v.loc.get_end() >= loc)
            .map(|(_, v)| v)
            .collect::<Vec<&'b VariantContext>>();

        let contains_insertion_at_loc = overlapping_events.par_iter().any(|v| {
            v.variant_type.as_ref().unwrap() == &VariantType::Indel
                && v.get_reference().length() == 1
        });

        let deletion_events_ending_at_loc = overlapping_events
            .iter()
            .filter(|v| {
                (v.variant_type.as_ref().unwrap() == &VariantType::Indel
                    && v.get_alternate_alleles()[0].length() == 1)
                    && (v.loc.get_end() == loc)
            })
            .map(|v| &**v)
            .collect::<Vec<&VariantContext>>();
        let contains_deletion_ending_at_loc = deletion_events_ending_at_loc.len() > 0;

        if contains_deletion_ending_at_loc && contains_insertion_at_loc {
            // We are at the end of a deletion and the start of an insertion;
            // only the insertion should be kept in this case.
            return overlapping_events
                .par_iter()
                .filter(|v| *v != &deletion_events_ending_at_loc[0])
                .map(|v| *v)
                .collect::<Vec<&VariantContext>>();
        } else {
            return overlapping_events;
        }
    }
}
