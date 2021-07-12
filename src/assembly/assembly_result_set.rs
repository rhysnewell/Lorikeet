use assembly::assembly_result::AssemblyResult;
use std::collections::{HashMap, HashSet, BTreeSet};
use haplotype::haplotype::Haplotype;
use utils::simple_interval::{Locatable, SimpleInterval};
use assembly::assembly_region::AssemblyRegion;
use model::variant_context::VariantContext;
use linked_hash_map::LinkedHashMap;
use linked_hash_set::LinkedHashSet;

/**
 * Collection of read assembly using several kmerSizes.
 *
 * <p>
 *     There could be a different assembly per each kmerSize. In turn, haplotypes are result of one of those
 *     assemblies.
 * </p>
 *
 * <p>
 *     Where there is more than one possible kmerSize that generates a haplotype we consider the smaller one.
 * </p>
 *
 * @original_author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 * @author Rhys Newell; rhys.newell@hdr.qut.edu.au; Rust re-implementation
 */
pub struct AssemblyResultSet<'a, L: Locatable> {
    assembly_result_by_kmer_size: LinkedHashMap<usize, AssemblyResult<'a>>,
    haplotypes: LinkedHashSet<Haplotype<L>>,
    assembly_result_by_haplotype: LinkedHashMap<Haplotype<L>, AssemblyResult<'a>>,
    region_for_genotyping: &'a AssemblyRegion,
    full_reference_with_padding: &'a [u8],
    padded_reference_loc: &'a SimpleInterval,
    variation_present: bool,
    ref_haplotype: Haplotype<L>,
    kmer_sizes: BTreeSet<usize>,
    variation_events: BTreeSet<VariantContext>,
    last_max_mnp_distance_used: Option<usize>,
}


impl<'a, L: Locatable> AssemblyResultSet<'a, L> {
    /**
     * Constructs a new empty assembly result set.
     */
    pub fn new(
        assembly_region: &'a AssemblyRegion,
        full_reference_with_padding: &'a [u8],
        ref_loc: &'a SimpleInterval,
        ref_haplotype: Haplotype<L>
    ) -> AssemblyResultSet<'a, L> {
        let mut haplotypes = LinkedHashSet::new();
        haplotypes.insert(ref_haplotype.clone());

        AssemblyResultSet {
            assembly_result_by_kmer_size: LinkedHashMap::new(),
            haplotypes: haplotypes,
            assembly_result_by_haplotype: LinkedHashMap::new(),
            region_for_genotyping: assembly_region,
            full_reference_with_padding: full_reference_with_padding,
            padded_reference_loc: ref_loc,
            variation_present: false,
            ref_haplotype: ref_haplotype,
            kmer_sizes: BTreeSet::new(),
            variation_events: BTreeSet::new(),
            last_max_mnp_distance_used: None
        }
    }

    /**
     * Adds a haplotype to the result set without indicating a generating assembly result.
     *
     * <p>
     *     It is possible to call this method with the same haplotype several times. In that the second and further
     *     calls won't have any effect (thus returning {@code false}).
     * </p>
     *
     * @param h the haplotype to add to the assembly result set.
     *
     * @throws NullPointerException if {@code h} is {@code null}
     * @throws IllegalArgumentException if {@code h} does not have a genome location.
     *
     * @return {@code true} if the assembly result set has been modified as a result of this call.
     */
    pub fn add(&mut self, h: Haplotype<'a, L>) -> bool {
        if self.haplotypes.contains(&h) {
            return false
        } else {
            self.haplotypes.insert(h);
            return true
        }
    }
}
