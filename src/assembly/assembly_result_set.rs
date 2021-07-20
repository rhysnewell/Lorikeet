use assembly::assembly_result::AssemblyResult;
use std::collections::BTreeSet;
use haplotype::haplotype::Haplotype;
use utils::simple_interval::{Locatable, SimpleInterval};
use assembly::assembly_region::AssemblyRegion;
use model::variant_context::VariantContext;
use linked_hash_map::LinkedHashMap;
use linked_hash_set::LinkedHashSet;
use graphs::base_edge::BaseEdge;
use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;

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
pub struct AssemblyResultSet<'a, E: BaseEdge, L: Locatable, A: AbstractReadThreadingGraph<'a>> {
    assembly_result_by_kmer_size: LinkedHashMap<usize, AssemblyResult<'a, E, L, A>>,
    haplotypes: LinkedHashSet<Haplotype<'a, L>>,
    assembly_result_by_haplotype: LinkedHashMap<Haplotype<'a, L>, AssemblyResult<'a, E, L, A>>,
    region_for_genotyping: &'a AssemblyRegion,
    full_reference_with_padding: &'a [u8],
    padded_reference_loc: &'a SimpleInterval,
    variation_present: bool,
    ref_haplotype: Haplotype<'a, L>,
    kmer_sizes: BTreeSet<usize>,
    variation_events: BTreeSet<VariantContext>,
    last_max_mnp_distance_used: Option<usize>,
}


impl<'a, E: BaseEdge, L: Locatable, A: AbstractReadThreadingGraph<'a>>
    AssemblyResultSet<'a, E, L, A> {
    /**
     * Constructs a new empty assembly result set.
     */
    pub fn new(
        assembly_region: &'a AssemblyRegion,
        full_reference_with_padding: &'a [u8],
        ref_loc: &'a SimpleInterval,
        ref_haplotype: Haplotype<'a, L>
    ) -> AssemblyResultSet<'a, E, L, A> {
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

    /**
     * Adds simultaneously a haplotype and the generating assembly-result.
     *
     * <p>
     *     Haplotypes and their assembly-result can be added multiple times although just the first call will have
     *     any effect (return value is {@code true}).
     * </p>
     *
     *
     * @param h haplotype to add.
     * @param ar assembly-result that is assumed to have given rise to that haplotype.
     *
     * @throws NullPointerException if {@code h} or {@code ar} is {@code null}.
     * @throws IllegalArgumentException if {@code h} has not defined genome location.
     *
     * @return {@code true} iff this called changes the assembly result set.
     */
    pub fn add_with_assembly_result(
        &mut self,
        h: Haplotype<'a, L>,
        ar: AssemblyResult<'a, E, L, A>
    ) -> bool {

    }

    /**
     * Add a assembly-result object.
     *
     * @param ar the assembly result to add.
     *
     * @throws NullPointerException if {@code ar} is {@code null}.
     * @throws IllegalStateException if there is an assembly result with the same kmerSize.
     * @return {@code true} iff this addition changed the assembly result set.
     */
    fn add_assembly_result(ar: AssemblyResult<'a, E, L, A>) {

    }
}
