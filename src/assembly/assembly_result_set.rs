use assembly::assembly_region::AssemblyRegion;
use assembly::assembly_result::AssemblyResult;
use graphs::base_edge::BaseEdge;
use haplotype::haplotype::Haplotype;
use linked_hash_map::LinkedHashMap;
use linked_hash_set::LinkedHashSet;
use model::variant_context::VariantContext;
use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use std::collections::BTreeSet;
use utils::simple_interval::{Locatable, SimpleInterval};

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
pub struct AssemblyResultSet<'a, L: Locatable, A: AbstractReadThreadingGraph> {
    // kmer size and assembly_results index hashmap
    pub(crate) assembly_result_by_kmer_size: LinkedHashMap<usize, usize>,
    pub(crate) haplotypes: LinkedHashSet<Haplotype<L>>,
    // haplotype and assembly_results index hashmap
    pub(crate) assembly_result_by_haplotype: LinkedHashMap<Haplotype<L>, usize>,
    pub(crate) region_for_genotyping: AssemblyRegion,
    pub(crate) full_reference_with_padding: &'a [u8],
    pub(crate) padded_reference_loc: SimpleInterval,
    pub(crate) variation_present: bool,
    pub(crate) ref_haplotype: Haplotype<L>,
    pub(crate) kmer_sizes: BTreeSet<usize>,
    pub(crate) variation_events: BTreeSet<VariantContext>,
    pub(crate) last_max_mnp_distance_used: Option<usize>,
    pub(crate) assembly_results: Vec<AssemblyResult<L, A>>,
}

impl<'a, L: Locatable, A: AbstractReadThreadingGraph> AssemblyResultSet<'a, L, A> {
    /**
     * Constructs a new empty assembly result set.
     */
    pub fn new(
        assembly_region: AssemblyRegion,
        full_reference_with_padding: &'a [u8],
        ref_loc: SimpleInterval,
        ref_haplotype: Haplotype<L>,
    ) -> AssemblyResultSet<'a, L, A> {
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
            last_max_mnp_distance_used: None,
            assembly_results: Vec::new(),
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
    pub fn add_haplotype(&mut self, h: Haplotype<L>) -> bool {
        if self.haplotypes.contains(&h) {
            return false;
        } else {
            self.haplotypes.insert(h);
            return true;
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
    pub fn add_haplotype_and_assembly_result(&mut self, h: Haplotype<L>, ar: usize) -> bool {
        let assembly_result_addition_return = (self.assembly_results.len() - 1) <= ar;
        if self.haplotypes.contains(&h) {
            let previous_ar = self.assembly_result_by_haplotype.get(&h);
            if previous_ar.is_none() {
                self.assembly_result_by_haplotype
                    .insert(h, self.assembly_results.len() - 1);
                return true;
            } else if previous_ar.unwrap() != &(self.assembly_results.len() - 1) {
                panic!("There is already a different assembly result for the input haplotype")
            } else {
                return assembly_result_addition_return;
            }
        } else {
            if !h.allele.is_ref {
                self.variation_present = true;
            };
            self.haplotypes.insert(h.clone());
            self.assembly_result_by_haplotype
                .insert(h, self.assembly_results.len() - 1);
            return true;
        }
    }

    /**
     * Add a assembly-result object.
     *
     * @param ar the assembly result to add.
     *
     * @throws NullPointerException if {@code ar} is {@code null}.
     * @throws IllegalStateException if there is an assembly result with the same kmerSize.
     * @return {@code usize} return index in assembly result array that this assembly result belongs
     */
    pub fn add_assembly_result(&mut self, ar: AssemblyResult<L, A>) -> usize {
        let kmer_size = ar.get_kmer_size();

        if self.assembly_result_by_kmer_size.contains_key(&kmer_size) {
            if ar
                .ne(&self.assembly_results
                    [*self.assembly_result_by_kmer_size.get(&kmer_size).unwrap()])
            {
                panic!("a different assembly result with the same kmerSize was already added");
            } else {
                return *self.assembly_result_by_kmer_size.get(&kmer_size).unwrap();
            }
        } else {
            self.assembly_results.push(ar);
            self.assembly_result_by_kmer_size
                .insert(kmer_size, self.assembly_results.len() - 1);
            self.kmer_sizes.insert(kmer_size);
            return self.assembly_results.len() - 1;
        }
    }

    // /**
    //  * Given whether a new haplotype that has been already added to {@link #haplotypes} collection is the
    //  * reference haplotype and updates {@link #refHaplotype} accordingly.
    //  *
    //  * <p>
    //  *     This method assumes that the colling code has verified that the haplotype was not already in {@link #haplotypes}
    //  *     I.e. that it is really a new one. Otherwise it will result in an exception if it happen to be a reference
    //  *     haplotype and this has already be set. This is the case even if the new haplotypes and the current reference
    //  *     are equal.
    //  * </p>
    //  *
    //  * @param newHaplotype the new haplotype.
    //  * @throws NullPointerException if {@code newHaplotype} is {@code null}.
    //  * @throws IllegalStateException if there is already a reference haplotype.
    //  */
    // fn update_reference_haplotype(&mut self, new_haplotype: Haplotype<L>) {
    //     if !new_haplotype.allele.is_ref {
    //         // pass
    //     } else {
    //         if self.ref_haplotype
    //     }
    // }
}
