use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;
use std::collections::HashSet;
use read_threading::abstract_read_threading_graph::SequenceForKmers;
use linked_hash_set::LinkedHashSet;

/**
 * Note: not final but only intended to be subclassed for testing.
 */
pub struct ReadThreadingGraph<'a> {
    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    non_unique_kmers: HashSet<Kmer<'a>>,
    min_matching_bases_to_dangling_end_recovery: i64,
    counter: usize,
    /**
     * Sequences added for read threading before we've actually built the graph
     */
    pending: LinkedHashMap<String, Vec<SequenceForKmers<'a>>>,
    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    kmer_to_vertex_map: LinkedHashMap<Kmer<'a>, MultiDeBruijnVertex<'a>>,
    debug_graph_transformations: bool,
    min_base_quality_to_use_in_assembly: u8,
    reference_path: Vec<MultiDeBruijnVertex<'a>>,
    already_built: bool,
    // --------------------------------------------------------------------------------
    // state variables, initialized in setToInitialState()
    // --------------------------------------------------------------------------------
    ref_source: Option<Kmer<'a>>,
    start_threading_only_at_existing_vertex: bool,
    max_mismatches_in_dangling_head: i32,
    increase_counts_through_branches: bool,
    num_pruning_samples: usize,
    pub base_graph: BaseGraph<MultiDeBruijnVertex<'a>, MultiSampleEdge>,
}

impl<'a> ReadThreadingGraph<'a> {
    pub fn new(
        kmer_size: usize,
        debug_graph_transformations: bool,
        min_base_quality_to_use_in_assembly: u8,
        min_pruning_samples: usize,
        min_matching_bases_to_dangling_end_recovery: i64
    ) -> Self {
        let base_graph = BaseGraph::new(kmer_size);
        Self {
            non_unique_kmers: HashSet::new(),
            min_matching_bases_to_dangling_end_recovery: min_matching_bases_to_dangling_end_recovery,
            counter: 0,
            pending: LinkedHashMap::new(),
            kmer_to_vertex_map: LinkedHashMap::new(),
            debug_graph_transformations: false,
            min_base_quality_to_use_in_assembly: 0,
            reference_path: Vec::new(),
            already_built: false,
            ref_source: None,
            start_threading_only_at_existing_vertex: false,
            max_mismatches_in_dangling_head: -1,
            increase_counts_through_branches: false,
            num_pruning_samples: min_pruning_samples,
            base_graph: base_graph,
        }
    }

    /**
     * Get the collection of non-unique kmers from sequence for kmer size kmerSize
     * @param seqForKmers a sequence to get kmers from
     * @param kmerSize the size of the kmers
     * @return a non-null collection of non-unique kmers in sequence
     */
    pub fn determine_non_unique_kmers(seq_for_kmers: SequenceForKmers, kmer_size: usize) -> Vec<Kmer<'_>> {
        // count up occurrences of kmers within each read
        let mut all_kmers = LinkedHashSet::new();
        let mut non_unique_kmers = Vec::new();
        let stop_position = seq_for_kmers.stop - kmer_size;
        for i in 0..stop_position {
            let kmer = Kmer::new_with_start_and_length(seq_for_kmers.sequence, i, kmer_size);
            if !all_kmers.insert(&kmer) {
                non_unique_kmers.push(kmer)
            }
        }

        return non_unique_kmers
    }
}