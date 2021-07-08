use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;

/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
pub struct AbstractReadThreadingGraph<'a> {
    min_matchin_bases_to_dangling_end_recovery: i64,
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

impl AbstractReadThreadingGraph<'_> {
    pub fn new(
        kmer_size: usize,
        debug_graph_transformations: bool,
        min_base_quality_to_use_in_assembly: u8,
        min_pruning_samples: usize,
    ) -> Self {
        let base_graph = BaseGraph::new(kmer_size);
        Self {
            min_matchin_bases_to_dangling_end_recovery: -1,
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
}

enum TraversalDirection {
    Downwards,
    Upwards,
}

/**
 * Keeps track of the information needed to add a sequence to the read threading assembly graph
 */
struct SequenceForKmers<'a> {
    name: String,
    sequence: &'a [u8],
    start: usize,
    stop: usize,
    count: usize,
    is_ref: bool
}

impl SequenceForKmers<'_> {

    /**
     * Create a new sequence for creating kmers
     */
    pub fn new(
        name: String,
        sequence: &'_ [u8],
        start: usize,
        stop: usize,
        count: usize,
        is_ref: bool
    ) -> SequenceForKmers {
        SequenceForKmers {
            name,
            sequence,
            start,
            stop,
            count,
            is_ref
        }
    }
}

