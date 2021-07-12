use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;

/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
pub trait AbstractReadThreadingGraph<'a> {
    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    fn is_threading_start(&self, kmer: &Kmer<'a>, start_threading_only_at_existing_vertex: bool) -> bool;

    // get the next kmerVertex for ChainExtension and validate if necessary.
    fn get_next_kmer_vertex_for_chain_extension(&self, kmer: &Kmer<'a>, is_ref: bool, prev_vertex: &MultiDeBruijnVertex) -> MultiDeBruijnVertex;
}

impl AbstractReadThreadingGraph {

}

enum TraversalDirection {
    Downwards,
    Upwards,
}

/**
 * Keeps track of the information needed to add a sequence to the read threading assembly graph
 */
pub struct SequenceForKmers<'a> {
    pub name: String,
    pub sequence: &'a [u8],
    pub start: usize,
    pub stop: usize,
    pub count: usize,
    pub is_ref: bool
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

