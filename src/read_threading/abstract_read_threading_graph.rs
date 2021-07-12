use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;
use reads::bird_tool_reads::BirdToolRead;

/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
pub trait AbstractReadThreadingGraph<'a> {
    const ANONYMOUS_SAMPLE: &'static str = "XXX_UNAMED_XXX";
    const WRITE_GRAPH: bool = false;
    const DEBUG_NON_UNIQUE_CALC: bool = false;
    const MAX_CIGAR_COMPLEXITY: usize = 3;
    const INCREASE_COUNTS_BACKWARDS: bool = true;

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

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    fn get_kmer_vertex(&self, kmer: &Kmer, allow_ref_source: bool) -> Option<&MultiDeBruijnVertex<'a>>;

    /**
     * Add bases in sequence to this graph
     *
     * @param seqName  a useful seqName for this read, for debugging purposes
     * @param sequence non-null sequence of bases
     * @param start    the first base offset in sequence that we should use for constructing the graph using this sequence, inclusive
     * @param stop     the last base offset in sequence that we should use for constructing the graph using this sequence, exclusive
     * @param count    the representative count of this sequence (to use as the weight)
     * @param isRef    is this the reference sequence.
     */
    fn add_sequence(
        &mut self,
        seq_name: String,
        sample_name: String,
        sequence: &'a [u8],
        start: usize,
        stop: usize,
        count: usize,
        is_ref: bool,
    );

    /**
     * Changes the threading start location policy.
     *
     * @param value {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *              it can start at any unique kmer.
     */
    fn set_threading_start_only_at_existing_vertex(&mut self, value: bool);

    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
    fn add_read(&mut self, read: BirdToolRead);

    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    fn base_is_usable_for_assembly(&self, base: u8, qual: u8) -> bool;

    /**
     * Since we want to duplicate non-unique kmers in the graph code we must determine what those kmers are
     */
    fn preprocess_reads(&mut self);

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    fn build_graph_if_necessary(&mut self);
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

