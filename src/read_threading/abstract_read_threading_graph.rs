use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;
use reads::bird_tool_reads::BirdToolRead;
use petgraph::graph::{NodeIndex, EdgeIndex};
use utils::smith_waterman_aligner::SmithWatermanAligner;
use rust_htslib::bam::record::Cigar;

/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
pub trait AbstractReadThreadingGraph<'a>: Sized + Send + Sync {
    const ANONYMOUS_SAMPLE: &'static str = "XXX_UNAMED_XXX";
    const WRITE_GRAPH: bool = false;
    const DEBUG_NON_UNIQUE_CALC: bool = false;
    const MAX_CIGAR_COMPLEXITY: usize = 3;
    const INCREASE_COUNTS_BACKWARDS: bool = true;

    fn has_cycles(&self) -> bool;

    // heuristic to decide if the graph should be thown away based on low complexity/poor assembly
    fn is_low_quality_graph(&self) -> bool;

    fn print_graph(&self, file_name: String, prune_factor: usize);

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    fn is_threading_start(&self, kmer: &Kmer<'a>, start_threading_only_at_existing_vertex: bool) -> bool;

    // get the next kmerVertex for ChainExtension and validate if necessary.
    fn get_next_kmer_vertex_for_chain_extension(
        &self, kmer: &Kmer<'a>, is_ref: bool, prev_vertex: NodeIndex
    ) -> Option<&NodeIndex>;

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

    // only add the new kmer to the map if it exists and isn't in our non-unique kmer list
    fn track_kmer(&mut self, kmer: Kmer, new_vertex: NodeIndex);

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

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     *
     * @param seqForKmers a non-null sequence
     */
    fn thread_sequence(&mut self, seq_for_kmers: &SequenceForKmers);

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    fn find_start(&self, seq_for_kmers: &SequenceForKmers) -> usize;

    // whether reads are needed after graph construction
    fn should_remove_reads_after_graph_construction(&self) -> bool;

    /**
     * Get the vertex for the kmer in sequence starting at start
     *
     * @param sequence the sequence
     * @param start    the position of the kmer start
     * @return a non-null vertex
     */
    fn get_or_create_kmer_vertex(&mut self, sequence: &[u8], start: usize) -> NodeIndex;

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     */
    fn recover_dangling_tails(
        &mut self,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool,
        aligner: &mut SmithWatermanAligner,
    );

    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     */
    fn recover_dangling_heads(
        &mut self,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool,
        aligner: &mut SmithWatermanAligner,
    );

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex                  the vertex to recover
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param aligner
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    fn recover_dangling_tail(
        &mut self,
        vertex: NodeIndex,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool,
        aligner: &mut SmithWatermanAligner,
    ) -> usize;

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex                  the vertex to recover
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    fn recover_dangling_head(
        &mut self,
        vertex: NodeIndex,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool,
        aligner: &mut SmithWatermanAligner,
    ) -> usize;

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param aligner
     * @param vertex      the sink of the dangling chain
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param recoverAll  recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    fn generate_cigar_against_downwards_reference_path(
        &self,
        vertex: NodeIndex,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool,
        aligner: &mut SmithWatermanAligner,
    ) -> DanglingChainMergeHelper;

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex         the original vertex
     * @param pruneFactor    the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     * has an ancestor with multiple incoming edges before hitting the reference path
     */
    fn find_path_upwards_to_lowest_common_ancestor(
        &self, vertex: NodeIndex, prune_factor: usize, give_up_at_branch: bool
    ) -> Option<Vec<NodeIndex>>;

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start           the reference vertex to start from
     * @param direction       describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdge edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    fn get_reference_path(
        &self,
        start: NodeIndex,
        direction: TraversalDirection,
        blacklisted_edge: Option<EdgeIndex>
    ) -> Vec<NodeIndex>;

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    fn get_kmer_vertex(&self, kmer: &Kmer, allow_ref_source: bool) -> Option<NodeIndex>;

    /**
     * Create a new vertex for kmer.  Add it to the kmerToVertexMap map if appropriate.
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    fn create_vertex(&mut self, kmer: Kmer) -> NodeIndex;

    fn increase_counts_in_matched_kmers(
        &mut self,
        seqs_for_kmers: &SequenceForKmers,
        vertex: NodeIndex,
        original_kmer: &[u8],
        offset: usize,
    );

    /**
     * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
     * the graph one bp according to the bases in sequence.
     *
     * @param prevVertex a non-null vertex where sequence was last anchored in the graph
     * @param sequence   the sequence we're threading through the graph
     * @param kmerStart  the start of the current kmer in graph we'd like to add
     * @param count      the number of observations of this kmer in graph (can be > 1 for GGA)
     * @param isRef      is this the reference sequence?
     * @return a non-null vertex connecting prevVertex to in the graph based on sequence
     */
    fn extend_chain_by_one(
        &mut self,
        prev_vertex: NodeIndex,
        sequence: &'a [u8],
        kmer_start: usize,
        count: usize,
        is_ref: bool
    ) -> NodeIndex;

    /**
     * Finds a path starting from a given vertex and satisfying various predicates
     *
     * @param vertex      the original vertex
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param done        test for whether a vertex is at the end of the path
     * @param returnPath  test for whether to return a found path based on its terminal vertex
     * @param nextEdge    function on vertices returning the next edge in the path
     * @param nextNode    function of edges returning the next vertex in the path
     * @return a path, if one satisfying all predicates is found, {@code null} otherwise
     */
    fn find_path(
        &self, vertex: NodeIndex, prune_factor: usize,
        done: &dyn Fn(NodeIndex) -> bool,
        return_path: &dyn Fn(NodeIndex) -> bool,
        next_edge: &dyn Fn(NodeIndex) -> EdgeIndex,
        next_node: &dyn Fn(EdgeIndex) -> NodeIndex,
    ) -> Option<Vec<NodeIndex>>;

    fn has_incident_ref_edge(&self, v: NodeIndex) -> bool;

    fn get_heaviest_incoming_edge(&self, v: NodeIndex) -> EdgeIndex;
}

impl AbstractReadThreadingGraph {

}

pub enum TraversalDirection {
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

/**
 * Class to keep track of the important dangling chain merging data
 */
pub struct DanglingChainMergeHelper<'a> {
    pub(crate) dangling_path: NodeIndex,
    pub(crate) reference_path: NodeIndex,
    pub(crate) dangling_path_string: &'a [u8],
    pub(crate) reference_path_string: &'a [u8],
    pub(crate) cigar: Vec<Cigar>
}

