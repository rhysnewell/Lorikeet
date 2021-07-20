use assembly::kmer::Kmer;
use graphs::base_graph::BaseGraph;
use reads::bird_tool_reads::BirdToolRead;
use petgraph::stable_graph::{NodeIndex, EdgeIndex};
use rust_htslib::bam::record::{Cigar, CigarString};
use petgraph::Direction;
use graphs::seq_graph::SeqGraph;
use graphs::base_edge::{BaseEdgeStruct, BaseEdge};
use graphs::base_vertex::BaseVertex;
use utils::simple_interval::Locatable;


/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
pub trait AbstractReadThreadingGraph<'a>: Sized + Send + Sync {
    const ANONYMOUS_SAMPLE: &'static str = "XXX_UNNAMED_XXX";
    const WRITE_GRAPH: bool = false;
    const DEBUG_NON_UNIQUE_CALC: bool = false;
    const MAX_CIGAR_COMPLEXITY: usize = 3;
    const INCREASE_COUNTS_BACKWARDS: bool = true;

    fn get_base_graph<V: BaseVertex, E: BaseEdge>(&self) -> &BaseGraph<V, E>;

    fn get_kmer_size(&self) -> usize;

    fn get_reference_source_vertex(&self) -> Option<NodeIndex>;

    fn get_reference_sink_vertex(&self) -> Option<NodeIndex>;

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
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar                the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    fn cigar_is_okay_to_merge(
        cigar: &CigarString, require_first_element_m: bool, require_last_element_m: bool
    ) -> bool {
        let num_elements = cigar.0.len();

        // don't allow more than a couple of different ops
        if num_elements == 0 || num_elements > Self::MAX_CIGAR_COMPLEXITY {
            return false
        };

        // the first element must be an M
        if require_first_element_m {
            match cigar.0[0] {
                Cigar::Match(_) => {
                    // correct operator but more checks to do
                },
                _ => return false
            }
        };

        // the last element must be an M
        if require_last_element_m {
            match cigar.0[num_elements - 1] {
                Cigar::Match(_) => {
                    // correct operator but more checks to do
                },
                _ => return false
            }
        };

        // note that there are checks for too many mismatches in the dangling branch later in the process
        return true
    }

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
     * calculates the longest suffix match between a sequence and a smaller kmer
     *
     * @param seq      the (reference) sequence
     * @param kmer     the smaller kmer sequence
     * @param seqStart the index (inclusive) on seq to start looking backwards from
     * @return the longest matching suffix
     */
    fn longest_suffix_match(seq: &[u8], kmer: &[u8], seq_start: i64) -> usize {
        for len in 1..(kmer.len() as i64 + 1) {
            let seq_i = seq_start - len + 1;
            let kmer_i = kmer.len() as i64 - len;
            if seq_i < 0 || seq_i[seq_i as usize] != kmer[kmer_i as usize] {
                return len as usize - 1;
            }
        }

        return kmer.len()
    }

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
    ) -> usize;

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_tail(&mut self, dangling_tail_merge_result: DanglingChainMergeHelper) -> usize;

    /**
     * Actually merge the dangling head if possible, this is the old codepath that does not handle indels
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_head_legacy(&mut self, dangling_head_merge_result: DanglingChainMergeHelper) -> usize;

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_head(&mut self, dangling_head_merge_result: DanglingChainMergeHelper) -> usize;

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Requires that at a minimum there are at least #getMinMatchingBases() matches between the reference and the read
     * at the end in order to emit an alignment offset.
     *
     * @param cigarElements cigar elements corresponding to the alignment between path1 and path2
     * @param path1  the first path
     * @param path2  the second path
     * @return an integer pair object where the key is the offset into path1 and the value is offset into path2 (both -1 if no path is found)
     */
    fn best_prefix_match(&self, cigar_elements: &Vec<Cigar>, path1: &[u8], path2: &[u8]) -> (i64, i64);

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1    the first path
     * @param path2    the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    fn best_prefix_match_legacy(&self, path1: &[u8], path2: &[u8], max_index: usize) -> Option<usize>;

    /**
     * NOTE: this method is only used for dangling heads and not tails.
     *
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch the length of the branch itself
     * @return positive integer
     */
    fn get_max_mismatches_legacy(&self, length_of_dangling_branch: usize) -> usize;

    /**
     * The minimum number of matches to be considered allowable for recovering dangling ends
     */
    fn get_min_matching_bases(&self) -> i64;

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
    ) -> Option<DanglingChainMergeHelper>;

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param aligner
     * @param vertex      the source of the dangling head
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param recoverAll  recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    fn generate_cigar_against_upwards_reference_path(
        &self,
        vertex: NodeIndex,
        prune_factor: usize,
        min_dangling_branch_length: usize,
        recover_all: bool
    ) -> Option<DanglingChainMergeHelper>;

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
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex         the original vertex
     * @param pruneFactor    the prune factor to use in ignoring chain pieces
     * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     * has a descendant with multiple outgoing edges before hitting the reference path
     */
    fn find_path_downwards_to_highest_common_desendant_of_reference(
        &self,
        vertex: NodeIndex,
        prune_factor: usize,
        give_up_at_branch: bool,
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

    fn get_heaviest_edge(&self, v: NodeIndex, direction: Direction) -> EdgeIndex;

    /**
     * The base sequence for the given path.
     *
     * @param path         the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return non-null sequence of bases corresponding to the given path
     */
    fn get_bases_for_path(&self, path: &Vec<NodeIndex>, expand_source: bool) -> &[u8];

    fn extend_dangling_path_against_reference(
        &mut self,
        dangling_head_merge_result: &DanglingChainMergeHelper,
        num_nodes_to_extend: usize,
    ) -> bool;

    fn remove_paths_not_connected_to_ref(&mut self);

    fn to_sequence_graph(&self) -> SeqGraph<BaseEdgeStruct>;

    // Method that will be called immediately before haplotype finding in the event there are
    // alteations that must be made to the graph based on implementation
    fn post_process_for_haplotype_finding<L: Locatable>(
        &mut self,
        debug_graph_output_path: String,
        ref_haplotype: &L
    );
}

pub enum TraversalDirection {
    Downwards,
    Upwards,
}

/**
 * Keeps track of the information needed to add a sequence to the read threading assembly graph
 */
#[derive(Debug, Clone)]
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
    pub(crate) dangling_path: Vec<NodeIndex>,
    pub(crate) reference_path: Vec<NodeIndex>,
    pub(crate) dangling_path_string: &'a [u8],
    pub(crate) reference_path_string: &'a [u8],
    pub(crate) cigar: CigarString
}

impl<'a> DanglingChainMergeHelper<'a> {
    pub fn new(
        dangling_path: Vec<NodeIndex>,
        reference_path: Vec<NodeIndex>,
        dangling_path_string: &'a [u8],
        reference_path_string: &'a [u8],
        cigar: CigarString
    ) -> DanglingChainMergeHelper<'a> {
        DanglingChainMergeHelper {
            dangling_path,
            reference_path,
            dangling_path_string,
            reference_path_string,
            cigar,
        }
    }
}

