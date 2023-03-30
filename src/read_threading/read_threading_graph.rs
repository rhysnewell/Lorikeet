use assembly::kmer::Kmer;
use compare::{Compare, Extract};
use gkl::smithwaterman::{OverhangStrategy, Parameters};
use graphs::base_edge::BaseEdge;
use graphs::base_edge::BaseEdgeStruct;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::multi_sample_edge::MultiSampleEdge;
use graphs::seq_graph::SeqGraph;
use hashlink::{LinkedHashMap, LinkedHashSet};
use pair_hmm::pair_hmm_likelihood_calculation_engine::AVXMode;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use rayon::prelude::*;
use read_threading::abstract_read_threading_graph::{
    AbstractReadThreadingGraph, DanglingChainMergeHelper, SequenceForKmers, TraversalDirection,
};
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use reads::alignment_utils::AlignmentUtils;
use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_utils::CigarUtils;
use rust_htslib::bam::record::Cigar;
use smith_waterman::smith_waterman_aligner::{SmithWatermanAligner, STANDARD_NGS};
use std::cmp::{max, min};
use std::collections::HashSet;
use utils::simple_interval::Locatable;

/**
 * Note: not final but only intended to be subclassed for testing.
 */
#[derive(Debug, Clone)]
pub struct ReadThreadingGraph {
    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    non_unique_kmers: HashSet<Kmer>,
    min_matching_bases_to_dangling_end_recovery: i32,
    counter: usize,
    /**
     * Sequences added for read threading before we've actually built the graph
     */
    // pending: LinkedHashMap<usize, Vec<SequenceForKmers>>,
    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    kmer_to_vertex_map: LinkedHashMap<Kmer, NodeIndex>,
    debug_graph_transformations: bool,
    min_base_quality_to_use_in_assembly: u8,
    pub reference_path: Vec<NodeIndex>,
    already_built: bool,
    // --------------------------------------------------------------------------------
    // state variables, initialized in setToInitialState()
    // --------------------------------------------------------------------------------
    ref_source: Option<Kmer>,
    start_threading_only_at_existing_vertex: bool,
    max_mismatches_in_dangling_head: i32,
    increase_counts_through_branches: bool,
    num_pruning_samples: usize,
    pub(crate) base_graph: BaseGraph<MultiDeBruijnVertex, MultiSampleEdge>,
    avx_mode: AVXMode,
}

impl ReadThreadingGraph {
    pub const ANONYMOUS_SAMPLE: &'static str = "XXX_UNNAMED_XXX";
    const WRITE_GRAPH: bool = false;
    const DEBUG_NON_UNIQUE_CALC: bool = false;
    pub const MAX_CIGAR_COMPLEXITY: usize = 3;
    const INCREASE_COUNTS_BACKWARDS: bool = true;

    pub fn new(
        kmer_size: usize,
        debug_graph_transformations: bool,
        min_base_quality_to_use_in_assembly: u8,
        min_pruning_samples: usize,
        min_matching_bases_to_dangling_end_recovery: i32,
        avx_mode: AVXMode,
    ) -> Self {
        let base_graph = BaseGraph::new(kmer_size);
        Self {
            non_unique_kmers: HashSet::new(),
            min_matching_bases_to_dangling_end_recovery,
            counter: 0,
            // pending: LinkedHashMap::new(),
            kmer_to_vertex_map: LinkedHashMap::new(),
            debug_graph_transformations: true,
            min_base_quality_to_use_in_assembly,
            reference_path: Vec::new(),
            already_built: false,
            ref_source: None,
            start_threading_only_at_existing_vertex: false,
            max_mismatches_in_dangling_head: -1,
            increase_counts_through_branches: false,
            num_pruning_samples: min_pruning_samples,
            base_graph,
            avx_mode,
        }
    }

    pub fn default_with_kmer_size(kmer_size: usize) -> Self {
        Self::new(kmer_size, false, 6, 1, -1, AVXMode::detect_mode())
    }

    /**
     * Get the collection of non-unique kmers from sequence for kmer size kmerSize
     * @param seqForKmers a sequence to get kmers from
     * @param kmerSize the size of the kmers
     * @return a non-null collection of non-unique kmers in sequence
     */
    pub fn determine_non_unique_kmers<'a>(
        seq_for_kmers: &SequenceForKmers<'a>,
        kmer_size: usize,
    ) -> Vec<Kmer> {
        // count up occurrences of kmers within each read
        let mut all_kmers = HashSet::new();
        let mut non_unique_kmers = Vec::new();

        let stop_position = seq_for_kmers.stop.checked_sub(kmer_size);
        match stop_position {
            None => {
                // pass
            }
            Some(stop_position) => {
                for i in 0..=stop_position {
                    let kmer =
                        Kmer::new_with_start_and_length(seq_for_kmers.sequence, i, kmer_size);
                    if all_kmers.contains(&kmer) {
                        non_unique_kmers.push(kmer);
                    } else {
                        all_kmers.insert(kmer);
                    };
                }
            }
        }

        return non_unique_kmers;
    }

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    fn get_all_pending_sequences<'a, 'b>(
        pending: &'b LinkedHashMap<usize, Vec<SequenceForKmers<'a>>>,
    ) -> Vec<&'b SequenceForKmers<'a>> {
        return pending
            .values()
            .flat_map(|one_sample_worth| one_sample_worth)
            .collect::<Vec<&SequenceForKmers>>();
    }

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param kmerSize the kmer size to check for non-unique kmers of
     * @return a non-null NonUniqueResult
     */
    fn determine_non_uniques(
        &self,
        kmer_size: usize,
        with_non_uniques: Vec<&SequenceForKmers>,
    ) -> HashSet<Kmer> {
        // loop over all sequences that have non-unique kmers in them from the previous iterator
        with_non_uniques
            .iter()
            .filter_map(|sequence_for_kmers| {
                let non_uniques_from_seq =
                    Self::determine_non_unique_kmers(sequence_for_kmers, kmer_size);
                if !non_uniques_from_seq.is_empty() {
                    // keep track of the non-uniques for this kmerSize, and keep it in the list of sequences that have non-uniques
                    Some(non_uniques_from_seq)
                } else {
                    None
                }
            })
            .flat_map(|non_uniques| non_uniques.into_iter())
            .collect::<HashSet<Kmer>>()
    }

    // /**
    //  * Returns the pending sequences, swapping out the old values for an empty hashmap. This
    //  * effectively clears out the result and gives ownership of pending the function call.
    //  * Used to avoid cloning in the build_graph_if_necessary()
    //  */
    // pub fn get_pending(&mut self) -> LinkedHashMap<usize, Vec<SequenceForKmers>> {
    //     std::mem::replace(&mut self.pending, LinkedHashMap::new())
    // }
    //
    // fn replace_pending(&mut self, used_pending: &mut LinkedHashMap<usize, Vec<SequenceForKmers>>) {
    //     std::mem::swap(used_pending, &mut self.pending)
    // }

    /**
     * Get the set of non-unique kmers in this graph.  For debugging purposes
     * @return a non-null set of kmers
     */
    pub fn get_non_uniques(&mut self) -> &HashSet<Kmer> {
        &self.non_unique_kmers
    }
}

impl AbstractReadThreadingGraph for ReadThreadingGraph {
    fn get_base_graph(&self) -> &BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {
        &self.base_graph
    }

    fn get_base_graph_mut(&mut self) -> &mut BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {
        &mut self.base_graph
    }

    fn find_kmer(&self, kmer: &Kmer) -> Option<&NodeIndex> {
        self.kmer_to_vertex_map.get(kmer)
    }

    fn set_increase_counts_through_branches(&mut self, increase_counts_through_branches: bool) {
        self.increase_counts_through_branches = increase_counts_through_branches;
    }

    fn set_min_matching_bases_to_dangling_end_recovery(&mut self, min_matching_bases: i32) {
        self.min_matching_bases_to_dangling_end_recovery = min_matching_bases;
    }

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    fn is_threading_start(
        &self,
        kmer: &Kmer,
        start_threading_only_at_existing_vertex: bool,
    ) -> bool {
        if self.start_threading_only_at_existing_vertex {
            self.kmer_to_vertex_map.contains_key(kmer)
        } else {
            !self.non_unique_kmers.contains(kmer)
        }
    }

    fn has_cycles(&self) -> bool {
        self.base_graph.has_cycles()
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity
     * as a situation where the number of non-unique kmers is more than 20%
     * of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    fn is_low_quality_graph(&self) -> bool {
        return self.non_unique_kmers.len() * 4 > self.kmer_to_vertex_map.len();
    }

    // get the next kmerVertex for ChainExtension and validate if necessary.
    fn get_next_kmer_vertex_for_chain_extension(
        &self,
        kmer: &Kmer,
        is_ref: bool,
        prev_vertex: NodeIndex,
    ) -> Option<&NodeIndex> {
        let unique_merge_vertex = self.get_kmer_vertex(kmer, false);
        assert!(
            !(is_ref && unique_merge_vertex.is_some()),
            "Did not find a unique vertex to merge into the reference path"
        );

        return unique_merge_vertex;
    }

    /**
     * Since we want to duplicate non-unique kmers in the graph code we must determine what those kmers are
     */
    fn preprocess_reads<'a>(&mut self, pending: &LinkedHashMap<usize, Vec<SequenceForKmers<'a>>>) {
        debug!(
            "All pending before preprocess {}",
            Self::get_all_pending_sequences(pending).len()
        );
        self.non_unique_kmers = self.determine_non_uniques(
            self.base_graph.get_kmer_size(),
            Self::get_all_pending_sequences(pending),
        );
        debug!(
            "All pending after preprocess {}",
            Self::get_all_pending_sequences(pending).len()
        );
        debug!("Non-uniques {}", self.non_unique_kmers.len());
    }

    // whether reads are needed after graph construction
    fn should_remove_reads_after_graph_construction(&self) -> bool {
        return true;
    }

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
    fn add_sequence<'a>(
        &self,
        pending: &mut LinkedHashMap<usize, Vec<SequenceForKmers<'a>>>,
        seq_name: String,
        sample_index: usize,
        sequence: &'a [u8],
        start: usize,
        stop: usize,
        count: usize,
        is_ref: bool,
    ) {
        // note that argument testing is taken care of in SequenceForKmers
        // get the list of sequences for this sample
        let sample_sequences = pending.entry(sample_index).or_insert(Vec::new());
        // add the new sequence to the list of sequences for sample
        sample_sequences.push(SequenceForKmers::new(
            seq_name, sequence, start, stop, count, is_ref,
        ))
    }

    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
    fn add_read<'a>(
        &mut self,
        read: &'a BirdToolRead,
        sample_names: &[String],
        count: &mut usize,
        pending: &mut LinkedHashMap<usize, Vec<SequenceForKmers<'a>>>,
    ) {
        let sequence = read.seq();
        let qualities = read.read.qual();

        let mut last_good = -1;
        for end in 0..=sequence.len() {
            if end as usize == sequence.len()
                || !self.base_is_usable_for_assembly(sequence[end], qualities[end])
            {
                // the first good base is at lastGood, can be -1 if last base was bad
                let start = last_good;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                let len = end as i32 - start;

                if start != -1 && len >= self.base_graph.get_kmer_size() as i32 {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    let name = format!(
                        "{}_{}_{}",
                        std::str::from_utf8(read.read.qname()).unwrap(),
                        start,
                        end
                    );
                    // debug!("Read {name}");
                    self.add_sequence(
                        pending,
                        name,
                        read.sample_index,
                        sequence,
                        start as usize,
                        end,
                        1,
                        false,
                    );
                    *count += 1;
                }
                last_good = -1;
            } else if last_good == -1 {
                last_good = end as i32;
            }
        }
    }

    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    fn base_is_usable_for_assembly(&self, base: u8, qual: u8) -> bool {
        return base.to_ascii_uppercase() != b'N'
            && qual >= self.min_base_quality_to_use_in_assembly;
    }

    fn set_threading_start_only_at_existing_vertex(&mut self, value: bool) {
        self.start_threading_only_at_existing_vertex = value
    }

    fn print_graph(&self, file_name: String, prune_factor: usize) {
        self.base_graph.print_graph(&file_name, true, 0);
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    fn build_graph_if_necessary<'a>(
        &mut self,
        pending: &mut LinkedHashMap<usize, Vec<SequenceForKmers<'a>>>,
    ) {
        if self.already_built {
            return;
        }

        // Capture the set of non-unique kmers for the given kmer size (if applicable)
        self.preprocess_reads(&pending);

        // let pending = self.get_pending();
        // go through the pending sequences, and add them to the graph
        for (name, sequences_for_samples) in pending.iter() {
            debug!("Sample {} reads {}", *name, sequences_for_samples.len());
            for sequence_for_kmers in sequences_for_samples.iter() {
                self.thread_sequence(sequence_for_kmers);
                if Self::WRITE_GRAPH {
                    self.base_graph.print_graph(
                        &format!(
                            "threading.{}.{}.dot",
                            self.counter,
                            sequence_for_kmers.name.replace(" ", "_")
                        ),
                        true,
                        0,
                    );
                }
                self.counter += 1;
            }
            debug!(
                "Threaded {} Nodes {} Edges {}",
                self.counter,
                self.base_graph.graph.node_count(),
                self.base_graph.graph.edge_count()
            );
            // flush the single sample edge values from the graph
            for e in self.base_graph.graph.edge_weights_mut() {
                e.flush_single_sample_multiplicity()
            }
            debug!(
                "Flushed {} Nodes {} Edges {}",
                self.counter,
                self.base_graph.graph.node_count(),
                self.base_graph.graph.edge_count()
            );
        }

        // self.replace_pending(&mut pending);

        // clear the pending reads pile to conserve memory
        // if self.should_remove_reads_after_graph_construction() {
        //     self.pending.clear();
        // }

        self.already_built = true;
        for (_, v) in self.kmer_to_vertex_map.iter_mut() {
            let mut node = self.base_graph.graph.node_weight_mut(*v).unwrap();
            node.set_additional_info(format!("{}+", node.get_additional_info()));
        }
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     *
     * @param seqForKmers a non-null sequence
     */
    fn thread_sequence(&mut self, seq_for_kmers: &SequenceForKmers) {
        let start_pos_option = self.find_start(seq_for_kmers);

        match start_pos_option {
            None => {
                // do nothing :)
            }
            Some(start_pos) => {
                let starting_vertex =
                    self.get_or_create_kmer_vertex(&seq_for_kmers.sequence, start_pos);
                // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
                if Self::INCREASE_COUNTS_BACKWARDS {
                    // let original_kmer = self.base_graph.graph.node_weight(starting_vertex).unwrap().get_sequence().to_vec();
                    self.increase_counts_in_matched_kmers(
                        seq_for_kmers,
                        starting_vertex,
                        &seq_for_kmers.sequence,
                        self.base_graph.get_kmer_size().checked_sub(2),
                    );
                }

                if self.debug_graph_transformations {
                    self.base_graph
                        .graph
                        .node_weight_mut(starting_vertex)
                        .unwrap()
                        .add_read(seq_for_kmers.name.clone());
                }

                // keep track of information about the reference source
                if seq_for_kmers.is_ref {
                    if self.ref_source.is_some() {
                        panic!(
                            "Found two ref_sources! prev: {:?}, new: {:?}",
                            self.ref_source, seq_for_kmers
                        )
                    }
                    self.reference_path = Vec::with_capacity(
                        seq_for_kmers.sequence.len() - self.base_graph.get_kmer_size(),
                    );
                    self.reference_path.push(starting_vertex);
                    self.ref_source = Some(Kmer::new_with_start_and_length(
                        seq_for_kmers.sequence,
                        seq_for_kmers.start,
                        self.base_graph.get_kmer_size(),
                    ));
                };

                let mut vertex = starting_vertex;

                for i in start_pos + 1..=(seq_for_kmers.stop - self.base_graph.get_kmer_size()) {
                    vertex = self.extend_chain_by_one(
                        vertex,
                        &seq_for_kmers.sequence,
                        i,
                        seq_for_kmers.count,
                        seq_for_kmers.is_ref,
                    );
                    if seq_for_kmers.is_ref {
                        self.reference_path.push(vertex);
                    }
                    if self.debug_graph_transformations {
                        self.base_graph
                            .graph
                            .node_weight_mut(starting_vertex)
                            .unwrap()
                            .add_read(seq_for_kmers.name.clone());
                    }
                }
                // TODO: Here GATK changes reference_path to an unmodifiable list, I'm not sure
                //      if you can change mutability like that in Rust?
            }
        }
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    fn find_start(&self, seq_for_kmers: &SequenceForKmers) -> Option<usize> {
        if seq_for_kmers.is_ref {
            return Some(0);
        } else {
            let stop = seq_for_kmers.stop - self.base_graph.get_kmer_size();

            for i in seq_for_kmers.start..stop {
                let kmer1 = Kmer::new_with_start_and_length(
                    seq_for_kmers.sequence,
                    i,
                    self.base_graph.get_kmer_size(),
                );
                if self.is_threading_start(&kmer1, self.start_threading_only_at_existing_vertex) {
                    return Some(i);
                }
            }

            return None;
        }
    }

    /**
     * Get the vertex for the kmer in sequence starting at start
     *
     * @param sequence the sequence
     * @param start    the position of the kmer start
     * @return a non-null vertex
     */
    fn get_or_create_kmer_vertex(&mut self, sequence: &[u8], start: usize) -> NodeIndex {
        let kmer =
            Kmer::new_with_start_and_length(sequence, start, self.base_graph.get_kmer_size());
        let vertex = self.get_kmer_vertex(&kmer, true);
        match vertex {
            None => self.create_vertex(sequence, kmer),
            Some(vertex) => *vertex,
        }
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    fn get_kmer_vertex(&self, kmer: &Kmer, allow_ref_source: bool) -> Option<&NodeIndex> {
        if !allow_ref_source {
            match &self.ref_source {
                None => {
                    return self.kmer_to_vertex_map.get(kmer);
                }
                Some(ref_source) => {
                    if kmer == ref_source {
                        return None;
                    } else {
                        return self.kmer_to_vertex_map.get(kmer);
                    }
                }
            }
        } else {
            return self.kmer_to_vertex_map.get(kmer);
        }
    }

    /**
     * Create a new vertex for kmer.  Add it to the kmerToVertexMap map if appropriate.
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    fn create_vertex(&mut self, sequence: &[u8], mut kmer: Kmer) -> NodeIndex {
        let new_vertex = MultiDeBruijnVertex::new(kmer.bases(sequence).to_vec(), false);
        let node_index = self.base_graph.add_node(&new_vertex);

        self.track_kmer(kmer, node_index);

        return node_index;
    }

    // only add the new kmer to the map if it exists and isn't in our non-unique kmer list
    fn track_kmer(&mut self, kmer: Kmer, new_vertex: NodeIndex) {
        if !self.non_unique_kmers.contains(&kmer) && !self.kmer_to_vertex_map.contains_key(&kmer) {
            self.kmer_to_vertex_map.insert(kmer, new_vertex);
        };
    }

    fn increase_counts_in_matched_kmers(
        &mut self,
        seqs_for_kmers: &SequenceForKmers,
        vertex: NodeIndex,
        original_kmer: &[u8],
        offset: Option<usize>,
    ) {
        match offset {
            None => {} // do nothing :|,
            Some(offset) => {
                let outgoing_edges = self
                    .base_graph
                    .graph
                    .edges_directed(vertex, Direction::Incoming)
                    .map(|e| e.id())
                    .collect::<Vec<EdgeIndex>>();

                for edge in outgoing_edges {
                    let prev = self.base_graph.get_edge_source(edge);
                    let suffix = self
                        .base_graph
                        .graph
                        .node_weight(prev)
                        .unwrap()
                        .get_suffix();
                    let seq_base = original_kmer[offset];

                    if suffix == seq_base
                        && (self.increase_counts_through_branches
                            || self.base_graph.in_degree_of(vertex) == 1)
                    {
                        self.base_graph
                            .graph
                            .edge_weight_mut(edge)
                            .unwrap()
                            .inc_multiplicity(seqs_for_kmers.count);
                        self.increase_counts_in_matched_kmers(
                            seqs_for_kmers,
                            prev,
                            original_kmer,
                            offset.checked_sub(1),
                        );
                    }
                }
            }
        }
    }

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
        sequence: &[u8],
        kmer_start: usize,
        count: usize,
        is_ref: bool,
    ) -> NodeIndex {
        let outgoing_edges = self
            .base_graph
            .graph
            .edges_directed(prev_vertex, Direction::Outgoing);

        let next_pos = kmer_start + self.base_graph.get_kmer_size() - 1;
        let mut found_match = None;
        let mut return_target = None;
        for outgoing_edge in outgoing_edges {
            let target = self
                .base_graph
                .graph
                .edge_endpoints(outgoing_edge.id())
                .unwrap()
                .1;
            if self
                .base_graph
                .graph
                .node_weight(target)
                .unwrap()
                .get_suffix()
                == sequence[next_pos]
            {
                found_match = Some(outgoing_edge.id());
                return_target = Some(target);
                break;
            }
        }

        if found_match.is_some() {
            // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
            self.base_graph
                .graph
                .edge_weight_mut(found_match.unwrap())
                .unwrap()
                .inc_multiplicity(count);

            return return_target.unwrap();
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        let kmer =
            Kmer::new_with_start_and_length(sequence, kmer_start, self.base_graph.get_kmer_size());
        let merge_vertex =
            self.get_next_kmer_vertex_for_chain_extension(&kmer, is_ref, prev_vertex);

        let next_vertex = if merge_vertex.is_none() {
            self.create_vertex(sequence, kmer)
        } else {
            *merge_vertex.unwrap()
        };

        let edge_index = self.base_graph.graph.add_edge(
            prev_vertex,
            next_vertex,
            MultiSampleEdge::new(is_ref, count, self.num_pruning_samples),
        );
        // debug!("adding edge {:?} : {:?} -> {:?}, ref {} count {}, prune {}", edge_index, prev_vertex, next_vertex, is_ref, count, self.num_pruning_samples);

        return next_vertex;
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_end_sw_parameters: &Parameters,
    ) {
        if !self.already_built {
            panic!("recover_dangling_tails requires that graph already be built")
        }

        // we need to build a list of dangling tails because that process can modify the graph
        // (and otherwise attempt to hold multiple mutable references if we do it while iterating over the vertexes)
        let dangling_tails = self
            .base_graph
            .graph
            .node_indices()
            .filter(|v| self.base_graph.out_degree_of(*v) == 0 && !self.base_graph.is_ref_sink(*v))
            .collect::<Vec<NodeIndex>>();

        let mut attempted = 0;
        let mut n_recovered = 0;
        for v in dangling_tails {
            attempted += 1;
            n_recovered += self.recover_dangling_tail(
                v,
                prune_factor,
                min_dangling_branch_length,
                recover_all,
                dangling_end_sw_parameters,
            )
        }

        debug!("Recovered {} of {} dangling tails", n_recovered, attempted);
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_head_sw_parameters: &Parameters,
    ) {
        if !self.already_built {
            panic!("recover_dangling_tails requires that graph already be built")
        }

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        let dangling_heads = self
            .base_graph
            .graph
            .node_indices()
            .filter(|v| self.base_graph.in_degree_of(*v) == 0 && !self.base_graph.is_ref_source(*v))
            .collect::<Vec<NodeIndex>>();

        // now we can try to recover the dangling heads
        let mut attempted = 0;
        let mut n_recovered = 0;
        for v in dangling_heads {
            attempted += 1;
            n_recovered += self.recover_dangling_head(
                v,
                prune_factor,
                min_dangling_branch_length,
                recover_all,
                dangling_head_sw_parameters,
            )
        }

        debug!("Recovered {} of {} dangling heads", n_recovered, attempted);
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_tail_sw_parameters: &Parameters,
    ) -> usize {
        if self.base_graph.out_degree_of(vertex) != 0 {
            panic!(
                "Attempting to recover a dangling tail for {:?} but it has out-degree > 0",
                vertex
            )
        };

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        let dangling_tail_merge_result = self.generate_cigar_against_downwards_reference_path(
            vertex,
            prune_factor,
            min_dangling_branch_length,
            recover_all,
            dangling_tail_sw_parameters,
        );

        match dangling_tail_merge_result {
            None => return 0,
            Some(dangling_tail_merge_result) => {
                if !Self::cigar_is_okay_to_merge(&dangling_tail_merge_result.cigar, false, true) {
                    return 0;
                } else {
                    self.merge_dangling_tail(dangling_tail_merge_result)
                }
            }
        }
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_head_sw_parameters: &Parameters,
    ) -> usize {
        if self.base_graph.in_degree_of(vertex) != 0 {
            panic!(
                "Attempting to recover a dangling head for {:?} but it has in-degree > 0",
                vertex
            )
        };

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        let dangling_head_merge_result = self.generate_cigar_against_upwards_reference_path(
            vertex,
            prune_factor,
            min_dangling_branch_length,
            recover_all,
            dangling_head_sw_parameters,
        );

        match dangling_head_merge_result {
            None => return 0,
            Some(dangling_head_merge_result) => {
                if !Self::cigar_is_okay_to_merge(&dangling_head_merge_result.cigar, true, false) {
                    return 0;
                } else if self.min_matching_bases_to_dangling_end_recovery >= 0 {
                    self.merge_dangling_head(dangling_head_merge_result)
                } else {
                    self.merge_dangling_head_legacy(dangling_head_merge_result)
                }
            }
        }
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_tail(
        &mut self,
        dangling_tail_merge_result: DanglingChainMergeHelper,
    ) -> usize {
        let elements = &dangling_tail_merge_result.cigar.0;
        let last_element = elements[elements.len() - 1];

        // the last element must be an M
        match last_element {
            Cigar::Match(_) => {
                // correct operator but more checks to do
            }
            _ => panic!("Last cigar element must be Match(_)"),
        };

        let last_ref_index =
            CigarUtils::get_reference_length(&dangling_tail_merge_result.cigar) - 1;
        let matching_suffix = min(
            Self::longest_suffix_match(
                dangling_tail_merge_result.reference_path_string.as_bytes(),
                dangling_tail_merge_result.dangling_path_string.as_bytes(),
                last_ref_index as i64,
            ),
            last_element.len() as usize,
        );

        if self.min_matching_bases_to_dangling_end_recovery >= 0 {
            if matching_suffix < self.min_matching_bases_to_dangling_end_recovery as usize {
                return 0;
            }
        } else if matching_suffix == 0 {
            return 0;
        }

        let alt_index_to_merge = (CigarUtils::get_read_length(&dangling_tail_merge_result.cigar)
            as usize)
            .checked_sub(matching_suffix)
            .unwrap_or(0)
            .checked_sub(1)
            .unwrap_or(0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        let first_element_is_deletion = CigarUtils::cigar_elements_are_same_type(
            &dangling_tail_merge_result.cigar.0[0],
            &Some(Cigar::Del(0)),
        );

        let must_handle_leading_deletion_case = first_element_is_deletion
            && ((&dangling_tail_merge_result.cigar.0[0].len() + matching_suffix as u32)
                == last_ref_index + 1);
        let ref_index_to_merge = last_ref_index as i32 - matching_suffix as i32
            + 1
            + if must_handle_leading_deletion_case {
                1
            } else {
                0
            };

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if ref_index_to_merge <= 0 {
            return 0;
        }

        // it's safe to merge now
        self.base_graph.graph.add_edge(
            dangling_tail_merge_result.dangling_path[alt_index_to_merge],
            dangling_tail_merge_result.reference_path[ref_index_to_merge as usize],
            MultiSampleEdge::new(false, 1, self.num_pruning_samples),
        );

        return 1;
    }

    /**
     * Actually merge the dangling head if possible, this is the old codepath that does not handle indels
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_head_legacy(
        &mut self,
        mut dangling_head_merge_result: DanglingChainMergeHelper,
    ) -> usize {
        let first_element = &dangling_head_merge_result.cigar.0[0];

        // the last element must be an M
        match first_element {
            Cigar::Match(_) => {
                // correct operator but more checks to do
            }
            _ => panic!("First cigar element must be Match(_)"),
        };

        let indexes_to_merge = self.best_prefix_match_legacy(
            dangling_head_merge_result.reference_path_string.as_bytes(),
            dangling_head_merge_result.dangling_path_string.as_bytes(),
            first_element.len() as usize,
        );

        match indexes_to_merge {
            None => return 0,
            Some(indexes_to_merge) => {
                // we can't push back the reference path
                if indexes_to_merge >= dangling_head_merge_result.reference_path.len() - 1 {
                    return 0;
                }

                // but we can manipulate the dangling path if we need to
                let num_nodes_to_extend =
                    indexes_to_merge - dangling_head_merge_result.dangling_path.len() + 2;

                if indexes_to_merge >= dangling_head_merge_result.dangling_path.len()
                    && !self.extend_dangling_path_against_reference(
                        &mut dangling_head_merge_result,
                        num_nodes_to_extend,
                    )
                {
                    return 0;
                }

                // it's safe to merge now
                self.base_graph.graph.add_edge(
                    dangling_head_merge_result.reference_path[indexes_to_merge + 1 as usize],
                    dangling_head_merge_result.dangling_path[indexes_to_merge as usize],
                    MultiSampleEdge::new(false, 1, self.num_pruning_samples),
                );

                return 1;
            }
        }
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1    the first path
     * @param path2    the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    fn best_prefix_match_legacy(
        &self,
        path1: &[u8],
        path2: &[u8],
        max_index: usize,
    ) -> Option<usize> {
        let max_mismatches = self.get_max_mismatches_legacy(max_index);

        let mut mismatches = 0;
        let mut index = 0;
        let mut last_good_index = None;

        while index < max_index {
            if path1[index] != path2[index] {
                mismatches += 1;
                if mismatches > max_mismatches {
                    return None;
                }
                last_good_index = Some(index);
            }
            index += 1;
        }

        // if we got here then we hit the max index
        return last_good_index;
    }

    /**
     * NOTE: this method is only used for dangling heads and not tails.
     *
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch the length of the branch itself
     * @return positive integer
     */
    fn get_max_mismatches_legacy(&self, length_of_dangling_branch: usize) -> usize {
        return if self.max_mismatches_in_dangling_head > 0 {
            self.max_mismatches_in_dangling_head as usize
        } else {
            max(
                1,
                length_of_dangling_branch / self.base_graph.get_kmer_size(),
            )
        };
    }

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    fn merge_dangling_head(
        &mut self,
        mut dangling_head_merge_result: DanglingChainMergeHelper,
    ) -> usize {
        let first_element = &dangling_head_merge_result.cigar.0[0];

        // the last element must be an M
        match first_element {
            Cigar::Match(_) => {
                // correct operator but more checks to do
            }
            _ => panic!("First cigar element must be Match(_)"),
        };

        let indexes_to_merge = self.best_prefix_match(
            &dangling_head_merge_result.cigar.0,
            dangling_head_merge_result.reference_path_string.as_bytes(),
            dangling_head_merge_result.dangling_path_string.as_bytes(),
        );

        if indexes_to_merge.0 <= 0 || indexes_to_merge.1 <= 0 {
            return 0;
        };

        // we can't push back the reference path
        if indexes_to_merge.0 as usize >= dangling_head_merge_result.reference_path.len() - 1 {
            return 0;
        };

        // but we can manipulate the dangling path if we need to
        let num_nodes_to_extend =
            indexes_to_merge.1 - dangling_head_merge_result.dangling_path.len() as i64 + 2;
        if indexes_to_merge.1 >= dangling_head_merge_result.dangling_path.len() as i64
            && !self.extend_dangling_path_against_reference(
                &mut dangling_head_merge_result,
                (num_nodes_to_extend) as usize,
            )
        {
            return 0;
        }

        self.base_graph.graph.add_edge(
            dangling_head_merge_result.reference_path[(indexes_to_merge.0 + 1) as usize],
            dangling_head_merge_result.dangling_path[(indexes_to_merge.1) as usize],
            MultiSampleEdge::new(false, 1, self.num_pruning_samples),
        );

        return 1;
    }

    fn extend_dangling_path_against_reference(
        &mut self,
        dangling_head_merge_result: &mut DanglingChainMergeHelper,
        num_nodes_to_extend: usize,
    ) -> bool {
        let index_of_last_dangling_node = dangling_head_merge_result.dangling_path.len() - 1;
        let offset_for_ref_end_to_dangling_end = &dangling_head_merge_result
            .cigar
            .0
            .iter()
            .map(|ce| {
                let a = if CigarUtils::cigar_consumes_reference_bases(ce) {
                    ce.len() as i64
                } else {
                    0
                };

                let b = if CigarUtils::cigar_consumes_read_bases(ce) {
                    ce.len() as i64
                } else {
                    0
                };

                a - b
            })
            .sum::<i64>();

        let index_of_ref_node_to_use = index_of_last_dangling_node as i64
            + offset_for_ref_end_to_dangling_end
            + num_nodes_to_extend as i64;

        if index_of_ref_node_to_use >= dangling_head_merge_result.reference_path.len() as i64
            || index_of_ref_node_to_use < 0
        {
            return false;
        }

        let dangling_source = dangling_head_merge_result
            .dangling_path
            .remove(index_of_last_dangling_node);
        let ref_source_sequence = self
            .base_graph
            .graph
            .node_weight(
                dangling_head_merge_result.reference_path[index_of_ref_node_to_use as usize],
            )
            .unwrap()
            .get_sequence();

        let mut sequence_to_extend = Vec::new();
        for i in 0..num_nodes_to_extend {
            sequence_to_extend.push(ref_source_sequence[i])
        }
        sequence_to_extend.extend_from_slice(
            self.base_graph
                .graph
                .node_weight(dangling_source)
                .unwrap()
                .get_sequence(),
        );
        // let mut sequence_to_extend = sb.as_bytes();

        // clean up the source and edge
        let source_edge = self.get_heaviest_edge(dangling_source, Direction::Outgoing);
        let mut prev_v = self.base_graph.get_edge_target(source_edge);
        let removed_edge = self.base_graph.graph.remove_edge(source_edge).unwrap();

        // extend the path
        // TODO: Make sure there are tests for this
        for i in (1..=num_nodes_to_extend).rev() {
            let new_v = MultiDeBruijnVertex::new_with_sequence(
                sequence_to_extend[i..i + self.base_graph.get_kmer_size()].to_vec(),
            );
            let new_v_index = self.base_graph.add_node(&new_v);
            let new_e =
                MultiSampleEdge::new(false, removed_edge.multiplicity, self.num_pruning_samples);
            let new_e_index = self.base_graph.graph.add_edge(new_v_index, prev_v, new_e);
            dangling_head_merge_result.dangling_path.push(new_v_index);
            prev_v = new_v_index;
        }

        return true;
    }

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
    fn best_prefix_match(
        &self,
        cigar_elements: &Vec<Cigar>,
        path1: &[u8],
        path2: &[u8],
    ) -> (i64, i64) {
        let min_matching_bases = self.get_min_matching_bases() as i64;

        let mut ref_idx = cigar_elements
            .iter()
            .map(|cigar| {
                if CigarUtils::cigar_consumes_reference_bases(cigar) {
                    cigar.len() as i64
                } else {
                    0
                }
            })
            .sum::<i64>()
            - 1;
        let mut read_idx = path2.len() as i64 - 1;

        // NOTE: this only works when the last cigar element has a sufficient number of M bases,
        // so no indels within min-mismatches of the edge.
        'cigarLoop: for ce in cigar_elements.iter().rev() {
            if !(CigarUtils::cigar_consumes_reference_bases(ce)
                && CigarUtils::cigar_consumes_read_bases(ce))
            {
                break;
            } else {
                for j in 0..(ce.len()) {
                    if path1[ref_idx as usize] != path2[read_idx as usize] {
                        break 'cigarLoop;
                    };
                    ref_idx -= 1;
                    read_idx -= 1;
                    if ref_idx < 0 || read_idx < 0 {
                        break 'cigarLoop;
                    }
                }
            }
        }

        let matches = path2.len() as i64 - 1 - read_idx;
        return if matches < min_matching_bases {
            (-1, -1)
        } else {
            (ref_idx, read_idx)
        };
    }

    /**
     * The minimum number of matches to be considered allowable for recovering dangling ends
     */
    fn get_min_matching_bases(&self) -> i32 {
        self.min_matching_bases_to_dangling_end_recovery
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_tail_sw_parameters: &Parameters,
    ) -> Option<DanglingChainMergeHelper> {
        let min_tail_path_length = max(1, min_dangling_branch_length); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        let alt_path =
            self.find_path_upwards_to_lowest_common_ancestor(vertex, prune_factor, !recover_all);
        match alt_path {
            None => return None,
            Some(ref alt_path) => {
                if self.base_graph.is_ref_source(alt_path[0])
                    || (alt_path.len() as i32) < min_tail_path_length + 1
                {
                    return None;
                }
            }
        };

        let alt_path = alt_path.unwrap();
        // now get the reference path from the LCA
        let ref_path = self.get_reference_path(
            alt_path[0],
            TraversalDirection::Downwards,
            Some(self.get_heaviest_edge(alt_path[1], Direction::Incoming)),
        );

        // create the Smith-Waterman strings to use
        let ref_bases = self.get_bases_for_path(&ref_path, false);
        let alt_bases = self.get_bases_for_path(&alt_path, false);

        // run Smith-Waterman to determine the best alignment
        // (and remove trailing deletions since they aren't interesting)
        let alignment = SmithWatermanAligner::align(
            ref_bases.as_bytes(),
            alt_bases.as_bytes(),
            dangling_tail_sw_parameters,
            OverhangStrategy::LeadingInDel,
            self.avx_mode,
        );

        return Some(DanglingChainMergeHelper::new(
            alt_path,
            ref_path,
            alt_bases,
            ref_bases,
            AlignmentUtils::remove_trailing_deletions(alignment.cigar),
        ));
    }

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
        min_dangling_branch_length: i32,
        recover_all: bool,
        dangling_head_sw_parameters: &Parameters,
    ) -> Option<DanglingChainMergeHelper> {
        // find the highest common descendant path between vertex and the reference source if available
        let alt_path = self.find_path_downwards_to_highest_common_desendant_of_reference(
            vertex,
            prune_factor,
            !recover_all,
        );

        match alt_path {
            None => return None,
            Some(ref alt_path) => {
                if self.base_graph.is_ref_sink(alt_path[0])
                    || (alt_path.len() as i32) < min_dangling_branch_length + 1
                {
                    return None;
                }
            }
        }
        let alt_path = alt_path.unwrap();

        // now get the reference path from the LCA
        let ref_path = self.get_reference_path(alt_path[0], TraversalDirection::Upwards, None);

        // create the Smith-Waterman strings to use
        let ref_bases = self.get_bases_for_path(&ref_path, true);
        let alt_bases = self.get_bases_for_path(&alt_path, true);

        // run Smith-Waterman to determine the best alignment
        // (and remove trailing deletions since they aren't interesting)
        let alignment = SmithWatermanAligner::align(
            ref_bases.as_bytes(),
            alt_bases.as_bytes(),
            dangling_head_sw_parameters,
            OverhangStrategy::LeadingInDel,
            self.avx_mode,
        );

        return Some(DanglingChainMergeHelper::new(
            alt_path,
            ref_path,
            alt_bases,
            ref_bases,
            AlignmentUtils::remove_trailing_deletions(alignment.cigar),
        ));
    }

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
    ) -> Option<Vec<NodeIndex>> {
        return if give_up_at_branch {
            let done = |v: NodeIndex| -> bool {
                self.base_graph.is_reference_node(v) || self.base_graph.out_degree_of(v) != 1
            };

            let return_path = |v: NodeIndex| -> bool { self.base_graph.is_reference_node(v) };

            let next_edge =
                |v: NodeIndex| -> EdgeIndex { self.base_graph.outgoing_edge_of(v).unwrap() };

            let next_node = |e: EdgeIndex| -> NodeIndex { self.base_graph.get_edge_target(e) };

            self.find_path(
                vertex,
                prune_factor,
                &done,
                &return_path,
                &next_edge,
                &next_node,
            )
        } else {
            let done = |v: NodeIndex| -> bool {
                self.base_graph.is_reference_node(v) || self.base_graph.out_degree_of(v) == 0
            };

            let return_path = |v: NodeIndex| -> bool { self.base_graph.is_reference_node(v) };

            let next_edge =
                |v: NodeIndex| -> EdgeIndex { self.get_heaviest_edge(v, Direction::Outgoing) };

            let next_node = |e: EdgeIndex| -> NodeIndex { self.base_graph.get_edge_target(e) };

            self.find_path(
                vertex,
                prune_factor,
                &done,
                &return_path,
                &next_edge,
                &next_node,
            )
        };
    }

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
        blacklisted_edge: Option<EdgeIndex>,
    ) -> Vec<NodeIndex> {
        let mut path = Vec::new();

        let mut v = Some(start);
        while v.is_some() {
            path.push(v.unwrap());
            v = match direction {
                TraversalDirection::Downwards => {
                    self.base_graph
                        .get_next_reference_vertex(v, true, blacklisted_edge)
                }
                TraversalDirection::Upwards => self.base_graph.get_prev_reference_vertex(v),
            };
        }

        return path;
    }

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
        &self,
        vertex: NodeIndex,
        prune_factor: usize,
        give_up_at_branch: bool,
    ) -> Option<Vec<NodeIndex>> {
        return if give_up_at_branch {
            let done = |v: NodeIndex| -> bool {
                self.base_graph.in_degree_of(v) != 1 || self.base_graph.out_degree_of(v) >= 2
            };

            let return_path = |v: NodeIndex| -> bool { self.base_graph.out_degree_of(v) > 1 };

            let next_edge =
                |v: NodeIndex| -> EdgeIndex { self.base_graph.incoming_edge_of(v).unwrap() };

            let next_node = |e: EdgeIndex| -> NodeIndex { self.base_graph.get_edge_source(e) };

            self.find_path(
                vertex,
                prune_factor,
                &done,
                &return_path,
                &next_edge,
                &next_node,
            )
        } else {
            let done = |v: NodeIndex| -> bool {
                self.has_incident_ref_edge(v) || self.base_graph.in_degree_of(v) == 0
            };

            let return_path = |v: NodeIndex| -> bool {
                self.base_graph.out_degree_of(v) > 1 && self.has_incident_ref_edge(v)
            };

            let next_edge =
                |v: NodeIndex| -> EdgeIndex { self.get_heaviest_edge(v, Direction::Incoming) };

            let next_node = |e: EdgeIndex| -> NodeIndex { self.base_graph.get_edge_source(e) };

            self.find_path(
                vertex,
                prune_factor,
                &done,
                &return_path,
                &next_edge,
                &next_node,
            )
        };
    }

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
        &self,
        vertex: NodeIndex,
        prune_factor: usize,
        done: &dyn Fn(NodeIndex) -> bool,
        return_path: &dyn Fn(NodeIndex) -> bool,
        next_edge: &dyn Fn(NodeIndex) -> EdgeIndex,
        next_node: &dyn Fn(EdgeIndex) -> NodeIndex,
    ) -> Option<Vec<NodeIndex>> {
        let mut path = Vec::new();
        let mut visited_nodes = HashSet::new();

        let mut v = vertex;
        while !done(v) {
            let edge = next_edge(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if self
                .base_graph
                .graph
                .edge_weight(edge)
                .unwrap()
                .get_pruning_multiplicity()
                < prune_factor
            {
                // save the previously visited notes to protect us from riding forever 'neath the streets of boston
                visited_nodes.par_extend(path);
                path = Vec::new();
            } else {
                path.insert(0, v);
            }

            v = next_node(edge);
            // Check that we aren't stuck in a loop
            if path.contains(&v) || visited_nodes.contains(&v) {
                error!("Dangling end recovery killed because of a loop (find_path)");
                return None;
            }
        }
        path.insert(0, v);

        return if return_path(v) { Some(path) } else { None };
    }

    fn has_incident_ref_edge(&self, v: NodeIndex) -> bool {
        for edge in self.base_graph.graph.edges_directed(v, Direction::Incoming) {
            if self.base_graph.edge_is_ref(edge.id()) {
                return true;
            }
        }
        return false;
    }

    fn get_heaviest_edge(&self, v: NodeIndex, direction: Direction) -> EdgeIndex {
        let edges = self
            .base_graph
            .graph
            .edges_directed(v, direction)
            .map(|e| e.id())
            .collect::<Vec<EdgeIndex>>();
        return if edges.len() == 1 {
            edges[0]
        } else {
            let comparator = Extract::new(|edge: &EdgeIndex| {
                self.base_graph
                    .graph
                    .edge_weight(*edge)
                    .unwrap()
                    .multiplicity
            });

            *edges
                .iter()
                .max_by(|l, r| comparator.compare(l, r))
                .unwrap()
        };
    }

    /**
     * The base sequence for the given path.
     *
     * @param path         the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return non-null sequence of bases corresponding to the given path
     */
    fn get_bases_for_path(&self, path: &Vec<NodeIndex>, expand_source: bool) -> String {
        assert!(!path.is_empty(), "Path cannot be empty");

        let sb = path
            .par_iter()
            .map(|v| {
                if expand_source && self.base_graph.is_source(*v) {
                    let mut seq = self
                        .base_graph
                        .graph
                        .node_weight(*v)
                        .unwrap()
                        .get_sequence()
                        .to_vec();
                    seq.reverse();
                    std::str::from_utf8(&seq).unwrap().to_string()
                } else {
                    std::str::from_utf8(
                        self.base_graph
                            .graph
                            .node_weight(*v)
                            .unwrap()
                            .get_suffix_as_array(),
                    )
                    .unwrap()
                    .to_string()
                }
            })
            .collect::<String>();

        return sb;
    }

    fn remove_paths_not_connected_to_ref(&mut self) {
        self.base_graph.remove_paths_not_connected_to_ref()
    }

    fn to_sequence_graph(&self) -> SeqGraph<BaseEdgeStruct> {
        self.base_graph.to_sequence_graph()
    }

    fn get_kmer_size(&self) -> usize {
        self.base_graph.get_kmer_size()
    }

    fn get_reference_source_vertex(&self) -> Option<NodeIndex> {
        self.base_graph.get_reference_source_vertex()
    }

    fn get_reference_sink_vertex(&self) -> Option<NodeIndex> {
        self.base_graph.get_reference_sink_vertex()
    }

    // Method that will be called immediately before haplotype finding in the event there are
    // alteations that must be made to the graph based on implementation
    fn post_process_for_haplotype_finding<L: Locatable>(
        &mut self,
        debug_graph_output_path: Option<&String>,
        ref_haplotype: &L,
    ) {
        // Do nothing There is no processing required for this graph so simply return
    }
}
