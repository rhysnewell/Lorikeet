use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;
use std::collections::HashSet;
use read_threading::abstract_read_threading_graph::{SequenceForKmers, AbstractReadThreadingGraph};
use linked_hash_set::LinkedHashSet;
use reads::bird_tool_reads::BirdToolRead;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use graphs::base_vertex::BaseVertex;
use petgraph::stable_graph::node_index;

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
    kmer_to_vertex_map: LinkedHashMap<Kmer<'a>, NodeIndex>,
    debug_graph_transformations: bool,
    min_base_quality_to_use_in_assembly: u8,
    reference_path: Vec<NodeIndex>,
    already_built: bool,
    // --------------------------------------------------------------------------------
    // state variables, initialized in setToInitialState()
    // --------------------------------------------------------------------------------
    ref_source: Option<Kmer<'a>>,
    start_threading_only_at_existing_vertex: bool,
    max_mismatches_in_dangling_head: i32,
    increase_counts_through_branches: bool,
    num_pruning_samples: usize,
    pub(crate) base_graph: BaseGraph<MultiDeBruijnVertex<'a>, MultiSampleEdge>,
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

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    fn get_all_pending_sequences(&self) -> Vec<&SequenceForKmers> {
        return self.pending.par_iter()
            .flat_map(|one_sample_worth| one_sample_worth.1.iter())
            .collect::<Vec<&SequenceForKmers>>()
    }

    // only add the new kmer to the map if it exists and isn't in our non-unique kmer list
    fn track_kmer(&mut self, kmer: Kmer, new_vertex: NodeIndex) {
        if !self.non_unique_kmers.contains(&kmer) && !self.kmer_to_vertex_map.contains_key(&kmer) {
            self.kmer_to_vertex_map.insert(kmer, new_vertex)
        }
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
        &self, kmer_size: usize, with_non_uniques: Vec<&SequenceForKmers>
    ) -> HashSet<Kmer> {
        // loop over all sequences that have non-unique kmers in them from the previous iterator
        with_non_uniques.par_iter().filter_map(|sequence_for_kmers| {
            let non_uniques_from_seq = Self::determine_non_unique_kmers(sequence_for_kmers, kmer_size);
            if !non_uniques_from_seq.is_empty() {
                // keep track of the non-uniques for this kmerSize, and keep it in the list of sequences that have non-uniques
                non_uniques_from_seq
            }
        }).flat_map(|non_uniques| non_uniques.into_iter()).collect::<HashSet<Kmer>>()
    }
}

impl<'a> AbstractReadThreadingGraph<'a> for ReadThreadingGraph<'a> {
    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    fn is_threading_start(&self, kmer: &Kmer<'a>, start_threading_only_at_existing_vertex: bool) -> bool {
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
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    fn is_low_quality_graph(&self) -> bool {
        return self.non_unique_kmers.len() * 4 > self.kmer_to_vertex_map.len()
    }

    // get the next kmerVertex for ChainExtension and validate if necessary.
    fn get_next_kmer_vertex_for_chain_extension(
        &self, kmer: &Kmer<'a>, is_ref: bool, prev_vertex: NodeIndex
    ) -> Option<&NodeIndex> {
        let unique_merge_vertex = self.get_kmer_vertex(kmer, false);
        assert!(!(is_ref && unique_merge_vertex.is_some()),
                "Did not find a unique vertex to merge into the reference path");

        return unique_merge_vertex
    }

    /**
     * Since we want to duplicate non-unique kmers in the graph code we must determine what those kmers are
     */
    fn preprocess_reads(&mut self) {
        self.non_unique_kmers = self.determine_non_uniques(
            self.base_graph.get_kmer_size(), self.get_all_pending_sequences())
    }

    // whether reads are needed after graph construction
    fn should_remove_reads_after_graph_construction(&self) -> bool {
        return true
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
    fn add_sequence(
        &mut self,
        seq_name: String,
        sample_name: String,
        sequence: &'a [u8],
        start: usize,
        stop: usize,
        count: usize,
        is_ref: bool,
    ) {
        // note that argument testing is taken care of in SequenceForKmers
        // get the list of sequences for this sample
        let sample_sequences = self.pending.entry(&sample_name).or_insert(Vec::new());
        // add the new sequence to the list of sequences for sample
        sample_sequences.push(SequenceForKmers::new(
            seq_name,
            sequence,
            start,
            stop,
            count,
            is_ref
        ))
    }

    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
    fn add_read(&mut self, read: BirdToolRead, sample_names: &Vec<String>) {
        let sequence = read.read.seq().encoded;
        let qualities = read.read.qual();

        let mut last_good = -1;
        for end in 0..(sequence.len() as i32 + 1) {
            if end == sequence.len() || !self.base_is_usable_for_assembly(sequence[end], qualities[end]) {
                // the first good base is at lastGood, can be -1 if last base was bad
                let start = last_good;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                let len = end - start;

                if start != -1 && len >= self.base_graph.get_kmer_size() {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    let name = format!("{}_{}_{}",
                                       std::str::from_utf8(read.read.qname()).unwrap(), start, end);
                    self.add_sequence(
                        name,
                        sample_names[read.sample_index],
                        sequence,
                        start as usize,
                        end as usize,
                        1,
                        false
                    )
                }
                last_good = -1;
            } else if last_good == -1 {
                last_good = end;
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
        return base.to_ascii_uppercase() != 'N' as u8 && qual >= self.min_base_quality_to_use_in_assembly
    }

    fn set_threading_start_only_at_existing_vertex(&mut self, value: bool) {
        self.start_threading_only_at_existing_vertex = value
    }

    fn print_graph(&self, file_name: String, prune_factor: usize) {
        self.base_graph.print_graph(
            &file_name,
            true, 0);
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    fn build_graph_if_necessary(&mut self) {
        if self.already_built {
            // pass
        } else {
            // Capture the set of non-unique kmers for the given kmer size (if applicable)
            self.preprocess_reads();

            debug!("using {} kmer size for this assembly with the following non uniques:", self.base_graph.get_kmer_size());

            // go through the pending sequences, and add them to the graph
            for sequences_for_samples in self.pending.values() {
                for sequence_for_kmers in sequences_for_samples.iter() {
                    self.thread_sequence(sequence_for_kmers);
                    if Self::WRITE_GRAPH {
                        self.base_graph.print_graph(
                            &format!("threading.{}.{}.dot",
                                     self.counter, sequence_for_kmers.name.replace(" ", "_")),
                            true, 0);
                    }
                }

                // flush the single sample edge values from the graph
                for e in self.base_graph.graph.edge_weights_mut() {
                    e.flush_single_sample_multiplicity()
                }
            }

            // clear the pending reads pile to conserve memory
            if self.should_remove_reads_after_graph_construction() {
                self.pending.clear();
            }

            self.already_built = true;
            for (_, v) in self.kmer_to_vertex_map.iter_mut() {
                let mut node = self.base_graph.graph.node_weight_mut(v).unwrap();
                node.set_additional_info(format!("{}+", node.get_additional_info()));
            }
        }
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     *
     * @param seqForKmers a non-null sequence
     */
    fn thread_sequence(&mut self, seq_for_kmers: &SequenceForKmers) {
        let start_pos = self.find_start(seq_for_kmers);

        match start_pos {
            None => {
                // do nothing :)
            },
            Some(start_pos) => {
                let starting_vertex = self.get_or_create_kmer_vertex(seq_for_kmers.sequence, start_pos);
                // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
                if Self::INCREASE_COUNTS_BACKWARDS {
                    self.increase_counts_in_matched_kmers(
                        seq_for_kmers,
                        starting_vertex,
                        self.base_graph.graph.node_weight(starting_vertex).unwrap().get_sequence(),
                        self.base_graph.get_kmer_size().checked_sub(2)
                    )
                }

                if self.debug_graph_transformations {
                    self.base_graph.graph.node_weight_mut(starting_vertex).unwrap().add_read(
                        seq_for_kmers.name
                    );
                }

                // keep track of information about the reference source
                if seq_for_kmers.is_ref {
                    if self.ref_source.is_some() {
                        panic!("Found two ref_sources! prev: {:?}, new: {:?}", self.ref_source, seq_for_kmers)
                    }
                    self.reference_path = Vec::with_capacity(seq_for_kmers.sequence.len() - self.base_graph.get_kmer_size());
                    self.reference_path.push(starting_vertex);
                    self.ref_source = Some(Kmer::new_with_start_and_length(
                        seq_for_kmers.sequence, seq_for_kmers.start, self.base_graph.get_kmer_size()
                    ));
                };

                let mut vertex = starting_vertex;
                for i in start_pos + 1..(seq_for_kmers.stop - self.base_graph.get_kmer_size() + 1) {
                    vertex = self.extend_chain_by_one(
                        vertex, seq_for_kmers.sequence, i, seq_for_kmers.count, seq_for_kmers.is_ref
                    );
                    if seq_for_kmers.is_ref {
                        self.reference_path.push(vertex);
                    }
                    if self.debug_graph_transformations {
                        self.base_graph.graph.node_weight_mut(starting_vertex).unwrap().add_read(
                            seq_for_kmers.name
                        );
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
            return 0
        } else {
            for i in seq_for_kmers.start..(seq_for_kmers.stop - self.base_graph.get_kmer_size()) {
                let kmer1 = Kmer::new_with_start_and_length(seq_for_kmers.sequence, i, self.base_graph.get_kmer_size());
                if self.is_threading_start(&kmer1, self.start_threading_only_at_existing_vertex) {
                    return i
                }
            }

            return None
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
        let kmer = Kmer::new_with_start_and_length(sequence, start, self.base_graph.get_kmer_size());
        let vertex = self.get_kmer_vertex(&kmer, true);
        match vertex {
            None => self.create_vertex(kmer),
            Some(vertex) => vertex
        }
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    fn get_kmer_vertex(&self, kmer: &Kmer, allow_ref_source: bool) -> Option<&NodeIndex> {
        if !allow_ref_source && kmer == &self.ref_source {
            return None
        } else {
            return self.kmer_to_vertex_map.get(kmer)
        }
    }

    /**
     * Create a new vertex for kmer.  Add it to the kmerToVertexMap map if appropriate.
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    fn create_vertex(&mut self, kmer: Kmer) -> NodeIndex {
        let new_vertex = MultiDeBruijnVertex::new(kmer.bases(), false);
        let prev_size = self.base_graph.graph.node_indices().len();
        let node_index = self.base_graph.graph.add_node(new_vertex);

        self.track_kmer(Kmer, node_index);

        return node_index
    }

    fn increase_counts_in_matched_kmers(
        &mut self,
        seqs_for_kmers: &SequenceForKmers,
        vertex: NodeIndex,
        original_kmer: &[u8],
        offset: Option<usize>,
    ) {
        match offset {
            None => {}, // do nothing :|,
            Some(offset) => {
                for edge in self.base_graph.graph.edges_directed(vertex, Direction::Incoming) {
                    let prev = self.base_graph.get_edge_source(edge.id());
                    let suffix = self.base_graph.graph.node_weight(prev).unwrap().get_suffix();
                    let seq_base = original_kmer[offset];

                    if suffix == seq_base && (self.increase_counts_through_branches || self.base_graph.in_degree_of(vertex) == 1) {
                        self.base_graph.graph.edge_weight_mut(edge.id()).unwrap().inc_multiplicity(seqs_for_kmers.count);
                        self.increase_counts_in_matched_kmers(seqs_for_kmers, prev, original_kmer, offset.checked_sub(1));
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
        sequence: &'a [u8],
        kmer_start: usize,
        count: usize,
        is_ref: bool
    ) -> NodeIndex {
        let outgoing_edges = self.base_graph.graph.edges_directed(prev_vertex, Direction::Outgoing);

        let next_pos = kmer_start + self.base_graph.get_kmer_size() - 1;
        for outgoing_edge in outgoing_edges {
            let target = self.base_graph.graph.edge_endpoints(outgoing_edge.id()).unwrap().1;
            if self.base_graph.graph.node_weight(target).unwrap().get_suffix() == sequence[next_pos] {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                self.base_graph.graph.edge_weight_mut(outgoing_edge.id()).unwrap().inc_multiplicity(count);
                return target
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        let kmer = Kmer::new_with_start_and_length(sequence, kmer_start, self.base_graph.get_kmer_size());
        let merge_vertex = self.get_next_kmer_vertex_for_chain_extension(&kmer, is_ref, prev_vertex);

        let next_vertex = if merge_vertex.is_none() { self.create_vertex(kmer) } else { *merge_vertex.unwrap() };

        self.base_graph.graph.add_edge(prev_vertex, next_vertex, MultiSampleEdge::new(is_ref, count, self.num_pruning_samples));

        return next_vertex
    }
}