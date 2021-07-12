use linked_hash_map::LinkedHashMap;
use assembly::kmer::Kmer;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::base_graph::BaseGraph;
use graphs::multi_sample_edge::MultiSampleEdge;
use std::collections::HashSet;
use read_threading::abstract_read_threading_graph::{SequenceForKmers, AbstractReadThreadingGraph};
use linked_hash_set::LinkedHashSet;
use reads::bird_tool_reads::BirdToolRead;

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

    // get the next kmerVertex for ChainExtension and validate if necessary.
    fn get_next_kmer_vertex_for_chain_extension(
        &self, kmer: &Kmer<'a>, is_ref: bool, prev_vertex: &MultiDeBruijnVertex
    ) -> MultiDeBruijnVertex {
        let unique_merge_vertex = self.get
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    fn get_kmer_vertex(&self, kmer: &Kmer, allow_ref_source: bool) -> Option<&MultiDeBruijnVertex<'a>> {
        if !allow_ref_source && kmer == &self.ref_source {
            return None
        } else {
            return self.kmer_to_vertex_map.get(kmer)
        }
    }

    /**
     * Since we want to duplicate non-unique kmers in the graph code we must determine what those kmers are
     */
    fn preprocess_reads(&mut self) {
        
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

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    fn build_graph_if_necessary(&mut self) {
        if self.already_built {
            // pass
        } else {
            // Capture the set of non-unique kmers for the given kmer size (if applicable)
            self.pre
        }
    }
}