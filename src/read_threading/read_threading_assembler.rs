use graphs::adaptive_chain_pruner::AdaptiveChainPruner;
use graphs::chain_pruner::ChainPruner;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use graphs::multi_sample_edge::MultiSampleEdge;
use assembly::assembly_region::AssemblyRegion;
use haplotype::haplotype::Haplotype;
use utils::simple_interval::{Locatable, SimpleInterval};
use read_error_corrector::read_error_corrector::ReadErrorCorrector;
use utils::smith_waterman_aligner::SmithWatermanAligner;
use assembly::assembly_result_set::AssemblyResultSet;
use reads::read_clipper::ReadClipper;
use reads::bird_tool_reads::BirdToolRead;
use read_threading::abstract_read_threading_graph::{AbstractReadThreadingGraph, SequenceForKmers};
use graphs::seq_graph::SeqGraph;
use graphs::base_edge::BaseEdge;
use assembly::assembly_result::{AssemblyResult, Status};
use read_threading::read_threading_graph::ReadThreadingGraph;

pub struct ReadThreadingAssembler<'a, C: ChainPruner<MultiDeBruijnVertex<'a>, MultiSampleEdge<'a>>> {
    kmer_sizes: Vec<usize>,
    dont_increase_kmer_sizes_for_cycles: bool,
    allow_non_unique_kmers_in_ref: bool,
    generate_seq_graph: bool,
    recover_haplotypes_from_edges_not_covered_in_junction_trees: bool,
    num_pruning_samples: i32,
    num_best_haplotypes_per_graph: i32,
    prune_before_cycle_counting: bool,
    remove_paths_not_connected_to_ref: bool,
    just_return_raw_graph: bool,
    recover_dangling_branches: bool,
    recover_all_dangling_branches: bool,
    min_dangling_branch_length: i32,
    min_base_quality_to_use_in_assembly: u8,
    prune_factor: i32,
    min_matching_bases_to_dangling_end_recovery: i32,
    chain_pruner: C
}

impl<'a, C: ChainPruner<MultiDeBruijnVertex<'a>, MultiSampleEdge<'a>>,
    L: Locatable, R: ReadErrorCorrector, E: BaseEdge> ReadThreadingAssembler<'a, C> {
    const DEFAULT_NUM_PATHS_PER_GRAPH: usize = 128;
    const KMER_SIZE_ITERATION_INCREASE: usize = 10;
    const MAX_KMER_ITERATIONS_TO_ATTEMPT: usize = 6;

    /**
     * If false, we will only write out a region around the reference source
     */
    const PRINT_FILL_GRAPH_FOR_DEBUGGING: bool = true;
    const DEFAULT_MIN_BASE_QUALITY_TO_USE: u8 = 10;
    const MIN_HAPLOTYPE_REFERENCE_LENGTH: usize = 30;

    pub fn new(
        max_allowed_paths_for_read_threading_assembler: i32,
        kmer_sizes: Vec<usize>,
        dont_increase_kmer_sizes_for_cycles: bool,
        allow_non_unique_kmers_in_ref: bool,
        num_pruning_samples: usize,
        prune_factor: usize,
        use_adaptive_pruning: bool,
        initial_error_rate_for_pruning: f64,
        pruning_log_odds_threshold: f64,
        pruning_seeding_log_odds_threshold: f64,
        max_unpruned_variants: usize,
        use_linked_debruijn_graphs: bool,
        enable_legacy_graph_cycle_detection: bool,
        min_matching_bases_to_dangle_end_recovery: i32,
    ) -> ReadThreadingAssembler {
        assert!(max_allowed_paths_for_read_threading_assembler >= 1, "num_best_haplotypes_per_graph should be >= 1 but got {}", max_allowed_paths_for_read_threading_assembler);
        kmer_sizes.sort();

        let chain_pruner = AdaptiveChainPruner::new(initial_error_rate_for_pruning, pruning_log_odds_threshold, max_unpruned_variants);

        ReadThreadingAssembler {
            kmer_sizes,
            dont_increase_kmer_sizes_for_cycles,
            allow_non_unique_kmers_in_ref,
            num_pruning_samples,
            prune_factor,
            chain_pruner,
            generate_seq_graph: !use_linked_debruijn_graphs,
            prune_before_cycle_counting: !enable_legacy_graph_cycle_detection,
            remove_paths_not_connected_to_ref: true,
            just_return_raw_graph: false,
            recover_dangling_branches: true,
            recover_all_dangling_branches: false,
            min_dangling_branch_length: 0,
            num_best_haplotypes_per_graph: max_allowed_paths_for_read_threading_assembler,
            min_matching_bases_to_dangling_end_recovery: min_matching_bases_to_dangle_end_recovery,
            recover_haplotypes_from_edges_not_covered_in_junction_trees: true,
            min_base_quality_to_use_in_assembly: Self::DEFAULT_MIN_BASE_QUALITY_TO_USE,
        }
    }

    /**
     * Main entry point into the assembly engine. Build a set of deBruijn graphs out of the provided reference sequence and list of reads
     * @param assemblyRegion              AssemblyRegion object holding the reads which are to be used during assembly
     * @param refHaplotype              reference haplotype object
     * @param fullReferenceWithPadding  byte array holding the reference sequence with padding
     * @param refLoc                    GenomeLoc object corresponding to the reference sequence with padding
     * @param readErrorCorrector        a ReadErrorCorrector object, if read are to be corrected before assembly. Can be null if no error corrector is to be used.
     * @param aligner                   {@link SmithWatermanAligner} used to align dangling ends in assembly graphs to the reference sequence
     * @return                          the resulting assembly-result-set
     */
    pub fn run_local_assembly(
        &mut self,
        assembly_region: &'a AssemblyRegion,
        mut ref_haplotype: Haplotype<'a, L>,
        full_reference_with_padding: &'a [u8],
        ref_loc: &'a SimpleInterval,
        read_error_corrector: Option<R>,
        aligner: &mut SmithWatermanAligner,
    ) -> AssemblyResultSet {
        assert!(full_reference_with_padding.len() == ref_loc.size(), "Reference bases and reference loc must be the same size.");

        // Note that error correction does not modify the original reads,
        // which are used for genotyping TODO this might come before error correction /
        let mut corrected_reads = match read_error_corrector {
            // TODO: Is it possible to perform this
            //      without cloning? Perhaps get_reads() should just return ownership of reads?
            None => assembly_region.get_reads().clone(),
            Some(read_error_corrector) => read_error_corrector.correct_reads(assembly_region.get_reads())
        };

        // Revert clipped bases if necessary (since we do not want to assemble them)
        corrected_reads = corrected_reads.par_iter().map(|read|
                                                         {
                                                             ReadClipper::new(read).hard_clip_soft_clipped_bases()
                                                         }).collect::<Vec<BirdToolRead>>();

        let non_ref_rt_graphs: Vec<AbstractReadThreadingGraph> = Vec::new();
        let non_ref_seq_graphs: Vec<SeqGraph<E>> = Vec::new();
        let active_region_extended_location = assembly_region.get_padded_span();
        ref_haplotype.set_genome_location(active_region_extended_location.clone());
        let mut result_set = AssemblyResultSet::new(assembly_region, full_reference_with_padding, ref_loc, ref_haplotype.clone());
        // either follow the old method for building graphs and then assembling or assemble and haplotype call before expanding kmers
        // NOTE: We are just going to implement the new method here using kmer expansion heuristics
    }

    /**
     * Follow the kmer expansion heurisics as {@link #assemble(List, Haplotype, SAMFileHeader, SmithWatermanAligner)}, but in this case
     * attempt to recover haplotypes from the kmer graph and use them to assess whether to expand the kmer size.
     */
    fn assemble_graphs_and_expand_kmers_given_haplotypes(
        &mut self,
        ref_haplotype: Haplotype<'a, L>,
        ref_loc: &'a SimpleInterval,
        aligner: &mut SmithWatermanAligner,
        corrected_reads: Vec<BirdToolRead>,
        non_ref_rt_graphs: Vec<AbstractReadThreadingGraph>,
        result_set: AssemblyResultSet<L>,
        active_region_extended_location: &SimpleInterval,
    ) {
        let mut assembly_result_by_rt_graph = HashMap::new();

        let mut saved_assembly_results = Vec::new();

        let mut has_adequately_assembled_graph = false;
        let kmers_to_try = self.get_expanded_kmer_list();
        // first, try using the requested kmer sizes
        for i in 0..kmers_to_try.len() {
            let kmer_size = kmers_to_try[i];
            let is_last_cycle = i == kmers_to_try.len() - 1;
            if !has_adequately_assembled_graph {
                assembled_result = self.cre
            }
        }
    }

    /**
     * Method for getting a list of all of the specified kmer sizes to test for the graph including kmer expansions
     * @return
     */
    fn get_expanded_kmer_list(&self) -> Vec<usize> {
        let mut return_list = Vec::new();
        return_list.par_extend(self.kmer_sizes);
        if !self.dont_increase_kmer_sizes_for_cycles {
            let mut kmer_size = self.kmer_sizes.par_iter().max().unwrap() + Self::KMER_SIZE_ITERATION_INCREASE;
            let mut num_iterations = 1;
            while num_iterations <= Self::MAX_KMER_ITERATIONS_TO_ATTEMPT {
                return_list.push(kmer_size);
                kmer_size += Self::KMER_SIZE_ITERATION_INCREASE;
                num_iterations += 1;
            }
        }

        return return_list
    }

    /**
     * Creates the sequence graph for the given kmerSize
     *
     * @param reads            reads to use
     * @param refHaplotype     reference haplotype
     * @param kmerSize         kmer size
     * @param allowLowComplexityGraphs if true, do not check for low-complexity graphs
     * @param allowNonUniqueKmersInRef if true, do not fail if the reference has non-unique kmers
     * @param aligner {@link SmithWatermanAligner} used to align dangling ends to the reference sequence
     * @return sequence graph or null if one could not be created (e.g. because it contains cycles or too many paths or is low complexity)
     */
    fn create_graph(
        &mut self,
        reads: Vec<BirdToolRead>,
        ref_haplotype: &Haplotype<L>,
        kmer_size: usize,
        allow_low_complexity_graphs: bool,
        allow_non_unique_kmers_in_ref: bool,
        aligner: &mut SmithWatermanAligner,
    ) -> Option<AssemblyResult> {
        if ref_haplotype.len() < kmer_size {
            // happens in cases where the assembled region is just too small
            return AssemblyResult::new(Status::Failed, None, None)
        } else if !self.allow_non_unique_kmers_in_ref && !ReadThreadingGraph::determine_non_unique_kmers(
            SequenceForKmers::new("ref".to_string(), ref_haplotype.get_bases(), 0, ref_haplotype.get_bases().len(), 1, true),
            kmer_size
        ).is_empty() {
            return None
        } else {
            let mut rt_graph: AbstractReadThreadingGraph = if self.generate_seq_graph {
                ReadThreadingGraph::new(kmer_size, false, self.min_base_quality_to_use_in_assembly, self.num_pruning_samples, self.min_matching_bases_to_dangling_end_recovery)
            } else {
                // This is where the junction tree debruijn graph would go but considering it is experimental
                // we will leave it out for now
                ReadThreadingGraph::new(kmer_size, false, self.min_base_quality_to_use_in_assembly, self.num_pruning_samples, self.min_matching_bases_to_dangling_end_recovery)
            }
        }
    }
}