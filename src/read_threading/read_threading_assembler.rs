pub struct ReadThreadingAssembler {
    kmer_sizes: Vec<usize>,
    dont_increase_kmer_sizes_for_cycles: bool,
    allow_non_unique_kmers_in_ref: bool,
    generate_seq_graph: bool,
    recover_haplotypes_from_edges_not_covered_in_junction_trees: bool,
    num_pruning_samples: i32,
    num_best_haplotypes_per_graph: i32,
    prune_before_cycle_counting: bool,
    remove_paths__not_connected_to_ref: bool,
    just_return_raw_graph: bool,
    recover_dangling_branches: bool,
    recover_all_dangling_branches: bool,
    min_dangling_branch_length: i32,
    min_base_quality_to_use_in_assembly: u8,
    prune_factor: i32,
    min_matching_bases_to_dangling_end_recovery: i32
}

impl ReadThreadingAssembler {
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
        max_allowed_paths_for_read_threading_assembler: usize,
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
        min_matching_bases_to_dangle_end_recovery: usize
    ) -> ReadThreadingAssembler {
        assert!(max_allowed_paths_for_read_threading_assembler >= 1, "num_best_haplotypes_per_graph should be >= 1 but got {}", max_allowed_paths_for_read_threading_assembler);
        kmer_sizes.sort();
        ReadThreadingAssembler {
            kmer_sizes,
            dont_increase_kmer_sizes_for_cycles,
            allow_non_unique_kmers_in_ref,
            num_pruning_samples,
            prune_factor
        }
    }
}