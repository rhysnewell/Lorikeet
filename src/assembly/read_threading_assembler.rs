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
    pub fn new
}