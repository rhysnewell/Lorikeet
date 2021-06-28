use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;

pub struct AssemblyRegionWalker {
    evaluator: HaplotypeCallerEngine,
}

impl AssemblyRegionWalker {
    pub fn start(
        args: &clap::ArgMatches,
        ref_idx: usize,
        samples: Vec<String>,
        do_allele_specific_calcs: bool,
        sample_ploidy: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        indexed_bam_readers: &Vec<String>,
        m: &clap::ArgMatches,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        tree: &Arc<Mutex<Vec<&Elem>>>,
    ) -> AssemblyRegionWalker {

        let mut hc_engine = HaplotypeCallerEngine::new(
            args,
            ref_idx,
            indexed_bam_readers.clone(),
            false,
            args.value_of("ploidy").unwrap().parse().unwrap()
        );

        hc_engine.apply(
            &indexed_bam_readers,
            short_read_bam_count,
            long_read_bam_count,
            n_threads,
            ref_idx,
            args,
            genomes_and_contigs,
            &concatenated_genomes,
            flag_filters,
            tree
        );
    }
}