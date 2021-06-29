use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use std::collections::HashMap;
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use coverm::genomes_and_contigs::GenomesAndContigs;
use utils::reference_reader_utils::ReferenceReader;
use assembly::assembly_region_iterator::AssemblyRegionIterator;

pub struct AssemblyRegionWalker {
    evaluator: HaplotypeCallerEngine,
    shards: HashMap<usize, BandPassActivityProfile>,
    indexed_bam_readers: Vec<String>,
    short_read_bam_count: usize,
    long_read_bam_count: usize,
    ref_idx: usize,
    reference_reader: ReferenceReader,
    assembly_region_padding: usize,
    min_assembly_region_size: usize,
    max_assembly_region_size: usize,
    n_threads: u32,
}

impl AssemblyRegionWalker {
    pub fn start(
        args: &clap::ArgMatches,
        ref_idx: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        indexed_bam_readers: &Vec<String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        tree: &Arc<Mutex<Vec<&Elem>>>,
        mut reference_reader: ReferenceReader,
        n_threads: usize,
    ) -> AssemblyRegionWalker {

        let mut hc_engine = HaplotypeCallerEngine::new(
            args,
            ref_idx,
            indexed_bam_readers.clone(),
            false,
            args.value_of("ploidy").unwrap().parse().unwrap()
        );

        let shards = hc_engine.collect_activity_profile(
            &indexed_bam_readers,
            short_read_bam_count,
            long_read_bam_count,
            n_threads,
            ref_idx,
            args,
            genomes_and_contigs,
            &concatenated_genomes,
            flag_filters,
            tree,
            &mut reference_reader,
        );

        let assembly_region_padding = args.value_of("assembly-region-padding")
            .unwrap().parse().unwrap();
        let min_assembly_region_size = args.value_of("min-assembly-region-size")
            .unwrap().parse().unwrap();
        let max_assembly_region_size = args.value_of("max-assembly-region-size")
            .unwrap().parse().unwrap();

        AssemblyRegionWalker {
            evaluator: hc_engine,
            shards: shards,
            indexed_bam_readers: *indexed_bam_readers.clone(),
            short_read_bam_count,
            long_read_bam_count,
            ref_idx,
            reference_reader,
            assembly_region_padding,
            min_assembly_region_size,
            max_assembly_region_size,
            n_threads: n_threads as u32,
        }
    }

    pub fn traverse(&mut self) {
        for (tid, activity_profile) in self.shards.iter_mut() {
            self.process_shard(activity_profile)
        }
    }

    fn process_shard(&mut self, shard: &mut BandPassActivityProfile, args: &clap::ArgMatches) {

        let mut assembly_region_iter = AssemblyRegionIterator::new(
            &mut self.reference_reader,
            shard,
            &self.indexed_bam_readers,
            self.assembly_region_padding,
            self.min_assembly_region_size,
            self.max_assembly_region_size,
            self.n_threads,
        );

        for assembly_region in assembly_region_iter.pending_regions.iter_mut() {
            AssemblyRegionIterator::fill_next_assembly_region_with_reads(&mut assembly_region, self.n_threads);
            self.evaluator.call_region(assembly_region, args) // TODO: Implement Call Region
        }
    }
}