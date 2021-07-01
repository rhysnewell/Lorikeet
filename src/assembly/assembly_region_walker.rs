use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use std::collections::HashMap;
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use coverm::genomes_and_contigs::GenomesAndContigs;
use assembly::assembly_region_iterator::AssemblyRegionIterator;
use rust_htslib::bcf::{IndexedReader, Read};
use model::variant_context::VariantContext;
use coverm::FlagFilter;
use reference::reference_reader::ReferenceReader;

pub struct AssemblyRegionWalker {
    evaluator: HaplotypeCallerEngine,
    features: Option<IndexedReader>,
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

        let features_vcf = match args.value_of("features-vcf") {
            Some(vcf_path) => Some(IndexedReader::from_path(vcf_path)),
            None => None,
        };

        AssemblyRegionWalker {
            evaluator: hc_engine,
            shards: shards,
            indexed_bam_readers: *indexed_bam_readers.clone(),
            features: features_vcf,
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

    /**
    * Iterates through activity profiles per contig, sending each activity profile to be processed
    */
    pub fn traverse(
        &mut self,
        flag_filters: &FlagFilter,
        args: &clap::ArgMatches
    ) {
        for (tid, activity_profile) in self.shards.iter_mut() {
            self.process_shard(activity_profile, flag_filters, args)
        }
    }

    fn process_shard(
        &mut self,
        shard: &mut BandPassActivityProfile,
        flag_filters: &FlagFilter,
        args: &clap::ArgMatches
    ) {

        let mut assembly_region_iter = AssemblyRegionIterator::new(
            &mut self.reference_reader,
            shard,
            &self.indexed_bam_readers,
            self.assembly_region_padding,
            self.min_assembly_region_size,
            self.max_assembly_region_size,
            self.n_threads,
        );

        match &mut self.features {
            Some(indexed_vcf_reader) => {
                // Indexed readers don't have sync so we cannot parallelize this
                indexed_vcf_reader.set_threads(self.n_threads as usize);
                for assembly_region in assembly_region_iter.pending_regions.iter_mut() {
                    AssemblyRegionIterator::fill_next_assembly_region_with_reads(
                        &mut assembly_region, flag_filters, self.n_threads
                    );

                    let vcf_rid = VariantContext::get_contig_vcf_tid(
                        indexed_vcf_reader.header(),
                        self.reference_reader.retrieve_contig_name_from_tid(
                            assembly_region.get_contig()
                        )
                    );

                    let feature_variants = match vcf_rid {
                        Some(rid) => {
                            VariantContext::process_vcf_in_region(
                                indexed_vcf_reader,
                                rid,
                                assembly_region.get_start() as u64,
                                assembly_region.get_end() as u64
                            )
                        },
                        None => {
                            Vec::new()
                        }
                    };

                    self.evaluator.call_region(
                        assembly_region,
                        feature_variants,
                        args
                    ) // TODO: Implement Call Region
                }
            },
            None => {
                assembly_region_iter.pending_regions.iter_mut()
                    .for_each(|assembly_region| {
                        AssemblyRegionIterator::fill_next_assembly_region_with_reads(
                            &mut assembly_region, flag_filters, self.n_threads
                        );

                        self.evaluator.call_region(
                            assembly_region,
                            Vec::new(),
                            args
                        ) // TODO: Implement Call Region
                });
            }
        }
    }
}