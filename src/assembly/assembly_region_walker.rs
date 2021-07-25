use activity_profile::activity_profile::Profile;
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use assembly::assembly_region_iterator::AssemblyRegionIterator;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::FlagFilter;
use estimation::lorikeet_engine::Elem;
use graphs::chain_pruner::ChainPruner;
use graphs::multi_sample_edge::MultiSampleEdge;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use model::variant_context::VariantContext;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use reference::reference_reader::ReferenceReader;
use rust_htslib::bcf::{IndexedReader, Read};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

pub struct AssemblyRegionWalker {
    evaluator: HaplotypeCallerEngine,
    features: Option<Result<IndexedReader, rust_htslib::errors::Error>>,
    short_read_bam_count: usize,
    long_read_bam_count: usize,
    ref_idx: usize,
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
        n_threads: usize,
    ) -> AssemblyRegionWalker {
        let mut hc_engine = HaplotypeCallerEngine::new(
            args,
            ref_idx,
            indexed_bam_readers.clone(),
            false,
            args.value_of("ploidy").unwrap().parse().unwrap(),
        );

        let assembly_region_padding = args
            .value_of("assembly-region-padding")
            .unwrap()
            .parse()
            .unwrap();
        let min_assembly_region_size = args
            .value_of("min-assembly-region-size")
            .unwrap()
            .parse()
            .unwrap();
        let max_assembly_region_size = args
            .value_of("max-assembly-region-size")
            .unwrap()
            .parse()
            .unwrap();

        let features_vcf = match args.value_of("features-vcf") {
            Some(vcf_path) => Some(IndexedReader::from_path(vcf_path)),
            None => None,
        };

        AssemblyRegionWalker {
            evaluator: hc_engine,
            features: features_vcf,
            short_read_bam_count,
            long_read_bam_count,
            ref_idx,
            assembly_region_padding,
            min_assembly_region_size,
            max_assembly_region_size,
            n_threads: n_threads as u32,
        }
    }

    pub fn collect_shards(
        &mut self,
        args: &clap::ArgMatches,
        indexed_bam_readers: &Vec<String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        n_threads: usize,
        tree: &Arc<Mutex<Vec<&Elem>>>,
        reference_reader: &mut ReferenceReader,
    ) -> HashMap<usize, BandPassActivityProfile> {
        self.evaluator.collect_activity_profile(
            indexed_bam_readers,
            self.short_read_bam_count,
            self.long_read_bam_count,
            n_threads,
            self.ref_idx,
            args,
            genomes_and_contigs,
            concatenated_genomes,
            flag_filters,
            tree,
            reference_reader,
        )
    }

    /**
     * Iterates through activity profiles per contig, sending each activity profile to be processed
     */
    pub fn traverse(
        &mut self,
        shards: &mut HashMap<usize, BandPassActivityProfile>,
        flag_filters: &FlagFilter,
        args: &clap::ArgMatches,
        sample_names: &Vec<String>,
        reference_reader: &mut ReferenceReader,
    ) {
        for (tid, activity_profile) in shards.iter_mut() {
            self.process_shard(
                activity_profile,
                flag_filters,
                args,
                sample_names,
                reference_reader,
            )
        }
    }

    fn process_shard(
        &mut self,
        shard: &mut BandPassActivityProfile,
        flag_filters: &FlagFilter,
        args: &clap::ArgMatches,
        sample_names: &Vec<String>,
        reference_reader: &mut ReferenceReader,
    ) {
        let mut assembly_region_iter = AssemblyRegionIterator::new(sample_names, self.n_threads);

        let mut pending_regions = shard.pop_ready_assembly_regions(
            self.assembly_region_padding,
            self.min_assembly_region_size,
            self.max_assembly_region_size,
            false,
        );

        match &mut self.features {
            Some(indexed_vcf_reader) => {
                // Indexed readers don't have sync so we cannot parallelize this
                let mut indexed_vcf_reader = indexed_vcf_reader.as_mut().unwrap();
                indexed_vcf_reader.set_threads(self.n_threads as usize);

                for mut assembly_region in pending_regions.into_iter() {
                    assembly_region_iter.fill_next_assembly_region_with_reads(
                        &mut assembly_region,
                        flag_filters,
                        self.n_threads,
                    );

                    let vcf_rid = VariantContext::get_contig_vcf_tid(
                        indexed_vcf_reader.header(),
                        reference_reader
                            .retrieve_contig_name_from_tid(assembly_region.get_contig())
                            .unwrap(),
                    );

                    let feature_variants = match vcf_rid {
                        Some(rid) => VariantContext::process_vcf_in_region(
                            &mut indexed_vcf_reader,
                            rid,
                            assembly_region.get_start() as u64,
                            assembly_region.get_end() as u64,
                        ),
                        None => Vec::new(),
                    };

                    self.evaluator.call_region(
                        assembly_region,
                        reference_reader,
                        &feature_variants,
                        args,
                        sample_names,
                    );
                }
            }
            None => {
                let placeholder = Vec::new();
                for mut assembly_region in pending_regions.into_iter() {
                    assembly_region_iter.fill_next_assembly_region_with_reads(
                        &mut assembly_region,
                        flag_filters,
                        self.n_threads,
                    );

                    self.evaluator.call_region(
                        assembly_region,
                        reference_reader,
                        &placeholder,
                        args,
                        sample_names,
                    );
                }
            }
        }
    }
}
