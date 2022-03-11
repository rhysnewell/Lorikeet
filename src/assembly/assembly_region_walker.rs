use activity_profile::activity_profile::Profile;
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use assembly::assembly_region_iterator::AssemblyRegionIterator;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::FlagFilter;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use itertools::Itertools;
use model::variant_context::VariantContext;
use processing::lorikeet_engine::Elem;
use rayon::prelude::*;
use reference::reference_reader::ReferenceReader;
use reference::reference_reader_utils::ReferenceReaderUtils;
use rust_htslib::bcf::{IndexedReader, Read};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::path::Path;
use assembly::assembly_region::AssemblyRegion;

pub struct AssemblyRegionWalker<'c> {
    pub(crate) evaluator: HaplotypeCallerEngine<'c>,
    short_read_bam_count: usize,
    long_read_bam_count: usize,
    ref_idx: usize,
    assembly_region_padding: usize,
    min_assembly_region_size: usize,
    max_assembly_region_size: usize,
    n_threads: u32,
}

impl<'c> AssemblyRegionWalker<'c> {
    pub fn start(
        args: &'c clap::ArgMatches,
        ref_idx: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        indexed_bam_readers: &'c [String],
        n_threads: usize,
    ) -> AssemblyRegionWalker<'c> {
        let mut hc_engine = HaplotypeCallerEngine::new(
            args,
            ref_idx,
            indexed_bam_readers.to_vec(),
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

        AssemblyRegionWalker {
            evaluator: hc_engine,
            short_read_bam_count,
            long_read_bam_count,
            ref_idx,
            assembly_region_padding,
            min_assembly_region_size,
            max_assembly_region_size,
            n_threads: n_threads as u32,
        }
    }

    pub fn collect_shards_low_mem(
        &mut self,
        args: &clap::ArgMatches,
        indexed_bam_readers: &[String],
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        n_threads: usize,
        reference_reader: &mut ReferenceReader,
        output_prefix: &str,
    ) -> Vec<VariantContext> {
        self.evaluator.collect_activity_profile_low_mem(
            indexed_bam_readers,
            self.short_read_bam_count,
            self.long_read_bam_count,
            n_threads,
            self.ref_idx,
            args,
            genomes_and_contigs,
            concatenated_genomes,
            flag_filters,
            reference_reader,
            self.assembly_region_padding,
            self.min_assembly_region_size,
            self.max_assembly_region_size,
            self.short_read_bam_count,
            self.long_read_bam_count,
            args.value_of("max-input-depth").unwrap().parse().unwrap(),
            output_prefix,
        )

    }

    pub fn collect_shards(
        &mut self,
        args: &clap::ArgMatches,
        indexed_bam_readers: &[String],
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        n_threads: usize,
        tree: &Arc<Mutex<Vec<&Elem>>>,
        reference_reader: &mut ReferenceReader,
    ) -> HashMap<usize, Vec<BandPassActivityProfile>> {

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
    pub fn traverse<'a, 'b>(
        &'b mut self,
        mut shards: HashMap<usize, Vec<BandPassActivityProfile>>,
        flag_filters: &'a FlagFilter,
        args: &'a clap::ArgMatches,
        sample_names: &'a [String],
        reference_reader: &'b mut ReferenceReader,
        output_prefix: &'a str
    ) -> Vec<VariantContext> {
        let max_input_depth = args.value_of("max-input-depth").unwrap().parse().unwrap();
        let contexts = shards
            .into_par_iter()
            .flat_map(|(tid, mut activity_profiles)| {

                let ref_idx = &self.ref_idx;
                let n_threads = &self.n_threads;
                let assembly_region_padding = &self.assembly_region_padding;
                let min_assembly_region_size = &self.min_assembly_region_size;
                let max_assembly_region_size = &self.max_assembly_region_size;
                let short_read_bam_count = &self.short_read_bam_count;
                let long_read_bam_count = &self.long_read_bam_count;
                let evaluator = &self.evaluator;
                let reference_reader = &reference_reader;

                activity_profiles
                    .into_par_iter()
                    .flat_map(move |mut activity_profile| {
                        let mut inner_reader = ReferenceReader::new_from_reader_with_tid_and_rid(
                            reference_reader,
                            *ref_idx,
                            tid,
                        );

                        Self::process_shard(
                            activity_profile,
                            flag_filters,
                            args,
                            sample_names,
                            &inner_reader,
                            *n_threads,
                            *assembly_region_padding,
                            *min_assembly_region_size,
                            *max_assembly_region_size,
                            *short_read_bam_count,
                            *long_read_bam_count,
                            &evaluator,
                            max_input_depth,
                            output_prefix,
                        )
                        .into_par_iter()
                    })
            })
            .collect::<Vec<VariantContext>>();
        return contexts;
    }

    pub fn process_shard<'a, 'b>(
        mut shard: BandPassActivityProfile,
        flag_filters: &'a FlagFilter,
        args: &clap::ArgMatches,
        sample_names: &'a [String],
        reference_reader: &ReferenceReader,
        n_threads: u32,
        assembly_region_padding: usize,
        min_assembly_region_size: usize,
        max_assembly_region_size: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        evaluator: &HaplotypeCallerEngine,
        max_input_depth: usize,
        output_prefix: &'a str,
    ) -> Vec<VariantContext> {
        let mut assembly_region_iter = AssemblyRegionIterator::new(sample_names, n_threads);

        let mut pending_regions = shard.pop_ready_assembly_regions(
            assembly_region_padding,
            min_assembly_region_size,
            max_assembly_region_size,
            true, // not used, calculated in function
        );

        let features = args.value_of("features-vcf");
        match features {
            Some(indexed_vcf_reader) => {
                debug!("Attempting to extract features...");

                let contexts = pending_regions
                    .into_par_iter()
                    .flat_map(|mut assembly_region| {
                        let mut reference_reader = reference_reader.clone();
                        let mut evaluator = evaluator.clone();

                        // read in feature variants across the assembly region location
                        let mut feature_variants = retrieve_feature_variants(
                            indexed_vcf_reader,
                            &reference_reader,
                            &assembly_region
                        );

                        if long_read_bam_count > 0 && !args.is_present("do-not-call-svs") {
                            let sv_path = format!("{}/structural_variants.vcf.gz", output_prefix);
                            if Path::new(&sv_path).exists() {
                                // structural variants present so we will add them to feature variants
                                let structural_variants = retrieve_feature_variants(
                                    &sv_path,
                                    &reference_reader,
                                    &assembly_region,
                                );

                                feature_variants.extend(structural_variants);
                            }
                        }

                        debug!("Feature variants {:?}", &feature_variants);

                        assembly_region_iter.fill_next_assembly_region_with_reads(
                            &mut assembly_region,
                            flag_filters,
                            n_threads,
                            short_read_bam_count,
                            long_read_bam_count,
                            max_input_depth,
                            args,
                        );


                        evaluator
                            .call_region(
                                assembly_region,
                                &mut reference_reader,
                                feature_variants,
                                args,
                                sample_names,
                                flag_filters,
                            )
                            .into_par_iter()
                    })
                    .collect::<Vec<VariantContext>>();

                return contexts;
            }
            None => {
                let contexts = pending_regions
                    .into_par_iter()
                    .flat_map(|mut assembly_region| {
                        let mut reference_reader = reference_reader.clone();
                        let mut evaluator = evaluator.clone();

                        let feature_variants = if long_read_bam_count > 0 && !args.is_present("do-not-call-svs") {
                            let sv_path = format!("{}/structural_variants.vcf.gz", output_prefix);
                            if Path::new(&sv_path).exists() {
                                // structural variants present so we will add them to feature variants
                                retrieve_feature_variants(
                                    &sv_path,
                                    &reference_reader,
                                    &assembly_region,
                                )
                            } else {
                                Vec::new()
                            }
                        } else { Vec::new() };

                        debug!("Filling with reads...");
                        assembly_region_iter.fill_next_assembly_region_with_reads(
                            &mut assembly_region,
                            flag_filters,
                            n_threads,
                            short_read_bam_count,
                            long_read_bam_count,
                            max_input_depth,
                            args,
                        );

                        evaluator
                            .call_region(
                                assembly_region,
                                &mut reference_reader,
                                feature_variants,
                                args,
                                sample_names,
                                flag_filters,
                            )
                            .into_par_iter()
                    })
                    .collect::<Vec<VariantContext>>();

                return contexts;
            }
        }
    }
}


fn retrieve_feature_variants(
    indexed_vcf_reader: &str,
    reference_reader: &ReferenceReader,
    assembly_region: &AssemblyRegion
) -> Vec<VariantContext> {

    let mut indexed_vcf_reader =
        VariantContext::retrieve_indexed_vcf_file(indexed_vcf_reader);
    debug!("Retrieved indexed VCF...");

    let vcf_rid = VariantContext::get_contig_vcf_tid(
        indexed_vcf_reader.header(),
        reference_reader
            .retrieve_contig_name_from_tid(assembly_region.get_contig())
            .unwrap(),
    );
    debug!("VCF Rid {:?}", &vcf_rid);

    match vcf_rid {
        Some(rid) => VariantContext::process_vcf_in_region(
            &mut indexed_vcf_reader,
            rid,
            assembly_region.get_start() as u64,
            assembly_region.get_end() as u64,
        ),
        None => Vec::new(),
    }
}