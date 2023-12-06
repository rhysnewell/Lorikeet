use bio_types::sequence::SequenceRead;
use indicatif::ProgressStyle;
use mathru::special::gamma::{digamma, ln_gamma};
use ndarray::Array2;
use rayon::prelude::*;
use rust_htslib::bam::{record::Cigar, Record};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Writer};
use std::cmp::{max, min};
use std::collections::{HashSet, HashMap};
use std::fs::{create_dir_all, File};
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::sync::{Arc, Mutex};

use crate::ani_calculator::ani_calculator::ANICalculator;
use crate::annotator::variant_annotation::VariantAnnotations;
use crate::model::allele_likelihoods::AlleleLikelihoods;
use crate::model::byte_array_allele::ByteArrayAllele;
use crate::model::variant_context::VariantContext;
use crate::activity_profile::activity_profile::Profile;
use crate::activity_profile::activity_profile_state::{ActivityProfileState, ActivityProfileDataType};
use crate::activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use crate::annotator::variant_annotator_engine::VariantAnnotationEngine;
use crate::assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use crate::assembly::assembly_region::AssemblyRegion;
use crate::assembly::assembly_region_trimmer::AssemblyRegionTrimmer;
use crate::assembly::assembly_region_walker::AssemblyRegionWalker;
use crate::assembly::assembly_result_set::AssemblyResultSet;
use crate::reference::reference_reader_utils::GenomesAndContigs;
use crate::bam_parsing::{FlagFilter, bam_generator::*};
use crate::genotype::genotype_builder::Genotype;
use crate::genotype::genotype_prior_calculator::GenotypePriorCalculator;
use crate::genotype::genotyping_engine::GenotypingEngine;
use crate::haplotype::haplotype::Haplotype;
use crate::haplotype::haplotype_caller_genotyping_engine::HaplotypeCallerGenotypingEngine;
use crate::haplotype::ref_vs_any_result::RefVsAnyResult;
use crate::processing::lorikeet_engine::{ReadType, Elem};
use crate::read_orientation::beta_distribution_shape::BetaDistributionShape;
use crate::read_threading::read_threading_assembler::ReadThreadingAssembler;
use crate::read_threading::read_threading_graph::ReadThreadingGraph;
use crate::reads::alignment_utils::AlignmentUtils;
use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::cigar_utils::CigarUtils;
use crate::reads::read_utils::ReadUtils;
use crate::reference::reference_reader::ReferenceReader;
use crate::utils::errors::BirdToolError;
use crate::utils::interval_utils::IntervalUtils;
use crate::utils::math_utils::{MathUtils, RunningAverage};
use crate::utils::natural_log_utils::NaturalLogUtils;
use crate::utils::quality_utils::QualityUtils;
use crate::utils::simple_interval::{Locatable, SimpleInterval};
use crate::pair_hmm::pair_hmm_likelihood_calculation_engine::{
    AVXMode, PairHMMLikelihoodCalculationEngine,
};

#[derive(Debug, Clone)]
pub struct HaplotypeCallerEngine {
    active_region_evaluation_genotyper_engine: GenotypingEngine,
    genotyping_engine: HaplotypeCallerGenotypingEngine,
    genotype_prior_calculator: GenotypePriorCalculator,
    assembly_region_trimmer: AssemblyRegionTrimmer,
    assembly_engine: ReadThreadingAssembler,
    likelihood_calculation_engine: PairHMMLikelihoodCalculationEngine,
    ref_idx: usize,
    stand_min_conf: f64,
    mapping_quality_threshold: u8,
}

impl HaplotypeCallerEngine {
    pub const MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION: u8 = 6;
    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    const AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD: f32 = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    // const MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY: f64 = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    const MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY: usize = 2;

    /**
     * Reads with length lower than this number, after clipping off overhangs outside the active region,
     * won't be considered for genotyping.
     */
    const READ_LENGTH_FILTER_THRESHOLD: usize = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    pub const DEFAULT_READ_QUALITY_FILTER_THRESHOLD: usize = 20;

    /**
     * Surrogate quality score for no base calls.
     * <p>
     * This is the quality assigned to deletion (so without its own base-call quality) pile-up elements,
     * when assessing the confidence on the hom-ref call at that site.
     * </p>
     */
    pub const REF_MODEL_DELETION_QUAL: u8 = 30;

    /**
     * Only base calls with quality strictly greater than this constant,
     * will be considered high quality if they are part of a soft-clip.
     */
    const HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD: u8 = 28;

    // const NO_CALLS: Vec<Allele> = Vec::new();

    pub fn new(
        args: &clap::ArgMatches,
        ref_idx: usize,
        samples: Vec<String>,
        do_allele_specific_calcs: bool,
        sample_ploidy: usize,
    ) -> HaplotypeCallerEngine {
        let mut kmer_sizes = match args.get_many::<usize>("kmer-sizes") {
            Some(vals) => vals
                .map(|k_size| *k_size)
                .collect::<Vec<usize>>(),
            _ => vec![21, 33],
        };

        let mut prune_factor = if args.get_flag("use-adaptive-pruning") {
            0
        } else {
            *args.get_one::<usize>("min-prune-factor").unwrap()
        };

        let mut recover_all_dangling_branches = args.get_flag("recover-all-dangling-branches");
        let mut allow_non_unique_kmers_in_ref = args.get_flag("allow-non-unique-kmers-in-ref");

        Self::set_assembly_profile(
            args, 
            &mut kmer_sizes, 
            &mut prune_factor, 
            &mut allow_non_unique_kmers_in_ref, 
            &mut recover_all_dangling_branches
        );

        let mut assembly_engine = ReadThreadingAssembler::new(
            *args.get_one::<i32>("max-allowed-path-for-read-threading-assembler")
                .unwrap(),
            kmer_sizes,
            args.get_flag("dont-increase-kmer-sizes-for-cycles"),
            args.get_flag("allow-non-unique-kmers-in-ref"),
            *args.get_one::<i32>("num-pruning-samples")
                .unwrap(),
            prune_factor,
            args.get_flag("use-adaptive-pruning"),
            *args.get_one::<f64>("initial-error-rate-for-pruning")
                .unwrap(),
            MathUtils::log10_to_log(
                *args.get_one::<f64>("pruning-log-odds-threshold")
                    .unwrap(),
            ),
            MathUtils::log10_to_log(
                *args.get_one::<f64>("pruning-seeding-log-odds-threshold")
                    .unwrap(),
            ),
            *args.get_one::<usize>("max-unpruned-variants")
                .unwrap(),
            args.get_flag("use-linked-debruijn-graph"),
            args.get_flag("enable-legacy-graph-cycle-detection"),
            *args.get_one::<i32>("min-matching-bases-to-dangling-end-recovery")
                .unwrap(),
            args.get_flag("disable-prune-factor-correction")
        );

        assembly_engine.debug_graph_transformations =
            args.get_flag("debug-graph-transformations");
        assembly_engine.recover_dangling_branches =
            !args.get_flag("do-not-recover-dangling-branches");
        assembly_engine.recover_all_dangling_branches =
            recover_all_dangling_branches;
        assembly_engine.min_dangling_branch_length = *args
            .get_one::<i32>("min-dangling-branch-length")
            .unwrap();
        assembly_engine.graph_output_path = match args.get_one::<String>("graph-output") {
            Some(path) => Some(path.to_string()),
            None => None,
        };
        assembly_engine.debug_graph_output_path = match args.get_one::<String>("debug-graph-output") {
            Some(path) => Some(path.to_string()),
            None => None,
        };
        assembly_engine.min_base_quality_to_use_in_assembly =
            *args.get_one::<u8>("min-base-quality").unwrap();

        HaplotypeCallerEngine {
            active_region_evaluation_genotyper_engine: GenotypingEngine::make(
                args,
                samples.clone(),
                do_allele_specific_calcs,
                max(
                    sample_ploidy,
                    Self::MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY,
                ),
            ),
            genotyping_engine: HaplotypeCallerGenotypingEngine::new(
                args,
                samples,
                !args.get_flag("do-not-run-physical-phasing"),
                sample_ploidy,
            ),
            genotype_prior_calculator: GenotypePriorCalculator::make(args),
            stand_min_conf: *args
                .get_one::<f64>("standard-min-confidence-threshold-for-calling")
                .unwrap(),
            ref_idx,
            assembly_engine,
            assembly_region_trimmer: AssemblyRegionTrimmer::new(
                // *args.get_one::<usize>("assembly-region-padding")
                //     .unwrap(),
                *args.get_one::<usize>("indel-padding-for-genotyping")
                    .unwrap(),
                *args.get_one::<usize>("snp-padding-for-genotyping")
                    .unwrap(),
                *args.get_one::<usize>("str-padding-for-genotyping")
                    .unwrap(),
                *args.get_one::<usize>("max-extension-into-region-padding")
                    .unwrap(),
            ),
            likelihood_calculation_engine:
                AssemblyBasedCallerUtils::create_likelihood_calculation_engine(
                    args,
                    !args.get_flag("soft-clip-low-quality-ends"),
                ),
            mapping_quality_threshold: *args
                .get_one::<u8>("mapping-quality-threshold-for-genotyping")
                .unwrap(),
        }
    }

    fn set_assembly_profile(
        args: &clap::ArgMatches, 
        kmer_sizes: &mut Vec<usize>, 
        prune_factor: &mut usize,
        allow_non_unique_kmers_in_ref: &mut bool,
        recover_all_dangling_branches: &mut bool
    ) {
        if !args.contains_id("profile") {
            return
        }

        let profile = args.get_one::<String>("profile").unwrap().to_lowercase();
        match profile.as_str() {
            "very-fast" => {
                *prune_factor = 2; // prunes low weight chains
                if kmer_sizes.len() > 1 {
                    *kmer_sizes = vec![33]; // only a single k-mer assembly
                };
                *allow_non_unique_kmers_in_ref = false; // ignores repeat regions
                *recover_all_dangling_branches = false; // does not spend too much time recovering branches
            }
            "fast" => {
                *prune_factor = 2; // prunes low weight chains
                if kmer_sizes.len() > 1 {
                    *kmer_sizes = vec![21, 33]; // multiple k-mer assemblies
                };
                *allow_non_unique_kmers_in_ref = false; // ignores repeat regions
                *recover_all_dangling_branches = false; // does not spend too much time recovering branches
            },
            "precise" => {
                *prune_factor = 2; // prunes low weight chains
                *kmer_sizes = vec![21, 33, 45]; // multiple k-mer assemblies
                *allow_non_unique_kmers_in_ref = false; // ignores repeat regions
                *recover_all_dangling_branches = false; // does not spend too much time recovering branches
            },
            "sensitive" => {
                *prune_factor = 0; // does not prune low weight chains
                *kmer_sizes = vec![21, 33, 45]; // multiple k-mer assemblies
                // *allow_non_unique_kmers_in_ref = true; // does not ignore repeat regions
                *recover_all_dangling_branches = false; // does not spend too much time recovering branches
            },
            "super-sensitive" => {
                *prune_factor = 0; // does not prune low weight chains
                *kmer_sizes = vec![21, 33, 45, 57]; // multiple k-mer assemblies
                // *allow_non_unique_kmers_in_ref = true; // does not ignore repeat regions
                *recover_all_dangling_branches = false; // spends a lot of time recovering branches
            }
            _ => {
                // just use defaults or what is already set
                warn!("Unknown assembly profile: {}", profile);
            }
        }
    }

    pub fn stand_min_conf(&self) -> f64 {
        self.stand_min_conf
    }

    pub fn collect_activity_profile(
        &mut self,
        indexed_bam_readers: &[String],
        short_sample_count: usize,
        long_sample_count: usize,
        n_threads: usize,
        ref_idx: usize,
        m: &clap::ArgMatches,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        reference_reader: &mut ReferenceReader,
        assembly_region_padding: usize,
        min_assembly_region_size: usize,
        max_assembly_region_size: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        max_input_depth: usize,
        output_prefix: &str,
        pb_index: usize,
        pb_tree: &Arc<Mutex<Vec<&Elem>>>
    ) -> (Vec<VariantContext>, Array2<f32>) {
        // minimum PHRED base quality
        let bq = *m
            .get_one::<u8>("min-base-quality")
            .unwrap();

        let max_prob_prop = *m
            .get_one::<usize>("max-prob-propagation-distance")
            .unwrap();

        let active_prob_thresh = *m
            .get_one::<f32>("active-probability-threshold")
            .unwrap();

        let min_contig_length = *m
            .get_one::<u64>("min-contig-size")
            .unwrap();

        let min_mapq = *m.get_one::<u8>("min-mapq").unwrap();
        let min_long_read_size = *m
            .get_one::<usize>("min-long-read-size")
            .unwrap();
        let min_long_read_average_base_qual = *m
            .get_one::<usize>("min-long-read-average-base-qual")
            .unwrap();

        let limiting_interval = IntervalUtils::parse_limiting_interval(m);
        // debug!("Limiting {:?}", &limiting_interval);

        let ploidy: usize = max(
            *m.get_one::<usize>("ploidy").unwrap(),
            Self::MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY,
        );

        // let mut tids: HashMap<Arc<[u8]>, HashMap<usize, usize>> = HashMap::new();
        let mut tids: HashSet<usize> = HashSet::new();
        let mut found_contigs = HashMap::new();
        let reference = reference_reader.retrieve_reference_stem(ref_idx);


        indexed_bam_readers
            .iter()
            .enumerate()
            .for_each(|(_sample_idx, bam_generator)| {
                // Get the appropriate sample index based on how many references we are using
                let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                    vec![&bam_generator],
                    n_threads as u32,
                )
                .into_iter()
                .next()
                .unwrap();
                // get reference stats
                // let bam_generator = bam_generator.start();
                let header = bam_generator.header(); // bam header
                let header_names = header.target_names();
                header_names
                    .into_iter()
                    .enumerate()
                    .for_each(|(tid, contig_name)| {
                        let target_name = std::str::from_utf8(contig_name).unwrap();
                        let target_match = if target_name.contains("~") {
                            target_name.split_once("~").unwrap().0 == reference.as_str()
                        } else {
                            target_name.contains(&reference)
                        };
                        if target_match
                        {
                            debug!("Found reference: {} matching reference {}", target_name, &reference);
                            debug!("Ref idx {} tid {}", ref_idx, tid);
                            // Get contig stats
                            let _n_tids = reference_reader.update_ref_index_tids(ref_idx, tid);

                            reference_reader.add_target(contig_name, tid);
                            let target_len = header.target_len(tid as u32).unwrap();
                            reference_reader.add_length(tid, target_len);
                            tids.insert(tid);

                            let previous_tid = found_contigs.entry(contig_name.to_vec()).or_insert(tid);
                            if *previous_tid != tid {
                                warn!("Contig {} found more than once with different BAM header index", std::str::from_utf8(contig_name).unwrap());
                                warn!("Ensure Contigs occur in same order in all BAM files");
                                panic!("Contigs out of order in BAM files.");
                            }

                            // let tids_for_contig = tids.entry(Arc::from(contig_name)).or_insert_with(HashMap::new);
                            // tids_for_contig.entry(sample_idx).or_insert(tid);
                        }
                    })
            });
        
        let total_sample_count = short_sample_count + long_sample_count;
        let chunk_size = max(250000 / total_sample_count, max_assembly_region_size * 5);
        let genome_size = reference_reader.target_lens.values().sum::<u64>();

        // how many chunks are there going to be? for each contig, divide it's length by chunk size and count the ceiling
        let mut n_chunks = 0;
        for t_length in reference_reader.target_lens.values() {
            if *t_length < min_contig_length {
                continue;
            }

            n_chunks += (t_length / chunk_size as u64) as usize;
            if t_length % chunk_size as u64 != 0 {
                n_chunks += 1;
            }
        }

        {
            let pb = pb_tree.lock().unwrap();
            // simple spinner with a counter
            let new_style = ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>7}/{len:7}").unwrap();
            pb[pb_index].progress_bar.set_style(new_style);
            pb[pb_index].progress_bar.set_length(n_chunks as u64);
            pb[pb_index].progress_bar.set_message(format!("{}: Generating activity profile...", &pb[pb_index].key));
        }

        let contexts = tids
            .into_par_iter()
            .fold(
                || (Vec::with_capacity(genome_size as usize), Array2::default((total_sample_count, total_sample_count))), 
                |mut consolidator: (Vec<VariantContext>, Array2<f32>), tid: usize| {
                let target_length = reference_reader.target_lens[&tid];
                let mut reference_reader = reference_reader.clone();
                reference_reader.update_current_sequence_capacity(target_length as usize);
                // Update all contig information
                let retrieved =
                    match reference_reader.fetch_contig_from_reference_by_tid(tid, ref_idx) {
                        Ok(_) => true,
                        Err(_) => false,
                    };

                let context_depth_tuples = if retrieved {
                    // let mut contexts = Vec::new();
                    reference_reader.read_sequence_to_vec();
                    if target_length >= min_contig_length {

                        let result = (0..target_length as usize)
                            .into_par_iter()
                            .chunks(chunk_size)
                            .enumerate()
                            .fold(
                                || (Vec::with_capacity(chunk_size), Array2::default((total_sample_count, total_sample_count))), 
                                | mut consolidator: (Vec<VariantContext>, Array2<f32>), chunk_vals: (usize, Vec<usize>)| {
                                let (chunk_idx, positions) = chunk_vals;
                                let within_bounds = match &limiting_interval {
                                    Some(limit) => {
                                        let position_limit = SimpleInterval::new(
                                            0,
                                            positions[0],
                                            *positions.last().unwrap(),
                                        );
                                        position_limit.overlaps(limit)
                                    }
                                    None => true,
                                };

                                if within_bounds {
                                    let mut genotype_likelihoods = Vec::new();
                                    let first = positions[0];
                                    let last = *positions.last().unwrap();
                                    let length = last - first + 1;
                                    let mut per_contig_per_base_hq_soft_clips =
                                        vec![RunningAverage::new(); length];
                                    let chunk_location = SimpleInterval::new(tid, first, last);

                                    indexed_bam_readers.iter().enumerate().for_each(
                                        |(sample_idx, bam_generator)| {
                                            // Get the appropriate sample index based on how many references we are using
                                            // let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                                            //     vec![&bam_generator],
                                            //     n_threads as u32,
                                            // )
                                            // .into_iter()
                                            // .next()
                                            // .unwrap();
                                            // Get the appropriate sample index based on how many references we are using
                                            let bam_generator =
                                                generate_indexed_named_bam_readers_from_bam_files(
                                                    vec![&bam_generator],
                                                    n_threads as u32,
                                                )
                                                .into_iter()
                                                .next()
                                                .unwrap();
                                            let mut bam_generated = bam_generator.start();

                                            let mut read_type = ReadType::Short;

                                            if (m.contains_id("longreads")
                                                || m.contains_id("longread-bam-files"))
                                                && sample_idx >= short_sample_count
                                                && sample_idx
                                                    < (short_sample_count + long_sample_count)
                                            {
                                                // Get the appropriate sample index based on how many references we are using by tracking
                                                // changes in references
                                                read_type = ReadType::Long;
                                            }

                                            HaplotypeCallerEngine::update_activity_profile(
                                                &mut bam_generated,
                                                n_threads,
                                                ref_idx,
                                                read_type,
                                                ploidy,
                                                bq,
                                                genomes_and_contigs,
                                                concatenated_genomes,
                                                flag_filters,
                                                &mut per_contig_per_base_hq_soft_clips,
                                                &reference_reader,
                                                &limiting_interval,
                                                &mut genotype_likelihoods,
                                                min_contig_length,
                                                tid,
                                                &chunk_location,
                                                chunk_idx,
                                                target_length,
                                                min_mapq,
                                                min_long_read_size,
                                                min_long_read_average_base_qual,
                                            );
                                        },
                                    );
                                    // chunk_idx += 1;
                                    debug!("Beginning calling on chunk {}", chunk_idx);
                                    match self.calculate_activity_probabilities(
                                        genotype_likelihoods,
                                        per_contig_per_base_hq_soft_clips,
                                        &limiting_interval,
                                        tid,
                                        reference_reader.target_lens[&tid],
                                        ploidy,
                                        max_prob_prop,
                                        active_prob_thresh,
                                        ref_idx,
                                        indexed_bam_readers,
                                        min_contig_length,
                                        flag_filters,
                                        m,
                                        &reference_reader,
                                        n_threads as u32,
                                        assembly_region_padding,
                                        min_assembly_region_size,
                                        max_assembly_region_size,
                                        short_read_bam_count,
                                        long_read_bam_count,
                                        max_input_depth,
                                        &chunk_location,
                                        chunk_idx,
                                        output_prefix,
                                    ) {
                                        Ok(val) => {
                                            debug!("Finished calling on chunk {} of size {}", chunk_idx, positions.len());
                                            debug!("N. variant contexts {}", val.0.len());
                                            debug!("N. depth counts {:?}", val.1.shape());
                                            {
                                                let pb = pb_tree.lock().unwrap();
                                                pb[pb_index].progress_bar.inc(1);
                                            }
                                            let (vc_vec, concatenated_array) = val;
                                            consolidator.0.extend(vc_vec);
                                            (consolidator.0, consolidator.1 + &concatenated_array)
                                        },
                                        Err(e) => {
                                            panic!("Activity profile generation failed: {:?}", e)
                                        }
                                    }
                                } else {
                                    consolidator
                                }
                            })
                            .reduce(|| (Vec::with_capacity(target_length as usize), Array2::default((total_sample_count, total_sample_count))), |mut a, b| {
                                a.0.extend(b.0);
                                (a.0, a.1 + &b.1)
                            });
                        result
                    } else {
                        (Vec::new(), Array2::default((total_sample_count, total_sample_count)))
                    }
                    // contexts
                } else {
                    (Vec::new(), Array2::default((total_sample_count, total_sample_count)))
                };

                consolidator.0.extend(context_depth_tuples.0);

                (consolidator.0, consolidator.1 + &context_depth_tuples.1)
            })
            .reduce(|| (Vec::with_capacity(genome_size as usize), Array2::default((total_sample_count, total_sample_count))), |mut a, b| {
                a.0.extend(b.0);
                (a.0, a.1 + &b.1)
            });
        {
            let pb = pb_tree.lock().unwrap();
            pb[pb_index].progress_bar.finish_with_message(format!("{}: Finished generating activity profile.", &pb[pb_index].key));
        }
        (contexts.0, contexts.1)
    }

    pub fn update_activity_profile<'b>(
        bam_generated: &mut IndexedBamFileNamedReader,
        _split_threads: usize,
        _ref_idx: usize,
        readtype: ReadType,
        ploidy: usize,
        bq: u8,
        _genomes_and_contigs: &'b GenomesAndContigs,
        _concatenated_genomes: &'b Option<String>,
        flag_filters: &'b FlagFilter,
        per_contig_per_base_hq_soft_clips: &mut Vec<RunningAverage>,
        reference_reader: &ReferenceReader,
        limiting_interval: &Option<SimpleInterval>,
        current_likelihoods: &mut Vec<Vec<RefVsAnyResult>>,
        _min_contig_length: u64,
        tid: usize,
        outer_chunk_location: &SimpleInterval,
        _outer_chunk_index: usize,
        target_len: u64,
        min_mapq: u8,
        min_long_read_size: usize,
        min_long_read_average_base_qual: usize,
    ) {
        let likelihoodcount = ploidy + 1;
        let log10ploidy = (ploidy as f64).log10();

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        match readtype {
            ReadType::Short | ReadType::Long => {
                // The raw activity profile.
                // Frequency of bases not matching reference compared
                // to depth
                let mut likelihoods = (0..outer_chunk_location.size())
                    .into_iter()
                    .map(|pos| RefVsAnyResult::new(likelihoodcount, pos, tid))
                    .collect::<Vec<RefVsAnyResult>>();

                let mut positions = per_contig_per_base_hq_soft_clips
                    .iter_mut()
                    .zip(likelihoods.iter_mut())
                    .map(|(r, l)| (r, l))
                    .collect::<Vec<(&mut RunningAverage, &mut RefVsAnyResult)>>();

                // multiplier to help us map between chunk position
                // and actual reference position. This value represents the
                // starting reference base index of this chunk.
                // let chunk_multiplier = chunk_idx * chunk_size;

                bam_generated
                    .fetch((
                        tid as i32,
                        outer_chunk_location.start as i64,
                        min(outer_chunk_location.end + 1, target_len as usize - 1) as i64,
                    ))
                    .unwrap_or_else(|_| {
                        warn!("BAM index potentially outdated. Fetching of interval failed.");
                        panic!(
                            "Failed to fetch interval {}:{}-{} contig. Try regenerating BAM indices, or deleting old BAI files.",
                            tid,
                            outer_chunk_location.start,
                            min(outer_chunk_location.end + 1, target_len as usize - 1),
                            // std::str::from_utf8(bam_generated.header().tid2name(tid as u32)).unwrap()
                        )
                    });

                let mut record = Record::new();
                // pileup does not provide all of the alignments at a pos
                // Unsure why, but there are positions where entire
                // alignments are being ignored. Solution: create our
                // own version of a pileup and hopefully records won't be
                // filtered
                while bam_generated.read(&mut record) {
                    if ReadUtils::read_is_filtered(
                        &record,
                        &flag_filters,
                        min_mapq,
                        readtype,
                        limiting_interval,
                        min_long_read_size,
                        min_long_read_average_base_qual,
                    ) {
                        continue;
                    }

                    Self::parse_record(
                        &record,
                        flag_filters,
                        readtype,
                        &mut positions,
                        min(outer_chunk_location.start, target_len as usize),
                        bq,
                        likelihoodcount,
                        reference_reader,
                        log10ploidy,
                        outer_chunk_location.start,
                        min(outer_chunk_location.end + 1, target_len as usize),
                        false,
                    );
                }
                Self::update_ref_vs_any_results(
                    &mut likelihoods,
                    likelihoodcount,
                    log10ploidy,
                    outer_chunk_location.start,
                    false,
                );
                current_likelihoods.push(likelihoods);
            }
        }
    }

    fn update_ref_vs_any_results(
        likelihoods: &mut Vec<RefVsAnyResult>,
        likelihoodcount: usize,
        log10ploidy: f64,
        _outer_chunk_start: usize,
        _debug: bool,
    ) {
        for (_pos, result) in likelihoods.iter_mut().enumerate() {
            let denominator = result.read_counts as f64 * log10ploidy;

            for i in 0..likelihoodcount {
                result.genotype_likelihoods[i] -= denominator
            }
        }
    }

    fn parse_record(
        record: &Record,
        _flag_filters: &FlagFilter,
        _readtype: ReadType,
        positions: &mut Vec<(&mut RunningAverage, &mut RefVsAnyResult)>,
        subtractor: usize,
        bq: u8,
        likelihoodcount: usize,
        reference_reader: &ReferenceReader,
        log10ploidy: f64,
        bound_start: usize,
        bound_end: usize,
        _debug: bool,
    ) -> usize {
        // we need to iterate through each read pos (qpos)
        // and generate a pseudo pileup
        let _read_length = record.len();
        // let mut cigar_cursor = 0; // changes when cigar is consumed
        let mut read_cursor = 0; // changes when read bases are consumed
        let mut pos = record.pos() as usize; // read start alignment on reference.
                                             // updated as ref consumed
        let count = 0;
        let mut cig_index = 0;
        for cig in &record.cigar() {
            match cig {
                Cigar::Del(len) => {
                    // reference bases consumed

                    for _ in 0..*len as usize {
                        // if pos == 97055 {
                        //     count += 1;
                        // }

                        if pos < bound_start {
                            // read not within bounds yet
                            pos += 1;
                            continue;
                        } else if pos >= bound_end {
                            // read pass the bounds now so break
                            break;
                        };
                        let alignment = PosAlignment::new(None, true, false);
                        let (hq_soft_clips, result) =
                            &mut positions[pos.saturating_sub(subtractor)];
                        let refr_base = reference_reader.current_sequence[pos];
                        Self::alignment_context_creation(
                            &alignment,
                            &record,
                            result,
                            hq_soft_clips,
                            log10ploidy,
                            likelihoodcount,
                            refr_base,
                            bq,
                            cig_index,
                        );
                        pos += 1;
                    }
                    // cigar_cursor += *len as usize;
                }
                Cigar::RefSkip(_len) => {
                    panic!(
                        "Read contains N operator, should have been filtered prior to this point."
                    );
                }
                Cigar::Ins(len) => {
                    // read bases consumed
                    // if pos == 97055 {
                    //     count += 1;
                    // }
                    if pos < bound_start {
                        // read not within bounds yet
                        // cigar_cursor += *len as usize;
                        read_cursor += *len as usize;
                        continue;
                    } else if pos >= bound_end {
                        // read pass the bounds now so break
                        break;
                    };
                    // insertion like event
                    // read is consumed
                    let alignment = PosAlignment::new(Some(read_cursor), false, false);
                    let (hq_soft_clips, result) = &mut positions[pos.saturating_sub(subtractor)];
                    let refr_base = reference_reader.current_sequence[pos];
                    Self::alignment_context_creation(
                        &alignment,
                        &record,
                        result,
                        hq_soft_clips,
                        log10ploidy,
                        likelihoodcount,
                        refr_base,
                        bq,
                        cig_index,
                    );
                    // cigar_cursor += *len as usize;
                    read_cursor += *len as usize;
                }
                Cigar::Diff(len) | Cigar::Match(len) | Cigar::Equal(len) => {
                    // we need to check each position in
                    // these cigars
                    for _ in 0..*len as usize {
                        
                        if pos < bound_start {
                            // read not within bounds yet
                            read_cursor += 1;
                            pos += 1;
                            continue;
                        } else if pos >= bound_end {
                            // read pass the bounds now so break
                            break;
                        }

                        let alignment = PosAlignment::new(Some(read_cursor), false, false);
                        let (hq_soft_clips, result) =
                            &mut positions[pos.saturating_sub(subtractor)];
                        let refr_base = reference_reader.current_sequence[pos];
                        Self::alignment_context_creation(
                            &alignment,
                            &record,
                            result,
                            hq_soft_clips,
                            log10ploidy,
                            likelihoodcount,
                            refr_base,
                            bq,
                            cig_index,
                        );
                        read_cursor += 1;
                        pos += 1;
                    }
                    // cigar_cursor += *len as usize;
                }
                Cigar::SoftClip(len) => {
                    // cigar_cursor += *len as usize;
                    read_cursor += *len as usize;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {
                    // ignore these
                }
            }
            cig_index += 1;
        }

        count
    }

    /**
     * Create a collection of ActivityProfileState for each position on each contig
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     *
     * @return HashMap<usize, BandPassActivityProfile> - contig id and per base activity
     */
    pub fn calculate_activity_probabilities(
        &self,
        genotype_likelihoods: Vec<Vec<RefVsAnyResult>>,
        per_contig_per_base_hq_soft_clips: Vec<RunningAverage>,
        limiting_interval: &Option<SimpleInterval>,
        tid: usize,
        length: u64,
        ploidy: usize,
        max_prob_propagation: usize,
        active_prob_threshold: f32,
        ref_idx: usize,
        sample_names: &[String],
        _min_contig_length: u64,
        flag_filters: &FlagFilter,
        args: &clap::ArgMatches,
        reference_reader: &ReferenceReader,
        n_threads: u32,
        assembly_region_padding: usize,
        min_assembly_region_size: usize,
        max_assembly_region_size: usize,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
        max_input_depth: usize,
        chunk_location: &SimpleInterval,
        _chunk_index: usize,
        output_prefix: &str,
    ) -> Result<(Vec<VariantContext>, Array2<f32>), BirdToolError> {
        // let mut per_contig_activity_profiles = HashMap::new();
        let placeholder_vec = Vec::new();
        let depth_per_sample_filter = *args
            .get_one::<i64>("depth-per-sample-filter")
            .unwrap() as i32;
        
        // the total sample count will increase the number of RAM we will be using
        // each sample adds a "Genotype" struct which is a large struct with many fields
        let total_sample_count = short_read_bam_count + long_read_bam_count;

        // so we need to limit the inner_chunk_size to a reasonable size
        // to avoid OOM errors whilst also keeping it above max_assembly_region_size
        let inner_chunk_size = max(50000 / max(total_sample_count / 2, 1), max_assembly_region_size * 2);

        let variant_contexts = (0..chunk_location.size())
            .into_par_iter()
            .chunks(inner_chunk_size)
            .enumerate()
            .fold(|| (Vec::with_capacity(inner_chunk_size), Array2::default((total_sample_count, total_sample_count))), 
            |mut consolidator, (i, positions)| {
                let within_bounds = match &limiting_interval {
                    Some(limit) => {
                        let position_limit = SimpleInterval::new(
                            0,
                            chunk_location.start + i * inner_chunk_size,
                            chunk_location.end + i * inner_chunk_size,
                        );
                        // debug!("Position limit {:?}", &position_limit);
                        position_limit.overlaps(limit)
                    }
                    None => true,
                };
                let n_positions = positions.len();

                if within_bounds {
                    let mut active_region_evaluation_genotyper_engine =
                        self.active_region_evaluation_genotyper_engine.clone();

                    // Create bandpass
                    // debug!("Created bandpass profile");
                    // debug!(
                    //     "Calculating activity on {} of {:?} {} {}",
                    //     tid,
                    //     chunk_location,
                    //     positions[0],
                    //     positions.last().unwrap()
                    // );
                    let mut activity_profile = BandPassActivityProfile::new(
                        max_prob_propagation,
                        active_prob_threshold,
                        BandPassActivityProfile::MAX_FILTER_SIZE,
                        BandPassActivityProfile::DEFAULT_SIGMA,
                        true,
                        ref_idx,
                        tid,
                        length as usize,
                    );

                    let mut depth_of_position =
                        vec![Vec::with_capacity(n_positions); sample_names.len()];
                    let mut depths_counters = vec![0; sample_names.len()];
                    for pos in positions {
                        match &limiting_interval {
                            Some(limit) => {
                                let position_limit = SimpleInterval::new(
                                    0,
                                    chunk_location.start + pos,
                                    chunk_location.start + pos,
                                );
                                if !position_limit.overlaps(limit) {
                                    continue;
                                }
                            }
                            None => {
                                //
                            }
                        }
                        let mut genotypes = Vec::new();
                        let hq_soft_clips = per_contig_per_base_hq_soft_clips[pos];

                        // ANI should only be performed on "compared bases", that is bases that were >= depth per sample filter in both sample.
                        // First we need to store the number of bases in a sample above the require depth, for a contig of length 10:
                        //
                        // ```
                        // [5, 5, 5, 5, 0, 0, 0, 5, 5, 5] # sample 1
                        //
                        // [5, 5, 5, 0, 0, 5, 5, 5, 5, 0] # sample 2
                        // ```
                        //
                        // we could store the information compressed as such if the depth per sample filter = 5:
                        // - `[4, -3, 3]`
                        // - `[3, -2, 4, -1]`
                        for (idx, sample_likelihoods) in genotype_likelihoods.iter().enumerate() {
                            let ref_v_any = &sample_likelihoods[pos];
                            
                            // create compressed array of bases passing the depth filter.
                            // Used during ANI calculations to determine number of comparable
                            // bases between two samples.
                            if ref_v_any.get_dp() >= depth_per_sample_filter {
                                if depths_counters[idx] >= 0 {
                                    // positively increment current sample counter
                                    depths_counters[idx] += 1;
                                } else {
                                    depth_of_position[idx].push(depths_counters[idx]); // push previous stretch of negative bases
                                    depths_counters[idx] = 0;
                                    depths_counters[idx] += 1;
                                }
                            } else {
                                if depths_counters[idx] <= 0 {
                                    // negatively increment when a base does not have suffcient coverage
                                    depths_counters[idx] -= 1;
                                } else {
                                    depth_of_position[idx].push(depths_counters[idx]); // push previosu stretch of positive bases
                                    depths_counters[idx] = 0;
                                    depths_counters[idx] -= 1;
                                }
                            };
                            let result = ref_v_any.genotype_likelihoods.clone();
                            genotypes.push(Genotype::build(
                                ploidy,
                                result,
                                idx,
                            ))
                        }

                        let fake_alleles = ByteArrayAllele::create_fake_alleles();

                        let contig_position = chunk_location.start + pos;
                        let mut variant_context = VariantContext::build(
                            tid,
                            contig_position,
                            contig_position,
                            fake_alleles,
                        );

                        variant_context.add_genotypes(genotypes);

                        let vc_out = active_region_evaluation_genotyper_engine.calculate_genotypes(
                            variant_context,
                            ploidy,
                            &self.genotype_prior_calculator,
                            &placeholder_vec,
                            self.stand_min_conf,
                        );

                        let is_active_prob = match vc_out {
                            Some(vc) => {
                                QualityUtils::qual_to_prob(vc.get_phred_scaled_qual() as u8)
                            }
                            None => 0.0,
                        };

                        // debug!(
                        //     "{}-{} Active Prob {}",
                        //     chunk_location.start + pos,
                        //     chunk_location.start + pos,
                        //     is_active_prob
                        // );

                        let activity_profile_state = ActivityProfileState::new(
                            SimpleInterval::new(
                                tid,
                                chunk_location.start + pos,
                                chunk_location.start + pos,
                            ),
                            is_active_prob as f32,
                            ActivityProfileDataType::new(
                                hq_soft_clips.mean() as f32,
                                HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                            ),
                        );
                        activity_profile.add(activity_profile_state);
                    }

                    for (idx, sample_depth) in depths_counters.into_iter().enumerate() {
                        depth_of_position[idx].push(sample_depth);
                    }
                    
                    let comparable_bases = ANICalculator::calculate_compared_bases(Some(depth_of_position), n_positions as u64, total_sample_count);

                    let inner_reader = ReferenceReader::new_from_reader_with_tid_and_rid(
                        reference_reader,
                        ref_idx,
                        tid,
                    );

                    let processed = AssemblyRegionWalker::process_shard(
                        activity_profile,
                        flag_filters,
                        args,
                        sample_names,
                        &inner_reader,
                        n_threads,
                        assembly_region_padding,
                        min_assembly_region_size,
                        max_assembly_region_size,
                        short_read_bam_count,
                        long_read_bam_count,
                        &self,
                        max_input_depth,
                        output_prefix,
                    );

                    consolidator.0.extend(processed);
                    (consolidator.0, consolidator.1 + &comparable_bases)
                } else {
                    (Vec::new(), Array2::default((total_sample_count, total_sample_count)))
                }
                // activity_profile
            })
            .reduce(|| (Vec::with_capacity(inner_chunk_size), Array2::default((total_sample_count, total_sample_count))),
                |mut a, b| {
                    a.0.extend(b.0);
                    (a.0, a.1 + &b.1)
                });


        Ok((variant_contexts.0, variant_contexts.1))
    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    pub fn call_region<'b>(
        &mut self,
        mut region: AssemblyRegion,
        reference_reader: &'b mut ReferenceReader,
        given_alleles: Vec<VariantContext>,
        args: &'b clap::ArgMatches,
        sample_names: &'b [String],
        flag_filters: &'b FlagFilter,
    ) -> Vec<VariantContext> {
        let vc_priors = Vec::new();

        if !region.is_active() && given_alleles.is_empty() {
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }

        if given_alleles.is_empty() && region.len() == 0 {
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }

        let region_without_reads = region.clone_without_reads();

        // run the local assembler, getting back a collection of information on how we should proceed
        let mut untrimmed_assembly_result = AssemblyBasedCallerUtils::assemble_reads(
            region,
            &given_alleles,
            args,
            reference_reader,
            &mut self.assembly_engine,
            true,
            sample_names,
        );

        let all_variation_events = match untrimmed_assembly_result
            .get_variation_events(*args.get_one::<usize>("max-mnp-distance").unwrap())
        {
            Ok(result) => result,
            Err(_) => {
                return self.reference_model_for_no_variation(
                    &mut untrimmed_assembly_result.region_for_genotyping,
                    true,
                    &vc_priors,
                )
            }
        };

        // debug!(
        //     "Region {:?} All variation events  {:?}",
        //     &untrimmed_assembly_result.padded_reference_loc,
        //     &all_variation_events.len()
        // );

        let mut trimming_result = self.assembly_region_trimmer.trim(
            self.ref_idx,
            region_without_reads,
            all_variation_events,
            // args.get_flag("enable-legacy-assembly-region-trimming"),
            false,
            &reference_reader,
            untrimmed_assembly_result
                .full_reference_with_padding
                .as_slice(),
        );

        // debug!("Trim complete!");

        if !trimming_result.is_variation_present() && !args.get_flag("disable-optimizations") {
            return self.reference_model_for_no_variation(
                &mut trimming_result.original_region,
                false,
                &vc_priors,
            );
        }

        // debug!("Moving reads....");
        trimming_result.original_region.reads =
            untrimmed_assembly_result.region_for_genotyping.move_reads();
        // debug!(
        //     "Move complete! {}",
        //     trimming_result.original_region.reads.len()
        // );
        let assembly_result =
            untrimmed_assembly_result.trim_to(trimming_result.get_variant_region());

        // debug!(
        //     "Assembly result allele order after region trimming {:?}",
        //     &assembly_result.haplotypes
        // );

        let read_stubs = assembly_result
            .region_for_genotyping
            .reads
            .iter()
            .filter(|r| {
                AlignmentUtils::unclipped_read_length(r)
                    < AssemblyBasedCallerUtils::MINIMUM_READ_LENGTH_AFTER_TRIMMING
            })
            .cloned()
            .collect::<Vec<BirdToolRead>>();
        let assembly_result = assembly_result.remove_all(&read_stubs);
        // debug!(
        //     "Assembly result allele order after stub filter {:?} -> read stubs {}",
        //     &assembly_result.haplotypes.len(),
        //     read_stubs.len()
        // );

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        let (mut assembly_result, filtered_reads) =
            self.filter_non_passing_reads(assembly_result, flag_filters);
        // let filtered_reads = Vec::new();
        // debug!("Filtered reads {}", filtered_reads.len());
        // debug!(
        //     "Assembly result allele order after read filter {:?}",
        //     &assembly_result.haplotypes.len()
        // );
        let per_sample_filtered_read_list =
            AssemblyBasedCallerUtils::split_reads_by_sample(filtered_reads);

        if !assembly_result.variation_present || assembly_result.region_for_genotyping.len() == 0 {
            return self.reference_model_for_no_variation(
                &mut assembly_result.region_for_genotyping,
                false,
                &vc_priors,
            );
        };

        // debug!("==========================================================================");
        // debug!(
        //     "              Assembly region {:?} \n\
        //                    Alleles {:?}",
        //     &assembly_result.region_for_genotyping, &assembly_result.haplotypes
        // );
        // debug!("==========================================================================");
        let reads = AssemblyBasedCallerUtils::split_reads_by_sample(
            assembly_result.region_for_genotyping.move_reads(),
        );

        // debug!(
        //     "Before change: {}",
        //     reads.values().map(|r| r.len()).sum::<usize>()
        // );
        let sample_indices = sample_names
            .iter()
            .enumerate()
            .map(|(idx, _)| idx)
            .collect::<Vec<usize>>();
        let mut read_likelihoods: AlleleLikelihoods<Haplotype<SimpleInterval>> = self
            .likelihood_calculation_engine
            .compute_read_likelihoods(&mut assembly_result, sample_indices, reads);

        // if debug {
        // debug!(
        //     "Read by sample after compute: {:?}",
        //     (0..sample_names.len())
        //         .map(|s| read_likelihoods.sample_evidence_count(s))
        //         .collect::<Vec<usize>>()
        // );
        // }

        // debug!(
        //     "Read likelihoods first {:?} {:?} {:?}",
        //     read_likelihoods.alleles.len(),
        //     read_likelihoods
        //         .alleles
        //         .list
        //         .iter()
        //         .map(|a| a.is_ref())
        //         .collect::<Vec<bool>>(),
        //     read_likelihoods
        //         .alleles
        //         .list
        //         .iter()
        //         .map(|a| std::str::from_utf8(a.get_bases()).unwrap())
        //         .collect::<Vec<&str>>()
        // );
        if read_likelihoods.alleles.len() == 1 {
            return self.reference_model_for_no_variation(
                &mut assembly_result.region_for_genotyping,
                false,
                &vc_priors,
            );
        };
        // Realign reads to their best haplotype.
        let read_alignments = AssemblyBasedCallerUtils::realign_reads_to_their_best_haplotype(
            &mut read_likelihoods,
            &assembly_result.ref_haplotype,
            &assembly_result.padded_reference_loc,
            if args.get_flag("disable-avx") {
                AVXMode::None
            } else {
                AVXMode::detect_mode()
            },
        );
        read_likelihoods.change_evidence(read_alignments);

        // if debug {
        // debug!(
        //     "After change {:?}",
        //     (0..sample_names.len())
        //         .map(|s| read_likelihoods.sample_evidence_count(s))
        //         .collect::<Vec<usize>>()
        // );
        // // }

        // debug!(
        //     "Read likelihoods after change evidence {:?}",
        //     read_likelihoods.alleles.len()
        // );

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]
        let called_haplotypes = match self.genotyping_engine.assign_genotype_likelihoods(
            assembly_result.haplotypes.clone(),
            read_likelihoods,
            per_sample_filtered_read_list,
            assembly_result.full_reference_with_padding.as_slice(),
            &assembly_result.padded_reference_loc,
            &assembly_result.region_for_genotyping.active_span,
            given_alleles,
            false,
            *args.get_one::<usize>("max-mnp-distance").unwrap(),
            sample_names,
            *args.get_one::<usize>("ploidy").unwrap(),
            args,
            &reference_reader,
            self.stand_min_conf,
        ) {
            Ok(result) => result,
            Err(_) => {
                return self.reference_model_for_no_variation(
                    &mut assembly_result.region_for_genotyping,
                    false,
                    &vc_priors,
                );
            }
        };

        // TODO: Bam writing? Don't think
        //       Emit reference confidence? Maybe
        //

        return called_haplotypes.calls;
    }

    fn filter_non_passing_reads(
        &self,
        assembly_result: AssemblyResultSet<ReadThreadingGraph>,
        _flag_filters: &FlagFilter,
    ) -> (AssemblyResultSet<ReadThreadingGraph>, Vec<BirdToolRead>) {
        let reads_to_remove = assembly_result
            .region_for_genotyping
            .reads
            .par_iter()
            .filter(|r| {
                if AlignmentUtils::unclipped_read_length(r) < Self::READ_LENGTH_FILTER_THRESHOLD
                    || r.read.mapq() < self.mapping_quality_threshold
                    || (r.read.is_paired()
                        && (!r.read.is_mate_unmapped()
                            && (!r.read.is_unmapped() && r.read.tid() != r.read.mtid())))
                // || (!flag_filters.include_secondary && r.read.is_secondary())
                {
                    // debug!("Removing read {:?}", r);
                    true
                } else {
                    false
                }
            })
            .cloned()
            .collect::<Vec<BirdToolRead>>();

        return (
            assembly_result.remove_all(&reads_to_remove),
            reads_to_remove,
        );
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    fn reference_model_for_no_variation<'a>(
        &'a self,
        _region: &'a mut AssemblyRegion,
        _needs_to_be_finalized: bool,
        _vc_priors: &Vec<VariantContext>,
    ) -> Vec<VariantContext> {
        Vec::new() // TODO: Implement this potentially?
    }

    /**
     * Populates reference and non reference depth vectors
     */
    fn alignment_context_creation(
        alignment: &PosAlignment,
        record: &Record,
        result: &mut RefVsAnyResult,
        hq_soft_clips: &mut RunningAverage,
        log10ploidy: f64,
        likelihoodcount: usize,
        refr_base: u8,
        bq: u8,
        cig_index: usize,
    ) {
        let ref_likelihood;
        let non_ref_likelihood;

        // query position in read
        let qpos = &alignment.qpos;
        let record_qual = if alignment.is_del || alignment.is_refskip {
            Self::REF_MODEL_DELETION_QUAL
        } else {
            record.qual()[qpos.unwrap()]
        };
        let mut is_alt = false;

        if record_qual >= bq || alignment.is_del {
            result.read_counts += 1;

            is_alt = Self::is_alt(
                &record, qpos,
                refr_base,
                // Self::HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD,
                // debug,
            );

            if is_alt {
                result.non_ref_depth += 1;
                non_ref_likelihood = QualityUtils::qual_to_prob_log10(record_qual);
                ref_likelihood =
                    QualityUtils::qual_to_error_prob_log10(record_qual) + (-(3.0_f64.log10()))
            } else {
                result.ref_depth += 1;
                ref_likelihood = QualityUtils::qual_to_prob_log10(record_qual);
                non_ref_likelihood =
                    QualityUtils::qual_to_error_prob_log10(record_qual) + (-(3.0_f64.log10()));
            }

            HaplotypeCallerEngine::update_heterozygous_likelihood(
                result,
                likelihoodcount,
                log10ploidy,
                ref_likelihood,
                non_ref_likelihood,
            );
        }

        
        // add hq soft clips if possible
        if is_alt && Self::next_to_soft_clip(record, cig_index, alignment.qpos) {
            // TODO: Counts are different to GATK. hq soft clips seems to cap out
            //       at 101, not sure if that is an issue with rolling average or
            //       this section of code. GATK also has average soft clips values greater
            //       than length of read, so we might be correct here.
            let soft_clips = HaplotypeCallerEngine::count_high_quality_soft_clips(
                // &cig,
                &record,
                // read_cursor as usize,
                Self::HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD,
            );
            // debug!("Adding {}", soft_clips);
            hq_soft_clips.add(soft_clips);
        }
    }

    fn next_to_soft_clip(record: &Record, cig_index: usize, qpos: Option<usize>) -> bool {
        match qpos {
            Some(qpos) => Self::next_to_soft_clip_or_indel(record, qpos, false),
            None => {
                CigarUtils::cigar_is_soft_clip(&record.cigar().0[cig_index.saturating_sub(1)])
                    || CigarUtils::cigar_is_soft_clip(
                        &record.cigar().0[min(cig_index + 1, record.cigar_len() - 1)],
                    )
                    || CigarUtils::cigar_is_soft_clip(&record.cigar().0[cig_index])
            }
        }
    }

    /// Determine whether a pileup alignment is considered alternative
    /// conditions for alt are one of:
    ///     - base does not match ref
    ///     - is deletion
    ///     - is before a deletion
    ///     - is after a deletion
    ///     - is before insertion
    ///     - is after insertion
    ///     - is next to softclip
    fn is_alt(
        record: &Record,
        qpos: &Option<usize>,
        refr_base: u8,
        // min_soft_clip_qual: u8,
        // debug: bool
    ) -> bool {
        match qpos {
            Some(qpos) => {
                // there is a matching position in the read
                let read_char = record.seq()[*qpos];
                // is the alignment next to indel or softclip?
                let next_to_sc_indel = Self::next_to_soft_clip_or_indel(
                    record, *qpos, // None,
                    true,
                    // min_soft_clip_qual,
                );

                if read_char.to_ascii_uppercase() != refr_base.to_ascii_uppercase()
                    || next_to_sc_indel
                {
                    return true;
                }

                return false;
            }
            None => {
                // is deletion or ref_skip

                return true;
            }
        }
    }

    /// Checks if an alignment position is next to (immediately before or after) a softclip
    /// If so, return true otherwise false
    /// If hq_soft_clips is not None, add high quality soft clips to mutable running average
    /// If check indels is true, also return true when next INDEL
    fn next_to_soft_clip_or_indel(
        record: &Record,
        qpos: usize,
        // mut hq_soft_clips: Option<&mut RunningAverage>,
        check_indels: bool,
        // min_soft_clip_qual: u8,
    ) -> bool {
        // we want to check both the front and end of a cigar so two variable cursors
        let mut read_cursor = 0; // position at start of cigar
        let mut end_of_cigar_read_cursor = 0; // position at end of cigar
        let mut next_to_soft_clip = false;
        let mut past_query_pos;
        let qpos_to_cigar_cursor = qpos as i32 + 1; // convert to one-based for better
                                                    // read cursor comparison
                                                    // let cigar = &record.cigar().0;

        for cig in record.cigar().iter() {
            if CigarUtils::cigar_consumes_read_bases(cig) {
                end_of_cigar_read_cursor = read_cursor + cig.len() as i32;
            }

            if qpos_to_cigar_cursor == read_cursor {
                // debug!("Adjacent {} {} {:?} {:?}", qpos, read_cursor, cig, record.cigar());
                next_to_soft_clip = Self::check_position_against_cigar(
                    // &mut hq_soft_clips,
                    record,
                    cig,
                    read_cursor - 1,
                    // min_soft_clip_qual,
                    check_indels,
                )
            } else if qpos_to_cigar_cursor - 1 == end_of_cigar_read_cursor {
                next_to_soft_clip = Self::check_position_against_cigar(
                    // &mut hq_soft_clips,
                    record,
                    cig,
                    read_cursor - 1,
                    // min_soft_clip_qual,
                    check_indels,
                )
            }

            past_query_pos = read_cursor >= qpos as i32;
            // Cigar immediately before current position in read
            if past_query_pos || next_to_soft_clip {
                // should catch events after current position
                break;
            }

            // Progress the cigar cursor
            if CigarUtils::cigar_consumes_read_bases(cig) {
                read_cursor += cig.len() as i32;
            }
        }

        return next_to_soft_clip;
    }

    fn check_position_against_cigar(
        // hq_soft_clips: &mut Option<&mut RunningAverage>,
        _record: &Record,
        cig: &Cigar,
        _read_cursor: i32,
        // min_soft_clip_qual: u8,
        check_indels: bool,
    ) -> bool {
        let mut next_to_soft_clip = false;
        match cig {
            Cigar::SoftClip(_) => {
                next_to_soft_clip = true;
                // if let Some(ref mut hq_soft_clips) = hq_soft_clips {
                //     hq_soft_clips.add(
                //         HaplotypeCallerEngine::count_high_quality_soft_clips(
                //             // &cig,
                //             &record,
                //             // read_cursor as usize,
                //             min_soft_clip_qual,
                //         ),
                //     )
                // }
            }
            Cigar::Ins(_) | Cigar::Del(_) => {
                if check_indels {
                    next_to_soft_clip = true
                }
            }
            _ => {
                // Not a soft clip
            }
        }
        return next_to_soft_clip;
    }

    fn count_high_quality_soft_clips(
        // cig: &Cigar,
        record: &Record,
        // cigar_cursor: usize,
        min_soft_clip_qual: u8,
    ) -> f64 {
        // https://gatk.broadinstitute.org/hc/en-us/articles/360036227652?id=4147
        // If we have high quality soft clips then we want to
        // track them
        // get mean base quality
        let mut num_high_quality_soft_clips = 0.0;
        let mut align_pos = 0;

        for cig in record.cigar().iter() {
            match cig {
                Cigar::SoftClip(len) => {
                    for _ in 0..(*len as usize) {
                        let qual_pos = record.qual()[align_pos];
                        align_pos += 1;
                        if qual_pos > min_soft_clip_qual {
                            num_high_quality_soft_clips += 1.0
                        }
                    }
                }
                _ => {
                    if CigarUtils::cigar_consumes_read_bases(cig) {
                        align_pos += cig.len() as usize;
                    }
                }
            };
        }

        return num_high_quality_soft_clips;
    }

    fn update_heterozygous_likelihood(
        result: &mut RefVsAnyResult,
        likelihoodcount: usize,
        log10ploidy: f64,
        ref_likelihood: f64,
        non_ref_likelihood: f64,
    ) {
        // https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReferenceConfidenceModel.java
        // applyPileupElementRefVsNonRefLikelihoodAndCount
        // Homozygous likelihoods don't need the logSum trick.

        result.genotype_likelihoods[0] += ref_likelihood + log10ploidy;
        result.genotype_likelihoods[likelihoodcount - 1] += non_ref_likelihood + log10ploidy;
        // Heterozygous likelihoods need the logSum trick:

        let mut i = 1;
        let mut j = likelihoodcount - 2;
        while i < (likelihoodcount - 1) {
            result.genotype_likelihoods[i] += MathUtils::approximate_log10_sum_log10(
                ref_likelihood + (j as f64).log10(),
                non_ref_likelihood + (i as f64).log10(),
            );
            i += 1;
            j -= 1;
        }
    }

    /**
     *  this implements the isActive() algorithm described in docs/mutect/mutect.pdf
     *  the multiplicative factor is for the special case where we pass a singleton list
     *  of alt quals and want to duplicate that alt qual over multiple reads
     * @param nRef          ref read count
     * @param altQuals      Phred-scaled qualities of alt-supporting reads
     * @param repeatFactor  Number of times each alt qual is duplicated
     * @param afPrior       Beta prior on alt allele fraction
     * @return
     */
    pub fn log_likelihood_ratio(
        n_ref: usize,
        alt_quals: Vec<u8>,
        repeat_factor: usize,
        af_prior: Option<BetaDistributionShape>,
    ) -> f64 {
        let n_alt = repeat_factor * alt_quals.len();
        let n = n_ref + n_alt;

        let f_tilde_ratio = (digamma(n_ref as f64 + 1.) - digamma(n_alt as f64 + 1.)).exp();

        let read_sum = alt_quals
            .par_iter()
            .map(|qual| {
                let epsilon = QualityUtils::qual_to_error_prob(*qual);
                let z_bar_alt = (1.0 - epsilon) / (1.0 - epsilon + epsilon * f_tilde_ratio);
                let log_epsilon = NaturalLogUtils::qual_to_log_error_prob(*qual);
                let log_one_minus_epsilon = NaturalLogUtils::qual_to_log_prob(*qual);
                z_bar_alt * (log_one_minus_epsilon - log_epsilon)
                    + MathUtils::fast_bernoulli_entropy(z_bar_alt)
            })
            .sum::<f64>();

        let beta_entropy;
        match af_prior {
            Some(af_prior) => {
                let alpha = af_prior.get_alpha();
                let beta = af_prior.get_beta();
                beta_entropy = ln_gamma(alpha + beta)
                    - ln_gamma(alpha)
                    - ln_gamma(beta)
                    - ln_gamma(alpha + beta + n as f64)
                    + ln_gamma(alpha + n_alt as f64)
                    + ln_gamma(beta + n_ref as f64);
            }
            None => {
                beta_entropy = MathUtils::log10_to_log(
                    -MathUtils::log10_factorial(n as f64 + 1.0)
                        + MathUtils::log10_factorial(n_alt as f64)
                        + MathUtils::log10_factorial(n_ref as f64),
                );
            }
        }

        return beta_entropy + read_sum * repeat_factor as f64;
    }

    pub fn log_likelihood_ratio_constant_error(
        ref_count: usize,
        alt_count: usize,
        error_probability: f64,
    ) -> f64 {
        let qual = QualityUtils::error_prob_to_qual(error_probability);
        return Self::log_likelihood_ratio(ref_count, vec![qual], alt_count, None);
    }

    fn write_empty_vcf(
        &self,
        output_prefix: &str,
        sample_names: &[&str],
        reference_reader: &ReferenceReader,
        strain_info: bool,
    ) {
        // touch empy output file and return
        // initiate header
        let mut header = Header::new();
        // Add program info
        self.populate_vcf_header(sample_names, reference_reader, &mut header, strain_info);

        // ensure path exists
        create_dir_all(output_prefix).expect("Unable to create output directory");

        // Initiate writer
        let mut bcf_writer = Writer::from_path(
            format!(
                "{}/{}.vcf",
                output_prefix, &reference_reader.genomes_and_contigs.genomes[self.ref_idx],
            )
            .as_str(),
            &header,
            true,
            Format::Vcf, // uncompressed. Bcf compression seems busted?
        )
        .unwrap_or_else(|_| panic!("Unable to create VCF output: {}.vcf", output_prefix));
        
        let mut record = bcf_writer.empty_record();

        let mut phases = Vec::new();
        let mut pls = Vec::new();
        let mut ads = Vec::new();
        let mut gqs = Vec::new();
        let mut dps = Vec::new();
        for _ in 0..sample_names.len() {
            phases.extend(vec![GenotypeAllele::UnphasedMissing; 2]);
            pls.push(".");
            ads.push(".");
            dps.push(0);
            gqs.push(0);
        }

        record
            .push_genotypes(phases.as_slice())
            .expect("Unable to push genotypes");

        record
            .push_format_string(
                VariantAnnotations::PhredLikelihoods.to_key().as_bytes(),
                &pls.iter().map(|p| p.as_bytes()).collect::<Vec<&[u8]>>(),
            )
            .expect("Unable to push format tag");

        record
            .push_format_string(
                VariantAnnotations::DepthPerAlleleBySample
                    .to_key()
                    .as_bytes(),
                &ads.iter().map(|a| a.as_bytes()).collect::<Vec<&[u8]>>(),
            )
            .expect("Unable to push format tag");
        record
            .push_format_integer(
                VariantAnnotations::GenotypeQuality.to_key().as_bytes(),
                &gqs,
            )
            .expect("Unable to push format tag");
        record
            .push_format_integer(VariantAnnotations::Depth.to_key().as_bytes(), &dps)
            .expect("Unable to push format tag");

        bcf_writer.write(&record).expect("Unable to write empty record");
    }

    /// Takes a vector of VariantContexts and writes them to a single VCF4 file
    pub fn write_vcf(
        &self,
        output_prefix: &str,
        variant_contexts: &Vec<VariantContext>,
        sample_names: &[&str],
        reference_reader: &ReferenceReader,
        strain_info: bool,
    ) {
        if variant_contexts.len() == 0 {

            self.write_empty_vcf(
                output_prefix,
                sample_names,
                reference_reader,
                strain_info,
            );
            let out_file_name = format!(
                "{}/{}.vcf",
                output_prefix, &reference_reader.genomes_and_contigs.genomes[self.ref_idx],
            );
            let out_file_name_tmp = format!(
                "{}/{}.vcf.tmp",
                output_prefix, &reference_reader.genomes_and_contigs.genomes[self.ref_idx],
            );

            {
                // remove last line of file
                let file = File::open(out_file_name.as_str()).unwrap();
                let out_file = File::create(out_file_name_tmp.as_str()).unwrap();

                let reader = BufReader::new(&file);
                let mut writer = BufWriter::new(&out_file);
            
                for (_, line) in reader.lines().enumerate() {
                    let line = line.as_ref().unwrap();
                    if !line.contains("GT:PL:AD:GQ:DP") {
                        writeln!(writer, "{}", line).expect("Unable to write data to empty VCF");
                    }
                }
            }
            std::fs::rename(out_file_name_tmp.as_str(), out_file_name.as_str()).expect("Unable to rename VCF file");

            return;
        }

        // initiate header
        let mut header = Header::new();
        // Add program info
        self.populate_vcf_header(sample_names, reference_reader, &mut header, strain_info);

        // ensure path exists
        create_dir_all(output_prefix).expect("Unable to create output directory");

        // Initiate writer
        let mut bcf_writer = Writer::from_path(
            format!(
                "{}/{}.vcf",
                output_prefix, &reference_reader.genomes_and_contigs.genomes[self.ref_idx],
            )
            .as_str(),
            &header,
            true,
            Format::Vcf, // uncompressed. Bcf compression seems busted?
        )
        .unwrap_or_else(|_| panic!("Unable to create VCF output: {}.vcf", output_prefix));


        for vc in variant_contexts {
            vc.write_as_vcf_record(&mut bcf_writer, reference_reader, sample_names.len());
        }
    }

    fn populate_vcf_header(
        &self,
        sample_names: &[&str],
        reference_reader: &ReferenceReader,
        header: &mut Header,
        strain_info: bool,
    ) {
        header.push_record(format!("##source=lorikeet-v{}", env!("CARGO_PKG_VERSION")).as_bytes());

        // debug!("samples {:?}", &sample_names);
        for sample_idx in 0..sample_names.len() {
            // remove tmp file name from sample id
            header.push_record(
                format!(
                    "##sample=<ID={}, name={}>",
                    sample_idx + 1,
                    sample_names[sample_idx]
                )
                .as_bytes(),
            );
            header.push_sample(format!("{}", sample_idx + 1).as_bytes());
        }

        // Add contig info
        for tid in reference_reader
            .retrieve_tids_for_ref_index(self.ref_idx)
            .unwrap()
            .iter()
        {
            header.push_record(
                format!(
                    "##contig=<ID={}, length={}>",
                    std::str::from_utf8(reference_reader.get_target_name(*tid)).unwrap(),
                    reference_reader.target_lens.get(&tid).unwrap()
                )
                .as_bytes(),
            );
        }

        VariantAnnotationEngine::populate_vcf_header(header, strain_info);
    }
}

struct PosAlignment {
    qpos: Option<usize>,
    is_del: bool,
    is_refskip: bool,
}

impl PosAlignment {
    fn new(qpos: Option<usize>, is_del: bool, is_refskip: bool) -> Self {
        Self {
            qpos,
            is_del,
            is_refskip,
        }
    }
}
