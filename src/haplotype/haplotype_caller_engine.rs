use activity_profile::activity_profile::Profile;
use activity_profile::activity_profile_state::{ActivityProfileState, Type};
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use annotator::variant_annotator_engine::VariantAnnotationEngine;
use assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use assembly::assembly_region::AssemblyRegion;
use assembly::assembly_region_trimmer::AssemblyRegionTrimmer;
use assembly::assembly_result_set::AssemblyResultSet;
use bio::stats::{LogProb, PHREDProb};
use clap::ArgMatches;
use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::FlagFilter;
use genotype::genotype_builder::{Genotype, GenotypesContext};
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use genotype::genotyping_engine::GenotypingEngine;
use haplotype::haplotype::Haplotype;
use haplotype::haplotype_caller_genotyping_engine::HaplotypeCallerGenotypingEngine;
use haplotype::haplotype_clustering_engine::HaplotypeClusteringEngine;
use haplotype::ref_vs_any_result::RefVsAnyResult;
use itertools::Itertools;
use mathru::special::gamma::{digamma, ln_gamma};
use model::allele_likelihoods::AlleleLikelihoods;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::VariantContext;
use model::variants::*;
use pair_hmm::pair_hmm_likelihood_calculation_engine::{
    AVXMode, PairHMMLikelihoodCalculationEngine,
};
use processing::lorikeet_engine::Elem;
use processing::lorikeet_engine::ReadType;
use rayon::prelude::*;
use read_orientation::beta_distribution_shape::BetaDistributionShape;
use read_threading::read_threading_assembler::ReadThreadingAssembler;
use read_threading::read_threading_graph::ReadThreadingGraph;
use reads::alignment_utils::AlignmentUtils;
use reads::bird_tool_reads::BirdToolRead;
use reference::reference_reader::ReferenceReader;
use rust_htslib::bam::{self, record::Cigar, Record, pileup::Alignment};
use rust_htslib::bcf::{Format, Header, Writer};
use std::cmp::min;
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, read};
use std::sync::{Arc, Mutex};
use utils::math_utils::{MathUtils, RunningAverage};
use utils::natural_log_utils::NaturalLogUtils;
use utils::quality_utils::QualityUtils;
use utils::simple_interval::{Locatable, SimpleInterval};
use utils::utils::clean_sample_name;
use assembly::assembly_region_walker::AssemblyRegionWalker;
use hashlink::LinkedHashSet;
use num::traits::SaturatingSub;
use reads::cigar_utils::CigarUtils;
use std::ops::Deref;
use rust_htslib::bam::pileup::Indel;
use reads::read_utils::ReadUtils;
use bio_types::sequence::SequenceRead;
use utils::interval_utils::IntervalUtils;
use rust_htslib::bam::ext::BamRecordExtensions;
use mathru::algebra::abstr::Sign;

#[derive(Debug, Clone)]
pub struct HaplotypeCallerEngine<'c> {
    active_region_evaluation_genotyper_engine: GenotypingEngine,
    genotyping_engine: HaplotypeCallerGenotypingEngine,
    genotype_prior_calculator: GenotypePriorCalculator,
    assembly_region_trimmer: AssemblyRegionTrimmer,
    assembly_engine: ReadThreadingAssembler,
    likelihood_calculation_engine: PairHMMLikelihoodCalculationEngine<'c>,
    ref_idx: usize,
    stand_min_conf: f64,
    mapping_quality_threshold: u8,
}

impl<'c> HaplotypeCallerEngine<'c> {
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
    const MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY: f64 = 4.0;

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
    ) -> HaplotypeCallerEngine<'c> {
        let kmer_sizes = match args.values_of("kmer-sizes") {
            Some(vals) => {
                vals.map(|k_size| k_size.parse().unwrap())
                    .collect::<Vec<usize>>()
            },
            _ => vec![10, 25]
        };

        let mut assembly_engine = ReadThreadingAssembler::new(
            args.value_of("max-allowed-path-for-read-threading-assembler")
                .unwrap()
                .parse()
                .unwrap(),
            kmer_sizes,
            args.is_present("dont-increase-kmer-sizes-for-cycles"),
            args.is_present("allow-non-unique-kmers-in-ref"),
            args.value_of("num-pruning-samples")
                .unwrap()
                .parse()
                .unwrap(),
            if args.is_present("use-adaptive-pruning") {
                0
            } else {
                args.value_of("min-prune-factor").unwrap().parse().unwrap()
            },
            args.is_present("use-adaptive-pruning"),
            args.value_of("initial-error-rate-for-pruning")
                .unwrap()
                .parse()
                .unwrap(),
            MathUtils::log10_to_log(
                args.value_of("pruning-log-odds-threshold")
                    .unwrap()
                    .parse()
                    .unwrap(),
            ),
            MathUtils::log10_to_log(
                args.value_of("pruning-seeding-log-odds-threshold")
                    .unwrap()
                    .parse()
                    .unwrap(),
            ),
            args.value_of("max-unpruned-variants")
                .unwrap()
                .parse()
                .unwrap(),
            args.is_present("use-linked-debruijn-graph"),
            args.is_present("enable-legacy-graph-cycle-detection"),
            args.value_of("min-matching-bases-to-dangling-end-recovery")
                .unwrap()
                .parse()
                .unwrap(),
        );

        assembly_engine.debug_graph_transformations =
            args.is_present("debug-graph-transformations");
        assembly_engine.recover_dangling_branches =
            !args.is_present("do-not-recover-dangling-branches");
        assembly_engine.recover_all_dangling_branches =
            args.is_present("recover-all-dangling-branches");
        assembly_engine.min_dangling_branch_length = args
            .value_of("min-dangling-branch-length")
            .unwrap()
            .parse()
            .unwrap();
        assembly_engine.graph_output_path = match args.value_of("graph-output") {
            Some(path) => Some(path.to_string()),
            None => None,
        };
        assembly_engine.debug_graph_output_path = match args.value_of("debug-graph-output") {
            Some(path) => Some(path.to_string()),
            None => None,
        };
        assembly_engine.min_base_quality_to_use_in_assembly =
            args.value_of("min-base-quality").unwrap().parse().unwrap();

        HaplotypeCallerEngine {
            active_region_evaluation_genotyper_engine: GenotypingEngine::make(
                args,
                samples.clone(),
                do_allele_specific_calcs,
                std::cmp::max(sample_ploidy, Self::MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY),
            ),
            genotyping_engine: HaplotypeCallerGenotypingEngine::new(
                args,
                samples,
                !args.is_present("do-not-run-physical-phasing"),
                sample_ploidy,
            ),
            genotype_prior_calculator: GenotypePriorCalculator::make(args),
            stand_min_conf: args
                .value_of("standard-min-confidence-threshold-for-calling")
                .unwrap()
                .parse()
                .unwrap(),
            ref_idx,
            assembly_engine,
            assembly_region_trimmer: AssemblyRegionTrimmer::new(
                args.value_of("assembly-region-padding")
                    .unwrap()
                    .parse()
                    .unwrap(),
                args.value_of("indel-padding-for-genotyping")
                    .unwrap()
                    .parse()
                    .unwrap(),
                args.value_of("snp-padding-for-genotyping")
                    .unwrap()
                    .parse()
                    .unwrap(),
                args.value_of("str-padding-for-genotyping")
                    .unwrap()
                    .parse()
                    .unwrap(),
                args.value_of("max-extension-into-region-padding")
                    .unwrap()
                    .parse()
                    .unwrap(),
            ),
            likelihood_calculation_engine:
                AssemblyBasedCallerUtils::create_likelihood_calculation_engine(
                    args,
                    !args.is_present("soft-clip-low-quality-ends"),
                ),
            mapping_quality_threshold: args
                .value_of("mapping-quality-threshold-for-genotyping")
                .unwrap()
                .parse()
                .unwrap(),
        }
    }

    pub fn stand_min_conf(&self) -> f64 {
        self.stand_min_conf
    }

    pub fn collect_activity_profile_low_mem(
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
    ) -> Vec<(Vec<VariantContext>, Vec<Vec<i32>>)> {
        // minimum PHRED base quality
        let bq = m
            .value_of("min-base-quality")
            .unwrap()
            .parse::<u8>()
            .unwrap();

        let max_prob_prop = m
            .value_of("max-prob-propagation-distance")
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let active_prob_thresh = m
            .value_of("active-probability-threshold")
            .unwrap()
            .parse::<f32>()
            .unwrap();

        let min_contig_length = m
            .value_of("min-contig-size")
            .unwrap()
            .parse::<u64>()
            .unwrap();

        let limiting_interval = IntervalUtils::parse_limiting_interval(m);
        debug!("Limiting {:?}", &limiting_interval);

        let ploidy: usize = std::cmp::max(m.value_of("ploidy").unwrap().parse().unwrap(), Self::MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY);

        let mut tids: HashSet<usize> = HashSet::new();
        let reference = reference_reader.retrieve_reference_stem(ref_idx);

        indexed_bam_readers
            .iter()
            .enumerate()
            .for_each(|(sample_idx, bam_generator)| {
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
                            || reference_reader.match_target_name_and_ref_idx(ref_idx, target_name)
                        {
                            // Get contig stats
                            reference_reader.update_ref_index_tids(ref_idx, tid);
                            reference_reader.add_target(contig_name, tid);
                            let target_len = header.target_len(tid as u32).unwrap();
                            reference_reader.add_length(tid, target_len);
                            tids.insert(tid);
                        }
                    })
            });

        let chunk_size = 250000;
        let contexts = tids.into_iter().map(|tid| {
            let target_length = reference_reader.target_lens[&tid];
            // let mut reference_reader = reference_reader.clone();
            reference_reader.update_current_sequence_capacity(target_length as usize);
            // Update all contig information
            let retrieved = match reference_reader.fetch_contig_from_reference_by_tid(
                tid,
                ref_idx,
            ) {
                Ok(_) => true,
                Err(_) => false,
            };

            let context_depth_tuples = if retrieved {
                // let mut contexts = Vec::new();
                reference_reader.read_sequence_to_vec();
                if target_length >= min_contig_length {
                    let chunk_idx = 0;

                    (0..target_length as usize).into_par_iter().chunks(chunk_size).flat_map(|mut positions| {
                        let within_bounds = match &limiting_interval {
                            Some(limit) => {
                                let position_limit = SimpleInterval::new(0, positions[0], *positions.last().unwrap());
                                position_limit.overlaps(limit)
                            },
                            None => {
                                true
                            }
                        };

                        if within_bounds {
                            let mut genotype_likelihoods = Vec::new();
                            let first = positions[0];
                            let last = *positions.last().unwrap();
                            let length = last - first + 1;
                            let mut per_contig_per_base_hq_soft_clips = vec![RunningAverage::new(); length];
                            let chunk_location = SimpleInterval::new(
                                tid,
                                first,
                                last,
                            );

                            indexed_bam_readers
                                .iter()
                                .enumerate()
                                .for_each(|(sample_idx, bam_generator)| {
                                    // Get the appropriate sample index based on how many references we are using
                                    // let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                                    //     vec![&bam_generator],
                                    //     n_threads as u32,
                                    // )
                                    // .into_iter()
                                    // .next()
                                    // .unwrap();
                                    // Get the appropriate sample index based on how many references we are using
                                    let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                                        vec![&bam_generator],
                                        n_threads as u32,
                                    )
                                        .into_iter()
                                        .next()
                                        .unwrap();
                                    let mut bam_generated = bam_generator.start();

                                    let mut read_type = ReadType::Short;

                                    if (m.is_present("longreads") || m.is_present("longread-bam-files"))
                                        && sample_idx >= short_sample_count
                                        && sample_idx < (short_sample_count + long_sample_count)
                                    {
                                        // Get the appropriate sample index based on how many references we are using by tracking
                                        // changes in references
                                        read_type = ReadType::Long;
                                    }

                                    HaplotypeCallerEngine::update_activity_profile_low_mem(
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
                                    );
                                });
                            // chunk_idx += 1;
                            self.calculate_activity_probabilities_low_mem(
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
                            )
                        } else {
                            Vec::new()
                        }
                    }).collect::<Vec<(Vec<VariantContext>, Vec<Vec<i32>>)>>()
                } else {
                    Vec::new()
                }
                // contexts
            } else {
                Vec::new()
            };

            let mut contexts = Vec::with_capacity(target_length as usize);
            let mut passing_sites = vec![Vec::with_capacity(target_length as usize); indexed_bam_readers.len()];
            context_depth_tuples.into_iter().for_each(|result: (Vec<VariantContext>, Vec<Vec<i32>>)| {
                contexts.extend(result.0);
                // println!("depths {:?}", &result.1);
                for (i, mut depths) in result.1.into_iter().enumerate() {
                    let sample_passing_site = &mut passing_sites[i];
                    let passing_len = sample_passing_site.len();
                    if sample_passing_site.len() > 0 {
                        let first = match depths.first() {
                            Some(f) => *f,
                            None => continue,
                        };

                        let last = match sample_passing_site.last() {
                            Some(res) => *res,
                            _ => 0,
                        };

                        if last.sign() == first.sign() {
                            sample_passing_site[passing_len - 1] += first;
                        } else {
                            sample_passing_site.push(first);
                        }

                    }

                    if depths.len() > 1 {
                        sample_passing_site.extend(&depths[1..]);
                    }
                }

            });

            (contexts, passing_sites)
        }).collect::<Vec<(Vec<VariantContext>, Vec<Vec<i32>>)>>();

        // return genotype_likelihoods for each contig in current genome across samples
        contexts
    }

    pub fn update_activity_profile_low_mem<
        'b,
    >(
        bam_generated: &mut IndexedBamFileNamedReader,
        split_threads: usize,
        ref_idx: usize,
        readtype: ReadType,
        ploidy: usize,
        bq: u8,
        genomes_and_contigs: &'b GenomesAndContigs,
        concatenated_genomes: &'b Option<String>,
        flag_filters: &'b FlagFilter,
        per_contig_per_base_hq_soft_clips: &mut Vec<RunningAverage>,
        reference_reader: &ReferenceReader,
        limiting_interval: &Option<SimpleInterval>,
        current_likelihoods: &mut Vec<Vec<RefVsAnyResult>>,
        min_contig_length: u64,
        tid: usize,
        outer_chunk_location: &SimpleInterval,
        outer_chunk_index: usize,
        target_len: u64
    ) {

        let likelihoodcount = ploidy + 1;
        let log10ploidy = (ploidy as f64).log10();
        let chunk_size = 1000;

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        match readtype {
            ReadType::Short | ReadType::Long => {
                // The raw activity profile.
                // Frequency of bases not matching reference compared
                // to depth
                let mut likelihoods = (0..outer_chunk_location.size()).into_iter()
                        .map(|pos| {
                            RefVsAnyResult::new(likelihoodcount, pos, tid)
                        }).collect::<Vec<RefVsAnyResult>>();

                let outer_chunk_size = outer_chunk_location.size();

                let mut positions = per_contig_per_base_hq_soft_clips.iter_mut().zip(likelihoods.iter_mut())
                    .map(|(r, l)| (r, l)).collect::<Vec<(&mut RunningAverage, &mut RefVsAnyResult)>>();

                // multiplier to help us map between chunk position
                // and actual reference position. This value represents the
                // starting reference base index of this chunk.
                // let chunk_multiplier = chunk_idx * chunk_size;

                bam_generated
                    .fetch((
                        tid as i32,
                        outer_chunk_location.start as i64,
                        min(outer_chunk_location.end + 1, target_len as usize - 1) as i64)
                    ).unwrap_or_else(|_|
                        panic!(
                            "Failed to fetch interval {}:{}-{}", tid,
                            outer_chunk_location.start,
                            min(outer_chunk_location.end + 1, target_len as usize - 1)
                        )
                    );

                let mut record = Record::new();
                // pileup does not provide all of the alignments at a pos
                // Unsure why, but there are positions where entire
                // alignments are being ignored. Solution: create our
                // own version of a pileup and hopefully records won't be
                // filtered
                while bam_generated.read(&mut record) {
                    if ReadUtils::read_is_filtered(&record, &flag_filters, 20, readtype, limiting_interval)
                    {
                        continue;
                    }
                    // debug!("start {} end {} chunk {}-{}",
                    //        record.reference_start(),
                    //        record.reference_end(),
                    //        min(outer_chunk_location.start, target_len as usize),
                    //        min(outer_chunk_location.end,
                    //            target_len as usize),
                    // );
                    Self::parse_record(
                        &record,
                        flag_filters,
                        readtype,
                        &mut positions,
                        min(outer_chunk_location.start,
                            target_len as usize),
                        bq,
                        likelihoodcount,
                        reference_reader,
                        log10ploidy,
                        outer_chunk_location.start,
                        min(outer_chunk_location.end + 1,
                            target_len as usize),
                        false
                    );
                }
                Self::update_ref_vs_any_results(&mut likelihoods, likelihoodcount, log10ploidy, outer_chunk_location.start, false);
                current_likelihoods.push(likelihoods);
            },
            _ => {}
        }
    }

    fn update_ref_vs_any_results(
        likelihoods: &mut Vec<RefVsAnyResult>,
        likelihoodcount: usize,
        log10ploidy: f64,
        outer_chunk_start: usize,
        debug: bool
    ) {
        for (pos, result) in likelihoods.iter_mut().enumerate() {
            let denominator = result.read_counts as f64 * log10ploidy;

            for i in 0..likelihoodcount {
                result.genotype_likelihoods[i] -= denominator
            }
            if result.read_counts > 0 {
                debug!("Pos {} counts {} likelihood {:?} ref count {} non ref {}",
                       pos + outer_chunk_start, result.read_counts, &result.genotype_likelihoods, result.ref_depth, result.non_ref_depth);
            }
        }
    }

    fn parse_record(
        record: &Record,
        flag_filters: &FlagFilter,
        readtype: ReadType,
        positions: &mut Vec<(&mut RunningAverage, &mut RefVsAnyResult)>,
        subtractor: usize,
        bq: u8,
        likelihoodcount: usize,
        reference_reader: &ReferenceReader,
        log10ploidy: f64,
        bound_start: usize,
        bound_end: usize,
        debug: bool,
    ) -> usize {
        // we need to iterate through each read pos (qpos)
        // and generate a pseudo pileup
        let read_length = record.len();
        let mut cigar_cursor = 0; // changes when cigar is consumed
        let mut read_cursor = 0; // changes when read bases are consumed
        let mut pos = record.pos() as usize; // read start alignment on reference.
        // updated as ref consumed
        let mut count = 0;
        let mut cig_index = 0;
        for cig in &record.cigar() {

            match cig {
                Cigar::Del(len) => { // reference bases consumed

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
                            break
                        };
                        let alignment = PosAlignment::new(None, true, false);
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
                        pos += 1;
                    }
                    cigar_cursor += *len as usize;
                },
                Cigar::RefSkip(len) => {
                    panic!("Read contains N operator, should have been filtered prior to this point.");
                    std::process::exit(1);
                },
                Cigar::Ins(len) => { // read bases consumed
                    // if pos == 97055 {
                    //     count += 1;
                    // }
                    if pos < bound_start {
                        // read not within bounds yet
                        cigar_cursor += *len as usize;
                        read_cursor += *len as usize;
                        continue;
                    } else if pos >= bound_end {
                        // read pass the bounds now so break
                        break
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
                    cigar_cursor += *len as usize;
                    read_cursor += *len as usize;
                },
                Cigar::Diff(len)
                | Cigar::Match(len)
                | Cigar::Equal(len) => {
                    // we need to check each position in
                    // these cigars
                    for _ in 0..*len as usize {
                        // if pos == 97055 {
                        //     count += 1;
                        // }
                        if pos < bound_start {
                            // read not within bounds yet
                            read_cursor += 1;
                            pos += 1;
                            continue
                        } else if pos >= bound_end {
                            // read pass the bounds now so break
                            break
                        }

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
                        read_cursor += 1;
                        pos += 1;
                    }
                    cigar_cursor += *len as usize;
                },
                Cigar::SoftClip(len) => {
                    cigar_cursor += *len as usize;
                    read_cursor += *len as usize;
                },
                Cigar::HardClip(_)
                | Cigar::Pad(_) => {
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
        tree: &Arc<Mutex<Vec<&Elem>>>,
        reference_reader: &mut ReferenceReader,
    ) -> HashMap<usize, Vec<BandPassActivityProfile>> {
        // minimum PHRED base quality
        let bq = m
            .value_of("min-base-quality")
            .unwrap()
            .parse::<u8>()
            .unwrap();

        let max_prob_prop = m
            .value_of("max-prob-propagation-distance")
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let active_prob_thresh = m
            .value_of("active-probability-threshold")
            .unwrap()
            .parse::<f32>()
            .unwrap();

        let min_contig_length = m
            .value_of("min-contig-size")
            .unwrap()
            .parse::<u64>()
            .unwrap();

        let limiting_interval = IntervalUtils::parse_limiting_interval(m);

        let ploidy: usize = std::cmp::max(m.value_of("ploidy").unwrap().parse().unwrap(), Self::MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY);

        let mut genotype_likelihoods = Vec::new();

        let mut per_contig_per_base_hq_soft_clips = HashMap::new();

        indexed_bam_readers
            .iter()
            .enumerate()
            .for_each(|(sample_idx, bam_generator)| {
                // Get the appropriate sample index based on how many references we are using
                // let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                //     vec![&bam_generator],
                //     n_threads as u32,
                // )
                // .into_iter()
                // .next()
                // .unwrap();
                if sample_idx < short_sample_count {
                    HaplotypeCallerEngine::update_activity_profile(
                        bam_generator,
                        n_threads,
                        ref_idx,
                        ReadType::Short,
                        ploidy,
                        bq,
                        genomes_and_contigs,
                        concatenated_genomes,
                        flag_filters,
                        &mut per_contig_per_base_hq_soft_clips,
                        reference_reader,
                        &limiting_interval,
                        &mut genotype_likelihoods,
                        min_contig_length,
                    );
                } else if (m.is_present("longreads") || m.is_present("longread-bam-files"))
                    && sample_idx >= short_sample_count
                    && sample_idx < (short_sample_count + long_sample_count)
                {
                    // Get the appropriate sample index based on how many references we are using by tracking
                    // changes in references
                    HaplotypeCallerEngine::update_activity_profile(
                        bam_generator,
                        n_threads,
                        ref_idx,
                        ReadType::Long,
                        ploidy,
                        bq,
                        genomes_and_contigs,
                        concatenated_genomes,
                        flag_filters,
                        &mut per_contig_per_base_hq_soft_clips,
                        reference_reader,
                        &limiting_interval,
                        &mut genotype_likelihoods,
                        min_contig_length,
                    );
                }
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];

                    pb.progress_bar.set_message(format!(
                        "{}: Finding active regions in sample: {:.50}...",
                        pb.key,
                        clean_sample_name(sample_idx, indexed_bam_readers),
                    ));
                    // pb.progress_bar.inc(1);
                }
            });

        {
            let pb = &tree.lock().unwrap()[ref_idx + 2];

            pb.progress_bar
                .set_message(format!("{}: Calculating activity probabilities...", pb.key,));
        }
        // return genotype_likelihoods for each contig in current genome across samples
        return self.calculate_activity_probabilities(
            genotype_likelihoods,
            per_contig_per_base_hq_soft_clips,
            &reference_reader.target_lens,
            ploidy,
            max_prob_prop,
            active_prob_thresh,
            ref_idx,
            indexed_bam_readers,
            min_contig_length,
        );
    }

    pub fn update_activity_profile<
        'b,
        // R: IndexedNamedBamReader + Send,
        // G: NamedBamReaderGenerator<R> + Send,
    >(
        bam_generator: &str,
        split_threads: usize,
        ref_idx: usize,
        readtype: ReadType,
        ploidy: usize,
        bq: u8,
        genomes_and_contigs: &'b GenomesAndContigs,
        concatenated_genomes: &'b Option<String>,
        flag_filters: &'b FlagFilter,
        per_contig_per_base_hq_soft_clips: &mut HashMap<usize, Vec<RunningAverage>>,
        reference_reader: &mut ReferenceReader,
        limiting_interval: &Option<SimpleInterval>,
        current_likelihoods: &mut Vec<HashMap<usize, Vec<RefVsAnyResult>>>,
        min_contig_length: u64,
    ) {
        let mut tids: Vec<usize> = Vec::new();
        let reference = reference_reader.retrieve_reference_stem(ref_idx);

        {
            // Get the appropriate sample index based on how many references we are using
            let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                vec![&bam_generator],
                split_threads as u32,
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
                        || reference_reader.match_target_name_and_ref_idx(ref_idx, target_name)
                    {
                        // Get contig stats
                        reference_reader.update_ref_index_tids(ref_idx, tid);
                        reference_reader.add_target(contig_name, tid);
                        let target_len = header.target_len(tid as u32).unwrap();
                        reference_reader.add_length(tid, target_len);
                        tids.push(tid);
                    }
                })
        }

        let likelihoodcount = ploidy + 1;
        let log10ploidy = (ploidy as f64).log10();
        let mut genotype_likelihoods = HashMap::new();
        let chunk_size = 100000;

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        match readtype {
            ReadType::Short | ReadType::Long => {
                tids
                    .into_iter()
                    .for_each(|tid| {
                        {
                            let target_len = reference_reader.target_lens[&tid];
                            if target_len >= min_contig_length {

                                reference_reader.update_current_sequence_capacity(target_len as usize);
                                // Update all contig information
                                let retrieved = match reference_reader.fetch_contig_from_reference_by_tid(
                                    tid,
                                    ref_idx as usize,
                                ) {
                                    Ok(_) => true,
                                    _ => false,
                                };

                                if retrieved {
                                    reference_reader.read_sequence_to_vec();

                                    let per_base_hq_soft_clips = per_contig_per_base_hq_soft_clips.entry(tid).or_insert_with(|| vec![RunningAverage::new(); target_len as usize]);
                                    // The raw activity profile.
                                    // Frequency of bases not matching reference compared
                                    // to depth
                                    let likelihoods = genotype_likelihoods.entry(tid)
                                        .or_insert_with(|| (0..target_len as usize).into_iter()
                                            .map(|pos| {
                                                RefVsAnyResult::new(likelihoodcount, pos, tid)
                                            }).collect::<Vec<RefVsAnyResult>>());

                                    per_base_hq_soft_clips.par_iter_mut().zip(likelihoods.par_iter_mut())
                                        .chunks(chunk_size)
                                        .enumerate()
                                        .for_each(|(chunk_idx, mut positions)| {

                                            // multiplier to help us map between chunk position
                                            // and actual reference position. This value represents the
                                            // starting reference base index of this chunk.
                                            let chunk_multiplier = chunk_idx * chunk_size;

                                            // Get the appropriate sample index based on how many references we are using
                                            let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                                                vec![&bam_generator],
                                                split_threads as u32,
                                            )
                                                .into_iter()
                                                .next()
                                                .unwrap();


                                            let mut bam_generated = bam_generator.start();
                                            bam_generated
                                                .fetch((
                                                    tid as i32,
                                                    chunk_multiplier as i64,
                                                    min(chunk_multiplier + chunk_size, target_len as usize) as i64)
                                                ).unwrap_or_else(|_|
                                                panic!(
                                                    "Failed to fetch interval {}:{}-{}", tid,
                                                    chunk_multiplier as i64,
                                                    min(
                                                        chunk_multiplier + chunk_size,
                                                        target_len as usize
                                                    ) as i64
                                                )
                                            );
                                            let mut record = Record::new();
                                            // pileup does not provide all of the alignments at a pos
                                            // Unsure why, but there are positions where entire
                                            // alignments are being ignored. Solution: create our
                                            // own version of a pileup and hopefully records won't be
                                            // filtered
                                            let mut debug = false;
                                            // if (chunk_multiplier <= 97200 && chunk_multiplier + chunk_size >= 97200) && target_len >= 1000000 {
                                            //     debug = true;
                                            // }

                                            let mut filtered_count = 0;
                                            let mut parsed = 0;
                                            let mut hit_position = 0;

                                            while bam_generated.read(&mut record) {
                                                if ReadUtils::read_is_filtered(&record, &flag_filters, 20, readtype, limiting_interval)
                                                {
                                                    filtered_count += 1;
                                                    continue;
                                                }

                                                hit_position += Self::parse_record(
                                                    &record,
                                                    flag_filters,
                                                    readtype,
                                                    &mut positions,
                                                    chunk_multiplier,
                                                    bq,
                                                    likelihoodcount,
                                                    reference_reader.deref(),
                                                    log10ploidy,
                                                    chunk_multiplier,
                                                    min(chunk_multiplier + chunk_size, target_len as usize),
                                                    debug
                                                );
                                                parsed += 1;
                                            }
                                            // if debug {
                                            //     println!("read count in region {chunk_multiplier} - {} {filtered_count} {parsed} {}: hit {hit_position}", chunk_multiplier + chunk_size, filtered_count + parsed);
                                            // }
                                        });

                                    Self::update_ref_vs_any_results(
                                        likelihoods, likelihoodcount, log10ploidy, 0, false
                                    );
                                };
                            };
                        }
                    });
            }
            _ => {}
        }
        current_likelihoods.push(genotype_likelihoods);
    }

    /**
     * Create a collection of ActivityProfileState for each position on each contig
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     *
     * @return HashMap<usize, BandPassActivityProfile> - contig id and per base activity
     */
    pub fn calculate_activity_probabilities(
        &mut self,
        genotype_likelihoods: Vec<HashMap<usize, Vec<RefVsAnyResult>>>,
        per_contig_per_base_hq_soft_clips: HashMap<usize, Vec<RunningAverage>>,
        target_ids_and_lens: &HashMap<usize, u64>,
        ploidy: usize,
        max_prob_propagation: usize,
        active_prob_threshold: f32,
        ref_idx: usize,
        sample_names: &[String],
        min_contig_length: u64,
    ) -> HashMap<usize, Vec<BandPassActivityProfile>> {
        if genotype_likelihoods.len() == 0 { // don't use this, it seems slower?
            // Faster implementation for single sample analysis
            let per_contig_activity_profiles: HashMap<usize, Vec<BandPassActivityProfile>> =
                genotype_likelihoods.into_iter().next().unwrap()
                    .into_par_iter()
                    .map(|(tid, vec_of_ref_vs_any_result)| {
                        // Create bandpass
                        let mut activity_profile = BandPassActivityProfile::new(
                            max_prob_propagation,
                            active_prob_threshold,
                            BandPassActivityProfile::MAX_FILTER_SIZE,
                            BandPassActivityProfile::DEFAULT_SIGMA,
                            true,
                            ref_idx,
                            tid,
                            *target_ids_and_lens.get(&tid).unwrap() as usize,
                        );

                        // for each position determine per locus activity and add to bandpass
                        vec_of_ref_vs_any_result.into_iter().enumerate().for_each(
                            |(pos, ref_vs_any_result)| {
                                let is_active_prob = self
                                    .active_region_evaluation_genotyper_engine
                                    .allele_frequency_calculator
                                    .calculate_single_sample_biallelic_non_ref_posterior(
                                        &ref_vs_any_result.genotype_likelihoods,
                                        true,
                                    );
                                let per_base_hq_soft_clips =
                                    per_contig_per_base_hq_soft_clips.get(&tid).unwrap();

                                let hq_soft_clips = per_base_hq_soft_clips[pos];

                                let activity_profile_state = ActivityProfileState::new(
                                    ref_vs_any_result.loc,
                                    is_active_prob as f32,
                                    Type::new(
                                        hq_soft_clips.mean() as f32,
                                        HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                                    )
                                );

                                activity_profile.add(activity_profile_state);
                            },
                        );

                        (tid, vec![activity_profile])
                    })
                    .collect::<HashMap<usize, Vec<BandPassActivityProfile>>>();

            return per_contig_activity_profiles;
        } else {
            // let mut per_contig_activity_profiles = HashMap::new();
            let placeholder_vec = Vec::new();
            let per_contig_activity_profiles = target_ids_and_lens
                .par_iter()
                .filter(|(_, length)| *length >= &min_contig_length)
                .map(|(tid, length)| {
                    debug!("Calculating activity on {} of length {}", tid, length);
                    let per_base_hq_soft_clips =
                        per_contig_per_base_hq_soft_clips.get(tid).unwrap();

                    // Create bandpass
                    debug!("Created bandpass profile");

                    let activity_profile =
                        (0..(*length as usize)).into_par_iter().chunks(25000).map(|positions| {
                            let mut active_region_evaluation_genotyper_engine =
                                self.active_region_evaluation_genotyper_engine.clone();
                            let mut activity_profile = BandPassActivityProfile::new(
                                max_prob_propagation,
                                active_prob_threshold,
                                BandPassActivityProfile::MAX_FILTER_SIZE,
                                BandPassActivityProfile::DEFAULT_SIGMA,
                                true,
                                ref_idx,
                                *tid,
                                *target_ids_and_lens.get(&tid).unwrap() as usize,
                            );
                            for pos in positions {
                                let mut genotypes = Vec::new();

                                let hq_soft_clips = per_base_hq_soft_clips[pos];
                                let debug = genotype_likelihoods.iter().any(|sample_likelihoods| {
                                    sample_likelihoods[tid][pos].read_counts > 0
                                });

                                if !debug {
                                    let activity_profile_state = ActivityProfileState::new(
                                        SimpleInterval::new(*tid, pos, pos),
                                        0.0,
                                        Type::new(
                                            hq_soft_clips.mean() as f32,
                                            HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                                        ),
                                    );
                                    activity_profile.add(activity_profile_state);
                                    continue;
                                }

                                for (idx, sample_likelihoods) in genotype_likelihoods.iter().enumerate() {
                                    let result = sample_likelihoods[tid][pos].genotype_likelihoods.clone();
                                    genotypes.push(Genotype::build(
                                        ploidy,
                                        result,
                                        sample_names[idx].clone(),
                                    ))
                                }

                                let fake_alleles = ByteArrayAllele::create_fake_alleles();

                                let mut variant_context =
                                    VariantContext::build(*tid, pos, pos, fake_alleles);

                                variant_context.add_genotypes(genotypes);

                                let vc_out = active_region_evaluation_genotyper_engine.calculate_genotypes(
                                    variant_context,
                                    ploidy,
                                    &self.genotype_prior_calculator,
                                    &placeholder_vec,
                                    self.stand_min_conf,
                                );

                                let vc_is_some = vc_out.is_some();
                                let is_active_prob = match vc_out {
                                    Some(vc) => {
                                        QualityUtils::qual_to_prob(vc.get_phred_scaled_qual() as u8)
                                    }
                                    None => 0.0,
                                };

                                if debug {
                                    debug!("pos {} active {} vc is some {}", pos, is_active_prob, vc_is_some);
                                }

                                let activity_profile_state = ActivityProfileState::new(
                                    SimpleInterval::new(*tid, pos, pos),
                                    is_active_prob as f32,
                                    Type::new(
                                        hq_soft_clips.mean() as f32,
                                        HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                                    ),
                                );
                                activity_profile.add(activity_profile_state);
                            }
                            activity_profile
                        }).collect::<Vec<BandPassActivityProfile>>();


                    debug!("Finished {} of length {}", tid, length);
                    (*tid, activity_profile)
                }).collect::<HashMap<usize, Vec<BandPassActivityProfile>>>();
            return per_contig_activity_profiles;
        }
    }

    /**
     * Create a collection of ActivityProfileState for each position on each contig
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     *
     * @return HashMap<usize, BandPassActivityProfile> - contig id and per base activity
     */
    pub fn calculate_activity_probabilities_low_mem(
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
        min_contig_length: u64,
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
        chunk_index: usize,
        output_prefix: &str,
    ) -> Vec<(Vec<VariantContext>, Vec<Vec<i32>>)> {

        // let mut per_contig_activity_profiles = HashMap::new();
        let placeholder_vec = Vec::new();
        let depth_per_sample_filter = args.value_of("depth-per-sample-filter").unwrap().parse().unwrap();

        let inner_chunk_size = 50000;
        let variant_contexts =
            (0..chunk_location.size()).into_par_iter().chunks(inner_chunk_size).enumerate().map(|(i, positions)| {
                let within_bounds = match &limiting_interval {
                    Some(limit) => {
                        let position_limit = SimpleInterval::new(0, chunk_location.start + i * inner_chunk_size, chunk_location.end + i * inner_chunk_size);
                        debug!("Position limit {:?}", &position_limit);
                        position_limit.overlaps(limit)
                    },
                    None => {
                        true
                    }
                };

                if within_bounds {
                    let mut active_region_evaluation_genotyper_engine =
                        self.active_region_evaluation_genotyper_engine.clone();

                    // Create bandpass
                    debug!("Created bandpass profile");
                    debug!("Calculating activity on {} of {:?} {} {}", tid, chunk_location, positions[0], positions.last().unwrap());
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

                    let mut depth_of_position = vec![Vec::with_capacity(inner_chunk_size); sample_names.len()];
                    let mut depths_counters = vec![0; sample_names.len()];
                    for pos in positions {
                        match &limiting_interval {
                            Some(limit) => {
                                let position_limit = SimpleInterval::new(0, chunk_location.start + pos, chunk_location.start + pos);
                                if !position_limit.overlaps(limit) {
                                    continue
                                }
                            },
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
                                if depths_counters[idx] >= 0 { // positively increment current sample counter
                                    depths_counters[idx] += 1;
                                } else {
                                    depth_of_position[idx].push(depths_counters[idx]); // push previous stretch of negative bases
                                    depths_counters[idx] = 0;
                                    depths_counters[idx] += 1;
                                }
                            } else {
                                if depths_counters[idx] <= 0 { // negatively increment when a base does not have suffcient coverage
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
                                sample_names[idx].clone(),
                            ))
                        }

                        let fake_alleles = ByteArrayAllele::create_fake_alleles();

                        let contig_position = chunk_location.start + pos;
                        let mut variant_context =
                            VariantContext::build(tid, contig_position, contig_position, fake_alleles);

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

                        debug!("{}-{} Active Prob {}", chunk_location.start + pos, chunk_location.start + pos, is_active_prob);

                        let activity_profile_state = ActivityProfileState::new(
                            SimpleInterval::new(tid, chunk_location.start + pos, chunk_location.start + pos),
                            is_active_prob as f32,
                            Type::new(
                                hq_soft_clips.mean() as f32,
                                HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                            ),
                        );
                        activity_profile.add(activity_profile_state);
                    }

                    for (idx, sample_depth) in depths_counters.into_iter().enumerate() {
                        depth_of_position[idx].push(sample_depth);
                    }
                    let mut inner_reader = ReferenceReader::new_from_reader_with_tid_and_rid(
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
                    (processed, depth_of_position)
                } else {
                    (Vec::new(), Vec::new())
                }
                // activity_profile
            }).collect::<Vec<(Vec<VariantContext>, Vec<Vec<i32>>)>>();


        debug!("Finished {} of length {}", tid, length);

        variant_contexts
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
        debug!(
            "Region of size {} with {} reads",
            region.padded_span.size(),
            region.reads.len()
        );
        if !region.is_active() && given_alleles.is_empty() {
            debug!("Region was not active");
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }


        if given_alleles.is_empty() && region.len() == 0 {
            debug!("Region was of length 0");
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }

        // let debug = region.padded_span.start <= 1709882 && region.padded_span.end >= 1709882;

        // if debug {
        debug!("loc {:?} padded Loc {:?} reads {}", &region.active_span, &region.padded_span, region.reads.len());
        // }
        let mut region_without_reads = region.clone_without_reads();

        // run the local assembler, getting back a collection of information on how we should proceed
        let mut untrimmed_assembly_result = AssemblyBasedCallerUtils::assemble_reads(
            region,
            &given_alleles,
            args,
            reference_reader,
            &mut self.assembly_engine,
            !args.is_present("do-not-correct-overlapping-base-qualities"),
            sample_names,
        );

        let all_variation_events = match untrimmed_assembly_result
            .get_variation_events(args.value_of("max-mnp-distance").unwrap().parse().unwrap()) {
            Ok(result) => result,
            Err(_) => {
                return self.reference_model_for_no_variation(&mut untrimmed_assembly_result.region_for_genotyping, true, &vc_priors)
            }
        };

        debug!("Region {:?} All variation events  {:?}", &untrimmed_assembly_result.padded_reference_loc, &all_variation_events.len());

        let mut trimming_result = self.assembly_region_trimmer.trim(
            self.ref_idx,
            region_without_reads,
            all_variation_events,
            args.is_present("enable-legacy-assembly-region-trimming"),
            &reference_reader,
            untrimmed_assembly_result
                .full_reference_with_padding
                .as_slice(),
        );

        debug!("Trim complete!");

        if !trimming_result.is_variation_present() && !args.is_present("disable-optimizations") {
            return self.reference_model_for_no_variation(
                &mut trimming_result.original_region,
                false,
                &vc_priors,
            );
        }

        debug!("Moving reads....");
        trimming_result.original_region.reads =
            untrimmed_assembly_result.region_for_genotyping.move_reads();
        debug!("Move complete! {}", trimming_result.original_region.reads.len());
        let mut assembly_result =
            untrimmed_assembly_result.trim_to(trimming_result.get_variant_region());

        debug!(
            "Assembly result allele order after region trimming {:?}",
            &assembly_result.haplotypes
        );

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
        let mut assembly_result = assembly_result.remove_all(&read_stubs);
        debug!(
            "Assembly result allele order after stub filter {:?} -> read stubs {}",
            &assembly_result.haplotypes.len(), read_stubs.len()
        );


        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        let (mut assembly_result, filtered_reads) =
            self.filter_non_passing_reads(assembly_result, flag_filters);
        // let filtered_reads = Vec::new();
        debug!("Filtered reads {}", filtered_reads.len());
        debug!(
            "Assembly result allele order after read filter {:?}",
            &assembly_result.haplotypes.len()
        );
        let per_sample_filtered_read_list =
            AssemblyBasedCallerUtils::split_reads_by_sample(filtered_reads);

        if !assembly_result.variation_present || assembly_result.region_for_genotyping.len() == 0 {
            return self.reference_model_for_no_variation(
                &mut assembly_result.region_for_genotyping,
                false,
                &vc_priors,
            );
        };

        debug!("==========================================================================");
        debug!(
            "              Assembly region {:?} \n\
                           Alleles {:?}",
            &assembly_result.region_for_genotyping, &assembly_result.haplotypes
        );
        debug!("==========================================================================");
        let reads = AssemblyBasedCallerUtils::split_reads_by_sample(
            assembly_result.region_for_genotyping.move_reads(),
        );

        debug!("Before change: {}", reads.values().map(|r| r.len()).sum::<usize>());
        let mut read_likelihoods: AlleleLikelihoods<Haplotype<SimpleInterval>> = self
            .likelihood_calculation_engine
            .compute_read_likelihoods(&mut assembly_result, sample_names.to_vec(), reads);

        // if debug {
        debug!("Read by sample after compute: {:?}", (0..sample_names.len()).map(|s| read_likelihoods.sample_evidence_count(s)).collect::<Vec<usize>>());
        // }

        debug!(
            "Read likelihoods first {:?} {:?} {:?}",
            read_likelihoods.alleles.len(),
            read_likelihoods.alleles.list.iter().map(|a| a.is_ref()).collect::<Vec<bool>>(),
            read_likelihoods.alleles.list.iter().map(|a| std::str::from_utf8(a.get_bases()).unwrap()).collect::<Vec<&str>>()
        );
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
            if args.is_present("disable-avx") {
                AVXMode::None
            } else {
                AVXMode::detect_mode()
            },
        );
        read_likelihoods.change_evidence(read_alignments);

        // if debug {
        debug!("After change {:?}", (0..sample_names.len()).map(|s| read_likelihoods.sample_evidence_count(s)).collect::<Vec<usize>>());
        // }

        debug!(
            "Read likelihoods after change evidence {:?}",
            read_likelihoods.alleles.len()
        );

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]
        let called_haplotypes = match self.genotyping_engine.assign_genotype_likelihoods(
            assembly_result
                .haplotypes
                .clone(),
            read_likelihoods,
            per_sample_filtered_read_list,
            assembly_result.full_reference_with_padding.as_slice(),
            &assembly_result.padded_reference_loc,
            &assembly_result.region_for_genotyping.active_span,
            given_alleles,
            false,
            args.value_of("max-mnp-distance").unwrap().parse().unwrap(),
            sample_names,
            args.value_of("ploidy").unwrap().parse().unwrap(),
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
        mut assembly_result: AssemblyResultSet<ReadThreadingGraph>,
        flag_filters: &FlagFilter,
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
        region: &'a mut AssemblyRegion,
        needs_to_be_finalized: bool,
        vc_priors: &Vec<VariantContext>,
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
        let record_qual = if alignment.is_del || alignment.is_refskip { Self::REF_MODEL_DELETION_QUAL } else { record.qual()[qpos.unwrap()] };
        let mut is_alt = false;

        if record_qual >= bq || alignment.is_del {
            result.read_counts += 1;

            is_alt = Self::is_alt(
                &record,
                qpos,
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
                non_ref_likelihood
            );
        }

        // if debug {
        //     println!(
        //         "{} {} {} {} {} {} {} {:?}",
        //         std::str::from_utf8(record.qname()).unwrap(),
        //         is_alt,
        //         if !(alignment.is_del || alignment.is_refskip) {
        //             record.seq()[*qpos.as_ref().unwrap()].to_ascii_uppercase() != refr_base.to_ascii_uppercase()
        //         } else { true },
        //         if !(alignment.is_del || alignment.is_refskip) {
        //             Self::next_to_soft_clip_or_indel(
        //                 record,
        //                 *qpos.as_ref().unwrap(),
        //                 None,
        //                 false,
        //                 Self::HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD,
        //             )
        //         } else { true },
        //         if !(alignment.is_del || alignment.is_refskip) {
        //             Self::next_to_soft_clip_or_indel(
        //                 record,
        //                 *qpos.as_ref().unwrap(),
        //                 None,
        //                 true,
        //                 Self::HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD,
        //             )
        //         } else { true },
        //         refr_base as char,
        //         record.seq()[*qpos.as_ref().unwrap()] as char,
        //         &qpos
        //     );
        // }
        // add hq soft clips if possible
        if is_alt
            && Self::next_to_soft_clip(record, cig_index, alignment.qpos) {
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
            hq_soft_clips.add(
                soft_clips
            );
        }
    }

    fn next_to_soft_clip(
        record: &Record,
        cig_index: usize,
        qpos: Option<usize>
    ) -> bool {
        match qpos {
            Some(qpos) => {
                Self::next_to_soft_clip_or_indel(
                    record,
                    qpos,
                    false,
                )
            },
            None => {
                CigarUtils::cigar_is_soft_clip(&record.cigar().0[cig_index.saturating_sub(1)])
                || CigarUtils::cigar_is_soft_clip(&record.cigar().0[min(cig_index + 1, record.cigar_len() - 1)])
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
                    record,
                    *qpos,
                    // None,
                    true,
                    // min_soft_clip_qual,
                );

                if read_char.to_ascii_uppercase() != refr_base.to_ascii_uppercase() || next_to_sc_indel {
                    return true;
                }

                return false;
            },
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
        let mut past_query_pos = false;
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

        return next_to_soft_clip
    }

    fn check_position_against_cigar(
        // hq_soft_clips: &mut Option<&mut RunningAverage>,
        record: &Record,
        cig: &Cigar,
        read_cursor: i32,
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
            },
            Cigar::Ins(_) | Cigar::Del(_) => {
                if check_indels {
                    next_to_soft_clip = true
                }
            },
            _ => {
                // Not a soft clip
            }
        }
        return next_to_soft_clip
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
                },
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
        non_ref_likelihood: f64
    ) {
        // https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReferenceConfidenceModel.java
        // applyPileupElementRefVsNonRefLikelihoodAndCount
        // Homozygous likelihoods don't need the logSum trick.

        result.genotype_likelihoods[0] += ref_likelihood + log10ploidy;
        result.genotype_likelihoods[(likelihoodcount - 1)] += non_ref_likelihood + log10ploidy;
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

        let mut beta_entropy;
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

    /// Takes a vector of VariantContexts and writes them to a single VCF4 file
    pub fn write_vcf(
        &self,
        output_prefix: &str,
        variant_contexts: &Vec<VariantContext>,
        sample_names: &[&str],
        reference_reader: &ReferenceReader,
        strain_info: bool,
    ) {
        if variant_contexts.len() > 0 {
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
    }

    fn populate_vcf_header(
        &self,
        sample_names: &[&str],
        reference_reader: &ReferenceReader,
        header: &mut Header,
        strain_info: bool,
    ) {
        header.push_record(format!("##source=lorikeet-v{}", env!("CARGO_PKG_VERSION")).as_bytes());

        debug!("samples {:?}", &sample_names);
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
            is_refskip
        }
    }
}
