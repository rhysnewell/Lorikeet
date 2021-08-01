use activity_profile::activity_profile::Profile;
use activity_profile::activity_profile_state::{ActivityProfileState, Type};
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use assembly::assembly_region::AssemblyRegion;
use assembly::assembly_region_trimmer::AssemblyRegionTrimmer;
use assembly::assembly_result_set::AssemblyResultSet;
use bio::stats::{LogProb, PHREDProb};
use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::FlagFilter;
use estimation::lorikeet_engine::Elem;
use estimation::lorikeet_engine::ReadType;
use genotype::genotype_builder::Genotype;
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use genotype::genotyping_engine::GenotypingEngine;
use graphs::chain_pruner::ChainPruner;
use graphs::multi_sample_edge::MultiSampleEdge;
use haplotype::ref_vs_any_result::RefVsAnyResult;
use haplotype::reference_confidence_model::ReferenceConfidenceModel;
use mathru::special::gamma::{digamma, ln_gamma};
use model::variant_context::VariantContext;
use model::variants::*;
use pair_hmm::pair_hmm_likelihood_calculation_engine::PairHMMLikelihoodCalculationEngine;
use rayon::prelude::*;
use read_orientation::beta_distribution_shape::BetaDistributionShape;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use read_threading::read_threading_assembler::ReadThreadingAssembler;
use read_threading::read_threading_graph::ReadThreadingGraph;
use reads::alignment_utils::AlignmentUtils;
use reads::bird_tool_reads::BirdToolRead;
use reference::reference_reader::ReferenceReader;
use rust_htslib::bam::{self, record::Cigar, Record};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use utils::math_utils::{MathUtils, RunningAverage};
use utils::natural_log_utils::NaturalLogUtils;
use utils::quality_utils::QualityUtils;
use utils::simple_interval::SimpleInterval;

pub struct HaplotypeCallerEngine {
    genotyping_engine: GenotypingEngine,
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
    const AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD: f64 = 6.0;

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

    // pub const MIN_SOFT_CLIP_QUAL: usize = 29;

    // const NO_CALLS: Vec<Allele> = Vec::new();

    pub fn new(
        args: &clap::ArgMatches,
        ref_idx: usize,
        samples: Vec<String>,
        do_allele_specific_calcs: bool,
        sample_ploidy: usize,
    ) -> HaplotypeCallerEngine {
        let mut assembly_engine = ReadThreadingAssembler::new(
            args.value_of("max-allowed-path-for-read-threading-assembler")
                .unwrap()
                .parse()
                .unwrap(),
            args.values_of("kmer-sizes")
                .unwrap()
                .map(|k_size| k_size.parse().unwrap())
                .collect::<Vec<usize>>(),
            args.is_present("dont-increase-kmer-sizes-for-cycles"),
            args.is_present("allow-non-unique-kners-in-ref"),
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
            args.value_of("pruning-log-odds-threshold")
                .unwrap()
                .parse()
                .unwrap(),
            args.value_of("pruning-seeding-log-odds-threshold")
                .unwrap()
                .parse()
                .unwrap(),
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
            genotyping_engine: GenotypingEngine::make(
                args,
                samples,
                do_allele_specific_calcs,
                sample_ploidy,
            ),
            genotype_prior_calculator: GenotypePriorCalculator::make(args),
            stand_min_conf: args
                .value_of("standard-min-confidence-threshold-for-calling")
                .unwrap()
                .parse()
                .unwrap(),
            ref_idx,
            assembly_engine: assembly_engine,
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
                .value_of("mapping-quality-threshold")
                .unwrap()
                .parse()
                .unwrap(),
        }
    }

    pub fn stand_min_conf(&self) -> f64 {
        self.stand_min_conf
    }

    pub fn collect_activity_profile(
        &mut self,
        indexed_bam_readers: &Vec<String>,
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
    ) -> HashMap<usize, BandPassActivityProfile> {
        // minimum PHRED base quality
        let bq = m
            .value_of("base-quality-threshold")
            .unwrap()
            .parse::<i32>()
            .unwrap();

        let max_prob_prop = m
            .value_of("max-prob-propagation-distance")
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let active_prob_thresh = m
            .value_of("active-probability-threshold")
            .unwrap()
            .parse::<f64>()
            .unwrap();

        // let min_variant_quality = m
        //     .value_of("min-variant-quality")
        //     .unwrap()
        //     .parse::<i32>()
        //     .unwrap();
        let min_soft_clip_qual = 29;

        let ploidy: usize = m.value_of("ploidy").unwrap().parse().unwrap();

        let mut genotype_likelihoods = Vec::with_capacity(short_sample_count + long_sample_count);

        let mut per_contig_per_base_hq_soft_clips = HashMap::new();

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
                if sample_idx < short_sample_count {
                    genotype_likelihoods.push(HaplotypeCallerEngine::update_activity_profile(
                        bam_generator,
                        n_threads,
                        ref_idx,
                        ReadType::Short,
                        ploidy,
                        bq,
                        min_soft_clip_qual,
                        genomes_and_contigs,
                        concatenated_genomes,
                        flag_filters,
                        &mut per_contig_per_base_hq_soft_clips,
                        reference_reader,
                    ));
                } else if (m.is_present("longreads") | m.is_present("longread-bam-files"))
                    && sample_idx >= short_sample_count
                    && sample_idx < (short_sample_count + long_sample_count)
                {
                    debug!("Running structural variant detection...");
                    // Get the appropriate sample index based on how many references we are using by tracking
                    // changes in references
                    genotype_likelihoods.push(HaplotypeCallerEngine::update_activity_profile(
                        bam_generator,
                        n_threads,
                        ref_idx,
                        ReadType::Long,
                        ploidy,
                        bq,
                        min_soft_clip_qual,
                        genomes_and_contigs,
                        concatenated_genomes,
                        flag_filters,
                        &mut per_contig_per_base_hq_soft_clips,
                        reference_reader,
                    ));
                }
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];

                    pb.progress_bar.set_message(format!(
                        "{}: Variant calling on sample: {}",
                        pb.key, &indexed_bam_readers[sample_idx],
                    ));
                    pb.progress_bar.inc(1);
                }
                {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(1);
                }
            });

        // return genotype_likelihoods for each contig in current genome across samples
        return self.calculate_activity_probabilities(
            genotype_likelihoods,
            per_contig_per_base_hq_soft_clips,
            &reference_reader.target_lens,
            ploidy,
            max_prob_prop,
            active_prob_thresh,
            ref_idx,
        );
    }

    pub fn update_activity_profile<
        'b,
        R: IndexedNamedBamReader + Send,
        G: NamedBamReaderGenerator<R> + Send,
    >(
        bam_generator: G,
        split_threads: usize,
        ref_idx: usize,
        readtype: ReadType,
        ploidy: usize,
        bq: i32,
        min_soft_clip_qual: i32,
        genomes_and_contigs: &'b GenomesAndContigs,
        concatenated_genomes: &'b Option<String>,
        flag_filters: &'b FlagFilter,
        per_contig_per_base_hq_soft_clips: &mut HashMap<usize, Vec<RunningAverage>>,
        reference_reader: &mut ReferenceReader,
    ) -> HashMap<usize, Vec<RefVsAnyResult>> {
        let mut bam_generated = bam_generator.start();

        bam_generated.set_threads(split_threads);

        let header = bam_generated.header().clone(); // bam header
        let target_lens: Vec<u64> = (0..header.target_count())
            .into_iter()
            .map(|tid| header.target_len(tid).unwrap())
            .collect();
        let target_names = header.target_names();
        let reference = reference_reader.retrieve_reference_stem(ref_idx);

        let likelihoodcount = ploidy + 1;
        let log10ploidy = (likelihoodcount as f64).log10();

        let mut ref_vs_any_container = HashMap::new();
        // for each genomic position, only has hashmap when variants are present. Includes read ids
        match readtype {
            ReadType::Short | ReadType::Long => {
                target_names
                    .iter()
                    .enumerate()
                    .for_each(|(tid, contig_name)| {
                        let target_name = String::from_utf8(contig_name.to_vec()).unwrap();
                        if target_name.contains(&reference)
                            || reference_reader.match_target_name_and_ref_idx(ref_idx, &target_name)
                        {
                            // Get contig stats
                            reference_reader.add_target(contig_name, tid);
                            let target_len = target_lens[tid];
                            reference_reader.add_length(tid, target_len);
                            let per_base_hq_soft_clips = per_contig_per_base_hq_soft_clips.entry(tid).or_insert(vec![RunningAverage::new(); target_len as usize]);
                            // The raw activity profile.
                            // Frequency of bases not matching reference compared
                            // to depth
                            let mut genotype_likelihoods = Vec::with_capacity(target_len as usize);

                            {
                                bam_generated.fetch(tid as u32);

                                // Position based - Loop through the positions in the genome
                                // Calculate likelihood of activity being present at this position
                                match bam_generated.pileup() {
                                    Some(pileups) => {
                                        reference_reader.update_current_sequence_capacity(target_len as usize);
                                        // Update all contig information
                                        reference_reader.fetch_contig_from_reference_by_contig_name(
                                            &contig_name.to_vec(),
                                            ref_idx as usize,
                                        );

                                        reference_reader.read_sequence_to_vec();

                                        for p in pileups {
                                            let pileup = p.unwrap();
                                            let pos = pileup.pos() as usize;
                                            let hq_soft_clips = &mut per_base_hq_soft_clips[pos];
                                            let refr_base = reference_reader.current_sequence[pos];
                                            let mut result = RefVsAnyResult::new(
                                                likelihoodcount,
                                                pos,
                                                tid
                                            );
                                            for alignment in pileup.alignments() {
                                                let record = alignment.record();
                                                if (!flag_filters.include_supplementary
                                                    && record.is_supplementary()
                                                    && readtype != ReadType::Long)
                                                    || (!flag_filters.include_secondary
                                                    && record.is_secondary())
                                                    || (!flag_filters.include_improper_pairs
                                                    && !record.is_proper_pair()
                                                    && readtype != ReadType::Long)
                                                {
                                                    continue;
                                                } else if !flag_filters.include_secondary
                                                    && record.is_secondary()
                                                    && readtype == ReadType::Long
                                                {
                                                    continue;
                                                } else {
                                                    HaplotypeCallerEngine::alignment_context_creation(
                                                        &alignment,
                                                        &mut result,
                                                        hq_soft_clips,
                                                        log10ploidy,
                                                        likelihoodcount,
                                                        min_soft_clip_qual,
                                                        refr_base,
                                                        bq,
                                                    )
                                                }
                                            }

                                            let denominator = result.read_counts as f64 * log10ploidy;
                                            for i in 0..likelihoodcount {
                                                result.genotype_likelihoods[i] -= denominator
                                            }

                                            genotype_likelihoods.push(result)
                                        }
                                    }
                                    None => println!("no bam for pileups"),
                                };
                            };

                            // Add in genotype_likelihoods to return result
                            ref_vs_any_container.insert(tid, genotype_likelihoods);

                        }
                    });
            }
            _ => {}
        }

        bam_generated.finish();

        return ref_vs_any_container;
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
        mut per_contig_per_base_hq_soft_clips: HashMap<usize, Vec<RunningAverage>>,
        target_ids_and_lens: &HashMap<usize, u64>,
        ploidy: usize,
        max_prob_propagation: usize,
        active_prob_threshold: f64,
        ref_idx: usize,
    ) -> HashMap<usize, BandPassActivityProfile> {
        if genotype_likelihoods.len() == 1 {
            // Faster implementation for single sample analysis
            let per_contig_activity_profiles: HashMap<usize, BandPassActivityProfile> =
                genotype_likelihoods[0]
                    .par_iter()
                    .map(|(tid, vec_of_ref_vs_any_result)| {
                        // Create bandpass
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

                        // for each position determine per locus activity and add to bandpass
                        vec_of_ref_vs_any_result.iter().enumerate().for_each(
                            |(pos, ref_vs_any_result)| {
                                let is_active_prob = self
                                    .genotyping_engine
                                    .allele_frequency_calculator
                                    .calculate_single_sample_biallelic_non_ref_posterior(
                                        &ref_vs_any_result.genotype_likelihoods,
                                        true,
                                    );

                                let per_base_hq_soft_clips =
                                    per_contig_per_base_hq_soft_clips.get(tid).unwrap();

                                let hq_soft_clips = per_base_hq_soft_clips[pos];
                                let activity_profile_state = ActivityProfileState::new(
                        ref_vs_any_result.loc.clone(),
                        is_active_prob,
                        Type::new(
                            hq_soft_clips.mean(),
                            HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                        )
                    );

                                activity_profile.add(activity_profile_state);
                            },
                        );

                        (*tid, activity_profile)
                    })
                    .collect::<HashMap<usize, BandPassActivityProfile>>();

            return per_contig_activity_profiles;
        } else {
            let mut per_contig_activity_profiles = HashMap::new();
            for (tid, length) in target_ids_and_lens.iter() {
                let per_base_hq_soft_clips =
                    per_contig_per_base_hq_soft_clips.entry(*tid).or_insert(
                        vec![RunningAverage::new(); target_ids_and_lens[tid] as usize],
                    );

                // Create bandpass
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

                for pos in 0..(*length as usize) {
                    let mut genotypes = Vec::new();

                    let hq_soft_clips = per_base_hq_soft_clips[pos];

                    for sample_likelihoods in genotype_likelihoods.iter() {
                        let result = sample_likelihoods[tid][pos].genotype_likelihoods.clone();
                        genotypes.push(Genotype::build(ploidy, result))
                    }

                    let fake_alleles = Allele::create_fake_alleles();

                    let mut variant_context = VariantContext::build(*tid, pos, pos, fake_alleles);

                    variant_context.add_genotypes(genotypes);

                    let vc_out = self.genotyping_engine.calculate_genotypes(
                        variant_context,
                        ploidy,
                        &self.genotype_prior_calculator,
                        Vec::new(),
                        self.stand_min_conf,
                    );

                    let is_active_prob = match vc_out {
                        Some(vc) => QualityUtils::qual_to_prob(vc.get_phred_scaled_qual() as u8),
                        None => 0.0,
                    };

                    let activity_profile_state = ActivityProfileState::new(
                        SimpleInterval::new(*tid, pos, pos),
                        is_active_prob,
                        Type::new(
                            hq_soft_clips.mean(),
                            HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                        ),
                    );
                    activity_profile.add(activity_profile_state);
                }
                per_contig_activity_profiles.insert(*tid, activity_profile);
            }
            return per_contig_activity_profiles;
        }
    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    pub fn call_region(
        &mut self,
        mut region: AssemblyRegion,
        reference_reader: &mut ReferenceReader,
        given_alleles: &Vec<VariantContext>,
        args: &clap::ArgMatches,
        sample_names: &Vec<String>,
    ) -> Vec<VariantContext> {
        let vc_priors = Vec::new();
        if !region.is_active() {
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }

        if given_alleles.is_empty() && region.len() == 0 {
            return self.reference_model_for_no_variation(&mut region, true, &vc_priors);
        }

        let mut region_without_reads = region.clone_without_reads();

        // run the local assembler, getting back a collection of information on how we should proceed
        let mut untrimmed_assembly_result = AssemblyBasedCallerUtils::assemble_reads(
            region,
            given_alleles,
            args,
            reference_reader,
            &mut self.assembly_engine,
            args.is_present("correct-overlapping-quality"),
            sample_names,
        );

        debug!(
            "There were {} haplotypes found",
            untrimmed_assembly_result.haplotypes.len()
        );
        let all_variation_events = untrimmed_assembly_result
            .get_variation_events(args.value_of("max-mnp=distance").unwrap().parse().unwrap());

        let mut trimming_result = self.assembly_region_trimmer.trim(
            self.ref_idx,
            region_without_reads,
            all_variation_events,
            args.is_present("enable-legacy-assembly-region-trimming"),
            &reference_reader,
            &untrimmed_assembly_result.full_reference_with_padding,
        );

        if !trimming_result.is_variation_present() && args.is_present("disable-optimizations") {
            return self.reference_model_for_no_variation(
                &mut trimming_result.original_region,
                false,
                &vc_priors,
            );
        }

        let mut assembly_result =
            untrimmed_assembly_result.trim_to(trimming_result.get_variant_region());

        let read_stubs = assembly_result
            .region_for_genotyping
            .reads
            .par_iter()
            .filter(|r| {
                AlignmentUtils::unclipped_read_length(r)
                    < AssemblyBasedCallerUtils::MINIMUM_READ_LENGTH_AFTER_TRIMMING
            })
            .cloned()
            .collect::<Vec<BirdToolRead>>();
        let mut assembly_result = assembly_result.remove_all(&read_stubs);

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        let (mut assembly_result, filtered_reads) = self.filter_non_passing_reads(assembly_result);
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
            "              Assembly region {:?}",
            &assembly_result.region_for_genotyping
        );
        debug!("==========================================================================");
        let reads = AssemblyBasedCallerUtils::split_reads_by_sample(
            assembly_result.region_for_genotyping.get_reads_cloned(),
        );

        let read_likelihoods = self.likelihood_calculation_engine.compute_read_likelihoods(
            assembly_result,
            sample_names.clone(),
            reads,
        );

        // Realign reads to their best haplotype.

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        return vc_priors;
    }

    fn filter_non_passing_reads(
        &self,
        mut assembly_result: AssemblyResultSet<ReadThreadingGraph>,
    ) -> (AssemblyResultSet<ReadThreadingGraph>, Vec<BirdToolRead>) {
        let reads_to_remove = assembly_result
            .region_for_genotyping
            .reads
            .par_iter()
            .filter(|r| {
                AlignmentUtils::unclipped_read_length(r) < Self::READ_LENGTH_FILTER_THRESHOLD
                    || r.read.mapq() < self.mapping_quality_threshold
                    || !r.read.is_proper_pair()
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
    fn reference_model_for_no_variation(
        &self,
        region: &mut AssemblyRegion,
        needs_to_be_finalized: bool,
        vc_priors: &Vec<VariantContext>,
    ) -> Vec<VariantContext> {
        Vec::new() // TODO: Implement this potentially?
    }

    /**
     * Populates reference and non reference depth vectors
     */
    fn alignment_context_creation(
        alignment: &bam::pileup::Alignment,
        result: &mut RefVsAnyResult,
        hq_soft_clips: &mut RunningAverage,
        log10ploidy: f64,
        likelihoodcount: usize,
        min_soft_clip_qual: i32,
        refr_base: u8,
        bq: i32,
    ) {
        let ref_likelihood;
        let non_ref_likelihood;
        let record = alignment.record();

        if !alignment.is_del() && !alignment.is_refskip() {
            // query position in read
            let qpos = alignment.qpos().unwrap();
            let record_qual = record.qual()[qpos] as i32;
            let mut is_alt = false;
            if record_qual >= bq {
                let read_char = alignment.record().seq()[qpos];
                if refr_base != read_char {
                    is_alt = true;
                    result.non_ref_depth += 1;
                    non_ref_likelihood = f64::from(LogProb::from(PHREDProb(record_qual as f64)));
                    ref_likelihood = f64::from(LogProb::from(PHREDProb(record_qual as f64)))
                        + (-(3.0_f64.log10()))
                } else {
                    result.ref_depth += 1;
                    ref_likelihood = f64::from(LogProb::from(PHREDProb(record_qual as f64)));
                    non_ref_likelihood = f64::from(LogProb::from(PHREDProb(record_qual as f64)))
                        + (-(3.0_f64.log10()));
                }

                HaplotypeCallerEngine::update_heterozygous_likelihood(
                    result,
                    likelihoodcount,
                    log10ploidy,
                    ref_likelihood,
                    non_ref_likelihood,
                );
            }

            if is_alt {
                let mut read_cursor = 0;
                for cig in record.cigar().iter() {
                    // Cigar immediately before current position
                    // in read
                    if qpos == read_cursor - 1 {
                        match cig {
                            Cigar::SoftClip(_) => hq_soft_clips.add(
                                HaplotypeCallerEngine::count_high_quality_soft_clips(
                                    &cig,
                                    &record,
                                    read_cursor,
                                    min_soft_clip_qual,
                                ),
                            ),
                            _ => {
                                // Not a soft clip
                            }
                        }
                    } else if qpos == read_cursor + 1 {
                        match cig {
                            Cigar::SoftClip(_) => hq_soft_clips.add(
                                HaplotypeCallerEngine::count_high_quality_soft_clips(
                                    &cig,
                                    &record,
                                    read_cursor,
                                    min_soft_clip_qual,
                                ),
                            ),
                            _ => {
                                // Not a soft clip
                            }
                        }
                    } else if read_cursor > qpos {
                        // break out of loop since we have passed
                        // the position
                        break;
                    }
                    match cig {
                        // Progress the cigar cursor
                        Cigar::Match(_)
                        | Cigar::Diff(_)
                        | Cigar::Equal(_)
                        | Cigar::Ins(_)
                        | Cigar::SoftClip(_) => {
                            read_cursor += cig.len() as usize;
                        }
                        _ => {}
                    }
                }
            }
        }

        result.read_counts += 1;
    }

    fn count_high_quality_soft_clips(
        cig: &Cigar,
        record: &Record,
        cigar_cursor: usize,
        min_soft_clip_qual: i32,
    ) -> f64 {
        // https://gatk.broadinstitute.org/hc/en-us/articles/360036227652?id=4147
        // If we have high quality soft clips then we want to
        // track them
        // get mean base quality
        let mut num_high_quality_soft_clips = 0.0;
        for rpos in cigar_cursor..(cigar_cursor + cig.len() as usize) {
            let qual_pos = record.qual()[rpos] as i32;
            if qual_pos >= min_soft_clip_qual {
                num_high_quality_soft_clips += 1.0
            }
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
}
