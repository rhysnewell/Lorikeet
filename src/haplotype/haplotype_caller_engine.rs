// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use regex::Regex;
use serde_json;
use serde_json::Value;
use ordered_float::OrderedFloat;
use rust_htslib::bam::{self, record::Cigar, Read, Record};

use bio::stats::{LogProb, PHREDProb};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;


use coverm::bam_generator::*;
use coverm::FlagFilter;
use model::variants::*;
use model::variant_context::VariantContext;

use crate::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use statrs::statistics::Variance;
use utils::utils::{ReadType, retrieve_reference, Elem, fetch_contig_from_reference, read_sequence_to_vec};
use haplotype::ref_vs_any_result::RefVsAnyResult;
use genotype::genotype_builder::Genotype;
use activity_profile::activity_profile_state::{ActivityProfileState, Type};
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use genotype::genotyping_engine::GenotypingEngine;
use utils::quality_utils::QualityUtils;
use utils::math_utils::RunningAverage;
use utils::simple_interval::SimpleInterval;

pub struct HaplotypeCallerEngine {
    genotyping_engine: GenotypingEngine,
    genotype_prior_calculator: GenotypePriorCalculator,
    ref_idx: usize,
    stand_min_conf: f64,
}

impl HaplotypeCallerEngine {
    pub const MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION: usize = 6;
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
        sample_ploidy: usize
    ) -> HaplotypeCallerEngine {
        HaplotypeCallerEngine {
            genotyping_engine: GenotypingEngine::make(args, samples, do_allele_specific_calcs, sample_ploidy),
            genotype_prior_calculator: GenotypePriorCalculator::make(args),
            stand_min_conf: args.value_of("standard-min-confidence-threshold-for-calling").unwrap().parse().unwrap(),
            ref_idx,
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
        per_reference_samples: usize,
        m: &clap::ArgMatches,
        genomes_and_contigs: &GenomesAndContigs,
        reference_map: &HashMap<usize, String>,
        concatenated_genomes: &Option<String>,
        flag_filters: &FlagFilter,
        tree: &Mutex<Arc<Vec<&Elem>>>,
    ) -> HashMap<usize, Vec<ActivityProfileState>> {


        // minimum PHRED base quality
        let bq = m
            .value_of("base-quality-threshold")
            .unwrap()
            .parse::<i32>()
            .unwrap();


        // let min_variant_quality = m
        //     .value_of("min-variant-quality")
        //     .unwrap()
        //     .parse::<i32>()
        //     .unwrap();
        let min_soft_clip_qual = 29;

        let ploidy: usize = m.value_of("ploidy").unwrap().parse().unwrap();

        let mut genotype_likelihoods = Vec::with_capacity(
            short_sample_count + long_sample_count
        );

        let mut target_ids_and_lengths = HashMap::new();

        let mut per_contig_per_base_hq_soft_clips = HashMap::new();

        indexed_bam_readers.iter().enumerate().for_each(
            |(sample_idx, bam_generator)| {
                // Get the appropriate sample index based on how many references we are using
                let mut bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                    vec![&bam_generator],
                    n_threads as u32,
                )
                    .into_iter()
                    .next()
                    .unwrap();
                if sample_idx < short_sample_count {
                    genotype_likelihoods.push(
                        HaplotypeCallerEngine::update_activity_profile(
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
                            &mut target_ids_and_lengths,
                            &mut per_contig_per_base_hq_soft_clips,
                        )
                    );
                } else if (m.is_present("longreads") | m.is_present("longread-bam-files"))
                    && sample_idx >= short_sample_count
                    && sample_idx < (short_sample_count + long_sample_count)
                {
                    debug!("Running structural variant detection...");
                    // Get the appropriate sample index based on how many references we are using by tracking
                    // changes in references
                    genotype_likelihoods.push(
                        HaplotypeCallerEngine::update_activity_profile(
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
                            &mut target_ids_and_lengths,
                            &mut per_contig_per_base_hq_soft_clips,
                        )
                    );
                }
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];

                    pb.progress_bar.set_message(&format!(
                        "{}: Variant calling on sample: {}",
                        pb.key,
                        &indexed_bam_readers[sample_idx],
                    ));
                    pb.progress_bar.inc(1);
                }
                {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(1);
                }
            },
        );

        // return genotype_likelihoods for each contig in current genome across samples
        self.calculate_activity_probabilities(
            genotype_likelihoods,
            per_contig_per_base_hq_soft_clips,
            target_ids_and_lengths,
            ploidy,
        )

    }

    pub fn update_activity_profile<'b,
        R: IndexedNamedBamReader + Send,
        G: NamedBamReaderGenerator<R> + Send>(
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
        target_ids_and_lengths: &mut HashMap<usize, u64>,
        per_contig_per_base_hq_soft_clips: &mut HashMap<usize, Vec<RunningAverage<f64>>>,
    ) -> HashMap<usize, Vec<RefVsAnyResult>>{
        let mut bam_generated = bam_generator.start();

        let reference = &genomes_and_contigs.genomes[ref_idx];
        let mut reference_file = retrieve_reference(concatenated_genomes);

        bam_generated.set_threads(split_threads);

        let header = bam_generated.header().clone(); // bam header
        let target_lens: Vec<u64> = (0..header.target_count())
            .into_iter()
            .map(|tid| header.target_len(tid).unwrap())
            .collect();
        let target_names = header.target_names();

        let likelihoodcount = ploidy + 1;
        let log10ploidy = OrderedFloat((likelihoodcount as f64).log10());

        let mut contig_stats: HashMap<usize, Vec<f64>> = HashMap::new();

        let mut ref_vs_any_container = HashMap::new();
        // for each genomic position, only has hashmap when variants are present. Includes read ids
        match readtype {
            ReadType::Short | ReadType::Long => {
                target_names
                    .iter()
                    .enumerate()
                    .for_each(|(tid, contig_name)| {
                        let target_name = String::from_utf8(contig_name.to_vec()).unwrap();
                        if target_name.contains(reference)
                            || match genomes_and_contigs.contig_to_genome.get(&target_name) {
                            Some(ref_id) => *ref_id == ref_idx,
                            None => false,
                        }
                        {
                            // Get contig stats
                            let target_len = target_lens[tid];
                            let mut per_base_hq_soft_clips = per_contig_per_base_hq_soft_clips.entry(tid).or_insert(vec![RunningAverage::new(); target_len as usize]);

                            if !target_ids_and_lengths.contains_key(&tid) {
                                target_ids_and_lengths.insert(tid, target_len);
                            }
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
                                        let mut ref_seq = Vec::with_capacity(target_len as usize);
                                        // Update all contig information
                                        fetch_contig_from_reference(
                                            &mut reference_file,
                                            &contig_name.to_vec(),
                                            genomes_and_contigs,
                                            ref_idx as usize,
                                        );

                                        read_sequence_to_vec(
                                            &mut ref_seq,
                                            &mut reference_file,
                                            &contig_name.to_vec(),
                                        );

                                        for p in pileups {
                                            let pileup = p.unwrap();
                                            let pos = pileup.pos() as usize;
                                            let hq_soft_clips = &mut per_base_hq_soft_clips[pos];
                                            let refr_base = ref_seq[pos];
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

                                            let denominator = OrderedFloat(result.read_counts as f64) * log10ploidy;
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

                            // // Get coverage mean and standard dev
                            // let contig_cov = depth_sum as f64 / target_len as f64;
                            // let std_dev = 10. * coverage.std_dev();
                            // contig_stats.entry(tid).or_insert(vec![contig_cov, std_dev]);
                            //
                            // for (index, count) in soft_clips.iter().enumerate() {
                            //     if count >= 5 {
                            //         activity[index] += count as f32; // Add in good soft clips
                            //     }
                            //     activity[index] /= depths[index] as f32; // divide by depth
                            // }
                        }
                    });
            }
            _ => {}
        }

        bam_generated.finish();

        return ref_vs_any_container
    }

    /**
    * Create a collection of ActivityProfileState for each position on each contig
    * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
    * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
    *
    * @return HashMap<usize, Vec<ActivityProfileState>> - contig id and per base activity
    */
    pub fn calculate_activity_probabilities(
        &mut self,
        genotype_likelihoods: Vec<HashMap<usize, Vec<RefVsAnyResult>>>,
        per_contig_per_base_hq_soft_clips: HashMap<usize, Vec<RunningAverage<f64>>>,
        target_ids_and_lens: HashMap<usize, u64>,
        ploidy: usize,
    ) -> HashMap<usize, Vec<ActivityProfileState>> {
        if genotype_likelihoods.len() == 1 {
            // Faster implementation for single sample analysis
            let per_contig_activity_profile_states: HashMap<usize, Vec<ActivityProfileState>> = genotype_likelihoods[0]
                .par_iter().map(|&(tid, vec_of_ref_vs_any_result)| {
                let result = vec_of_ref_vs_any_result.iter().enumerate().map(|(pos, ref_vs_any_result)| {
                    let is_active_prob =
                        self.genotyping_engine.allele_frequency_calculator
                            .calculate_single_sample_biallelic_non_ref_posterior(
                                &ref_vs_any_result.genotype_likelihoods,
                                true
                            );

                    let mut per_base_hq_soft_clips = per_contig_per_base_hq_soft_clips.entry(tid)
                        .or_insert(vec![RunningAverage::new(); target_ids_and_lens[tid]]);

                    let hq_soft_clips = &mut per_base_hq_soft_clips[pos];
                    ActivityProfileState::new(
                        ref_vs_any_result.loc,
                        is_active_prob,
                        Type::new(
                            hq_soft_clips.mean(),
                            HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                        )
                    )
                }).collect::<Vec<ActivityProfileState>>();

                (tid, result)
            }).collect::<HashMap<usize, Vec<ActivityProfileState>>>();

            return per_contig_activity_profile_states
        } else {

            let mut per_contig_activity_profile_states = HashMap::new();
            for (tid, length) in target_ids_and_lens.iter() {

                let mut activity_profile_states = Vec::with_capacity(length);
                for pos in 0..length.iter() {
                    let mut activity_probability = 0.;
                    let mut genotypes = Vec::new();

                    let mut per_base_hq_soft_clips = per_contig_per_base_hq_soft_clips.entry(tid)
                        .or_insert(vec![RunningAverage::new(); target_ids_and_lens[tid]]);
                    let hq_soft_clips = &mut per_base_hq_soft_clips[pos];

                    for sample_likelihoods in genotype_likelihoods.iter() {
                        let result = sample_likelihoods[tid][pos].genotype_likelihoods.clone();
                        genotypes.push(Genotype::build(ploidy, result))
                    }

                    let fake_alleles = Allele::create_fake_alleles();

                    let mut variant_context = VariantContext::build(
                        *tid, pos, pos, fake_alleles
                    );

                    variant_context.add_genotypes(
                        genotypes
                    );


                    let vc_out = self.genotyping_engine.calculate_genotypes(
                        variant_context,
                        ploidy,
                        self.genotype_prior_calculator,
                        Vec::new(),
                        self.stand_min_conf
                    );

                    let is_active_prob = match vc_out {
                        Some(vc) => {
                            QualityUtils::qual_to_prob(vc.get_phred_scaled_qual())
                        },
                        None => 0.0
                    };

                    let activity_profile_state = ActivityProfileState::new(
                        SimpleInterval::new(tid, pos, pos),
                        is_active_prob,
                        Type::new(
                            hq_soft_clips.mean(),
                            HaplotypeCallerEngine::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
                        )
                    );
                    activity_profile_states.push(activity_profile_state);
                }
                per_contig_activity_profile_states.insert(tid, activity_profile_states)
            }
            return per_contig_activity_profile_states
        }
    }


    /**
    * Populates reference and non reference depth vectors
    */
    fn alignment_context_creation(
        alignment: &bam::pileup::Alignment,
        result: &mut RefVsAnyResult,
        hq_soft_clips: &mut RunningAverage<f64>,
        log10ploidy: OrderedFloat<f64>,
        likelihoodcount: usize,
        min_soft_clip_qual: i32,
        refr_base: u8,
        bq: i32,
    ) {
        let mut ref_likelihood = OrderedFloat(0.0);
        let mut non_ref_likelihood = OrderedFloat(0.0);
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
                    non_ref_likelihood = OrderedFloat(f64::from(LogProb::from(PHREDProb(record_qual as f64))));
                    ref_likelihood = OrderedFloat(f64::from(LogProb::from(PHREDProb(record_qual as f64))) + (-(3.0.log10())))
                } else {
                    result.ref_depth += 1;
                    ref_likelihood = OrderedFloat(f64::from(LogProb::from(PHREDProb(record_qual as f64))));
                    non_ref_likelihood = OrderedFloat(f64::from(LogProb::from(PHREDProb(record_qual as f64))) + (-(3.0.log10())));
                }

                HaplotypeCallerEngine::update_heterozygous_likelihood(
                    result,
                    likelihoodcount,
                    log10ploidy,
                    ref_likelihood,
                    non_ref_likelihood
                );
            }

            if is_alt {
                let mut read_cursor = 0;
                for cig in record.cigar().iter() {
                    // Cigar immediately before current position
                    // in read
                    if qpos == read_cursor - 1 {
                        match cig {
                            Cigar::SoftClip(_) => {
                                hq_soft_clips.add(
                                    HaplotypeCallerEngine::count_high_quality_soft_clips(
                                        &cig,
                                        &record,
                                        &mut read_cursor,
                                        min_soft_clip_qual
                                    )
                                )
                            },
                            _ => {
                                // Not a soft clip
                            }
                        }
                    } else if qpos == read_cursor + 1 {
                        match cig {
                            Cigar::SoftClip(_) => {
                                hq_soft_clips.add(
                                    HaplotypeCallerEngine::count_high_quality_soft_clips(
                                        &cig,
                                        &record,
                                        &mut read_cursor,
                                        min_soft_clip_qual
                                    )
                                )

                            },
                            _ => {
                                // Not a soft clip
                            }
                        }
                    } else if read_cursor > qpos {
                        // break out of loop since we have passed
                        // the position
                        break
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
                        _ => {},
                    }

                }
            }

        }

        result.read_counts += 1;
    }


    fn count_high_quality_soft_clips(
        cig: &Cigar,
        record: &Record,
        cigar_cursor: &mut usize,
        min_soft_clip_qual: i32
    ) -> usize {
        // https://gatk.broadinstitute.org/hc/en-us/articles/360036227652?id=4147
        // If we have high quality soft clips then we want to
        // track them
        // get mean base quality
        let mut num_high_quality_soft_clips = 0;
        for rpos in cigar_cursor..(cigar_cursor + cig.len() as usize) {
            let qual_pos = record.qual()[rpos];
            if qual_pos >= min_soft_clip_qual {
                num_high_quality_soft_clips += 1
            }
        }
        return num_high_quality_soft_clips
    }

    fn update_heterozygous_likelihood(
        result: &mut RefVsAnyResult,
        likelihoodcount: usize,
        log10ploidy: OrderedFloat<f64>,
        ref_likelihood: OrderedFloat<f64>,
        non_ref_likelihood: OrderedFloat<f64>
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
            result.genotype_likelihoods[i] += (
                ref_likelihood + OrderedFloat((j as f64).log10()),
                non_ref_likelihood + OrderedFloat((i as f64).log10())
            );
            i += 1;
            j -= 1;
        }
    }
}



fn compute_moving_average(window: &[(u32, u32)]) -> (u32, f32) {
    let window_size = window.len();

    let current_year = window[window_size / 2].0;

    let sum: u32 = window.iter().map(|&(_, val)| val).sum();
    let sum = sum as f32 / window_size as f32;

    (current_year, sum)
}

fn extract_moving_average_for_year(year: u32, moving_average: &[(u32, f32)]) -> Option<f32> {
    moving_average
        .iter()
        .find(|(yr, _)| yr == year)
        .map(|&(_, val)| val)
}