// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::str::FromStr;
use regex::Regex;
use std::fmt;
use serde_json;
use serde_json::Value;
use std::cmp::{max, min};
use std::fs::File;
use ordered_float::OrderedFloat;
use rust_htslib::bam::{self, record::Cigar, Read, IndexedReader, Record};

use bird_tool_utils::command;
use rust_htslib::errors::Error;
use bio::stats::{LogProb, PHREDProb};
use std;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};


use coverm::bam_generator::*;
use coverm::FlagFilter;
use external_command_checker;
use model::variants::*;
use model::variant_context::VariantContext;
use model::genotype_builder::Genotype;
use utils::*;

use crate::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use scoped_threadpool::Pool;
use statrs::statistics::Variance;
use std::path::Path;
use std::str;

#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub activity: f32,
}

impl Locus {
    pub fn new(chrom: String, start: u32, end: u32) -> Locus {
        Locus {
            chrom: chrom,
            start: start,
            end: end,
            activity: 0.0,
        }
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 */
pub struct RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    pub genotype_likelihoods: Vec<OrderedFloat<f32>>,
    pub final_phred_scaled_genotype_likelihoods: Vec<i32>,
    pub ref_depth: i32,
    pub non_ref_depth: i32,
    pub soft_clips: i32,
    pub read_counts: i32,
}

impl RefVsAnyResult {
    pub fn new(likelihood_capacity: usize) -> RefVsAnyResult {
        RefVsAnyResult {
            genotype_likelihoods: vec![OrderedFloat(0.0); likelihood_capacity],
            final_phred_scaled_genotype_likelihoods: vec![0; likelihood_capacity],
            ref_depth: 0,
            non_ref_depth: 0,
            soft_clips: 0,
            read_counts: 0,
        }
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    pub fn get_dp(&self) -> i32 {
        return self.ref_depth + self.non_ref_depth
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    pub fn get_ad(&self) -> Vec<i32> {
        return vec![*self.ref_depth, *self.non_ref_depth]
    }

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihood_capacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihood_capacity} is negative.
     */
    pub fn ref_vs_any_result(&mut self, likelihood_capacity: usize) {
        self.genotype_likelihoods = vec![OrderedFloat(0.0); likelihood_capacity];
        self.final_phred_scaled_genotype_likelihoods = vec![0; likelihood_capacity];
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     */
    pub fn get_genotype_likelihoods_capped_by_hom_ref_likelihood(&self) -> Vec<OrderedFloat<f32>> {
        let mut output = vec![OrderedFloat(0.0); self.genotype_likelihoods.len()];

        for i in 0..self.genotype_likelihoods.len() {
            output[i] = std::cmp::min(self.genotype_likelihoods[i], self.genotype_likelihoods[0])
        }

        return output
    }
}

pub fn collect_activity_profile(
    indexed_bam_readers: &Vec<String>,
    n_threads: usize,
    ref_idx: usize,
    per_reference_samples: usize,
    m: &clap::ArgMatches,
    genomes_and_contigs: &GenomesAndContigs,
    reference_map: &HashMap<usize, String>,
    concatenated_genomes: &Option<String>,
    flag_filters: &FlagFilter,
    tree: &Mutex<Arc<Vec<&Elem>>>,
) {


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

    let ploidy: f32 = m.value_of("ploidy").unwrap().parse().unwrap();

    let mut genotype_likelihoods = Vec::with_capacity(
        short_sample_count + long_sample_count
    );

    let mut target_ids_and_lengths = HashMap::new();

    indexed_bam_readers.into_iter().enumerate().for_each(
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
                genotype_likelihood.push(
                    update_activity_profile(
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
                        &mut target_ids,
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
                    update_activity_profile(
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
                        &mut target_ids,
                    )
                );
            }
            {
                let pb = &tree.lock().unwrap()[ref_idx + 2];

                pb.progress_bar.set_message(&format!(
                    "{}: Variant calling on sample: {}",
                    pb.key,
                    variant_matrix.get_sample_name(sample_idx),
                ));
                pb.progress_bar.inc(1);
            }
            {
                let pb = &tree.lock().unwrap()[0];
                pb.progress_bar.inc(1);
            }
        },
    );

    // return genotype_likelihoods
    calculate_activity_probabilities(
        genotype_likelihoods,
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
    ploidy: f32,
    bq: i32,
    min_soft_clip_qual: i32,
    genomes_and_contigs: &'b GenomesAndContigs,
    concatenated_genomes: &'b Option<String>,
    flag_filters: &'b FlagFilter,
    target_ids_and_lengths: &mut HashMap<usize, u64>
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

    let likelihoodcount = (ploidy + 1.) as usize;
    let log10ploidy = OrderedFloat((likelihoodcount as f32).log10());

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
                        if !target_ids_and_lengths.contains_key(&tid) {
                            target_ids_and_lengths.insert(tid, target_len);
                        }
                        // The raw activity profile.
                        // Frequency of bases not matching reference compared
                        // to depth
                        let mut genotype_likelihoods = Vec::with_capacity(target_len as usize);

                        {
                            bam_generated.fetch((tid as u32));

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
                                        let depth = pileup.depth();
                                        coverage.push(depth as f64);
                                        depth_sum += depth;
                                        let refr_base = ref_seq[pos];
                                        let mut result = RefVsAnyResult::new(likelihoodcount);
                                        for alignment in pileup.alignments() {
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
                                                alignment_context_creation(
                                                    &alignment,
                                                    &mut result,
                                                    log10ploidy,
                                                    likelihoodcount,
                                                    min_soft_clip_qual,
                                                    refr_base,
                                                    bq,
                                                )
                                            }
                                        }

                                        let denominator = OrderedFloat(result.read_counts as f32) * log10ploidy;
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

                        // Get coverage mean and standard dev
                        let contig_cov = depth_sum as f64 / target_len as f64;
                        let std_dev = 10. * coverage.std_dev();
                        contig_stats.entry(tid).or_insert(vec![contig_cov, std_dev]);

                        for (index, count) in soft_clips.iter().enumerate() {
                            if count >= 5 {
                                activity[index] += count as f32; // Add in good soft clips
                            }
                            activity[index] /= depths[index] as f32; // divide by depth
                        }
                    }
                });
        }
        _ => {}
    }

    bam_generated.finish();

    return ref_vs_any_container
}

pub fn calculate_activity_probabilities(
    genotype_likelihoods: Vec<HashMap<usize, Vec<RefVsAnyResult>>>,
    target_ids_and_lens: HashMap<usize, u64>,
    ploidy: f32,
) {
    if genotype_likelihoods.len() == 1 {
        // Faster implementation for single sample analysis
        for tid in target_ids.iter() {

        }
    } else {
        for (tid, length) in target_ids_and_lens.iter() {
            for pos in 0..length.iter() {
                let mut activity_probability = 0.;
                let mut genotypes = Vec::new();
                for sample_likelihoods in genotype_likelihoods.iter() {
                    let result = sample_likelihoods[tid][pos].genotype_likelihoods.clone();
                    genotypes.push(Genotype::build(ploidy as i32, result))
                }

                let mut fake_alleles = Allele::create_fake_alleles();

                let mut variant_context = VariantContext::build(
                    *tid, pos, pos, fake_alleles
                );

                variant_context.add_genotypes(
                    genotypes
                );
                variant_context.calculate_genotypes(
                    activity_probability,
                );

            }
        }
    }
}

pub fn calculate_genotypes(
    variant_context: VariantContext,
    activity_probability: f32,
) {
    if 180 < alleles.len() - 1 {}
}

pub calculate_most_likely_alleles(

)

/**
 * Calculate the posterior probability that a single biallelic genotype is non-ref
 *
 * The nth genotype (n runs from 0 to the sample ploidy, inclusive) contains n copies of the alt allele
 * @param log10GenotypeLikelihoods
 */
fn calculate_single_sample_biallelic_non_ref_posterior(
    log10_genotype_likelihoods: &Vec<OrderedFloat<f32>>,
    return_zero_if_ref_is_max: Bool,
) -> f32 {
    if return_zero_if_ref_is_max
        && log10_genotype_likelihoods.iter()
        .position(|&item| item == *(log10_genotype_likelihoods.iter().max().unwrap())).unwrap() == 0
        && (log10_genotype_likelihoods[0] != OrderedFloat(0.5) && log10_genotype_likelihoods.len() == 2){
        return 0.
    }

    let ploidy = log10_genotype_likelihoods.len() - 1;
}

/**
* Populates reference and non reference depth vectors
*/
fn alignment_context_creation(
    alignment: &bam::pileup::Alignment,
    result: &mut RefVsAnyResult,
    log10ploidy: OrderedFloat<f32>,
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
                non_ref_likelihood = OrderedFloat(f32::from(LogProb::from(PHREDProb(record_qual as f64))));
                ref_likelihood = OrderedFloat(f32::from(LogProb::from(PHREDProb(record_qual as f64))) + (-(3.0.log10())))
            } else {
                result.ref_depth += 1;
                ref_likelihood = OrderedFloat(f32::from(LogProb::from(PHREDProb(record_qual as f64))));
                non_ref_likelihood = OrderedFloat(f32::from(LogProb::from(PHREDProb(record_qual as f64))) + (-(3.0.log10())));
            }

            update_heterozygous_likelihood(
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
                            count_high_quality_soft_clips(
                                &cig,
                                &record,
                                result,
                                &mut read_cursor,
                                min_soft_clip_qual
                            );
                        },
                        _ => {
                            // Not a soft clip
                        }
                    }
                } else if qpos == read_cursor + 1 {
                    match cig {
                        Cigar::SoftClip(_) => {
                            count_high_quality_soft_clips(
                                &cig,
                                &record,
                                result,
                                &mut read_cursor,
                                min_soft_clip_qual
                            );

                        },
                        _ => {
                            // Not a soft clip
                        }
                    }
                } else if cigar_cursor > qpos {
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
    result: &mut RefVsAnyResult,
    cigar_cursor: &mut usize,
    min_soft_clip_qual: i32
) {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036227652?id=4147
    // If we have high quality soft clips then we want to
    // track them
    // get mean base quality
    for rpos in cigar_cursor..(cigar_cursor + cig.len() as usize) {
        let qual_pos = record.qual()[rpos];
        if qual_pos >= min_soft_clip_qual {
            result.soft_clips += 1
        }
    }
}

fn update_heterozygous_likelihood(
    result: &mut RefVsAnyResult,
    likelihoodcount: usize,
    log10ploidy: OrderedFloat<f32>,
    ref_likelihood: OrderedFloat<f32>,
    non_ref_likelihood: OrderedFloat<f32>
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
            ref_likelihood + OrderedFloat((j as f32).log10()),
            non_ref_likelihood + OrderedFloat((i as f32).log10())
        );
        i += 1;
        j -= 1;
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