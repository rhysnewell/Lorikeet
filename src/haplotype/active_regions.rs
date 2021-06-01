// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::str::FromStr;
use regex::Regex;
use std::fmt;
use serde_json;
use serde_json::Value;
use std::cmp::{max, min};
use std::fs::File;
use std::path::Path;
use ordered_float:OrderedFloat;
use rust_htslib::bam::{Read, IndexedReader};

use bird_tool_utils::command;
use rust_htslib::errors::Error;
use rust_htslib::{bcf, bcf::Read};
use std;
use std::collections::{HashMap, HashSet};

use coverm::bam_generator::*;
use coverm::FlagFilter;
use external_command_checker;
use model::variants::*;
use utils::*;

use crate::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use scoped_threadpool::Pool;
use statrs::statistics::Variance;
use std::io::Write;
use std::path::Path;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::Builder;

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
            activity: 0.0
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
    pub genotypeLikelihoods: Vec<OrderedFloat<f32>>,
    pub finalPhredScaledGenotypeLikelihoods: Vec<i32>,
    pub refDepth: i32,
    pub nonRefDepth: i32,
}

impl RefVsAnyResult {
    pub fn new(likelihoodCapacity: i32) -> RefVsAnyResult {
        RefVsAnyResult {
            genotypeLikelihoods: vec![0.0; likelihoodCapacity],
            finalPhredScaledGenotypeLikelihoods: vec![0; likelihoodCapacity],
            refDepth: 0,
            nonRefDepth: 0,
        }
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    pub fn getDP(&self) -> &i32 {
        return self.refDepth + self.nonRefDepth
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    pub fn getAD(&self) -> Vec<i32> {
        return vec![*self.refDepth, *self.nonRefDepth]
    }

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    pub fn RefVsAnyResult(&mut self, likelihoodCapacity: i32) {
        self.genotypeLikelihoods = vec![OrderedFloat(0.0); likelihoodCapacity];
        self.finalPhredScaledGenotypeLikelihoods = vec![0; likelihoodCapacity];
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     */
    pub fn getGenotypeLikelihoodsCappedByHomRefLikelihood(&self) -> Vec<OrderedFloat<f32>> {
        let mut output = vec![0.0; self.genotypeLikelihoods.len()];

        for i in 0..self.genotypeLikelihoods.len() {
            output[i] = std::cmp::min(self.genotypeLikelihoods[i], self.genotypeLikelihoods[0])
        }

        return output
    }

}

#[allow(unused)]
pub fn retrieve_active_loci<'b, R: IndexedNamedBamReader + Send, G: NamedBamReaderGenerator<R> + Send>(
    bam_generator: G,
    split_threads: usize,
    ref_idx: usize,
    sample_idx: usize,
    mut sample_count: usize,
    readtype: ReadType,
    m: &'b clap::ArgMatches,
    genomes_and_contigs: &'b GenomesAndContigs,
    reference_map: &'b HashMap<usize, String>,
    mut short_sample_count: usize,
    concatenated_genomes: &'b Option<String>,
    flag_filters: &'b FlagFilter,
) {
    let mut bam_generated = bam_generator.start();
    let mut stoit_name = bam_generated.name().to_string();

    let reference = &genomes_and_contigs.genomes[ref_idx];
    let mut reference_file = retrieve_reference(concatenated_genomes);

    bam_generated.set_threads(split_threads);

    let header = bam_generated.header().clone(); // bam header
    let target_lens: Vec<u64> = (0..header.target_count())
        .into_iter()
        .map(|tid| header.target_len(tid).unwrap())
        .collect();
    let target_names = header.target_names();

    let bam_path = bam_generated.path().to_string();

    // minimum PHRED base quality
    let bq = m
        .value_of("base-quality-threshold")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    // Minimum MAPQ value
    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());
    // Minimum count for a SNP to be considered
    let min_variant_depth: i32 = m.value_of("min-variant-depth").unwrap().parse().unwrap();
    let min_variant_quality = m
        .value_of("min-variant-quality")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let ploidy: i32 = m.value_of("ploidy").unwrap().parse().unwrap();
    let likelihoodcount = ploidy + 1;
    let


    let mut total_records = 0;
    let mut ref_target_names = Vec::new();
    let mut ref_target_lengths = Vec::new();
    let mut contig_stats: HashMap<usize, Vec<f64>> = HashMap::new();
    let valid_chars: HashSet<u8> = vec![
        "A".as_bytes()[0],
        "T".as_bytes()[0],
        "C".as_bytes()[0],
        "G".as_bytes()[0],
    ]
        .into_iter()
        .collect::<HashSet<u8>>();
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
                        // use pileups to call SNPs for low quality variants
                        // That are usually skipped by GATK
                        ref_target_names.push(target_name.clone());
                        let target_len = target_lens[tid];
                        let mut ref_seq = Vec::new();
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

                        // The raw activity profile.
                        // Frequency of bases not matching reference compared
                        // to depth
                        let mut activity = vec![0.0; target_len as usize];
                        let mut refdepths = vec![0; target_len as usize];
                        let mut nonrefdepths = vec![0; target_len as usize];


                        // Number of soft clips at target location
                        let mut soft_clips = vec![0; target_len as usize];
                        let mut coverage = Vec::with_capacity(target_len as usize);
                        let mut depth_sum = 0;
                        ref_target_lengths.push(target_len);

                        {
                            bam_generated.fetch((tid as u32));

                            while bam_generated.read(&mut record) == true
                            // Read records into the empty recorde
                            {
                                if (!flag_filters.include_supplementary
                                    && record.is_supplementary()
                                    && readtype != ReadType::Long) // We want supp alignments for longreads
                                    || (!flag_filters.include_secondary
                                    && record.is_secondary())
                                    || (!flag_filters.include_improper_pairs
                                    && !record.is_proper_pair()
                                    && readtype != ReadType::Long)
                                // Check against filter flags and current sample type
                                {
                                    continue;
                                }

                                // if reference has changed, print the last record
                                let tid = record.tid();
                                if !record.is_unmapped() {
                                    // if mapped
                                    if record.seq().len() == 0 {
                                        continue;
                                    } else if record.mapq() < mapq_threshold {
                                        continue;
                                    }

                                    if !record.is_supplementary() {
                                        num_mapped_reads_in_current_contig += 1;
                                    }

                                    // for each chunk of the cigar string
                                    let mut cursor: usize = record.pos() as usize;
                                    let quals = record.qual();
                                    let mut read_cursor: usize = 0;
                                    let read_len = record.seq().len();
                                    for cig in record.cigar().iter() {
                                        match cig {
                                            Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                                                // if M, X, or = check if bases match
                                                for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                                                    // See if read is match MNV
                                                    let qual_pos = record.qual()[qpos] as f64;
                                                    let read_char = record.seq()[qpos];
                                                    let refr_base = ref_seq[cursor + qpos];
                                                    depths[cursor + qpos] += 1;
                                                    depth_sum += 1;

                                                    if refr_base != "N".as_bytes()[0] {
                                                        if refr_base != read_char {
                                                            let record_qual = quals[qpos] as f64;
                                                            if record_qual >= bq {
                                                                activity[cursor + qpos] += 1.
                                                            }
                                                        }
                                                    }
                                                }
                                                read_cursor += cig.len() as usize;
                                            }
                                            Cigar::Del(del) => {
                                                let record_qual = quals[read_cursor] as f64;
                                                if record_qual >= bq {
                                                    activity[cursor] += 1.
                                                }
                                                cursor += cig.len() as usize;
                                            }
                                            Cigar::RefSkip(del) => {
                                                let record_qual = quals[read_cursor] as f64;
                                                if record_qual >= bq {
                                                    activity[cursor] += 1.
                                                }
                                                cursor += cig.len() as usize;
                                            }
                                            Cigar::Ins(ins) => {
                                                let qual_sum = quals[read_cursor..cig.len()].sum() as f64;
                                                if (qual_sum / cig.len()) >= bq {
                                                    activity[cursor] += 1.
                                                }
                                                read_cursor += cig.len() as usize;
                                                total_indels_in_current_contig += cig.len() as u64;
                                            }
                                            Cigar::SoftClip(_) => {
                                                // https://gatk.broadinstitute.org/hc/en-us/articles/360036227652?id=4147
                                                // If we have high quality soft clips then we want to
                                                // track them
                                                for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                                                    // See if read is match MNV
                                                    let qual_pos = record.qual()[qpos] as f64;
                                                    if qual_pos >= std::cmp::min(29, 2 * bq) {
                                                        soft_clips[cursor + qpos] += 1
                                                    }
                                                }

                                                read_cursor += cig.len() as usize;
                                            }
                                            Cigar::HardClip(_) | Cigar::Pad(_) => {}
                                        }
                                    }
                                }
                            }
                        };

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