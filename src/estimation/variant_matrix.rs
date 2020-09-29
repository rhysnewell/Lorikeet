use bio::io::gff;
use bio::stats::{
    bayesian,
    probs::{LogProb, PHREDProb, Prob},
};
use bird_tool_utils::command;
use coverm::genomes_and_contigs::GenomesAndContigs;
use dbscan::fuzzy;
use estimation::codon_structs::CodonTable;
use estimation::codon_structs::*;
use estimation::contig_variants::*;
use estimation::genotype_abundances;
use estimation::linkage::*;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::izip;
use itertools::Itertools;
use model::variants::*;
use ordered_float::NotNan;
use rayon::current_num_threads;
use rayon::prelude::*;
use rust_htslib::bcf::{self};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::process::Stdio;
use std::str;
use std::sync::{Arc, Mutex};
use tempfile;
use utils::generate_faidx;

#[derive(Debug, Clone)]
/// Container for all variants within a genome and associated clusters
pub enum VariantMatrix<'b> {
    VariantContigMatrix {
        coverages: HashMap<i32, Vec<f64>>,
        average_genotypes: HashMap<i32, Vec<f64>>,
        variances: HashMap<i32, Vec<f64>>,
        // Ref idx, TID, Position, Variant Type, Base
        all_variants: HashMap<usize, HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>>,
        // Placeholder hashmap for the depths of each contig for a sample
        // Deleted after use
        depths: HashMap<i32, Vec<i32>>,
        target_names: HashMap<usize, BTreeMap<i32, String>>,
        target_lengths: HashMap<usize, HashMap<i32, f64>>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        variant_counts: HashMap<usize, HashMap<i32, usize>>,
        variant_sums: HashMap<usize, HashMap<i32, Vec<Vec<f64>>>>,
        variant_info: HashMap<usize, Vec<fuzzy::Var>>,
        geom_mean_var: HashMap<usize, Vec<f64>>,
        geom_mean_dep: HashMap<usize, Vec<f64>>,
        geom_mean_frq: HashMap<usize, Vec<f64>>,
        pred_variants: HashMap<
            usize,
            HashMap<usize, HashMap<i32, HashMap<i64, HashMap<fuzzy::Category, HashSet<Variant>>>>>,
        >,
        _m: std::marker::PhantomData<&'b ()>,
        //        pred_variants_all: HashMap<usize, HashMap<i32, HashMap<i32, HashSet<String>>>>,
    },
}

impl<'b> VariantMatrix<'b> {
    pub fn new_matrix(sample_count: usize) -> VariantMatrix<'b> {
        VariantMatrix::VariantContigMatrix {
            coverages: HashMap::new(),
            variances: HashMap::new(),
            average_genotypes: HashMap::new(),
            all_variants: HashMap::new(),
            depths: HashMap::new(),
            target_names: HashMap::new(),
            target_lengths: HashMap::new(),
            sample_names: vec!["".to_string(); sample_count],
            kfrequencies: BTreeMap::new(),
            variant_counts: HashMap::new(),
            variant_sums: HashMap::new(),
            variant_info: HashMap::new(),
            geom_mean_var: HashMap::new(),
            geom_mean_dep: HashMap::new(),
            geom_mean_frq: HashMap::new(),
            pred_variants: HashMap::new(),
            _m: std::marker::PhantomData::default(),
        }
    }
}

pub trait VariantMatrixFunctions {
    fn setup(&mut self);

    /// returns sample name by index
    fn get_sample_name(&self, sample_idx: usize) -> &str;

    /// Returns the total amount of variant alleles
    fn get_variant_count(&self, ref_idx: usize) -> i64;

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize);

    fn add_info(&mut self, ref_idx: usize, tid: usize, target_name: Vec<u8>, target_len: u64);

    /// Adds the variant information retrieved from VCF records on a per reference basis
    fn add_variant_to_matrix(&mut self, sample_idx: usize, base: &Base, tid: usize, ref_idx: usize);

    /// Per sample FDR calculation
    fn remove_false_discoveries(&mut self, alpha: f64, reference: &str);

    /// Returns the alleles at the current position
    /// as a mutable reference
    fn variants(
        &mut self,
        ref_idx: usize,
        tid: i32,
        pos: i64,
    ) -> Option<&mut HashMap<Variant, Base>>;

    /// Returns all variants found within a contig
    fn variants_of_contig(
        &mut self,
        ref_idx: usize,
        tid: i32,
    ) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>>;

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(
        &mut self,
        variant_stats: VariantStats,
        sample_count: usize,
        sample_idx: usize,
        reference_index: usize,
    );

    /// Converts all variants into fuzzy::Var format
    fn generate_distances(&mut self);

    /// Perform fuzzy DBSCAN clustering using proportionality
    fn run_fuzzy_scan(
        &mut self,
        e_min: f64,
        e_max: f64,
        pts_min: f64,
        pts_max: f64,
        phi: f64,
        anchor_size: usize,
        anchor_similarity: f64,
        minimum_reads_in_link: usize,
        reference_map: &HashMap<usize, String>,
        multi: &Arc<MultiProgress>,
    );

    /// Takes clusters from DBSCAN and linkage method and writes variants to file as genotype
    fn generate_genotypes(
        &mut self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
    );

    /// EM Algorithm for caclulating the abundances of each genotype in each sample
    fn calculate_strain_abundances(
        &mut self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
    );

    /// Prints the per reference and per contig variant information e.g. How many SNPs were seen
    /// along a contig over the given window size on average
    fn print_variant_stats(
        &self,
        window_size: f64,
        output_prefix: &str,
        genomes_and_contigs: &GenomesAndContigs,
    );

    /// Calculates the dN/dS values for each gene
    fn calc_gene_mutation(
        &self,
        gff_map: &mut HashMap<usize, HashMap<String, Vec<bio::io::gff::Record>>>,
        genomes_and_contigs: &GenomesAndContigs,
        reference_map: &HashMap<usize, String>,
        codon_table: &CodonTable,
        output_prefix: &str,
        concatenated_genomes: &Option<String>,
    );

    fn polish_genomes(
        &self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
    );

    fn calculate_sample_distances(
        &self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
    );

    fn write_vcf(&self, output_prefix: &str, genomes_and_contigs: &GenomesAndContigs);
}

#[allow(unused)]
impl VariantMatrixFunctions for VariantMatrix<'_> {
    fn setup(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ref mut kfrequencies,
                ..
            } => {
                *coverages = HashMap::new();
                *average_genotypes = HashMap::new();
                *all_variants = HashMap::new();
                *target_names = HashMap::new();
                *target_lengths = HashMap::new();
                *sample_names = vec![];
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn get_variant_count(&self, ref_idx: usize) -> i64 {
        match self {
            VariantMatrix::VariantContigMatrix { variant_counts, .. } => {
                match variant_counts.get(&ref_idx) {
                    Some(contigs_map) => {
                        let mut total = 0;
                        for (_, count) in contigs_map {
                            total += count;
                        }
                        total as i64
                    }
                    None => 0,
                }
            }
        }
    }

    fn get_sample_name(&self, sample_idx: usize) -> &str {
        match self {
            VariantMatrix::VariantContigMatrix { sample_names, .. } => {
                match &sample_names[sample_idx][..4] {
                    ".tmp" => &sample_names[sample_idx][15..],
                    _ => &sample_names[sample_idx],
                }
            }
        }
    }

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ..
            } => {
                debug!("adding sample {} at index {}", &sample_name, &sample_idx);
                sample_names[sample_idx] = sample_name;
            }
        }
    }

    fn add_info(&mut self, ref_idx: usize, tid: usize, target_name: Vec<u8>, target_len: u64) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                let ref_target_names = target_names.entry(ref_idx).or_insert(BTreeMap::new());

                ref_target_names
                    .entry(tid as i32)
                    .or_insert(String::from_utf8(target_name).unwrap());

                // let target_len = header.target_len(tid).unwrap();
                let ref_target_lengths = target_lengths.entry(ref_idx).or_insert(HashMap::new());

                ref_target_lengths
                    .entry(tid as i32)
                    .or_insert(target_len as f64);
            }
        }
    }

    fn add_variant_to_matrix(
        &mut self,
        sample_idx: usize,
        base: &Base,
        tid: usize,
        ref_idx: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ref mut variant_counts,
                ..
            } => {
                // Initialize contig id in variant hashmap
                let reference_variants = all_variants.entry(ref_idx).or_insert(HashMap::new());

                let contig_variants = reference_variants
                    .entry(tid as i32)
                    .or_insert(HashMap::new());

                let mut ref_var_counts = variant_counts.entry(ref_idx).or_insert(HashMap::new());
                let mut con_var_counts = ref_var_counts.entry(tid as i32).or_insert(0);
                // Apppend the sample index to each variant abundance
                // Initialize the variant position index
                // Also turns out to be the total number of variant positions
                let position_variants = contig_variants.entry(base.pos).or_insert(HashMap::new());
                let allele = position_variants
                    .entry(base.variant.clone())
                    .or_insert(base.clone());
                let allele = position_variants
                    .entry(base.variant.clone())
                    .or_insert(base.clone());
                allele.combine_sample(base, sample_idx, 0);
                *con_var_counts += 1;

                // if position_variants.contains_key(&base.variant) {
                //     let allele = position_variants.entry(base.variant.clone()).or_insert(base.clone());
                //     allele.combine_sample(base, sample_idx, 0);
                // } else {
                //     let allele = position_variants.entry(base.variant.clone()).or_insert(base.clone());
                //     *con_var_counts += 1;
                // }
            }
        }
    }

    fn remove_false_discoveries(&mut self, alpha: f64, reference: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                // collect variant probabilities e.g. PHRED Sums
                let mut prob_dist = Vec::new();
                for (ref_idx, contig_map) in all_variants.iter() {
                    for (tid, position_map) in contig_map.iter() {
                        position_map.iter().for_each(|(pos, alleles)| {
                            for (variant, allele) in alleles.iter() {
                                if allele.variant != Variant::None {
                                    // get the mean PHRED value
                                    // This is probably not the best method
                                    let mut allele_probs = allele.quals.par_iter().sum();
                                    // allele_probs.par_sort();

                                    // let mut allele_probs = allele_probs
                                    //     .iter()
                                    //     .map(|p| {
                                    //         LogProb::from(Prob(10_f64.powf(-(**p as f64) / 10)))
                                    //     })
                                    //     .collect_vec();

                                    prob_dist.push(
                                        NotNan::new(allele_probs)
                                            .expect("Unable to convert to NotNan"),
                                    );
                                }
                            }
                        });
                    }
                }
                // sort prob_dist
                prob_dist.par_sort();
                // turn into logprob values
                let prob_dist = prob_dist
                    .into_iter()
                    .rev()
                    .map(|p| LogProb::from(PHREDProb(*p)))
                    .collect_vec();

                // info!("Prob {:?}", &prob_dist);

                // turn prob dist into posterior exact probabilities
                debug!("Calculating PEP dist...");
                let pep_dist = prob_dist
                    .iter()
                    .map(|p| p.ln_one_minus_exp())
                    .collect::<Vec<LogProb>>();

                // info!("PEP {:?}", &pep_dist);

                // calculate fdr values
                let fdrs = bayesian::expected_fdr(&pep_dist);
                // info!("FDRs {:?}", &fdrs);

                let alpha = LogProb::from(Prob(alpha));
                let mut threshold = None;
                if fdrs.is_empty() {
                    debug!("Genome {}: FDR calculations are empty...", reference);
                    threshold = None;
                } else if fdrs[0] > alpha {
                    debug!(
                        "Genome {}: FDR threshold for alpha of {:?} calculated as {:?}",
                        reference,
                        alpha,
                        LogProb::ln_one()
                    );
                    threshold = Some(LogProb::ln_one());
                } else {
                    // find the largest pep for which fdr <= alpha
                    // do not let peps with the same value cross the boundary
                    for i in (0..fdrs.len()).rev() {
                        if fdrs[i] <= alpha && (i == 0 || pep_dist[i] != pep_dist[i - 1]) {
                            let prob = prob_dist[i];
                            debug!(
                                "Genome {}: FDR threshold for alpha of {:?} calculated as {:?}",
                                reference, alpha, &prob
                            );
                            threshold = Some(prob);
                            break;
                        }
                    }
                }
                match threshold {
                    Some(prob) => {
                        let mut total_removed = 0;
                        let mut total_kept = 0;
                        for (ref_idx, contig_map) in all_variants.iter_mut() {
                            for (tid, position_map) in contig_map.iter_mut() {
                                let mut positions_to_remove = Vec::new();
                                for (pos, alleles) in position_map.iter_mut() {
                                    let mut to_remove = Vec::new();
                                    let mut variant_alleles = 0;
                                    for (variant, allele) in alleles.iter_mut() {
                                        if allele.variant != Variant::None {
                                            let mut allele_probs = allele.quals.par_iter().sum();
                                            let sum_prob = LogProb::from(PHREDProb(allele_probs));
                                            if sum_prob > prob {
                                                to_remove.push(variant.clone());
                                                total_removed += 1;
                                            } else {
                                                total_kept += 1;
                                            }
                                            variant_alleles += 1;
                                        }
                                    }
                                    if to_remove.len() > 1 {
                                        if to_remove.len() == variant_alleles {
                                            positions_to_remove.push(*pos);
                                        } else {
                                            for removing in to_remove.iter() {
                                                alleles
                                                    .remove_entry(removing)
                                                    .expect("Unable to remove variant entry");
                                            }
                                        }
                                    }
                                }
                                if positions_to_remove.len() > 0 {
                                    for position in positions_to_remove.iter() {
                                        position_map
                                            .remove_entry(position)
                                            .expect("Unable to remove position");
                                    }
                                }
                            }
                        }
                        debug!("Genome {}: {} variants passed FDR threshold, {} did not pass FDR threshold", reference, total_kept, total_removed);
                        *all_variants = all_variants.clone();
                    }
                    None => {}
                }
            }
        }
    }

    fn variants(
        &mut self,
        ref_idx: usize,
        tid: i32,
        pos: i64,
    ) -> Option<&mut HashMap<Variant, Base>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => match all_variants.get_mut(&ref_idx) {
                Some(ref_variants) => match ref_variants.get_mut(&tid) {
                    Some(contig_variants) => contig_variants.get_mut(&pos),
                    None => None,
                },
                _ => None,
            },
        }
    }

    fn variants_of_contig(
        &mut self,
        ref_idx: usize,
        tid: i32,
    ) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => match all_variants.get_mut(&ref_idx) {
                Some(ref_variants) => ref_variants.get_mut(&tid),
                _ => None,
            },
        }
    }

    fn add_contig(
        &mut self,
        variant_stats: VariantStats,
        sample_count: usize,
        sample_idx: usize,
        reference_index: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut all_variants,
                target_names,
                //                ref mut target_lengths,
                ref mut variances,
                ref mut depths,
                ..
            } => {
                match variant_stats {
                    VariantStats::VariantContigStats {
                        tid,
                        //                        target_name,
                        //                        target_len,
                        coverage,
                        variance,
                        depth,
                        //                        variations_per_n,
                        ..
                    } => {
                        let var = variances
                            .entry(tid)
                            .or_insert(vec![0.0 as f64; sample_count]);
                        var[sample_idx] = variance;
                        let cov = coverages
                            .entry(tid)
                            .or_insert(vec![0.0 as f64; sample_count]);
                        cov[sample_idx] = coverage;

                        // // reference index from contig index
                        // let reference_index = genomes_and_contigs
                        //     .genome_index_of_contig(
                        //         &String::from_utf8(
                        //             target_names[][tid as usize].to_vec()).unwrap()).unwrap();
                        // copy across depths
                        let ref_variants = all_variants
                            .entry(reference_index)
                            .or_insert(HashMap::new());
                        let contig_variants = ref_variants.entry(tid).or_insert(HashMap::new());
                        for (pos, d) in depth.iter().enumerate() {
                            let position_variants =
                                contig_variants.entry(pos as i64).or_insert(HashMap::new());
                            for (_variant, base_info) in position_variants.iter_mut() {
                                base_info.add_depth(sample_idx, *d);
                            }
                        }
                        depths.entry(tid).or_insert(depth);
                    }
                }
            }
        }
    }

    fn generate_distances(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                all_variants,
                sample_names,
                coverages,
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut geom_mean_frq,
                ..
            } => {
                let sample_count = sample_names.len();

                let variant_info_all = Arc::new(Mutex::new(HashMap::new()));

                // A vector the length of the number of samples
                // Cumulatively calculates the product of variant depths
                let geom_mean_v_all = Arc::new(Mutex::new(HashMap::new()));

                // product of total depth
                let geom_mean_d_all = Arc::new(Mutex::new(HashMap::new()));

                // product of reference frequency
                let geom_mean_f_all = Arc::new(Mutex::new(HashMap::new()));

                all_variants
                    .par_iter_mut()
                    .for_each(|(ref_idx, ref_variants)| {
                        let variant_info_ref = Arc::new(Mutex::new(Vec::new()));

                        // A vector the length of the number of samples
                        // Cumulatively calculates the product of variant depths
                        let geom_mean_v =
                            Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                        // product of total depth
                        let geom_mean_d =
                            Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                        // product of reference frequency
                        let geom_mean_f =
                            Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                        // get basic variant info and store as fuzzy::Var

                        ref_variants
                            .par_iter_mut()
                            .for_each(|(tid, variant_abundances)| {
                                variant_abundances
                                    .par_iter_mut()
                                    .for_each(|(position, hash)| {
                                        // loop through each position that has variants ignoring positions that
                                        // only contained the reference in all samples
                                        if hash.keys().len() > 0 {
                                            for (variant, base_info) in hash.iter_mut() {
                                                match variant {
                                                    Variant::SNV(_) | Variant::None => {
                                                        let _abundance: f64 = 0.;
                                                        let _mean_var: f64 = 0.;
                                                        let mut rel_abund =
                                                            vec![0.0; sample_count as usize];

                                                        // Get the mean abundance across samples
                                                        for index in (0..sample_count).into_iter() {
                                                            let mut geom_mean_v =
                                                                geom_mean_v.lock().unwrap();
                                                            let mut geom_mean_d =
                                                                geom_mean_d.lock().unwrap();
                                                            let mut geom_mean_f =
                                                                geom_mean_f.lock().unwrap();

                                                            let mut var_depth =
                                                                base_info.truedepth[index] as f64;

                                                            let total_depth =
                                                                base_info.totaldepth[index] as f64;
                                                            //                                                println!("var_depth {} tot {}", var_depth, total_depth);
                                                            //                                                base_info.freq[index] = ;
                                                            if total_depth <= 0. {
                                                                rel_abund[index] = var_depth / 1.;
                                                            } else {
                                                                rel_abund[index] =
                                                                    var_depth / total_depth;
                                                            }

                                                            geom_mean_v[index] +=
                                                                (var_depth + 1.).ln();
                                                            geom_mean_d[index] +=
                                                                (total_depth + 1.).ln();
                                                            geom_mean_f[index] += ((var_depth
                                                                + 1.)
                                                                / (total_depth + 1.))
                                                                .ln();
                                                        }

                                                        let mut variant_info_ref =
                                                            variant_info_ref.lock().unwrap();
                                                        //                                            base_info.rel_abunds = rel_abund;
                                                        let point = fuzzy::Var {
                                                            pos: *position,
                                                            var: variant.clone(),
                                                            deps: base_info.totaldepth.clone(),
                                                            vars: base_info.depth.clone(),
                                                            //                                                rel_abunds: rel_abund,
                                                            tid: *tid,
                                                            reads: base_info.reads.clone(),
                                                        };
                                                        variant_info_ref.push(point);
                                                    }
                                                    _ => {
                                                        let mut rel_abund =
                                                            vec![0.0; sample_count as usize];

                                                        // Get the mean abundance across samples
                                                        for index in (0..sample_count).into_iter() {
                                                            let mut geom_mean_v =
                                                                geom_mean_v.lock().unwrap();
                                                            let mut geom_mean_d =
                                                                geom_mean_d.lock().unwrap();
                                                            let mut geom_mean_f =
                                                                geom_mean_f.lock().unwrap();

                                                            let mut var_depth =
                                                                base_info.depth[index] as f64;
                                                            if var_depth <= 0. {
                                                                var_depth = base_info.truedepth
                                                                    [index]
                                                                    as f64;
                                                                base_info.depth[index] =
                                                                    base_info.truedepth[index];
                                                            }
                                                            let total_depth =
                                                                base_info.totaldepth[index] as f64;
                                                            //                                                base_info.freq[index] = ;
                                                            if total_depth <= 0. {
                                                                rel_abund[index] = var_depth / (1.);
                                                            } else {
                                                                rel_abund[index] =
                                                                    var_depth / total_depth;
                                                            }

                                                            geom_mean_v[index] +=
                                                                (var_depth + 1.).ln();
                                                            geom_mean_d[index] +=
                                                                (total_depth + 1.).ln();
                                                            geom_mean_f[index] += ((var_depth
                                                                + 1.)
                                                                / (total_depth + 1.))
                                                                .ln();
                                                        }

                                                        let mut variant_info_ref =
                                                            variant_info_ref.lock().unwrap();
                                                        //                                            base_info.rel_abunds = rel_abund;
                                                        let point = fuzzy::Var {
                                                            pos: *position,
                                                            var: variant.clone(),
                                                            deps: base_info.totaldepth.clone(),
                                                            vars: base_info.depth.clone(),
                                                            //                                                rel_abunds: rel_abund,
                                                            tid: *tid,
                                                            reads: base_info.reads.clone(),
                                                        };
                                                        variant_info_ref.push(point);
                                                    }
                                                }
                                            }
                                        }
                                    });
                            });

                        // Add variant info for reference to hashmap
                        let variant_info_ref = variant_info_ref.lock().unwrap().clone();
                        let mut variant_info_all = variant_info_all.lock().unwrap();
                        variant_info_all.insert(*ref_idx, variant_info_ref);

                        // Helper fn to calculate geom mean from sum of log values
                        let geom_mean = |input: &Vec<f64>| -> Vec<f64> {
                            let output = input
                                .iter()
                                .map(|sum| (sum / variant_info_all[ref_idx].len() as f64).exp())
                                .collect::<Vec<f64>>();
                            return output;
                        };

                        debug!(
                            "geoms {:?} {:?} {:?}",
                            geom_mean_d, geom_mean_v, geom_mean_f
                        );

                        let geom_mean_v = geom_mean_v.lock().unwrap().clone();
                        let geom_mean_v = geom_mean(&geom_mean_v);
                        debug!("Geom Mean Var {:?}", geom_mean_v);

                        let geom_mean_d = geom_mean_d.lock().unwrap().clone();
                        let geom_mean_d = geom_mean(&geom_mean_d);
                        debug!("Geom Mean Dep {:?}", geom_mean_d);

                        let geom_mean_f = geom_mean_f.lock().unwrap().clone();
                        let geom_mean_f = geom_mean(&geom_mean_f);
                        debug!("Geom Mean Frq {:?}", geom_mean_f);
                        debug!(
                            "geoms {:?} {:?} {:?}",
                            geom_mean_d, geom_mean_v, geom_mean_f
                        );

                        let mut geom_mean_v_all = geom_mean_v_all.lock().unwrap();
                        let mut geom_mean_d_all = geom_mean_d_all.lock().unwrap();
                        let mut geom_mean_f_all = geom_mean_f_all.lock().unwrap();

                        // Add geom means per ref to geom mean maps
                        geom_mean_v_all.insert(*ref_idx, geom_mean_v);
                        geom_mean_d_all.insert(*ref_idx, geom_mean_d);
                        geom_mean_f_all.insert(*ref_idx, geom_mean_f);
                    });

                let variant_info_all = variant_info_all.lock().unwrap().clone();
                let geom_mean_v_all = geom_mean_v_all.lock().unwrap().clone();
                let geom_mean_d_all = geom_mean_d_all.lock().unwrap().clone();
                let geom_mean_f_all = geom_mean_f_all.lock().unwrap().clone();

                *variant_info = variant_info_all;
                *geom_mean_var = geom_mean_v_all;
                *geom_mean_dep = geom_mean_d_all;
                *geom_mean_frq = geom_mean_f_all;
            }
        }
    }

    fn run_fuzzy_scan(
        &mut self,
        e_min: f64,
        e_max: f64,
        pts_min: f64,
        pts_max: f64,
        phi: f64,
        anchor_size: usize,
        anchor_similarity: f64,
        minimum_reads_in_link: usize,
        reference_map: &HashMap<usize, String>,
        multi: &Arc<MultiProgress>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut geom_mean_frq,
                ref mut pred_variants,
                ref mut all_variants,
                target_lengths,
                ..
            } => {
                let prediction_variants_all = Arc::new(Mutex::new(HashMap::new()));

                let all_variants_genotyped = Arc::new(Mutex::new(all_variants.clone()));
                // For each reference genome we will perform the DBSCAN clustering
                variant_info.par_iter().for_each(|(ref_idx, variant_info_vec)| {
                    if variant_info_vec.len() > 1 {
                        let mut fuzzy_scanner = fuzzy::FuzzyDBSCAN {
                            eps_min: e_min,
                            eps_max: e_max,
                            pts_min: match pts_min {
                                _ if pts_min > 1. => pts_min,
                                _ => pts_min * variant_info_vec.len() as f64,
                            },
                            pts_max: match pts_max {
                                _ if pts_max > 1. => pts_max,
                                _ => pts_max * variant_info_vec.len() as f64,
                            },
                            phi,
                            geom_var: geom_mean_var[ref_idx].to_owned(),
                            geom_dep: geom_mean_dep[ref_idx].to_owned(),
                            geom_frq: geom_mean_frq[ref_idx].to_owned(),
                        };


                        let links = Vec::new();
                        // run fuzzy DBSCAN
                        let reference_path = Path::new(
                            reference_map
                                .get(ref_idx)
                                .expect("reference index not found"),
                        );
                        if links.len() > 0 {
                            debug!(
                                "Genome {}: Running Seeded fuzzyDBSCAN with {} initial clusters",
                                reference_path.file_stem().unwrap().to_str().unwrap(),
                                links.len(),
                            );
                        } else {
                            debug!(
                                "Genome {}: No initial clusters formed, running fuzzyDBSCAN with no seeds.",
                                reference_path.file_stem().unwrap().to_str().unwrap(),
                            )
                        }
                        let clusters = fuzzy_scanner.cluster(
                            &variant_info_vec[..],
                            links,
                            reference_path.file_stem().unwrap().to_str().unwrap(),
                            multi,
                            *ref_idx,
                        );

                        // Perform read phasing clustering and return expanded clusters
                        let clusters: Vec<Vec<fuzzy::Assignment>> = linkage_clustering_of_variants(
                            clusters,
                            &variant_info_vec,
                            anchor_size,
                            anchor_similarity,
                            minimum_reads_in_link,
                            multi,
                            *ref_idx,
                        );

                        // Since these are hashmaps, I'm using Arc and Mutex here since not sure how
                        // keep hashmaps updated using channel()
                        let prediction_variants =
                            Arc::new(
                            Mutex::new(
                                HashMap::new()));

                        // Organize clusters into genotypes by recollecting full variant information
                        clusters.par_iter().enumerate().for_each(|(rank, cluster)| {
                            // Sets for each cluster keeping track of which variant types are present in
                            // a cluster
                            let prediction_set = Arc::new(Mutex::new(
                                HashSet::new()));
                            cluster.par_iter().for_each(|assignment| {
                                let variant: &fuzzy::Var = &variant_info_vec[assignment.index];

                                let mut prediction_set = prediction_set.lock().unwrap();
                                prediction_set.insert(variant.var.to_owned());

                                let mut prediction_variants = prediction_variants
                                    .lock()
                                    .unwrap();

                                let variant_tid = prediction_variants
                                    .entry(rank + 1)
                                    .or_insert(HashMap::new());

                                let variant_pos = variant_tid
                                    .entry(variant.tid)
                                    .or_insert(HashMap::new());

                                let variant_cat = variant_pos
                                    .entry(variant.pos)
                                    .or_insert(HashMap::new());

                                let variant_set = variant_cat
                                    .entry(assignment.category)
                                    .or_insert(HashSet::new());

                                variant_set.insert(variant.var.to_owned());

                                // Add genotypes to base
                                let mut all_variants_genotyped =
                                    all_variants_genotyped.lock().unwrap();

                                let ref_variants =
                                    all_variants_genotyped
                                        .entry(*ref_idx)
                                        .or_insert(HashMap::new());

                                let base_tid =
                                    ref_variants
                                        .entry(variant.tid)
                                        .or_insert(HashMap::new());

                                let base_pos = base_tid
                                    .entry(variant.pos)
                                    .or_insert(HashMap::new());

                                match base_pos.get_mut(&variant.var) {
                                    Some(base_base) => {
                                        base_base.genotypes.insert(rank as i32 + 1);
                                    },
                                    None => {},
                                };
                            });
                        });
                        let prediction_variants = prediction_variants.lock().unwrap().clone();
                        let mut prediction_variants_all = prediction_variants_all.lock().unwrap();
                        prediction_variants_all.insert(*ref_idx, prediction_variants);
                    }
                });

                let all_variants_genotyped = all_variants_genotyped.lock().unwrap();

                let prediction_variants_all = prediction_variants_all.lock().unwrap().clone();

                *pred_variants = prediction_variants_all;
                *all_variants = all_variants_genotyped.clone();
            }
        }
    }

    #[allow(unused)]
    fn generate_genotypes(
        &mut self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                ref mut pred_variants,
                ..
            } => {
                pred_variants
                    .par_iter()
                    .for_each(|(ref_index, strain_map)| {
                        strain_map.par_iter().for_each(|(strain_index, genotype)| {
                            let reference_path = Path::new(
                                reference_map
                                    .get(ref_index)
                                    .expect("reference index not found"),
                            );

                            // info!("{:?}", reference_path);

                            let mut reference = match concatenated_genomes {
                                Some(temp_file) => {
                                    match bio::io::fasta::IndexedReader::from_file(&temp_file) {
                                        Ok(reader) => reader,
                                        Err(_e) => generate_faidx(&temp_file),
                                    }
                                },
                                None => {
                                    match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                                        Ok(reader) => reader,
                                        Err(_e) => generate_faidx(&reference_path.to_str().unwrap()),
                                    }
                                }
                            };


                            let file_name = format!(
                                "{}/{}_strain_{}.fna",
                                output_prefix.to_string(),
                                reference_path.file_stem().unwrap().to_str().unwrap(),
                                strain_index,
                            );

                            debug!("Genotype file name {}", &file_name);

                            let file_path = Path::new(&file_name);
                            // Open haplotype file or create one
                            let mut file_open = File::create(file_path)
                                .expect("No Read or Write Permission in current directory");

                            let mut original_contig = Vec::new();

                            let mut genotype = genotype.clone();

                            // Generate the variant genome
                            debug!(
                                "genotyping Reference index {} target names {:?}",
                                &ref_index, &target_names
                            );
                            let mut total_variant_alleles = 0;
                            let mut total_reference_alleles = 0;

                            // If we find multiple variants at one site, then we
                            // have overfitted this gentoype
                            let mut overclustered = false;
                            let mut multivariant_sites = 0;

                            for (tid, target_name) in target_names[ref_index].iter() {
                                let mut contig = String::new();
                                original_contig = Vec::new();
                                {
                                    match reference.fetch_all(
                                        std::str::from_utf8(target_name.as_bytes()).unwrap(),
                                    ) {
                                        Ok(reference) => reference,
                                        Err(e) => {
                                            warn!("Cannot read sequence from reference {:?}", e);
                                            std::process::exit(1)
                                        }
                                    };
                                    match reference.read(&mut original_contig) {
                                        Ok(reference) => reference,
                                        Err(e) => {
                                            warn!("Cannot read sequence from reference {:?}", e);
                                            std::process::exit(1)
                                        }
                                    };
                                }
                                let mut skip_n = 0;
                                let mut skip_cnt = 0;
                                //                            let char_cnt = 0;
                                let mut variations = 0;
                                let mut ref_alleles = 0;

                                if genotype.contains_key(&tid) {
                                    for (pos, base) in original_contig.iter().enumerate() {
                                        if skip_cnt < skip_n {
                                            skip_cnt += 1;
                                        } else {
                                            skip_n = 0;
                                            skip_cnt = 0;
                                            //                                if genotype.contains_key(&tid) {
                                            let tid_genotype = genotype.get_mut(&tid).unwrap();

                                            if tid_genotype.contains_key(&(pos as i64)) {
                                                let categories = &genotype[&tid][&(pos as i64)];
                                                let hash;
                                                if categories.contains_key(&fuzzy::Category::Core) {
                                                    hash =
                                                        categories[&fuzzy::Category::Core].clone();
                                                } else if categories
                                                    .contains_key(&fuzzy::Category::Border)
                                                {
                                                    hash = categories[&fuzzy::Category::Border]
                                                        .clone();
                                                } else {
                                                    hash =
                                                        categories[&fuzzy::Category::Noise].clone();
                                                }

                                                let mut max_var = Variant::None;

                                                for var in hash.iter() {
                                                    // If there are two variants possible for
                                                    // a single site and one is the reference
                                                    // we will choose the reference
                                                    if max_var == Variant::None {
                                                        max_var = var.clone();
                                                    }
                                                }
                                                if hash.len() > 1 {
                                                    multivariant_sites += 1;
                                                    overclustered = true;
                                                    debug!("Multi hash {:?} {:?}", hash, max_var);
                                                }
                                                match max_var {
                                                    Variant::Deletion(size) => {
                                                        // Skip the next n bases but rescue the reference prefix
                                                        skip_n = size;
                                                        skip_cnt = 0;
                                                        contig = contig
                                                            + str::from_utf8(&[*base]).unwrap();
                                                        // If we had the sequence we would rescue first base like this
                                                        //                                                let first_byte = max_var.as_bytes()[0];
                                                        //                                                contig = contig + str::from_utf8(
                                                        //                                                    &[first_byte]).unwrap();
                                                        variations += 1;
                                                    }
                                                    Variant::Insertion(alt) => {
                                                        // Remove prefix from variant
                                                        let removed_first_base =
                                                            str::from_utf8(&alt[1..]).unwrap();
                                                        contig = contig + removed_first_base;
                                                        variations += 1;
                                                    }
                                                    Variant::Inversion(alt) | Variant::MNV(alt) => {
                                                        // Skip the next n bases
                                                        skip_n = alt.len() as u32 - 1;
                                                        skip_cnt = 0;
                                                        // Inversions and MNVs don't have a first base prefix, so take
                                                        // wholes tring
                                                        let inversion =
                                                            str::from_utf8(&alt).unwrap();
                                                        contig = contig + inversion;
                                                        variations += 1;
                                                    }
                                                    Variant::None => {
                                                        contig = contig
                                                            + str::from_utf8(&[*base]).unwrap();
                                                        ref_alleles += 1;
                                                    }
                                                    Variant::SNV(alt) => {
                                                        contig = contig
                                                            + str::from_utf8(&[alt]).unwrap();
                                                        variations += 1;
                                                    }
                                                    _ => {
                                                        contig = contig
                                                            + str::from_utf8(&[*base]).unwrap();
                                                        ref_alleles += 1;
                                                    }
                                                }
                                            } else {
                                                contig = contig + str::from_utf8(&[*base]).unwrap();
                                            }
                                            //                                } else {
                                            //                                    contig = str::from_utf8(&original_contig)
                                            //                                        .expect("Can't convert to str").to_string();
                                            //                                }
                                        }
                                    }
                                } else {
                                    contig = str::from_utf8(&original_contig)
                                        .expect("Can't convert to str")
                                        .to_string();
                                }
                                writeln!(
                                    file_open,
                                    // ">{}_strain_{}_alt_alleles_{}_ref_alleles_{}",
                                    ">{}_strain_{}",
                                    target_names[ref_index][&tid],
                                    strain_index,
                                    // variations,
                                    // ref_alleles
                                )
                                    .expect("Unable to write to file");

                                for line in contig
                                    .as_bytes()
                                    .to_vec()[..]
                                    .chunks(60)
                                    .into_iter() {
                                    file_open.write(line).unwrap();
                                    file_open.write(b"\n").unwrap();
                                }
                                total_variant_alleles += variations;
                                total_reference_alleles += ref_alleles;
                            }
                            if overclustered {
                                debug!(
                                    "Genome {}: Cluster {} contains {} multi-variant sites. It may be overclustered.",
                                    &genomes_and_contigs.genomes[*ref_index],
                                    strain_index,
                                    multivariant_sites,
                                )
                            }
                            debug!(
                                "Genome {}: Cluster {} contains {} variant alleles and {} reference alleles",
                                &genomes_and_contigs.genomes[*ref_index],
                                strain_index,
                                total_variant_alleles,
                                total_reference_alleles,
                            );
                        });
                    });
            }
        }
    }

    fn calculate_strain_abundances(
        &mut self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                pred_variants,
                all_variants,
                sample_names,
                ..
            } => {
                let number_of_samples = sample_names.len();
                all_variants
                    .par_iter()
                    .for_each(|(ref_index, ref_variants)| {
                        match pred_variants.get(&ref_index) {
                            Some(genotype_map) => {
                                debug!(
                                    "Calculating genotype abundances for reference: {}",
                                    reference_map[&ref_index],
                                );
                                let number_of_genotypes = genotype_map.keys().len();

                                // The initialization vector for the EM algorithm
                                let mut genotype_vectors = vec![
                                    vec![
                                            genotype_abundances::Genotype::new(0);
                                            number_of_genotypes
                                        ];
                                    number_of_samples
                                ];

                                // A key tracking the genotype_index and strain_id values
                                let mut genotype_key: HashMap<usize, usize> = HashMap::new();
                                for (genotype_idx, (strain_id, _)) in
                                    genotype_map.iter().enumerate()
                                {
                                    genotype_key.insert(*strain_id, genotype_idx);
                                }

                                debug!(
                                    "Populating sample_vec with genotypes... {:?}",
                                    &genotype_key,
                                );
                                debug!(
                                    "Number of genotypes {} and genotype vectors {:?}",
                                    number_of_genotypes, genotype_vectors
                                );
                                genotype_vectors.par_iter_mut().for_each(|sample_vec| {
                                    for (strain_id, _) in genotype_map.iter() {
                                        let genotype_idx = genotype_key.get(strain_id).unwrap();
                                        sample_vec[*genotype_idx] =
                                            genotype_abundances::Genotype::new(*strain_id);
                                    }
                                });

                                debug!("Collecting variants...");
                                genotype_vectors.par_iter_mut().enumerate().for_each(
                                    |(sample_index, sample_vec)| {
                                        for (tid, contig_variants) in ref_variants.iter() {
                                            for (pos, pos_variants) in contig_variants.iter() {
                                                if pos_variants.len() > 1 {
                                                    for (variant, base_info) in pos_variants.iter()
                                                    {
                                                        match variant {
                                                            Variant::SNV(_)
                                                            | Variant::MNV(_)
                                                            | Variant::SV(_)
                                                            | Variant::Insertion(_)
                                                            | Variant::Inversion(_)
                                                            | Variant::Deletion(_) => {
                                                                // We divide the total depth of variant here
                                                                // by the total amount of genotypes that
                                                                // variant occurs in.
                                                                // E.g. if a variant had a depth of 6
                                                                // and occurred in 3 genotypes, then for each
                                                                // genotype its initialization value would be 2
                                                                let lorikeet_depth = base_info
                                                                    .truedepth[sample_index]
                                                                    as f64;
                                                                // let relative_abundance = lorikeet_depth / base_info.totaldepth[sample_index] as f64;
                                                                let weight = lorikeet_depth
                                                                    / base_info.genotypes.len()
                                                                        as f64;
                                                                for strain_id in
                                                                    base_info.genotypes.iter()
                                                                {
                                                                    let genotype_idx = genotype_key
                                                                        .get(&(*strain_id as usize))
                                                                        .unwrap();
                                                                    sample_vec[*genotype_idx]
                                                                        .variant_weights
                                                                        .push(weight);
                                                                    let genotype_indices =
                                                                        base_info
                                                                            .genotypes
                                                                            .par_iter()
                                                                            .map(|idx| {
                                                                                *genotype_key
                                                                            .get(&(*idx as usize))
                                                                            .unwrap()
                                                                            })
                                                                            .collect::<Vec<usize>>(
                                                                            );
                                                                    sample_vec[*genotype_idx]
                                                                        .variant_genotype_ids
                                                                        .push(genotype_indices);
                                                                }
                                                            }
                                                            Variant::None => {
                                                                // debug!("Reference variant {:?}", base_info.truedepth);
                                                                let lorikeet_depth = base_info
                                                                    .truedepth[sample_index]
                                                                    as f64;
                                                                // let relative_abundance = lorikeet_depth / base_info.totaldepth[sample_index] as f64;
                                                                let weight = lorikeet_depth
                                                                    / base_info.genotypes.len()
                                                                        as f64;
                                                                for strain_id in
                                                                    base_info.genotypes.iter()
                                                                {
                                                                    let genotype_idx = genotype_key
                                                                        .get(&(*strain_id as usize))
                                                                        .unwrap();
                                                                    sample_vec[*genotype_idx]
                                                                        .variant_weights
                                                                        .push(weight);
                                                                    let genotype_indices =
                                                                        base_info
                                                                            .genotypes
                                                                            .par_iter()
                                                                            .map(|idx| {
                                                                                *genotype_key
                                                                            .get(&(*idx as usize))
                                                                            .unwrap()
                                                                            })
                                                                            .collect::<Vec<usize>>(
                                                                            );
                                                                    sample_vec[*genotype_idx]
                                                                        .variant_genotype_ids
                                                                        .push(genotype_indices);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    },
                                );

                                debug!("Calculating abundances...");
                                // Set up multi progress bars
                                // let multi = MultiProgress::new();
                                // let sty = ProgressStyle::default_bar()
                                //     .template("[{elapsed_precise}] {bar:40.green/blue} {pos:>7}/{len:7} {msg}")
                                //     .progress_chars("##-");
                                //
                                //
                                // let pb1 = multi.insert(0, ProgressBar::new(genotype_vectors.len() as u64));
                                // pb1.set_style(sty.clone());
                                // pb1.set_message("Calculating genotype abundances...");
                                //
                                // let _ = std::thread::spawn(move || {
                                //     multi.join_and_clear().unwrap();
                                // });
                                genotype_vectors.par_iter_mut().enumerate().for_each(
                                    |(idx, sample_genotypes)| {
                                        debug!(
                                            "Genotype Vector before EM {} {:?}",
                                            idx, sample_genotypes
                                        );
                                        genotype_abundances::calculate_abundances(sample_genotypes);
                                        debug!(
                                            "Genotype Vector after EM {} {:?}",
                                            idx, sample_genotypes
                                        );
                                        // pb1.inc(1);
                                    },
                                );
                                // pb1.finish_with_message("Done!");

                                let reference_stem = &genomes_and_contigs.genomes[*ref_index];
                                debug!("Printing strain coverages {}", &reference_stem);
                                let file_name = format!(
                                    "{}/{}_strain_coverages.tsv",
                                    &output_prefix, &reference_stem,
                                );

                                let file_path = Path::new(&file_name);

                                let mut file_open = match File::create(file_path) {
                                    Ok(fasta) => fasta,
                                    Err(e) => {
                                        println!("Cannot create file {:?}", e);
                                        std::process::exit(1)
                                    }
                                };

                                // Snp density summary start

                                write!(file_open, "Sample").unwrap();
                                let mut printing_order = vec![];
                                for (strain_id, _) in genotype_key.iter() {
                                    write!(file_open, "\tstrain_{}", strain_id).unwrap();
                                    printing_order.push(*strain_id);
                                }

                                write!(file_open, "\n").unwrap();

                                for (sample_idx, genotype) in genotype_vectors.iter().enumerate() {
                                    let mut sample_name = &sample_names[sample_idx];
                                    // remove tmp file name from sample id
                                    let sample_name = match &sample_name[..4] {
                                        ".tmp" => &sample_name[15..],
                                        _ => &sample_name,
                                    };
                                    write!(file_open, "{} Expected Coverage", &sample_name,)
                                        .unwrap();
                                    for strain_id in printing_order.iter() {
                                        let genotype_idx = genotype_key.get(strain_id).unwrap();
                                        write!(
                                            file_open,
                                            "\t{}",
                                            genotype[*genotype_idx].abundance_weight
                                        )
                                        .unwrap();
                                    }
                                    write!(file_open, "\n").unwrap();
                                }
                            }
                            None => {
                                // info!("No genotypes created for: {}", reference_map[&ref_index],);
                            }
                        }
                    });
            }
        }
    }

    fn print_variant_stats(
        &self,
        window_size: f64,
        output_prefix: &str,
        genomes_and_contigs: &GenomesAndContigs,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                sample_names,
                all_variants,
                variant_sums,
                variant_info,
                ..
            } => {
                for (ref_idx, ref_variants) in all_variants.iter() {
                    let reference_stem = &genomes_and_contigs.genomes[*ref_idx];
                    debug!("Profiling {}", &reference_stem);
                    let file_name = format!("{}/{}_summary.tsv", &output_prefix, &reference_stem,);
                    let snp_locs =
                        format!("{}/{}_snp_locations.tsv", &output_prefix, &reference_stem,);

                    let file_path = Path::new(&file_name);

                    let mut file_open = match File::create(file_path) {
                        Ok(fasta) => fasta,
                        Err(e) => {
                            println!("Cannot create file {:?}", e);
                            std::process::exit(1)
                        }
                    };

                    let snp_loc_path = Path::new(&snp_locs);
                    let mut snp_loc_open = match File::create(snp_loc_path) {
                        Ok(tsv) => tsv,
                        Err(e) => {
                            println!("Cannot create file {:?}", e);
                            std::process::exit(1)
                        }
                    };

                    // Snp density summary start
                    write!(file_open, "contigName\tcontigLen").unwrap();

                    // Snp location start
                    write!(snp_loc_open, "SNP\tchr\tpos").unwrap();
                    debug!("Sample Names {:?}", sample_names);
                    for sample_name in sample_names.iter() {
                        write!(
                            file_open,
                            "\t{}.snvsPer{}kb\t{}.svsPer{}kb\t{}.snvCount\t{}.svCount",
                            &sample_name,
                            &window_size,
                            &sample_name,
                            &window_size,
                            &sample_name,
                            &sample_name
                        )
                        .unwrap();
                    }
                    write!(file_open, "\n").unwrap();
                    write!(snp_loc_open, "\n").unwrap();

                    let mut snps = 0;
                    debug!(
                        "Reference index {} target names {:?}",
                        &ref_idx, &target_names
                    );

                    match target_names.get(ref_idx) {
                        Some(target_name_set) => {
                            for (tid, contig_name) in target_name_set.iter() {
                                let contig_len = target_lengths[ref_idx][&tid];
                                write!(file_open, "{}\t{}", contig_name, contig_len).unwrap();
                                match ref_variants.get(tid) {
                                    Some(variants_in_contig) => {
                                        // Set up channels that receive a vector of values
                                        // for each sample
                                        let mut snps_cnt_vec = vec![0; sample_names.len()];
                                        let svs_cnt_vec =
                                            Arc::new(Mutex::new(vec![0; sample_names.len()]));

                                        let window = contig_len / window_size;
                                        variants_in_contig.iter().enumerate()
                                            .for_each(
                                                |(index,
                                                     (position, variants))| {
                                                // Get how many alleles are present at loci
                                                let alleles = variants.len();
                                                for (var, base) in variants {
                                                    match var {
                                                        Variant::SNV(_) => {
                                                            write!(snp_loc_open, "SNP{}\t{}\t{}\n",
                                                                   index, contig_name, position).unwrap();
                                                            base.truedepth
                                                                .iter()
                                                                .enumerate()
                                                                .zip(base.depth.iter().enumerate())
                                                                .for_each(|((index, count_1), (_index_2, count_2))| {
                                                                    if count_1 > &0 || count_2 > &0 {
                                                                        snps_cnt_vec[index] += 1;
                                                                    }
                                                                });
                                                            snps += 1;
                                                        },
                                                        Variant::None => {
                                                            // If biallelic or multiallelic then
                                                            // I don't think we want the ref depth?
                                                        },
                                                        _ => {
                                                            base.truedepth
                                                                .par_iter()
                                                                .enumerate()
                                                                .zip(base.depth.par_iter().enumerate())
                                                                .for_each(|((index, count_1), (_, count_2))| {
                                                                    if count_1 > &0 || count_2 > &0 {
                                                                        let mut svs_cnt_vec
                                                                            = svs_cnt_vec.lock().unwrap();
                                                                        svs_cnt_vec[index] += 1;
                                                                    }
                                                                })
                                                        }
                                                    };
                                                };
                                            });
                                        let svs_cnt_vec = svs_cnt_vec.lock().unwrap().clone();
                                        let snps_per_win: Vec<_> = snps_cnt_vec
                                            .iter()
                                            .map(|count| *count as f64 / window)
                                            .collect();
                                        let svs_per_win: Vec<_> = svs_cnt_vec
                                            .iter()
                                            .map(|count| *count as f64 / window)
                                            .collect();

                                        for (snp_w, svs_w, snp_c, svs_c) in izip!(
                                            &snps_per_win,
                                            &svs_per_win,
                                            &snps_cnt_vec,
                                            &svs_cnt_vec
                                        ) {
                                            write!(
                                                file_open,
                                                "\t{}\t{}\t{}\t{}",
                                                snp_w, svs_w, snp_c, svs_c
                                            )
                                            .unwrap();
                                        }
                                        write!(file_open, "\n").unwrap();
                                    }
                                    None => {
                                        // Write out zeros for contigs with no variants
                                        for (sample_idx, _sample_name) in
                                            sample_names.iter().enumerate()
                                        {
                                            write!(
                                                file_open,
                                                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                                0., 0., 0., 0., 0., 0., 0.,
                                            )
                                            .unwrap();
                                        }
                                    }
                                }
                            }
                        }
                        None => {}
                    }

                    if snps > 1 {
                        debug!("Creating SNP density plot for {}...", &reference_stem);
                        let plot_command = format!(
                            "set -eou pipefail; snp_density_plots.R {} {} && \
                                     mv SNP-Density*.pdf {}_snp_density_plot.pdf",
                            format!("{}/{}", &output_prefix, &reference_stem,),
                            window_size,
                            format!("{}/{}", &output_prefix, &reference_stem,),
                        );
                        command::finish_command_safely(
                            std::process::Command::new("bash")
                                .arg("-c")
                                .arg(&plot_command)
                                .stderr(Stdio::piped())
                                .stdout(Stdio::piped())
                                .spawn()
                                .expect("Unable to execute Rscript"),
                            "CMplot",
                        );
                    } else {
                        std::fs::remove_file(snp_loc_path).expect("Unable to remove file");
                    }
                }
            }
        }
    }

    fn calc_gene_mutation(
        &self,
        gff_map: &mut HashMap<usize, HashMap<String, Vec<bio::io::gff::Record>>>,
        genomes_and_contigs: &GenomesAndContigs,
        reference_map: &HashMap<usize, String>,
        codon_table: &CodonTable,
        output_prefix: &str,
        concatenated_genomes: &Option<String>,
    ) {
        // Calculates parse the per reference information to functions in codon_structs
        // Writes the resulting gff records
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                all_variants,
                ..
            } => {
                for (ref_idx, ref_variants) in all_variants.iter() {
                    // Retrieve path to reference
                    let reference_path = Path::new(
                        reference_map
                            .get(ref_idx)
                            .expect("reference index not found"),
                    );
                    let mut reference = match concatenated_genomes {
                        Some(temp_file) => {
                            match bio::io::fasta::IndexedReader::from_file(&temp_file) {
                                Ok(reader) => reader,
                                Err(_e) => generate_faidx(&temp_file),
                            }
                        }
                        None => match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                            Ok(reader) => reader,
                            Err(_e) => generate_faidx(&reference_path.to_str().unwrap()),
                        },
                    };

                    // Reference specific GFF records
                    let mut gff_ref = gff_map.get_mut(&ref_idx).expect(&format!(
                        "No GFF records for reference {:?}",
                        reference_path
                    ));

                    let file_name = format!(
                        "{}/{}_dnds_values.gff",
                        output_prefix.to_string(),
                        reference_path.file_stem().unwrap().to_str().unwrap(),
                    );

                    let file_path = Path::new(&file_name);

                    let mut gff_writer =
                        gff::Writer::to_file(file_path, bio::io::gff::GffType::GFF3)
                            .expect("unable to create GFF file");
                    for (tid, contig_name) in target_names[ref_idx].iter() {
                        let mut ref_sequence = Vec::new();
                        match reference
                            .fetch_all(std::str::from_utf8(contig_name.as_bytes()).unwrap())
                        {
                            Ok(reference) => reference,
                            Err(e) => {
                                println!("Cannot read sequence from reference {:?}", e);
                                std::process::exit(1)
                            }
                        };
                        match reference.read(&mut ref_sequence) {
                            Ok(reference) => reference,
                            Err(e) => {
                                println!("Cannot read sequence from reference {:?}", e);
                                std::process::exit(1)
                            }
                        };
                        let mut placeholder = Vec::new();
                        let gff_records = match gff_ref.get_mut(contig_name) {
                            Some(records) => records,
                            None => &mut placeholder,
                        };
                        debug!(
                            "Calculating population dN/dS from reads for {} genes on contig {}",
                            gff_records.len(),
                            contig_name
                        );

                        let placeholder_map = HashMap::new();
                        let variants = match ref_variants.get(tid) {
                            Some(map) => map,
                            None => &placeholder_map,
                        };

                        // Set up multi progress bars
                        let multi = MultiProgress::new();
                        let sty = ProgressStyle::default_bar()
                            .template(
                                "[{elapsed_precise}] {bar:40.green/blue} {pos:>7}/{len:7} {msg}",
                            )
                            .progress_chars("##-");

                        let pb1 = multi.insert(0, ProgressBar::new(gff_records.len() as u64));
                        pb1.set_style(sty.clone());
                        pb1.set_message("Calculating dN/dS for each ORF...");

                        let _ = std::thread::spawn(move || {
                            multi.join_and_clear().unwrap();
                        });

                        gff_records.iter_mut().enumerate().for_each(|(_id, gene)| {
                            let dnds = codon_table.find_mutations(gene, variants, &ref_sequence);
                            gene.attributes_mut()
                                .insert(format!("dNdS"), format!("{}", dnds));
                            debug!("gene {:?} attributes {:?}", gene.seqname(), gene);
                            gff_writer.write(gene);
                            pb1.inc(1);
                        });
                    }
                }
            }
        }
    }

    fn polish_genomes(
        &self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
        concatenated_genomes: &Option<String>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                sample_names,
                all_variants,
                variant_sums,
                variant_info,
                ..
            } => {
                all_variants
                    .par_iter()
                    .for_each(|(ref_index, ref_variants)| {


                        let reference_path =
                            Path::new(reference_map
                                .get(ref_index)
                                .expect("reference index not found"));


                        let mut reference = match concatenated_genomes {
                            Some(temp_file) => {
                                match bio::io::fasta::IndexedReader::from_file(&temp_file) {
                                    Ok(reader) => reader,
                                    Err(_e) => generate_faidx(&temp_file),
                                }
                            },
                            None => {
                                match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                                    Ok(reader) => reader,
                                    Err(_e) => generate_faidx(&reference_path.to_str().unwrap()),
                                }
                            }
                        };

                        debug!("Reference index {} target names {:?}", &ref_index, &target_names);

                        match target_names.get(ref_index) {
                            Some(target_name_set) => {
                                sample_names.iter().enumerate().for_each(|(sample_idx, sample_name)| {

                                    // remove tmp file name from sample id
                                    let sample_name = match &sample_name[..4] {
                                        ".tmp" => &sample_name[15..],
                                        _ => &sample_name,
                                    };

                                    let file_name =
                                        format!(
                                            "{}/{}_polished_{}.fna",
                                            output_prefix.to_string(),
                                            reference_path.file_stem().unwrap().to_str().unwrap(),
                                            &sample_name,
                                        );

                                    let file_path = Path::new(&file_name);
                                    // Open haplotype file or create one
                                    let mut file_open = File::create(file_path)
                                        .expect("No Read or Write Permission in current directory");

                                    let mut original_contig = Vec::new();

                                    let mut total_variant_alleles = 0;
                                    let mut total_reference_alleles = 0;

                                    for (tid, target_name) in target_name_set.iter() {
                                        // let contig_len = target_lengths[ref_idx][&tid];
                                        let mut contig = String::new();
                                        original_contig = Vec::new();
                                        {
                                            match reference.fetch_all(
                                                std::str::from_utf8(target_name.as_bytes()).unwrap()) {
                                                Ok(reference) => reference,
                                                Err(e) => {
                                                    warn!("Cannot read sequence from reference {:?}", e);
                                                    std::process::exit(1)
                                                },
                                            };
                                            match reference.read(&mut original_contig) {
                                                Ok(reference) => reference,
                                                Err(e) => {
                                                    warn!("Cannot read sequence from reference {:?}", e);
                                                    std::process::exit(1)
                                                },
                                            };
                                        }
                                        let mut skip_n = 0;
                                        let mut skip_cnt = 0;
                                        let mut variations = 0;
                                        let mut ref_alleles = 0;

                                        match ref_variants.get(tid) {
                                            Some(variants_in_contig) => {

                                                for (pos, base) in original_contig.iter().enumerate() {
                                                    if skip_cnt < skip_n {
                                                        skip_cnt += 1;
                                                    } else {
                                                        skip_n = 0;
                                                        skip_cnt = 0;

                                                        if variants_in_contig.contains_key(&(pos as i64)) {

                                                            let hash = variants_in_contig
                                                                .get(&(pos as i64))
                                                                .unwrap();

                                                            let mut max_var = Variant::None;
                                                            let mut max_depth = 0;

                                                            for (var, base) in hash.iter() {
                                                                // If there are two variants possible for
                                                                // a single site and one is the reference
                                                                // we will choose the reference
                                                                if max_depth == 0 {
                                                                    max_var = var.clone();
                                                                    max_depth = base.truedepth[sample_idx]
                                                                } else if base.truedepth[sample_idx] > max_depth {
                                                                    max_var = var.clone();
                                                                    max_depth = base.truedepth[sample_idx]
                                                                }
                                                            }
                                                            if hash.len() > 1 {
                                                                //                                            multivariant_sites += 1;
                                                                debug!("Multi hash {:?} {:?}", hash, max_var)
                                                            }
                                                            match max_var {
                                                                Variant::Deletion(size) => {
                                                                    // Skip the next n bases but rescue the reference prefix
                                                                    skip_n = size;
                                                                    skip_cnt = 0;
                                                                    contig = contig + str::from_utf8(&[*base]).unwrap();
                                                                    // If we had the sequence we would rescue first base like this
                                                                    //                                                let first_byte = max_var.as_bytes()[0];
                                                                    //                                                contig = contig + str::from_utf8(
                                                                    //                                                    &[first_byte]).unwrap();
                                                                    variations += 1;
                                                                },
                                                                Variant::Insertion(alt) => {

                                                                    // Remove prefix from variant
                                                                    let removed_first_base = str::from_utf8(
                                                                        &alt[1..]).unwrap();
                                                                    contig = contig + removed_first_base;
                                                                    variations += 1;
                                                                },
                                                                Variant::Inversion(alt) | Variant::MNV(alt) => {
                                                                    // Skip the next n bases
                                                                    skip_n = alt.len() as u32 - 1;
                                                                    skip_cnt = 0;
                                                                    // Inversions and MNVs don't have a first base prefix, so take
                                                                    // wholes tring
                                                                    let inversion = str::from_utf8(
                                                                        &alt).unwrap();
                                                                    contig = contig + inversion;
                                                                    variations += 1;
                                                                },
                                                                Variant::None => {
                                                                    contig = contig + str::from_utf8(&[*base]).unwrap();
                                                                    ref_alleles += 1;
                                                                },
                                                                Variant::SNV(alt) => {
                                                                    contig = contig + str::from_utf8(&[alt]).unwrap();
                                                                    variations += 1;
                                                                },
                                                                _ => {
                                                                    contig = contig + str::from_utf8(&[*base]).unwrap();
                                                                    ref_alleles += 1;
                                                                }
                                                            }
                                                        } else {
                                                            contig = contig + str::from_utf8(&[*base]).unwrap();
                                                        }

                                                    }
                                                };
                                            },
                                            None => {
                                                // Write out zeros for contigs with no variants
                                                contig = str::from_utf8(&original_contig)
                                                    .expect("Can't convert to str").to_string();
                                            }
                                        };

                                        writeln!(
                                            file_open,
                                            ">{}_polished_{}_alt_alleles_{}_ref_alleles_{}",
                                            target_names[ref_index][&tid],
                                            sample_name,
                                            variations,
                                            ref_alleles
                                        ).expect("Unable to write to file");

                                        for line in contig.as_bytes().to_vec()[..].chunks(60).into_iter() {
                                            file_open.write(line).unwrap();
                                            file_open.write(b"\n").unwrap();
                                        };
                                    };
                                });
                            },
                            None => {},
                        }
                });
            }
        }
    }

    fn calculate_sample_distances(
        &self,
        output_prefix: &str,
        reference_map: &HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                pred_variants,
                all_variants,
                sample_names,
                ..
            } => {
                let number_of_samples = sample_names.len();
                all_variants
                    .par_iter()
                    .for_each(|(ref_index, ref_variants)| {
                        if ref_variants.len() > 0 {
                            debug!(
                                "Generating variant adjacency matrix for {}",
                                reference_map[&ref_index],
                            );

                            // The initialization of adjacency matrix n*n n=samples
                            let adjacency_matrix =
                                Arc::new(Mutex::new(vec![
                                    vec![0; number_of_samples];
                                    number_of_samples
                                ]));

                            let sample_totals = Arc::new(Mutex::new(vec![0; number_of_samples]));

                            debug!("Collecting variants...");

                            for (tid, contig_variants) in ref_variants.iter() {
                                for (pos, pos_variants) in contig_variants.iter() {
                                    if pos_variants.len() > 1 {
                                        for (variant, base_info) in pos_variants.iter() {
                                            match variant {
                                                Variant::SNV(_)
                                                | Variant::MNV(_)
                                                | Variant::SV(_)
                                                | Variant::Insertion(_)
                                                | Variant::Inversion(_)
                                                | Variant::Deletion(_) => {
                                                    // Update adjacency matrix for each variant
                                                    // Evaluate each sample pairwise
                                                    (0..number_of_samples)
                                                        .into_iter()
                                                        .combinations(2)
                                                        .collect::<Vec<Vec<usize>>>()
                                                        .into_par_iter()
                                                        .for_each(|indices| {
                                                            let d1 =
                                                                base_info.truedepth[indices[0]];
                                                            let d2 =
                                                                base_info.truedepth[indices[1]];

                                                            if d1 > 0 && d2 > 0 {
                                                                // udpate both
                                                                let mut adjacency_matrix =
                                                                    adjacency_matrix
                                                                        .lock()
                                                                        .unwrap();
                                                                let mut sample_totals =
                                                                    sample_totals.lock().unwrap();
                                                                sample_totals[indices[0]] += 1;
                                                                sample_totals[indices[1]] += 1;

                                                                adjacency_matrix[indices[0]]
                                                                    [indices[1]] += 1;
                                                                adjacency_matrix[indices[1]]
                                                                    [indices[0]] += 1;
                                                            } else if d1 > 0 && d2 <= 0 {
                                                                // update only one
                                                                let mut sample_totals =
                                                                    sample_totals.lock().unwrap();
                                                                sample_totals[indices[0]] += 1;
                                                            } else if d1 <= 0 && d2 > 0 {
                                                                let mut sample_totals =
                                                                    sample_totals.lock().unwrap();
                                                                sample_totals[indices[1]] += 1;
                                                            }
                                                        });
                                                }
                                                Variant::None => {
                                                    // do nothing
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            let adjacency_matrix = adjacency_matrix.lock().unwrap();
                            let sample_totals = sample_totals.lock().unwrap();

                            let reference_stem = &genomes_and_contigs.genomes[*ref_index];
                            debug!("Printing strain coverages {}", &reference_stem);
                            let file_name = format!(
                                "{}/{}_adjacency_matrix.tsv",
                                &output_prefix, &reference_stem,
                            );

                            let file_path = Path::new(&file_name);

                            let mut file_open = match File::create(file_path) {
                                Ok(fasta) => fasta,
                                Err(e) => {
                                    println!("Cannot create file {:?}", e);
                                    std::process::exit(1)
                                }
                            };

                            // Adjacency summary start
                            write!(file_open, "{: <20}", "").unwrap();
                            for sample_name in sample_names.iter() {
                                // remove tmp file name from sample id
                                let sample_name = match &sample_name[..4] {
                                    ".tmp" => &sample_name[15..],
                                    _ => &sample_name,
                                };
                                write!(file_open, "\t{: <20}", sample_name).unwrap();
                            }

                            write!(file_open, "\n").unwrap();

                            for (sample_idx_1, distance_vec) in adjacency_matrix.iter().enumerate()
                            {
                                let mut sample_name = &sample_names[sample_idx_1];
                                // remove tmp file name from sample id
                                let sample_name = match &sample_name[..4] {
                                    ".tmp" => &sample_name[15..],
                                    _ => &sample_name,
                                };
                                write!(file_open, "{: <20}", &sample_name,).unwrap();
                                for (sample_idx_2, count) in distance_vec.iter().enumerate() {
                                    // Return the jaccard's similarity between sets of variants in
                                    // samples
                                    // let jaccards_sim = if sample_idx_1 == sample_idx_2 {
                                    //     1.
                                    // } else if sample_totals[sample_idx_1] > 0
                                    //     || sample_totals[sample_idx_2] > 0
                                    // {
                                    //     *count as f32
                                    //         / (sample_totals[sample_idx_1] as f32
                                    //             + sample_totals[sample_idx_2] as f32
                                    //             - *count as f32)
                                    // } else {
                                    //     0.
                                    // };

                                    // Return shared variant count
                                    let count = if sample_idx_1 == sample_idx_2 {
                                        sample_totals[sample_idx_1] as f32
                                            / (number_of_samples - 1) as f32
                                    } else {
                                        *count as f32
                                    };

                                    write!(file_open, "\t{: <20}", count,).unwrap();
                                }
                                write!(file_open, "\n").unwrap();
                            }
                        }
                    });
            }
        }
    }

    fn write_vcf(&self, output_prefix: &str, genomes_and_contigs: &GenomesAndContigs) {
        match self {
            VariantMatrix::VariantContigMatrix {
                all_variants,
                target_names,
                target_lengths,
                sample_names,
                variant_info,
                ..
            } => {
                for (ref_idx, variant_info_ref) in variant_info.iter() {
                    if variant_info_ref.len() > 0 {
                        // initiate header
                        let mut header = bcf::Header::new();
                        // Add program info
                        header.push_record(
                            format!("##source=lorikeet-v{}", env!("CARGO_PKG_VERSION")).as_bytes(),
                        );

                        debug!("samples {:?}", &sample_names);
                        for sample_name in sample_names.iter() {
                            // remove tmp file name from sample id
                            let sample_name = match &sample_name[..4] {
                                ".tmp" => &sample_name[15..],
                                _ => &sample_name,
                            };
                            header.push_sample(&sample_name.to_string().into_bytes()[..]);
                        }

                        // Add contig info
                        for (tid, contig_name) in target_names[ref_idx].iter() {
                            header.push_record(
                                format!(
                                    "##contig=<ID={}, length={}>",
                                    contig_name, target_lengths[ref_idx][&tid]
                                )
                                .as_bytes(),
                            );
                        }

                        // Add INFO flags
                        header.push_record(
                            format!(
                                "##INFO=<ID=TYPE,Number=A,Type=String,\
                    Description=\"The type of allele, either SNV, MNV, INS, DEL, or INV.\">"
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            format!(
                                "##INFO=<ID=SVLEN,Number=1,Type=Integer,\
                    Description=\"Length of structural variant\">"
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            b"##INFO=<ID=TDP,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=TAD,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth of alternative allele\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=TRD,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth of reference allele\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=ST,Number=.,Type=Integer,\
                            Description=\"The strain IDs assigned to this variant\">",
                        );

                        // Add FORMAT flags
                        header.push_record(
                            format!(
                                "##FORMAT=<ID=DP,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth in each sample\">",
                                sample_names.len()
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            format!("##FORMAT=<ID=AD,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth of alternative allele in each sample\">",
                                    sample_names.len()).as_bytes(),
                        );

                        header.push_record(
                            format!("##FORMAT=<ID=RD,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth of reference allele in each sample\">",
                                    sample_names.len()).as_bytes(),
                        );

                        header.push_record(
                            format!(
                                "##FORMAT=<ID=QA,Number={},Type=Float,\
                            Description=\"Quality scores for allele in each sample\">",
                                sample_names.len()
                            )
                            .as_bytes(),
                        );

                        let vcf_presort = tempfile::Builder::new()
                            .prefix("lorikeet-fifo")
                            .tempfile()
                            .expect("Unable to create VCF tempfile");

                        // Initiate writer
                        let mut bcf_writer = bcf::Writer::from_path(
                            format!(
                                "{}/{}.vcf",
                                &output_prefix, &genomes_and_contigs.genomes[*ref_idx],
                            )
                            .as_str(),
                            &header,
                            true,
                            bcf::Format::VCF,
                        )
                        .expect(
                            format!("Unable to create VCF output: {}.vcf", output_prefix).as_str(),
                        );

                        bcf_writer.set_threads(current_num_threads()).unwrap();

                        for (tid, position_variants) in all_variants[ref_idx].iter() {
                            let contig_name = &target_names[ref_idx][&tid];
                            for (pos, variants) in position_variants.iter() {
                                for (variant, base) in variants.iter() {
                                    // Create empty record for this variant
                                    let mut record = bcf_writer.empty_record();
                                    record.set_rid(Some(
                                        bcf_writer
                                            .header()
                                            .name2rid(contig_name.as_bytes())
                                            .unwrap(),
                                    ));
                                    record.set_pos(*pos);

                                    // Sum the quality scores across samples
                                    let qual_sum = base.quals.par_iter().sum::<f64>();
                                    record.set_qual(qual_sum as f32);

                                    // Collect strain information
                                    let mut strains =
                                        base.genotypes.iter().cloned().collect::<Vec<i32>>();
                                    strains.sort();

                                    // Push info tags to record
                                    record
                                        .push_info_integer(b"TDP", &[base.totaldepth.iter().sum()]);
                                    record
                                        .push_info_integer(b"TAD", &[base.truedepth.iter().sum()]);
                                    record.push_info_integer(
                                        b"TRD",
                                        &[base.referencedepth.iter().sum()],
                                    );
                                    record.push_info_integer(b"ST", &strains[..]);

                                    // Push format flags to record
                                    record.push_format_integer(b"DP", &base.totaldepth[..]);
                                    record.push_format_integer(b"AD", &base.truedepth[..]);
                                    record.push_format_integer(b"RD", &base.referencedepth[..]);
                                    record.push_format_float(
                                        b"QA",
                                        &base.quals.iter().map(|p| *p as f32).collect_vec()[..],
                                    );

                                    let refr = &base.refr[..];

                                    match variant {
                                        Variant::SNV(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = &[*alt];
                                            let mut collect_alleles = vec![refr, alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"SNV"]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::MNV(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"MNV"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Inversion(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"INV"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Insertion(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"INS"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Deletion(alt) => {
                                            // Collect and set the alleles to record
                                            // Create deletion variant, attach refr to head of N
                                            record
                                                .set_alleles(&vec![
                                                    refr,
                                                    format!("{}N", refr[0] as char).as_bytes(),
                                                ])
                                                .unwrap();

                                            record.push_info_string(b"TYPE", &[b"DEL"]);
                                            record.push_info_integer(b"SVLEN", &[*alt as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        _ => {}
                                    }
                                }
                            }
                        }
                        debug!("Finished writing VCF file for {}", &output_prefix);
                    }
                }
            }
        }
    }
}

/// Add read count entry to cluster hashmap
pub fn add_entry(
    shared_read_counts: &mut HashMap<usize, HashMap<usize, usize>>,
    clust1: usize,
    clust2: usize,
    count: usize,
) {
    let entry = shared_read_counts.entry(clust1).or_insert(HashMap::new());
    entry.insert(clust2, count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::variants;
    use std::collections::HashSet;

    fn create_base(ref_sequence: &Vec<u8>, var_char: u8, pos: i64, sample_count: usize) -> Base {
        Base {
            tid: 0,
            pos,
            refr: ref_sequence[pos as usize..(pos as usize + 1)].to_vec(),
            variant: Variant::SNV(var_char),
            depth: vec![0, 5],
            truedepth: vec![0, 5],
            totaldepth: vec![5, 5],
            genotypes: HashSet::new(),
            quals: vec![0.; sample_count],
            referencedepth: vec![5, 0],
            freq: vec![0.; sample_count],
            rel_abunds: vec![0.; sample_count],
            reads: HashSet::new(),
        }
    }

    #[test]
    fn test_stats_and_mat() {
        let ref_sequence = "ATGAAACCCGGGTTTTAA".as_bytes().to_vec();
        let sample_count = 2;

        let mut variant_abundances: HashMap<i64, HashMap<Variant, Base>> = HashMap::new();

        // Setup variant matrix
        let mut var_mat = VariantMatrix::new_matrix(2);

        // Create fake variants
        let var_1 = create_base(&ref_sequence, "G".bytes().nth(0).unwrap(), 7, 2);
        let var_2 = create_base(&ref_sequence, "C".bytes().nth(0).unwrap(), 11, 2);
        let var_3 = create_base(&ref_sequence, "A".bytes().nth(0).unwrap(), 13, 2);
        let var_4 = create_base(&ref_sequence, "C".bytes().nth(0).unwrap(), 14, 2);

        let mut ups_and_downs = vec![0; ref_sequence.len()];
        ups_and_downs[0] = 5;
        ups_and_downs[ref_sequence.len() - 1] = -5;

        // Setup sample 1
        let mut var_stats = VariantStats::new_contig_stats(0., 5., 0);

        var_stats.add_contig(
            Some(&mut variant_abundances),
            0,
            0,
            b"test".to_vec(),
            ref_sequence.len(),
            0,
            vec![10., 10., 0.],
            ups_and_downs,
        );
        var_mat.add_contig(var_stats, 2, 0, 0);

        {
            // Add variants in
            let hash = variant_abundances.entry(7).or_insert(HashMap::new());
            hash.insert(var_1.variant.clone(), var_1);

            let hash = variant_abundances.entry(11).or_insert(HashMap::new());
            hash.insert(var_2.variant.clone(), var_2);

            let hash = variant_abundances.entry(13).or_insert(HashMap::new());
            hash.insert(var_3.variant.clone(), var_3);

            let hash = variant_abundances.entry(14).or_insert(HashMap::new());
            hash.insert(var_4.variant.clone(), var_4);
        }
        // Add sample 2
        let mut ups_and_downs = vec![0; ref_sequence.len()];
        ups_and_downs[0] = 5;
        ups_and_downs[ref_sequence.len() - 1] = -5;
        let mut var_stats = VariantStats::new_contig_stats(0., 5., 0);
        var_stats.add_contig(
            Some(&mut variant_abundances),
            0,
            0,
            b"test".to_vec(),
            ref_sequence.len(),
            1,
            vec![10., 10., 0.],
            ups_and_downs,
        );

        var_mat.add_contig(var_stats, 2, 1, 0);

        var_mat.generate_distances();
        let mut ref_map = HashMap::new();
        ref_map.insert(0, "test".to_string());
        let multi = Arc::new(MultiProgress::new());
        var_mat.run_fuzzy_scan(0.01, 0.05, 0.01, 0.01, 0., 0, 0., 0, &ref_map, &multi)
    }
}
