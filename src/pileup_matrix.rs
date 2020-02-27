use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use pileup_structs::*;
use matrix_handling::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use ndarray::{Array2, prelude::*};
use ndarray_linalg::{norm::*};
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::process;
use std::cmp;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use crate::{tempfile, finish_command_safely};
use crate::factorization::{nmf::*, seeding::*};
use crate::dbscan::fuzzy;
use itertools::{Itertools};
use ordered_float::NotNan;


#[derive(Debug)]
pub enum PileupMatrix {
    PileupContigMatrix {
        coverages: HashMap<i32, Vec<f64>>,
        average_genotypes: HashMap<i32, Vec<f64>>,
        variances: HashMap<i32, Vec<f64>>,
        variants: HashMap<i32, HashMap<i32, BTreeMap<String, Vec<(f64, f64)>>>>,
        snps_map: HashMap<i32, HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>,
        indels_map: HashMap<i32, HashMap<i32, BTreeMap<String, BTreeSet<i64>>>>,
        contigs: HashMap<i32, Vec<u8>>,
        target_names: HashMap<i32, String>,
        target_lengths: HashMap<i32, f64>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        clusters: HashMap<i32, HashMap<i32, BTreeMap<String, (i32, usize)>>>,
        clusters_mean: HashMap<i32, f64>,
        variant_counts: HashMap<usize, HashMap<i32, usize>>,
        variant_sums: HashMap<usize, HashMap<i32, Vec<Vec<f64>>>>,
        variant_info: Vec<(i32, String, (Vec<f64>, Vec<f64>), i32)>,
        geom_mean_var: Vec<f64>,
        geom_mean_dep: Vec<f64>,
        pred_variants: HashMap<usize, HashMap<i32, HashMap<i32, HashMap<fuzzy::Category, HashSet<String>>>>>,
        pred_variants_all: HashMap<usize, HashMap<i32, HashMap<i32, HashSet<String>>>>,
    }
}

impl PileupMatrix {
    pub fn new_matrix() -> PileupMatrix {
        PileupMatrix::PileupContigMatrix {
            coverages: HashMap::new(),
            variances: HashMap::new(),
            average_genotypes: HashMap::new(),
            variants: HashMap::new(),
            snps_map: HashMap::new(),
            indels_map: HashMap::new(),
            contigs: HashMap::new(),
            target_names: HashMap::new(),
            target_lengths: HashMap::new(),
            sample_names: Vec::new(),
            kfrequencies: BTreeMap::new(),
            clusters: HashMap::new(),
            clusters_mean: HashMap::new(),
            variant_counts: HashMap::new(),
            variant_sums: HashMap::new(),
            variant_info: Vec::new(),
            geom_mean_var: Vec::new(),
            geom_mean_dep: Vec::new(),
            pred_variants: HashMap::new(),
            pred_variants_all: HashMap::new(),
        }
    }
}

pub trait PileupMatrixFunctions {
    fn setup(&mut self);

    fn add_sample(&mut self, sample_name: String);

    fn add_contig(&mut self,
                  pileup_stats: PileupStats,
                  sample_count: usize,
                  sample_idx: usize,
                  contig: Vec<u8>);

    fn generate_distances(&mut self, threads: usize, output_prefix: &str);

    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64);

    fn run_nmf(&mut self);

    fn generate_genotypes(&mut self,
                          output_prefix: &str);

    fn print_variant_stats(&self, output_prefix: &str);


}

impl PileupMatrixFunctions for PileupMatrix{
    fn setup(&mut self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut variants,
                ref mut snps_map,
                ref mut indels_map,
                ref mut contigs,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ref mut kfrequencies,
                ..
            } => {
                *coverages = HashMap::new();
                *average_genotypes = HashMap::new();
                *variants = HashMap::new();
                *snps_map = HashMap::new();
                *indels_map = HashMap::new();
                *contigs = HashMap::new();
                *target_names = HashMap::new();
                *target_lengths = HashMap::new();
                *sample_names = vec!();
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn add_sample(&mut self, sample_name: String) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut sample_names,
                ..
            } => {
                sample_names.push(sample_name);
            }
        }
    }

    fn add_contig(&mut self,
                  mut pileup_stats: PileupStats,
                  sample_count: usize,
                  sample_idx: usize,
                  contig: Vec<u8>) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut variants,
                ref mut snps_map,
                ref mut indels_map,
                ref mut contigs,
                ref mut target_names,
                ref mut target_lengths,
                ref mut variances,
                ref mut variant_counts,
                ref mut variant_sums,
                ..
            } => {
                match pileup_stats {
                    PileupStats::PileupContigStats {
                        tid,
                        target_name,
                        target_len,
                        coverage,
                        variance,
                        mean_genotypes,
                        variant_abundances,
                        indels,
                        nucfrequency,
                        total_variants,
                        depth,
//                        variations_per_n,
                        ..
                    } => {
                        let ag = average_genotypes.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]);
                        ag[sample_idx] = mean_genotypes;
                        let var = variances.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]
                        );
                        var[sample_idx] = variance;
                        let cov = coverages.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]
                        );
                        cov[sample_idx] = coverage;
                        target_names.entry(tid)
                            .or_insert(str::from_utf8(&target_name).unwrap().to_string());
                        target_lengths.entry(tid).or_insert(target_len);

                        // Initialize contig id in variant hashmap
                        let mut contig_variants = variants.entry(tid)
                            .or_insert(HashMap::new());


                        // Apppend the sample index to each variant abundance... so many loops >:(
                        // Initialize the variant position index
                        // Also turns out to be the total number of variant positions
                        for (pos, total_depth) in depth.iter().enumerate() {
                            let position_variants = contig_variants.entry(pos as i32)
                                .or_insert(BTreeMap::new());
                            let mut variant_depth: f64 = 0.;

                            if variant_abundances.contains_key(&(pos as i32)) {
                                let abundance_map = variant_abundances.get(&(pos as i32)).unwrap();
                                for (variant, abundance) in abundance_map.iter() {
                                    let sample_map = position_variants.entry(variant.clone())
                                        .or_insert(vec![(0., 0.); sample_count]);
                                    variant_depth += (abundance.0) as f64;
                                    sample_map[sample_idx] = (abundance.0 as f64, *total_depth);

                                }
                            }
                            // add pseudocounts
                            let ref_depth = *total_depth
                                                    - variant_depth;

                            //Add Reference as variant
                            let sample_map = position_variants.entry("R".to_string())
                                .or_insert(vec![(0., 0.); sample_count]);
                            sample_map[sample_idx] = (ref_depth, *total_depth);

                        }

                        let contig_indels = indels_map.entry(tid)
                            .or_insert(HashMap::new());
                        if contig_indels.len() == 0 {
                            *contig_indels = indels;
                        } else {
                            for (pos, indel_map) in indels.iter(){
                                let position_indels = contig_indels.entry(*pos)
                                    .or_insert(BTreeMap::new());
                                if position_indels.len() == 0 {
                                    *position_indels = indel_map.clone();
                                } else {
                                    for (indel, read_set) in indel_map.iter() {
                                        let read_map = position_indels.entry(indel.clone())
                                            .or_insert(BTreeSet::new());
                                        if read_map.len() == 0 {
                                            *read_map = read_set.clone();
                                        } else {
                                            let new_read_set = read_map.union(read_set).cloned().collect();
                                            *read_map = new_read_set;
                                        }
                                    }
                                }
                            }
                        }

                        let contig_snps = snps_map.entry(tid)
                            .or_insert(HashMap::new());
                        if contig_snps.len() == 0 {
                            *contig_snps = nucfrequency;
                        } else {
                            for (pos, snp_map) in nucfrequency.iter(){
                                let mut position_snps = contig_snps.entry(*pos)
                                    .or_insert(BTreeMap::new());
                                if position_snps.len() == 0 {
                                    *position_snps = snp_map.clone();
                                } else {
                                    for (snp, read_set) in snp_map.iter() {
                                        let read_map = position_snps.entry(*snp)
                                            .or_insert(BTreeSet::new());
                                        if read_map.len() == 0 {
                                            *read_map = read_set.clone();
                                        } else {
                                            let new_read_set = read_map.union(read_set).cloned().collect();
                                            *read_map = new_read_set;
                                        }
                                    }
                                }
                            }
                        }

                        contigs.entry(tid).or_insert(contig);
                        
                    }
                }
            }
        }
    }

    fn generate_distances(&mut self, threads: usize, output_prefix: &str) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                target_names,
                sample_names,
                coverages,
                contigs,
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut pred_variants_all,
                ..
            } => {

                let sample_count = sample_names.len() as f64;


                let variant_info_all =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                // A vector the length of the number of samples
                // Cumulatively calculates the product of abundaces from each sample
                let geom_mean_v =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                let geom_mean_d =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                // get basic variant info
                variants.par_iter().for_each(|(tid, variant_abundances)| {
                    let contig_coverages = coverages.get(tid)
                        .expect("Unable to retrieve contig coverage");

                    let max_coverage = contig_coverages.iter().cloned().fold1(f64::max)
                        .expect("Unable to retrieve max coverage");

                    variant_abundances.par_iter().for_each(
                        |(position, hash)| {
                            // loop through each position that has variants ignoring positions that
                            // only contained the reference in all samples
                            if hash.keys().len() > 1 {
                                let mut depths = Vec::new();
                                let a_vec = &hash[&"R".to_string()];
                                for (_a, d) in a_vec.iter() {
                                    depths.push(*d);
                                };
                                for (variant, abundances_vector) in hash.iter() {
                                    let mut abundance: f64 = 0.;
                                    let mut mean_var: f64 = 0.;
//                                let mut total_d: f64 = 0.;
                                    let mut freqs = Vec::new();
                                    // Get the mean abundance across samples
                                    let mut sample_idx: usize = 0;
                                    abundances_vector.iter().for_each(|(var, d)| {
//                                        mean_var += *var;
                                        // Total depth of location
//                                        total_d += *d;
//                                            let freq = (*var + 1.) / (*d + 1.);
                                        let mut geom_mean_v =
                                            geom_mean_v.lock().unwrap();
                                        let mut geom_mean_d =
                                            geom_mean_d.lock().unwrap();

                                        let sample_coverage = contig_coverages[sample_idx];

//                                            freqs.push(freq * (sample_coverage / max_coverage));
                                        freqs.push(*var);
//                                        if variant == &"R".to_string() {
                                            geom_mean_v[sample_idx] += ((*var + 1.) as f64).ln();
                                            geom_mean_d[sample_idx] += ((depths[sample_idx] + 1.) as f64).ln();
//                                        }
                                        sample_idx += 1;
                                    });


                                    let mut variant_info_all = variant_info_all
                                        .lock().unwrap();


                                    variant_info_all.push(
                                        (*position, variant.to_string(),
                                         (depths.clone(), freqs), *tid));
                                }
                            }
                        });
                });
                let mut variant_info_all = variant_info_all.lock().unwrap().clone();
                let geom_mean = |input: &Vec<f64>| -> Vec<f64> {
                    let output = input.iter()
                        .map(|sum| {
                            (sum / variant_info_all.len() as f64).exp()
                        }).collect::<Vec<f64>>();
                    return output
                };

                let mut geom_mean_v = geom_mean_v.lock().unwrap().clone();
                let geom_mean_v = geom_mean(&geom_mean_v);
                debug!("Geom Mean Var {:?}", geom_mean_v);

                let mut geom_mean_d = geom_mean_d.lock().unwrap().clone();
                let geom_mean_d = geom_mean(&geom_mean_d);
                debug!("Geom Mean Dep {:?}", geom_mean_d);

                *variant_info = variant_info_all;
                *geom_mean_var = geom_mean_v;
                *geom_mean_dep = geom_mean_d;

            }
        }
    }

    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut pred_variants,
                ..
            } => {

                info!("Running fuzzyDBSCAN with {} Variants", variant_info.len());

                let points = Arc::new(Mutex::new(Vec::new()));
                variant_info.into_par_iter().for_each(|info| {
                    let point = fuzzy::Point{
                       pos: info.0,
                       var: info.1.clone(),
                       geom_var: geom_mean_var.clone(),
                       geom_dep: geom_mean_dep.clone(),
                       deps: (info.2).0.clone(),
                       vars: (info.2).1.clone(),
                       tid: info.3
                    };
                    let mut points = points.lock().unwrap();
                    points.push(point);
                });
                let mut fuzzy_scanner = fuzzy::FuzzyDBSCAN {
                    eps_min: e_min,
                    eps_max: e_max,
                    pts_min: pts_min*variant_info.len() as f64,
                    pts_max: pts_max*variant_info.len() as f64,
                };
                let points = points.lock().unwrap();
                let clusters = fuzzy_scanner.cluster(&points[..]);
                let prediction_variants = Arc::new(
                    Mutex::new(
                        HashMap::new()));
                let prediction_count = Arc::new(
                    Mutex::new(
                        HashMap::new()));
                let prediction_features = Arc::new(
                    Mutex::new(
                        HashMap::new()));



                clusters.par_iter().enumerate().for_each(|(rank, cluster)|{
                    cluster.par_iter().for_each(|assignment|{
                        let variant = &points[assignment.index];

                        let mut prediction_variants = prediction_variants.lock().unwrap();
                        let variant_tid = prediction_variants.entry(rank + 1)
                            .or_insert(HashMap::new());

                        let variant_pos = variant_tid.entry(variant.tid).or_insert(HashMap::new());

                        let variant_cat = variant_pos.entry(variant.pos).or_insert(HashMap::new());
                        let variant_set = variant_cat.entry(assignment.category).or_insert(HashSet::new());
                        variant_set.insert(variant.var.clone());

                        let mut prediction_count = prediction_count.lock().unwrap();
                        let count = prediction_count.entry(rank + 1).or_insert(HashSet::new());
                        count.insert(assignment.index);

                        let mut prediction_features = prediction_features.lock().unwrap();
                        let feature = prediction_features.entry(rank+1).or_insert(HashMap::new());
                        let category = feature.entry(assignment.category).or_insert(0);
                        *category += 1;

                    });
                });
                let prediction_variants = prediction_variants.lock().unwrap().clone();
                let prediction_count = prediction_count.lock().unwrap();
                let prediction_features = prediction_features.lock().unwrap();
                // Pairs of clusters that shared border points
//                let mut to_combine = HashMap::new();
//                // Clusters that were completely contained within another cluster
//                let mut to_remove = HashSet::new();
//                for combination in prediction_count.iter().combinations(2) {
//                    let intersect: HashSet<usize> = combination[0].1.intersection(combination[1].1).collect();
//                    if intersect.len() == combination[0].1.len() {
//                        to_remove.insert(combination[0].0);
//                    } else if intersect.len() == combination[1].1.len() {
//                        to_remove.insert(combination[1].0);
//                    } else if intersect.len() > 0 {
//                        if combination[0].1.len() >= combination[1].1.len() {
//
//                        }
//                    }
//                }

                debug!("Predictions {:?}", prediction_variants);
                debug!("Prediction count {:?}", prediction_count);
                debug!("Prediction categories {:?}", prediction_features);
                *pred_variants = prediction_variants;
            }
        }
    }

    fn run_nmf(&mut self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut pred_variants,
                sample_names,
                indels_map,
                snps_map,
                ..
            } => {
                if variant_info.len() > 1 {
                    let sample_count = sample_names.len() as f64;

                    info!("Generating Variant Distances with {} Variants", variant_info.len());



                    // Generate distance or covariance matrix and pass that NMF functions
                    let v = get_condensed_distances(&variant_info[..],
                                                    indels_map,
                                                    snps_map,
                                                    &geom_mean_var[..],
                                                    &geom_mean_dep[..],
                                                    sample_count as i32);

                    let v = v.lock().unwrap();

                    let mut v = v.get_array2();

                    let v = v.clone() / v.norm();
                    let mut nmf = Factorization::new_nmf(
                        v,
                        20,
                        1,
                        "euclidean",
                        "conn",
                        Seed::new_random_vcol(),
                        30,
                        500,
                        1e-5);

                    nmf.estimate_rank();

                    info!("EVAR: {} RSS: {}", nmf.evar(), nmf.rss());

                    let mut predictions = nmf.probabilities("samples");

                    let mut prediction_count = Arc::new(Mutex::new(HashMap::new()));
                    let mut prediction_features = Arc::new(Mutex::new(HashMap::new()));
                    let mut prediction_variants = Arc::new(Mutex::new(HashMap::new()));
    //                    let mut prediction_variants_all = Arc::new(Mutex::new(HashMap::new()));

                    let mut max_cnt = 0;
                    let mut max_strain = 0;
                    let thresh = 1. / predictions.shape()[1] as f64;
                    // check if feature score is greater than geom mean score / SD
                    // if so then place into that rank
                    // If not, then prediction could realistically be any available rank
                    variant_info.par_iter().enumerate().for_each(|(row, info)| {
                        let scores = predictions.slice(s![row, ..]).clone();
                        let notnan_scores: Vec<NotNan<f64>> = scores.into_par_iter().cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();
                        let max = notnan_scores.par_iter().max().expect("No maximum found");

                        scores.iter().enumerate().for_each(|(rank, score)|{
                            if score >= &(thresh) || NotNan::from(*score) == *max {
                                let mut prediction_count = prediction_count.lock().unwrap();
                                let mut prediction_features = prediction_features.lock().unwrap();
                                let mut prediction_variants = prediction_variants.lock().unwrap();

                                let rank: usize = rank + 1;
                                let count = prediction_count.entry(rank)
                                    .or_insert(0);
                                let variant_tid = prediction_variants.entry(rank)
                                    .or_insert(HashMap::new());

                                // variant_info_all.push((position, var.to_string(), (depths, freqs), tid));
                                let variant_pos = variant_tid.entry(info.3).or_insert(HashMap::new());

                                let variant_cat = variant_pos.entry(info.0).or_insert(HashMap::new());
                                let variant = variant_cat.entry(fuzzy::Category::Core).or_insert(HashSet::new());
                                variant.insert(info.1.clone());

                                let feature = prediction_features.entry(rank)
                                    .or_insert(0.);
                                *feature += score.ln();

                                *count += 1;
                            }
                        });
                    });
                    let mut prediction_count = prediction_count.lock().unwrap();
                    let mut prediction_features = prediction_features.lock().unwrap();
                    let mut prediction_variants = prediction_variants.lock().unwrap().clone();
                    // get the strain with maximum members
                    let mut max_cnt = 0;
                    let mut max_strain = 0;
                    prediction_count.iter().for_each(|(strain, cnt)| {
                        if &max_cnt <= cnt {
                            max_strain = *strain;
                            max_cnt = *cnt;
                        }
                    });


                    let mut prediction_geom = HashMap::new();
                    prediction_features.iter()
                        .for_each(|(pred, sum)| {
                            prediction_geom.insert(pred, (*sum / prediction_count[pred] as f64).exp());
                        });

                    debug!("Prediction Geom Means {:?}", prediction_geom);
                    debug!("Prediction Counts {:?}", prediction_count);

                    *pred_variants = prediction_variants;
                } else {
                    debug!("Not enough variants found for NMF calculation: {:?}", variant_info);
                    process::exit(1)
                }
            }
        }
    }

    fn generate_genotypes(&mut self, output_prefix: &str) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                indels_map,
                snps_map,
                target_names,
                sample_names,
                coverages,
                contigs,
                ref mut pred_variants,
                ref mut pred_variants_all,
                ..
            } => {

                pred_variants.par_iter().for_each(|(strain_index, genotype)|{
                    if *strain_index != 0 {

                        let file_name = format!("{}_strain_{}.fna", output_prefix.to_string(), strain_index);

                        let file_path = Path::new(&file_name);

                        // Open haplotype file or create one
                        let mut file_open = File::create(file_path)
                            .expect("No Read or Write Permission in current directory");

                        let mut genotype = genotype.clone();

                        // Generate the variant genome
                        for (tid, original_contig) in contigs.iter() {
                            let mut contig = String::new();

                            let mut skip_n = 0;
                            let mut skip_cnt = 0;
                            let mut char_cnt = 0;
                            let mut variations = 0;

                            for (pos, base) in original_contig.iter().enumerate() {
                                if skip_cnt < skip_n {
                                    skip_cnt += 1;
                                } else {

                                    skip_n = 0;
                                    skip_cnt = 0;
                                    if genotype.contains_key(&tid) {
                                        let mut tid_genotype = genotype.get_mut(&tid).unwrap();

//                                        if pred_variants_all.contains_key(&0) {
//                                            if pred_variants_all[&0].contains_key(&tid) {
//                                                tid_genotype
//                                                    .extend(pred_variants_all[&0][&tid].clone());
//                                            }
//                                        };

                                        if tid_genotype.contains_key(&(pos as i32)) {

                                            let categories = &genotype[tid][&(pos as i32)];
                                            let mut hash= HashSet::new();
                                            if categories.contains_key(&fuzzy::Category::Core) {
                                                hash = categories[&fuzzy::Category::Core].clone();
                                            } else if categories.contains_key(&fuzzy::Category::Border) {
                                                hash = categories[&fuzzy::Category::Border].clone();
                                            }

                                            let mut max_var = "";
                                            for var in hash.iter() {
                                                if max_var == "" {
                                                    max_var = var;
                                                } else if max_var == "R" {
                                                    max_var = var;
                                                }
                                                variations += 1;
                                            }
                                            if max_var.contains("N") {
                                                // Skip the next n bases but rescue the reference prefix
                                                skip_n = max_var.len() - 1;
                                                skip_cnt = 0;
                                                let first_byte = max_var.as_bytes()[0];
                                                contig = contig + str::from_utf8(
                                                    &[first_byte]).unwrap()
                                            } else if max_var.len() > 1 {
                                                // Insertions have a reference prefix that needs to be removed
                                                let removed_first_base = str::from_utf8(
                                                    &max_var.as_bytes()[1..]).unwrap();
                                                contig = contig + removed_first_base;
                                            } else if max_var.contains("R") {
                                                contig = contig + str::from_utf8(&[*base]).unwrap();
                                            } else {
                                                contig = contig + max_var;
                                            }
                                        } else {
                                            contig = contig + str::from_utf8(&[*base]).unwrap();
                                        }
                                    } else {
                                        contig = str::from_utf8(&original_contig)
                                            .expect("Can't convert to str").to_string();
                                    }
                                }
                            };
                            writeln!(file_open, ">{}_strain_{}\t#variants_{}",
                                     target_names[tid],
                                     strain_index,
                                     variations).expect("Unable to write to file");


                            for line in contig.as_bytes().to_vec()[..].chunks(60).into_iter() {
                                file_open.write(line).unwrap();
                                file_open.write(b"\n").unwrap();
                            };
                        }
                    }
                });
            }
        }
    }

    fn print_variant_stats(&self, output_prefix: &str) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                target_names,
                target_lengths,
                sample_names,
                variances,
                variant_counts,
                variant_sums,
                ..
            } => {
                let file_name = output_prefix.to_string()
                    + &".tsv".to_owned();
                let file_path = Path::new(&file_name);
                let mut file_open = match File::create(file_path) {
                    Ok(fasta) => fasta,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    },
                };
                write!(file_open, "contigName\tcontigLen").unwrap();
                for sample_name in sample_names.iter(){
                    write!(file_open,
                           "\t{}.subsPer10kb\t{}.variants\t{}.meanRefAbd\
                            \t{}.refStdDev\t{}.meanVarAbd\t{}.varStdDev",
                           &sample_name, &sample_name, &sample_name,
                           &sample_name, &sample_name, &sample_name).unwrap();
                }
                write!(file_open, "\n").unwrap();
                for (tid, contig_name) in target_names.iter() {
                    let contig_len =  target_lengths[tid];
                    write!(file_open, "{}\t{}", contig_name, contig_len).unwrap();
                    for (sample_idx, _sample_name) in sample_names.iter().enumerate() {
                        let ten_kbs = contig_len / 10000.;
                        let total_variants = variant_counts[&sample_idx][tid] as f64;
                        if total_variants > 0. {
                            let var_ten_kbs = total_variants / ten_kbs;
                            let sample_sums = &variant_sums[&sample_idx][tid];

//                            let var_ratios = sample_sums[0]
//                                .iter().zip(&sample_sums[1])
//                                .map(|(var, dep)| { var / dep }).collect::<Vec<f64>>();
//
//                            let refr_ratios = sample_sums[2]
//                                .iter().zip(&sample_sums[1])
//                                .map(|(refr, dep)| { refr / dep }).collect::<Vec<f64>>();

                            let var_ratios_mean: f64 = sample_sums[0].iter().sum::<f64>()
                                / sample_sums[1].len() as f64;

                            let refr_ratios_mean: f64 = sample_sums[2].iter().sum::<f64>()
                                / sample_sums[1].len() as f64;

                            let mut var_std: f64 = sample_sums[0].iter().map(|x|
                                {(*x - var_ratios_mean).powf(2.)}).collect::<Vec<f64>>().iter().sum::<f64>();
                            var_std = (var_std / (sample_sums[1].len()) as f64).powf(1./2.);

                            let mut ref_std: f64 = sample_sums[2].iter().map(|x|
                                {(*x - refr_ratios_mean).powf(2.)}).collect::<Vec<f64>>().iter().sum::<f64>();
                            ref_std = (ref_std / (sample_sums[1].len()) as f64).powf(1./2.);

                            writeln!(file_open,
                                     "\t{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                                     var_ten_kbs, total_variants,
                                     refr_ratios_mean, ref_std,
                                     var_ratios_mean, var_std).unwrap();
                        } else {
                            writeln!(file_open,
                                     "\t{}\t{}\t{}\t{}\t{}\t{}",
                                     0., 0., 0., 0., 0., 0.,).unwrap();
                        }
                    }
                }
            }
        }
    }
}