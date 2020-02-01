use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use pileup_structs::*;
use matrix_handling::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use ndarray::{Array2};
use ndarray_npy::read_npy;
use kodama::{Dendrogram};
use std::sync::{Arc, Mutex};
use haplotypes_and_genotypes::*;
use std::fs::File;
use std::process;
use std::cmp;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use ::{tempfile, finish_command_safely};
use itertools::Itertools;


#[derive(Debug)]
pub enum PileupMatrix {
    PileupContigMatrix {
        coverages: HashMap<i32, Vec<f32>>,
        average_genotypes: HashMap<i32, Vec<f32>>,
        variances: HashMap<i32, Vec<f32>>,
        variants: HashMap<i32, HashMap<i32, BTreeMap<String, Vec<(f32, f32)>>>>,
        snps_map: HashMap<i32, HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>,
        indels_map: HashMap<i32, HashMap<i32, BTreeMap<String, BTreeSet<i64>>>>,
        contigs: HashMap<i32, Vec<u8>>,
        target_names: HashMap<i32, String>,
        target_lengths: HashMap<i32, f32>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        dendrogram: Dendrogram<f32>,
        clusters: HashMap<i32, HashMap<i32, BTreeMap<String, (i32, usize)>>>,
        clusters_mean: HashMap<i32, f32>,
        variant_counts: HashMap<usize, HashMap<i32, usize>>,
        variant_sums: HashMap<usize, HashMap<i32, Vec<Vec<f32>>>>,
        pred_variants: HashMap<i32, HashMap<i32, HashMap<i32, HashSet<String>>>>,
        pred_variants_all: HashMap<i32, HashMap<i32, HashMap<i32, HashSet<String>>>>,
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
            sample_names: vec!(),
            kfrequencies: BTreeMap::new(),
            dendrogram: Dendrogram::new(0),
            clusters: HashMap::new(),
            clusters_mean: HashMap::new(),
            variant_counts: HashMap::new(),
            variant_sums: HashMap::new(),
            pred_variants: HashMap::new(),
            pred_variants_all: HashMap::new(),
        }
    }
}

pub trait PileupMatrixFunctions {
    fn setup(&mut self);

    fn add_sample(&mut self, sample_name: String);

    fn add_kmers(&mut self,
                 tid: i32,
                 number_of_contigs: usize,
                 k_freq: BTreeMap<Vec<u8>, usize>);

    fn add_contig(&mut self,
                  pileup_stats: PileupStats,
                  sample_count: usize,
                  sample_idx: usize,
                  contig: Vec<u8>);

    fn generate_distances(&mut self, threads: usize, output_prefix: &str);

    fn generate_genotypes(&mut self,
                          output_prefix: &str);

    fn print_matrix(&self);

    fn print_variant_stats(&self, output_prefix: &str);

    fn print_kmers(&self, output_prefix: &str, kmer_size: &usize);

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

    fn add_kmers(&mut self, tid: i32, number_of_contigs: usize, k_freq: BTreeMap<Vec<u8>, usize>) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut kfrequencies,
                target_lengths,
                ..
            } => {
                if !target_lengths.contains_key(&tid) {
                    let contig_order_id = tid as usize;
                    for (tet, count) in k_freq.iter() {
                        let count_vec = kfrequencies.entry(tet.to_vec())
                                                    .or_insert(vec![0; number_of_contigs]);
                        count_vec[contig_order_id] = *count;
                    }
                }
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
//                        variations_per_n,
                        ..
                    } => {
                        let ag = average_genotypes.entry(tid).or_insert(
                            vec![0.0 as f32; sample_count]);
                        ag[sample_idx] = mean_genotypes;
                        let var = variances.entry(tid).or_insert(
                            vec![0.0 as f32; sample_count]
                        );
                        var[sample_idx] = variance;
                        let cov = coverages.entry(tid).or_insert(
                            vec![0.0 as f32; sample_count]
                        );
                        cov[sample_idx] = coverage;
                        target_names.entry(tid)
                            .or_insert(str::from_utf8(&target_name).unwrap().to_string());
                        target_lengths.entry(tid).or_insert(target_len);

                        // Initialize contig id in variant hashmap
                        let mut contig_variants = variants.entry(tid)
                            .or_insert(HashMap::new());

                        let mut sample_sums = variant_sums.entry(sample_idx)
                            .or_insert(HashMap::new());
                        if total_variants > 0 {
                            let mut contig_sums = sample_sums.entry(tid)
                                .or_insert(vec![vec![0.; total_variants as usize]; 3]);

                            // Apppend the sample index to each variant abundance... so many loops >:(
                            // Initialize the variant position index
                            // Also turns out to be the total number of variant positions
                            let mut variant_index = 0;
                            for (pos, abundance_map) in variant_abundances.iter() {
                                let position_variants = contig_variants.entry(*pos)
                                    .or_insert(BTreeMap::new());
                                let mut variant_depth: f32 = 0.;
                                let mut total_depth: f32 = 0.;
                                for (variant, abundance) in abundance_map.iter() {
                                    let sample_map = position_variants.entry(variant.clone())
                                        .or_insert(vec![(0., 0.); sample_count]);
                                    variant_depth += abundance.0 as f32;
                                    total_depth = abundance.1 as f32 + 1 as f32;
                                    sample_map[sample_idx] = (abundance.0 as f32, abundance.1 as f32);
                                }
                                // add pseudocounts
                                let ref_depth = total_depth
                                                        - variant_depth;
                                variant_depth += 1.;

                                let geom_mean = ((variant_depth / total_depth)
                                    * (ref_depth / total_depth)).powf(1./2.);

                                contig_sums[0][variant_index] = variant_depth / total_depth;
                                contig_sums[2][variant_index] = ref_depth / total_depth;
                                contig_sums[1][variant_index] = total_depth;

                                variant_index += 1;
                            }

                            let contig_variant_counts = variant_counts.entry(sample_idx)
                                .or_insert(HashMap::new());
                            contig_variant_counts.insert(tid, variant_index);
                        } else {
                            let mut contig_sums = sample_sums.entry(tid)
                                .or_insert(vec![vec![0.]; 3]);

                            let contig_variant_counts = variant_counts.entry(sample_idx)
                                .or_insert(HashMap::new());
                            contig_variant_counts.insert(tid, 0);
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

                let sample_count = sample_names.len() as f32;


                let variant_info_all =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                // A vector the length of the number of samples
                // Cumulatively calculates the product of abundaces from each sample
                let geom_mean_var =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                let geom_mean_dep =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                // get basic variant info
                variants.par_iter().for_each(|(tid, variant_abundances)| {
                    let contig_coverages = coverages.get(tid)
                        .expect("Unable to retrieve contig coverage");

                    let max_coverage = contig_coverages.iter().cloned().fold1(f32::max)
                        .expect("Unable to retrieve max coverage");

                    variant_abundances.par_iter().for_each(
                        |(position, hash)| {
                            // loop through each position that has variants
                            for (var, abundances_vector) in hash.iter() {
                                let mut abundance: f32 = 0.;
                                let mut mean_var: f32 = 0.;
//                                let mut total_d: f32 = 0.;
                                let mut freqs = Vec::new();
                                let mut depths = Vec::new();
                                if !var.contains("R") {
                                    // Get the mean abundance across samples
                                    let mut sample_idx: usize = 0;
                                    abundances_vector.iter().for_each(|(var, d)| {
//                                        mean_var += *var;
                                        // Total depth of location
//                                        total_d += *d;
                                        if var > &0. {
//                                            let freq = (*var + 1.) / (*d + 1.);
                                            let mut geom_mean_var =
                                                geom_mean_var.lock().unwrap();
                                            let mut geom_mean_dep =
                                                geom_mean_dep.lock().unwrap();

                                            let sample_coverage = contig_coverages[sample_idx];

//                                            freqs.push(freq * (sample_coverage / max_coverage));
                                            freqs.push(*var + 1.);
                                            geom_mean_var[sample_idx] += ((*var + 1.) as f64).ln();
                                            geom_mean_dep[sample_idx] += ((*d + 1.) as f64).ln();

                                            depths.push(*d + 1.);
                                            abundance += *var / *d;
                                        } else {
                                            freqs.push(1.);
                                        }
                                        sample_idx += 1;
                                    });

//                                    mean_var = mean_var / sample_count;
//                                    abundance = abundance / sample_count;

                                    let mut variant_info_all = variant_info_all
                                        .lock().unwrap();


                                    variant_info_all.push(
                                        (position, var.to_string(),
                                         (depths, freqs), tid));
                                }
                            }
                        });
                });
                let mut variant_info_all = variant_info_all.lock().unwrap();
                if variant_info_all.len() > 1 {

                    info!("Generating Variant Distances with {} Variants", variant_info_all.len());
                    let geom_mean = |input: &Vec<f64>| -> Vec<f64> {
                        let output = input.iter()
                            .map(|sum| {
                                (sum / variant_info_all.len() as f64).exp()
                            }).collect::<Vec<f64>>();
                        return output
                    };
                    debug!("Geom Mean Var {:?}", geom_mean_var);

                    let mut geom_mean_var = geom_mean_var.lock().unwrap();
                    let geom_mean_var = geom_mean(&geom_mean_var);
                    debug!("Geom Mean Var {:?}", geom_mean_var);

                    let mut geom_mean_dep = geom_mean_dep.lock().unwrap();
                    let geom_mean_dep = geom_mean(&geom_mean_dep);
                    debug!("Geom Mean Dep {:?}", geom_mean_dep);

//                    geom_means.iter().enumerate().for_each(|(sample, geom_means)| {
//                        variant_info_all.iter_mut().for_each(|var| {
//
//                        });
//                    });

                    let tmp_dir = TempDir::new("lorikeet_fifo")
                        .expect("Unable to create temporary directory");
                    let fifo_path = tmp_dir.path().join("foo.pipe");

                    // create new fifo and give read, write and execute rights to the owner.
                    // This is required because we cannot open a Rust stream as a BAM file with
                    // rust-htslib.
                    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
                        .expect(&format!("Error creating named pipe {:?}", fifo_path));

                    let mut distances_file = tempfile::Builder::new()
                        .prefix("lorikeet-distances")
                        .tempfile_in(tmp_dir.path())
                        .expect(&format!("Failed to create distances tempfile"));

                    let mut constraints_file = tempfile::Builder::new()
                        .prefix("lorikeet-constraints")
                        .tempfile_in(tmp_dir.path())
                        .expect(&format!("Failed to create constraints tempfile"));

//                writeln!(distances_file, "{:?}", variant_distances).expect("Unable to write to tempfile");
                    let tmp_path_dist = distances_file.path().to_str()
                        .expect("Failed to convert tempfile path to str").to_string();

                    let tmp_path_cons = constraints_file.path().to_str()
                        .expect("Failed to convert tempfile path to str").to_string();

                    // TODO: Move NMF calculation to be within rust. Need extern crate
                    //       None are availble currently 29/01/2020

                    get_condensed_distances(&variant_info_all[..],
                                            indels_map,
                                            snps_map,
                                            &geom_mean_var[..],
                                            &geom_mean_dep[..],
                                            sample_count as i32,
                                            &tmp_path_dist,
                                            &tmp_path_cons);

                    run_nmf(&tmp_path_dist,
                                 &tmp_dir,
                                 threads,
                                 output_prefix,
                                 variant_info_all.len(),
                                sample_names);

                    let mut predictions: Array2<f32> = read_npy(tmp_path_dist + ".npy")
                        .expect("Unable to read predictions");

                    tmp_dir.close().expect("Unable to close temp directory");
                    debug!("Predictions {:?}", predictions);
                    let mut unique_ranks = HashSet::new();

                    predictions
                        .outer_iter().for_each(|row| {
                        unique_ranks.insert(row[0] as i32);
                    });

                    debug!("Unique ranks {:?}", unique_ranks);

                    let mut prediction_map = HashMap::new();
                    let mut prediction_count = HashMap::new();
                    let mut prediction_features = HashMap::new();
                    let mut prediction_variants = HashMap::new();
                    let mut prediction_variants_all = HashMap::new();

                    let mut max_cnt = 0;
                    let mut max_strain = 0;
                    let thresh = 1. / unique_ranks.len() as f32;
                    // check if prediction probability is greater than certain amount
                    // if so then place into that rank
                    // If not, then prediction could realistically be any available rank
                    for (row, variant_info) in variant_info_all.iter().enumerate() {
                        let prob = predictions[[row, 1]];
                        if prob >= thresh {
                            let rank = predictions[[row, 0]] as i32;

                            let prediction = prediction_map.entry(rank + 1)
                                .or_insert(0.);
                            let count = prediction_count.entry(rank + 1)
                                .or_insert(0.);
                            let variant_tid = prediction_variants.entry(rank + 1)
                                .or_insert(HashMap::new());

                            // variant_info_all.push((position, var.to_string(), (depths, freqs), tid));
                            let variant_pos = variant_tid.entry(*variant_info.3).or_insert(HashMap::new());

                            let variant = variant_pos.entry(*variant_info.0).or_insert(HashSet::new());
                            variant.insert(variant_info.1.to_owned());


                            let feature = prediction_features.entry(rank + 1).or_insert(0.);
                            *feature += predictions[[row, 2]];

                            *count += 1.;
                            *prediction += predictions[[row, 1]].ln();
                        } else {
                            // we add the variant to all strains
                            let rank = 0;

                            let prediction = prediction_map.entry(rank)
                                .or_insert(0.);
                            let count = prediction_count.entry(rank)
                                .or_insert(0.);
                            let variant_tid = prediction_variants_all.entry(rank)
                                .or_insert(HashMap::new());

                            // variant_info_all.push((position, var.to_string(), (depths, freqs), tid));
                            let variant_pos = variant_tid
                                .entry(*variant_info.3).or_insert(HashMap::new());

                            let variant = variant_pos
                                .entry(*variant_info.0).or_insert(HashSet::new());
                            variant.insert(variant_info.1.to_owned());


                            let feature = prediction_features.entry(rank).or_insert(0.);
                            *feature += predictions[[row, 2]];

                            *count += 1.;
                            *prediction += predictions[[row, 1]].ln();
                        }
                    }

                    // get the strain with maximum members
                    let mut max_cnt = 0.;
                    let mut max_strain = 0;
                    prediction_count.iter().map(|(strain, cnt)| {
                        if &max_cnt <= cnt {
                            max_strain = *strain;
                            max_cnt = *cnt;
                        }
                    });


                    let mut prediction_geom = HashMap::new();
                    prediction_map.iter()
                        .for_each(|(pred, sum)| {
                            prediction_geom.insert(pred, (sum / prediction_count[pred]).exp());
                        });

                    println!("Prediction Geom Means {:?}", prediction_geom);
                    println!("Prediction Counts {:?}", prediction_count);
                    println!("Prediction Features {:?}", prediction_features);

                    *pred_variants = prediction_variants;
                    *pred_variants_all = prediction_variants_all;


                } else {
                    debug!("Not enough variants found {:?}, Non heterogeneous population", variant_info_all);
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

                for (strain_index, genotype) in pred_variants.iter_mut() {
                    if strain_index != &0 {

                        let file_name = format!("{}_strain_{}.fna", output_prefix.to_string(), strain_index);

                        let file_path = Path::new(&file_name);

                        // Open haplotype file or create one
                        let mut file_open = File::create(file_path)
                            .expect("No Read or Write Permission in current directory");

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
                                    let mut max_var = "";

                                    skip_n = 0;
                                    skip_cnt = 0;
                                    if genotype.contains_key(&tid) {
                                        let mut tid_genotype = genotype.get_mut(&tid).unwrap();

                                        if pred_variants_all.contains_key(&0) {
                                            if pred_variants_all[&0].contains_key(&tid) {
                                                tid_genotype
                                                    .extend(pred_variants_all[&0][&tid].clone());
                                            }
                                        };

                                        if tid_genotype.contains_key(&(pos as i32)) {

                                            let hash = &genotype[tid][&(pos as i32)];

                                            for var in hash.iter() {
                                                max_var = var;
                                                variations += 1;
                                                break
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
                                     variations);


                            for line in contig.as_bytes().to_vec()[..].chunks(60).into_iter() {
                                file_open.write(line).unwrap();
                                file_open.write(b"\n").unwrap();
                            };
                        }
                    }
                }
            }
        }
    }

    fn print_matrix(&self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                target_names,
                sample_names,
                ..
            } => {

                for (tid, contig_name) in target_names.iter() {
                    for (pos, variant_map) in variants[tid].iter() {
                        for (variant, variant_depths) in variant_map.iter() {
                            print!("{}\t{}", contig_name, pos);
                            for sample_depths in variant_depths.iter(){
                                print!("\t{}", sample_depths.0/sample_depths.1);
                            }
                            print!("\n");
                        }
                    }
                }
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
                        let total_variants = variant_counts[&sample_idx][tid] as f32;
                        if total_variants > 0. {
                            let var_ten_kbs = total_variants / ten_kbs;
                            let sample_sums = &variant_sums[&sample_idx][tid];

//                            let var_ratios = sample_sums[0]
//                                .iter().zip(&sample_sums[1])
//                                .map(|(var, dep)| { var / dep }).collect::<Vec<f32>>();
//
//                            let refr_ratios = sample_sums[2]
//                                .iter().zip(&sample_sums[1])
//                                .map(|(refr, dep)| { refr / dep }).collect::<Vec<f32>>();

                            let var_ratios_mean: f32 = sample_sums[0].iter().sum::<f32>()
                                / sample_sums[1].len() as f32;

                            let refr_ratios_mean: f32 = sample_sums[2].iter().sum::<f32>()
                                / sample_sums[1].len() as f32;

                            let mut var_std: f32 = sample_sums[0].iter().map(|x|
                                {(*x - var_ratios_mean).powf(2.)}).collect::<Vec<f32>>().iter().sum::<f32>();
                            var_std = (var_std / (sample_sums[1].len()) as f32).powf(1./2.);

                            let mut ref_std: f32 = sample_sums[2].iter().map(|x|
                                {(*x - refr_ratios_mean).powf(2.)}).collect::<Vec<f32>>().iter().sum::<f32>();
                            ref_std = (ref_std / (sample_sums[1].len()) as f32).powf(1./2.);

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

    fn print_kmers(&self, output_prefix: &str, kmer_size: &usize) {
        match self {
            PileupMatrix::PileupContigMatrix {
                kfrequencies,
                target_names,
                ..
            } => {
                let file_name = output_prefix.to_string() + &"_".to_owned()
                                + &kmer_size.clone().to_string() + &"mer_counts".to_owned()
                                + &".tsv".to_owned();
                let file_path = Path::new(&file_name);
                let mut file_open = match File::create(file_path) {
                    Ok(fasta) => fasta,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    },
                };
                for (tid, name) in target_names.iter() {
                    write!(file_open, "{}\t",
                           name).unwrap();
                    for (_kmer, counts) in kfrequencies.iter(){
                        write!(file_open, "{}\t", counts[*tid as usize]).unwrap();
                    }
                    write!(file_open, "\n").unwrap();
                }
            }
        }
    }
}


fn run_nmf(dist_file_path: &str,
           tmp_dir: &TempDir,
           threads: usize,
           output_prefix: &str,
           n_variants: usize,
           sample_names: &Vec<String>) {
    let sample_count = sample_names.len();
    let mut converged = false;
    let mut max_rank = cmp::min(15, n_variants);
    let mut min_rank = cmp::min(4, n_variants);
    let mut best_rank = 0;
    let mut best_rss = 0.;
    while !converged {
        let mut ranks_rss = Arc::new(
            Mutex::new(vec![0.; max_rank - min_rank]));

        let mut in_threads = threads / (max_rank - min_rank - 1);
        if in_threads < 1 {
            in_threads = 1;
        }

        (min_rank..max_rank).into_par_iter().for_each(|rank| {
            let cmd_string = format!(
                "set -e -o pipefail; \
                     nice nmf.py {} True {} {} {} {} {}",
                // NMF
                rank + 1,
                10,
                dist_file_path,
                "tmp_path_cons",
                sample_count as i32,
                in_threads, );
            info!("Queuing cmd_string: {}", cmd_string);
            let mut python = std::process::Command::new("bash")
                .arg("-c")
                .arg(&cmd_string)
                .stderr(process::Stdio::piped())
                .stdout(process::Stdio::piped())
                .spawn()
                .expect("Unable to execute bash");

            python = finish_command_safely(python, "run_nmf");

            let mut out = String::new();
            python.stdout.expect("Failed to grab stdout from NMF").read_to_string(&mut out)
                .expect("Failed to read stdout to string");
            let mut ranks_rss = ranks_rss.lock().expect("Unable to lock RSS vec");
            let rss: f32 = match out.trim().parse() {
                Ok(value) => value,
                Err(error) => {
                    debug!("Unable to parse RSS {}", error);
                    0.
                }
            };
            ranks_rss[rank as usize - min_rank] = rss;
        });

        let ranks_rss = ranks_rss.lock().expect("unable to lock rss vec");

        debug!("RSS Values {:?}", ranks_rss);

        for (rank, rss) in ranks_rss.iter().enumerate() {
            if best_rank == 0 && best_rss == 0. && rank == 0 {
                best_rank = rank + min_rank + 1;
                best_rss = *rss;
            } else if &best_rss >= rss {
                best_rss = *rss;
                best_rank = rank + min_rank + 1;
            } else if rss > &best_rss {
                break
            }
        }
        if best_rank < max_rank {
            converged = true;
        } else {
            max_rank += 11;
            min_rank += 11;
        }
    }

    let cmd_string = format!(
        "set -e -o pipefail; \
                     nmf.py {} False {} {} {} {} {}",
        // NMF
        best_rank,
        30,
        dist_file_path,
        "tmp_path_cons",
        sample_count as i32,
        threads);
    info!("Queuing cmd_string: {}", cmd_string);
    let mut python = std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd_string)
        .stderr(process::Stdio::piped())
        .stdout(process::Stdio::piped())
        .spawn()
        .expect("Unable to execute bash");

    python = finish_command_safely(python, "run_nmf");
    let mut out = String::new();
    python.stdout.expect("Failed to grab stdout from NMF").read_to_string(&mut out)
        .expect("Failed to read stdout to string");
    println!("{}", sample_names[0]);
    println!("{}", out);
}