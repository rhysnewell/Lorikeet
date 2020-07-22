use std::collections::{HashMap, HashSet, BTreeMap};
use estimation::contig_variants::*;
use estimation::codon_structs::*;
use estimation::linkage::*;
use model::variants::*;
use std::{str};
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::fs::File;
use dbscan::fuzzy;
use itertools::{izip};
use rust_htslib::{bcf::{self}, bam::HeaderView};
use bio::io::gff;
use bird_tool_utils::command;
use std::process::{Stdio};
use tempfile;
use external_command_checker;
use rayon::current_num_threads;
use estimation::codon_structs::CodonTable;

#[derive(Debug)]
/// Container for all variants within a genome and associated clusters
pub enum VariantMatrix {
    VariantContigMatrix {
        coverages: HashMap<i32, Vec<f64>>,
        average_genotypes: HashMap<i32, Vec<f64>>,
        variances: HashMap<i32, Vec<f64>>,
        // TID, Position, Base, Var Depth, Total Depth
        all_variants: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
        // Placeholder hashmap for the depths of each contig for a sample
        // Deleted after use
        depths: HashMap<i32, Vec<i32>>,
        target_names: BTreeMap<i32, String>,
        target_lengths: HashMap<i32, f64>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        variant_counts: HashMap<usize, HashMap<i32, usize>>,
        variant_sums: HashMap<usize, HashMap<i32, Vec<Vec<f64>>>>,
        variant_info: Vec<fuzzy::Var>,
        geom_mean_var: Vec<f64>,
        geom_mean_dep: Vec<f64>,
        geom_mean_frq: Vec<f64>,
        pred_variants: HashMap<usize, HashMap<i32, HashMap<i64, HashMap<fuzzy::Category, HashSet<Variant>>>>>,
//        pred_variants_all: HashMap<usize, HashMap<i32, HashMap<i32, HashSet<String>>>>,
    }
}

impl VariantMatrix {
    pub fn new_matrix(sample_count: usize) -> VariantMatrix {
        VariantMatrix::VariantContigMatrix {
            coverages: HashMap::new(),
            variances: HashMap::new(),
            average_genotypes: HashMap::new(),
            all_variants: HashMap::new(),
            depths: HashMap::new(),
            target_names: BTreeMap::new(),
            target_lengths: HashMap::new(),
            sample_names: vec!["".to_string(); sample_count],
            kfrequencies: BTreeMap::new(),
            variant_counts: HashMap::new(),
            variant_sums: HashMap::new(),
            variant_info: Vec::new(),
            geom_mean_var: Vec::new(),
            geom_mean_dep: Vec::new(),
            geom_mean_frq: Vec::new(),
            pred_variants: HashMap::new(),
        }
    }
}

pub trait VariantMatrixFunctions {
    fn setup(&mut self);

    fn add_sample(&mut self, sample_name: String, sample_idx: usize,
                  variant_records: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
                  header: &HeaderView);

    /// Returns the alleles at the current position
    /// as a mutable reference
    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>>;

    fn variants_of_contig(&mut self, tid: i32) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>>;

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(&mut self,
                  variant_stats: VariantStats,
                  sample_count: usize,
                  sample_idx: usize);

    /// Converts all variants into fuzzy::Var format
    fn generate_distances(&mut self, threads: usize, output_prefix: &str);

    /// Perform fuzzy DBSCAN clustering using proportionality
    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64, phi: f64,
                      anchor_size: usize, anchor_similarity: f64, minimum_reads_in_link: usize);

    /// Takes clusters from DBSCAN and linkage method and writes variants to file as genotype
    fn generate_genotypes(&mut self,
                          output_prefix: &str,
                          reference: &mut bio::io::fasta::IndexedReader<File>);

    fn print_variant_stats(&self, output_prefix: &str, window_size: f64);

    fn calc_gene_mutation(&self,
                          gff_map: &mut HashMap<String, Vec<bio::io::gff::Record>>,
                          reference: &mut bio::io::fasta::IndexedReader<File>,
                          codon_table: &CodonTable,
                          ref_name: &str,
                          output_prefix: &str);

    fn write_vcf(&self, output_prefix: &str);
}

#[allow(unused)]
impl VariantMatrixFunctions for VariantMatrix {
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
                *target_names = BTreeMap::new();
                *target_lengths = HashMap::new();
                *sample_names = vec!();
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn add_sample(&mut self, sample_name: String, sample_idx: usize,
                  variant_records: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
                   header: &HeaderView) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                info!("adding sample {} at index {}", &sample_name, &sample_idx);
                sample_names[sample_idx] = sample_name;
                let target_count = header.target_count();
                let tid_names = header.target_names();

                for target_name in (tid_names).into_iter() {
                    let tid = header.tid(target_name).unwrap();
                    target_names.entry(tid as i32)
                        .or_insert(str::from_utf8(&target_name).unwrap().to_string());
                    let target_len = header.target_len(tid).unwrap();
                    target_lengths.entry(tid as i32).or_insert(target_len as f64);
                    // Initialize contig id in variant hashmap
                    let contig_variants = all_variants.entry(tid as i32)
                        .or_insert(HashMap::new());
                    let placeholder = HashMap::new();
                    let mut variants = match variant_records.get(&(tid as i32)) {
                        Some(map) => map,
                        _ => &placeholder,
                    };
                    let target_len = header.target_len(tid).unwrap();
                    // Apppend the sample index to each variant abundance
                    // Initialize the variant position index
                    // Also turns out to be the total number of variant positions
                    for (pos, abundance_map) in variants.iter() {
                        let position_variants = contig_variants.entry(*pos as i64)
                            .or_insert(HashMap::new());
                        for (variant, base_info) in abundance_map {
                            let sample_map = position_variants.entry(variant.clone())
                                .or_insert(base_info.clone());
                            sample_map.combine_sample(base_info, sample_idx, 0);
                        }
                    }
                }
            }
        }
    }

    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => {
                match all_variants.get_mut(&tid) {
                    Some(contig_variants) => {
                        contig_variants.get_mut(&pos.clone())
                    },
                    _ => {
                        None
                    },
                }
            }
        }
    }

    fn variants_of_contig(&mut self, tid: i32) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => {
                all_variants.get_mut(&tid)
            }
        }
    }

    fn add_contig(&mut self,
                  variant_stats: VariantStats,
                  sample_count: usize,
                  sample_idx: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut all_variants,
//                ref mut target_names,
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
                        let var = variances.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]
                        );
                        var[sample_idx] = variance;
                        let cov = coverages.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]
                        );
                        cov[sample_idx] = coverage;

                        // copy across depths
                        let contig_variants = all_variants.entry(tid)
                            .or_insert(HashMap::new());
                        for (pos, d) in depth.iter().enumerate() {
                            let position_variants = contig_variants.entry(pos as i64)
                                .or_insert(HashMap::new());
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

    fn generate_distances(&mut self, _threads: usize, output_prefix: &str) {
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

                let variant_info_all =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                // A vector the length of the number of samples
                // Cumulatively calculates the product of variant depths
                let geom_mean_v =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                // product of total depth
                let geom_mean_d =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                // product of reference frequency
                let geom_mean_f =
                    Arc::new(
                        Mutex::new(
                            vec![1. as f64; sample_count as usize]));

                // get basic variant info and store as fuzzy::Var
                all_variants.par_iter_mut().for_each(|(tid, variant_abundances)| {

                    variant_abundances.par_iter_mut().for_each(
                        |(position, hash)| {
                            // loop through each position that has variants ignoring positions that
                            // only contained the reference in all samples
                            if hash.keys().len() > 0 {

                                for (variant, base_info) in hash.iter_mut() {
                                    match variant {

                                        Variant::SNV(_) | Variant::None => {
                                            let _abundance: f64 = 0.;
                                            let _mean_var: f64 = 0.;
                                            let mut rel_abund = vec![0.0; sample_count as usize];

                                            // Get the mean abundance across samples
                                            for index in (0..sample_count).into_iter() {
                                                let mut geom_mean_v =
                                                    geom_mean_v.lock().unwrap();
                                                let mut geom_mean_d =
                                                    geom_mean_d.lock().unwrap();
                                                let mut geom_mean_f =
                                                    geom_mean_f.lock().unwrap();

                                                let mut var_depth
                                                    = base_info.truedepth[index] as f64;

                                                let total_depth
                                                    = base_info.totaldepth[index] as f64;
//                                                println!("var_depth {} tot {}", var_depth, total_depth);
//                                                base_info.freq[index] = ;
                                                if total_depth <= 0. {
                                                    rel_abund[index] =
                                                        var_depth / 1.;
                                                } else {
                                                    rel_abund[index] =
                                                        var_depth / total_depth;
                                                }

                                                geom_mean_v[index] += (var_depth + 1.).ln();
                                                geom_mean_d[index] += (total_depth + 1.).ln();
                                                geom_mean_f[index] += ((var_depth + 1.)
                                                    / (total_depth + 1.)).ln();
                                            };


                                            let mut variant_info_all = variant_info_all
                                                .lock().unwrap();
//                                            base_info.rel_abunds = rel_abund;
                                            let point = fuzzy::Var {
                                                pos: *position,
                                                var: variant.clone(),
                                                deps: base_info.totaldepth.clone(),
                                                vars: base_info.depth.clone(),
//                                                rel_abunds: rel_abund,
                                                tid: *tid,
                                                reads: base_info.reads.clone()
                                            };
                                            variant_info_all.push(point);
                                        },
                                        _ => {
                                            let mut rel_abund = vec![0.0; sample_count as usize];

                                            // Get the mean abundance across samples
                                            for index in (0..sample_count).into_iter() {
                                                let mut geom_mean_v =
                                                    geom_mean_v.lock().unwrap();
                                                let mut geom_mean_d =
                                                    geom_mean_d.lock().unwrap();
                                                let mut geom_mean_f =
                                                    geom_mean_f.lock().unwrap();

                                                let mut var_depth
                                                    = base_info.depth[index] as f64;
                                                if var_depth <= 0. {

                                                    var_depth
                                                        = base_info.truedepth[index] as f64;
                                                    base_info.depth[index] = base_info.truedepth[index];
                                                }
                                                let total_depth
                                                    = base_info.totaldepth[index] as f64;
//                                                base_info.freq[index] = ;
                                                if total_depth <= 0. {
                                                    rel_abund[index] =
                                                        var_depth / (1.);
                                                } else {
                                                    rel_abund[index] =
                                                        var_depth / total_depth;
                                                }

                                                geom_mean_v[index] += (var_depth + 1.).ln();
                                                geom_mean_d[index] += (total_depth + 1.).ln();
                                                geom_mean_f[index] += ((var_depth + 1.)
                                                    / (total_depth + 1.)).ln();
                                            };


                                            let mut variant_info_all = variant_info_all
                                                .lock().unwrap();
//                                            base_info.rel_abunds = rel_abund;
                                            let point = fuzzy::Var {
                                                pos: *position,
                                                var: variant.clone(),
                                                deps: base_info.totaldepth.clone(),
                                                vars: base_info.depth.clone(),
//                                                rel_abunds: rel_abund,
                                                tid: *tid,
                                                reads: base_info.reads.clone()
                                            };
                                            variant_info_all.push(point);
                                        },
                                    }
                                }
                            }
                        });
                });

                let variant_info_all = variant_info_all.lock().unwrap().clone();
                let geom_mean = |input: &Vec<f64>| -> Vec<f64> {
                    let output = input.iter()
                        .map(|sum| {
                            (sum / variant_info_all.len() as f64).exp()
                        }).collect::<Vec<f64>>();
                    return output
                };

                debug!("geoms {:?} {:?} {:?}", geom_mean_d, geom_mean_v, geom_mean_f);

                let geom_mean_v = geom_mean_v.lock().unwrap().clone();
                let geom_mean_v = geom_mean(&geom_mean_v);
                debug!("Geom Mean Var {:?}", geom_mean_v);

                let geom_mean_d = geom_mean_d.lock().unwrap().clone();
                let geom_mean_d = geom_mean(&geom_mean_d);
                debug!("Geom Mean Dep {:?}", geom_mean_d);

                let geom_mean_f = geom_mean_f.lock().unwrap().clone();
                let geom_mean_f = geom_mean(&geom_mean_f);
                debug!("Geom Mean Frq {:?}", geom_mean_f);
                debug!("geoms {:?} {:?} {:?}", geom_mean_d, geom_mean_v, geom_mean_f);

                info!("Beginning analysis of {} with {} variants", &output_prefix, &variant_info_all.len());

                *variant_info = variant_info_all;
                *geom_mean_var = geom_mean_v;
                *geom_mean_dep = geom_mean_d;
                *geom_mean_frq = geom_mean_f;

            }
        }
    }

    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64, phi: f64,
                      anchor_size: usize, anchor_similarity: f64, minimum_reads_in_link: usize) {
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

                if variant_info.len() == 1 {
                    warn!("Where did the variants go? {:?}", variant_info);
                }
                let fuzzy_scanner = fuzzy::FuzzyDBSCAN {
                    eps_min: e_min,
                    eps_max: e_max,
                    pts_min: match pts_min {
                        _ if pts_min > 1. => pts_min,
                        _ => pts_min*variant_info.len() as f64,
                        },
                    pts_max: match pts_max {
                        _ if pts_max > 1. => pts_max,
                        _ => pts_max*variant_info.len() as f64,
                    },
                    phi,
                    geom_var: geom_mean_var.to_owned(),
                    geom_dep: geom_mean_dep.to_owned(),
                    geom_frq: geom_mean_frq.to_owned(),
                };

                // Perform read phasing clustering and return initial clusters
                let links = linkage_clustering_of_variants(&variant_info,
                                                           anchor_size, anchor_similarity,
                                                           minimum_reads_in_link);

                // run fuzzy DBSCAN
                if links.len() > 0 {
                    info!("Running Seeded fuzzyDBSCAN with {} initial clusters", links.len());
                } else {
                    warn!("No initial clusters formed, running fuzzyDBSCAN with no seeds. Perhaps lower linkage thresholds?")
                }
                let clusters = fuzzy_scanner.cluster(
                    &variant_info[..],
                    links);

                // Since these are hashmaps, I'm using Arc and Mutex here since not sure how
                // keep hashmaps updated using channel()
                let prediction_variants =
                    Mutex::new(
                        HashMap::new());


                let all_variants =
                    Mutex::new(
                        all_variants);

                // Organize clusters into genotypes by recollecting full variant information
                clusters.par_iter().enumerate().for_each(|(rank, cluster)|{
                    // Sets for each cluster keeping track of which variant types are present in
                    // a cluster
                    let prediction_set = Mutex::new(
                        HashSet::new());
                    cluster.par_iter().for_each(|assignment|{

                        let variant: &fuzzy::Var = &variant_info[assignment.index];

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
                        let mut all_variants = all_variants.lock().unwrap();
                        let base_tid = all_variants.entry(variant.tid)
                            .or_insert(HashMap::new());
                        let base_pos = base_tid.entry(variant.pos)
                            .or_insert(HashMap::new());
                        match base_pos.get_mut(&variant.var) {
                            Some(base_base) => {
                                base_base.genotypes.insert(rank as i32 + 1);
                            },
                            None => {},
                        };

                    });
                    let mut prediction_set = prediction_set.lock().unwrap();
                    if !(prediction_set.len() == 1 && prediction_set.contains(&Variant::None)) {

                    } else {
                        let mut prediction_variants = prediction_variants
                            .lock()
                            .unwrap();

                        prediction_variants
                            .remove_entry(&(rank + 1)).expect("Unable to remove cluster");
                    }

                });

                let prediction_variants = prediction_variants.lock().unwrap();
                let mut all_variants = all_variants.lock().unwrap();
//                debug!("Prediction count {:?}", prediction_count);
//                debug!("Prediction categories {:?}", prediction_features);
                *pred_variants = prediction_variants.clone();
                **all_variants = all_variants.clone();
            }
        }
    }

    #[allow(unused)]
    fn generate_genotypes(&mut self, output_prefix: &str,
                          reference: &mut bio::io::fasta::IndexedReader<File>) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                ref mut pred_variants,
                ..
            } => {
                let reference = Mutex::new(reference);
                pred_variants.par_iter().for_each(|(strain_index, genotype)|{
                    let file_name = format!("{}_strain_{}.fna", output_prefix.to_string(), strain_index);

                    let file_path = Path::new(&file_name);
                    // Open haplotype file or create one
                    let mut file_open = File::create(file_path)
                        .expect("No Read or Write Permission in current directory");

                    let mut original_contig = Vec::new();

                    let mut genotype = genotype.clone();

                    let mut total_variant_alleles = 0;
                    let mut total_reference_alleles = 0;

                    // Generate the variant genome
                    for (tid, target_name) in target_names.iter() {
                        let mut contig = String::new();
                        original_contig = Vec::new();
                        {
                            let mut reference = reference.lock().unwrap();
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
                                            let categories = &genotype[tid][&(pos as i64)];
                                            let hash;
                                            if categories.contains_key(&fuzzy::Category::Core) {
                                                hash = categories[&fuzzy::Category::Core].clone();
                                            } else if categories.contains_key(&fuzzy::Category::Border) {
                                                hash = categories[&fuzzy::Category::Border].clone();
                                            } else {
                                                hash = categories[&fuzzy::Category::Noise].clone();
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
    //                                } else {
    //                                    contig = str::from_utf8(&original_contig)
    //                                        .expect("Can't convert to str").to_string();
    //                                }
                                }
                            };
                        } else {
                            contig = str::from_utf8(&original_contig)
                                .expect("Can't convert to str").to_string();
                        }

                        writeln!(file_open, ">{}_strain_{}_alt_alleles_{}_ref_alleles_{}",
                                 target_names[tid],
                                 strain_index,
                                 variations, ref_alleles).expect("Unable to write to file");

                        for line in contig.as_bytes().to_vec()[..].chunks(60).into_iter() {
                            file_open.write(line).unwrap();
                            file_open.write(b"\n").unwrap();
                        };
                        total_variant_alleles += variations;
                        total_reference_alleles += ref_alleles;
                    }
                    info!("Cluster {} contains {} variant alleles and {} reference alleles",
                          strain_index, total_variant_alleles, total_reference_alleles);
                });
            }
        }
    }

    fn print_variant_stats(&self, output_prefix: &str, window_size: f64) {
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

                let file_name = output_prefix.to_string()
                    + &".tsv".to_owned();
                let snp_locs = format!("{}_snp_locations.tsv", output_prefix);
                let file_path = Path::new(&file_name);
                let mut file_open = match File::create(file_path) {
                    Ok(fasta) => fasta,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    },
                };

                let snp_loc_path = Path::new(&snp_locs);
                let mut snp_loc_open = match File::create(snp_loc_path) {
                    Ok(tsv) => tsv,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    },
                };

                // Snp density summary start
                write!(file_open, "contigName\tcontigLen").unwrap();

                // Snp location start
                write!(snp_loc_open, "SNP\tchr\tpos").unwrap();
                debug!("Sample Names {:?}", sample_names);
                for sample_name in sample_names.iter(){
                    write!(file_open,
                           "\t{}.snvsPer{}kb\t{}.svsPer{}kb\t{}.snvCount\t{}.svCount",
                           &sample_name, &window_size, &sample_name, &window_size,
                           &sample_name, &sample_name).unwrap();
                    write!(snp_loc_open, "\t{}_depth", &sample_name).unwrap();
                }
                write!(file_open, "\n").unwrap();
                write!(snp_loc_open, "\n").unwrap();

                let mut snps = 0;

                for (tid, contig_name) in target_names.iter() {
                    let contig_len =  target_lengths[tid];
                    write!(file_open, "{}\t{}", contig_name, contig_len).unwrap();
                    match all_variants.get(tid) {
                        Some(variants_in_contig) => {
                            // Set up channels that receive a vector of values for each sample
                            let mut snps_cnt_vec = vec![0; sample_names.len()];
                            let svs_cnt_vec = Mutex::new(vec![0; sample_names.len()]);
//                            let (snp_freq_s, snp_freq_r) = channel();
//                            let (sv_freq_s, sv_freq_r) = channel();
//                            let (snp_std_s, snp_std_r) = channel();
//                            let (sv_std_s, sv_std_r) = channel();
                            let window = contig_len / window_size;
                            variants_in_contig.iter().enumerate()
                                .for_each(|(index, (position, variants))| {
                                    // Get how many alleles are present at loci
                                    let alleles = variants.len();
                                    for (var, base) in variants {
                                        match var {
                                            Variant::SNV(_) => {
                                                write!(snp_loc_open, "SNP{}\t{}\t{}",
                                                       index, contig_name, position).unwrap();
                                                base.truedepth
                                                    .iter()
                                                    .enumerate()
                                                    .zip(base.depth.iter().enumerate())
                                                    .for_each(|((index, count_1), (_index_2, count_2))| {
                                                        write!(snp_loc_open, "\t{}", count_1);
                                                        if count_1 > &0 || count_2 > &0 {
                                                            snps_cnt_vec[index] += 1;
                                                        }
                                                    });
                                                write!(snp_loc_open, "\n").unwrap();
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
                            let snps_per_win: Vec<_> = snps_cnt_vec.iter().map(|count| *count as f64 / window).collect();
                            let svs_per_win: Vec<_> = svs_cnt_vec.iter().map(|count| *count as f64 / window).collect();

                            for (snp_w, svs_w, snp_c, svs_c) in izip!(&snps_per_win, &svs_per_win, &snps_cnt_vec, &svs_cnt_vec) {
                                write!(file_open, "\t{}\t{}\t{}\t{}", snp_w, svs_w, snp_c, svs_c).unwrap();
                            }
                            write!(file_open, "\n").unwrap();

                        },
                        None => {
                            // Write out zeros for contigs with no variants
                            for (sample_idx, _sample_name) in sample_names.iter().enumerate() {
                                write!(file_open,
                                         "\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                         0., 0., 0., 0., 0., 0., 0.,).unwrap();
                            }
                        }
                    }
                };
                if snps > 0 {
                    let plot_command = format!("set -eou pipefail; snp_density_plots.R {} {} && \
                                                    mv SNP-Density*.pdf {}_snp_density_plot.pdf",
                                               output_prefix, window_size, output_prefix);
                    command::finish_command_safely(
                        std::process::Command::new("bash")
                            .arg("-c")
                            .arg(&plot_command)
                            .stderr(Stdio::piped())
                            .stdout(Stdio::piped())
                            .spawn()
                            .expect("Unable to execute Rscript"), "CMplot");
                } else {
                    std::fs::remove_file(snp_loc_path).expect("Unable to remove file");
                }
            }
        }
    }

    fn calc_gene_mutation(&self,
                          gff_map: &mut HashMap<String, Vec<bio::io::gff::Record>>,
                          reference: &mut bio::io::fasta::IndexedReader<File>,
                          codon_table: &CodonTable,
                          ref_name: &str,
                          output_prefix: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                all_variants,
                ..
            } => {
                let file_name = format!("{}_{}_dnds_values.gff",
                                        output_prefix.to_string(), ref_name);

                let file_path = Path::new(&file_name);
//                let mut file_open = File::create(file_path)
//                    .expect("Unable to create GFF file");
                let mut gff_writer = gff::Writer::to_file(file_path,
                                                          bio::io::gff::GffType::GFF3)
                                                          .expect("unable to create GFF file");
                for (tid, contig_name) in target_names.iter() {


                    let mut ref_sequence = Vec::new();
                    match reference.fetch_all(
                        std::str::from_utf8(contig_name.as_bytes()).unwrap()) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)
                        },
                    };
                    match reference.read(&mut ref_sequence) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)
                        },
                    };
                    let mut placeholder = Vec::new();
                    let gff_records = match gff_map.get_mut(contig_name){
                        Some(records) => records,
                        None => &mut placeholder,
                    };
                    debug!("Calculating population dN/dS from reads for {} genes on contig {}",
                           gff_records.len(), contig_name);

                    let placeholder_map = HashMap::new();
                    let variants = match all_variants.get(tid) {
                        Some(map) => map,
                        None => &placeholder_map
                    };

                    gff_records.iter_mut().enumerate().for_each(|(_id, gene)| {
                        let dnds = codon_table.find_mutations(gene, variants, &ref_sequence);
                        gene.attributes_mut().insert(format!("dNdS"), format!("{}", dnds));
                        debug!("gene {:?} attributes {:?}", gene.seqname(), gene);
                        gff_writer.write(gene);
//                        let attributes = if !gene.attributes().is_empty() {
//                            gene
//                                .attributes()
//                                .iter()
//                                .map(|(a, b)| format!("{}{}{}", a, "=", b))
//                                .join(";")
//                        } else {
//                            "".to_owned()
//                        };
//
//                        let gff_line = format!(
//                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
//                                            gene.seqname(),
//                                            gene.source(),
//                                            gene.feature_type(),
//                                            gene.start(),
//                                            gene.end(),
//                                            gene.score().unwrap_or(0),
//                                            gene.strand().unwrap_or(Strand::Unknown),
//                                            gene.frame(),
//                                            attributes);
//                        println!("{}", gff_line);
//                        writeln!(file_open, "{}", gff_line);
                    });

                }
            }
        }
    }

    fn write_vcf(&self, output_prefix: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                all_variants,
                target_names,
                target_lengths,
                sample_names,
                variant_info,
                ..
            } => {
                if variant_info.len() > 0 {
                    // initiate header
                    let mut header = bcf::Header::new();
                    // Add program info
                    header.push_record(format!("##source=lorikeet-v{}",
                                               env!("CARGO_PKG_VERSION")).as_bytes());

                    debug!("samples {:?}", &sample_names);
                    for sample in sample_names.iter() {
                        header.push_sample(&sample.clone().into_bytes()[..]);
                    }

                    // Add contig info
                    for (tid, contig_name) in target_names {
                        header.push_record(
                            format!("##contig=<ID={}, length={}>",
                                    contig_name, target_lengths[tid]).as_bytes()
                        );
                    }

                    // Add INFO flags
                    header.push_record(
                        format!("##INFO=<ID=TYPE,Number=A,Type=String,\
                    Description=\"The type of allele, either SNV, MNV, INS, DEL, or INV.\">").as_bytes());

                    header.push_record(
                        format!("##INFO=<ID=SVLEN,Number=1,Type=Integer,\
                    Description=\"Length of structural variant\">").as_bytes());

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
                        format!("##FORMAT=<ID=DP,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth in each sample\">",
                                sample_names.len()).as_bytes(),
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
                        format!("##FORMAT=<ID=QA,Number={},Type=Float,\
                            Description=\"Quality scores for allele in each sample\">",
                                sample_names.len()).as_bytes(),
                    );

                    let vcf_presort = tempfile::Builder::new()
                        .prefix("lorikeet-fifo")
                        .tempfile().expect("Unable to create VCF tempfile");

                    // Initiate writer
                    let mut bcf_writer = bcf::Writer::from_path(
                        format!("{}.vcf",
                                output_prefix).as_str(),
                        &header,
                        true,
                        bcf::Format::VCF).expect(
                        format!("Unable to create VCF output: {}.vcf",
                                output_prefix).as_str());

                    bcf_writer.set_threads(current_num_threads()).unwrap();


                    for (tid, position_variants) in all_variants.into_iter() {
                        let contig_name = &target_names[tid];
                        for (pos, variants) in position_variants.iter() {
                            for (variant, base) in variants.iter() {
                                // Create empty record for this variant
                                let mut record = bcf_writer.empty_record();
                                record.set_rid(
                                    Some(bcf_writer.header()
                                        .name2rid(contig_name.as_bytes()).unwrap()));
                                record.set_pos(*pos);

                                // Sum the quality scores across samples
                                let qual_sum = base.quals.iter().sum();
                                record.set_qual(qual_sum);

                                // Collect strain information
                                let mut strains = base.genotypes.iter().cloned()
                                    .collect::<Vec<i32>>();
                                strains.sort();

                                // Push info tags to record
                                record.push_info_integer(b"TDP", &[base.totaldepth.iter().sum()]);
                                record.push_info_integer(b"TAD", &[base.truedepth.iter().sum()]);
                                record.push_info_integer(b"TRD", &[base.referencedepth.iter().sum()]);
                                record.push_info_integer(b"ST", &strains[..]);

                                // Push format flags to record
                                record.push_format_integer(b"DP", &base.totaldepth[..]);
                                record.push_format_integer(b"AD", &base.truedepth[..]);
                                record.push_format_integer(b"RD", &base.referencedepth[..]);
                                record.push_format_float(b"QA", &base.quals[..]);

                                let refr = &base.refr[..];


                                match variant {
                                    Variant::SNV(alt) => {

                                        // Collect and set the alleles to record
                                        let alt = &[*alt];
                                        let mut collect_alleles = vec![refr, alt];
                                        record.set_alleles(&collect_alleles).unwrap();

                                        record.push_info_string(b"TYPE", &[b"SNV"]);

                                        bcf_writer.write(&record).expect("Unable to write record");
                                    },
                                    Variant::MNV(alt) => {

                                        // Collect and set the alleles to record
                                        let alt = [&refr[..], &alt[..]].concat();
                                        let mut collect_alleles = vec![refr, &alt];
                                        record.set_alleles(&collect_alleles).unwrap();

                                        record.push_info_string(b"TYPE", &[b"MNV"]);
                                        record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                        bcf_writer.write(&record).expect("Unable to write record");
                                    },
                                    Variant::Inversion(alt) => {
                                        // Collect and set the alleles to record
                                        let alt = [&refr[..], &alt[..]].concat();
                                        let mut collect_alleles = vec![refr, &alt];
                                        record.set_alleles(&collect_alleles).unwrap();

                                        record.push_info_string(b"TYPE", &[b"INV"]);
                                        record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                        bcf_writer.write(&record).expect("Unable to write record");
                                    },
                                    Variant::Insertion(alt) => {
                                        // Collect and set the alleles to record
                                        let alt = [&refr[..], &alt[..]].concat();
                                        let mut collect_alleles = vec![refr, &alt];
                                        record.set_alleles(&collect_alleles).unwrap();

                                        record.push_info_string(b"TYPE", &[b"INS"]);
                                        record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                        bcf_writer.write(&record).expect("Unable to write record");
                                    },
                                    Variant::Deletion(alt) => {
                                        // Collect and set the alleles to record
                                        // Create deletion variant, attach refr to head of N
                                        record.set_alleles(&vec![refr, format!(
                                            "{}N", refr[0] as char)
                                            .as_bytes()]).unwrap();

                                        record.push_info_string(b"TYPE", &[b"DEL"]);
                                        record.push_info_integer(b"SVLEN", &[*alt as i32]);

                                        bcf_writer.write(&record).expect("Unable to write record");
                                    },
                                    _ => {},
                                }
                            }
                        }
                    }

                    // Initiate sorting command
                    // Since variants are in HashMap, they are unsorted.
                    // Sorting using bcftools is fast so just do it for the user.
                    external_command_checker::check_for_bcftools();

                    // let command_string = format!(
                    //     "mv {} {}.vcf",
                    //     &vcf_presort
                    //         .path()
                    //         .to_str()
                    //         .expect("Failed to convert tempfile to path"),
                    //     output_prefix
                    // );
                    // This fails randomly, not sure what causes it so will keep unsorted for now
//                    let command_string = format!(
//                        "bcftools sort {} > {}.vcf",
//                        &vcf_presort
//                            .path()
//                            .to_str()
//                            .expect("Failed to convert tempfile to path"),
//                        output_prefix
//                    );
//                    if log_enabled!(Level::Debug) {
//                        vcf_presort.keep().unwrap();
//                    }

                    // command::finish_command_safely(
                    //     Command::new("bash")
                    //         .arg("-c")
                    //         .arg(&command_string)
                    //         .stderr(Stdio::piped())
                    //         .spawn()
                    //         .expect("Unable to execute bash"),
                    //     "bcftools"
                    // );
                    debug!("Finished writing VCF file for {}", &output_prefix);
                }
            }
        }
    }
}

/// Add read count entry to cluster hashmap
pub fn add_entry(shared_read_counts: &mut HashMap<usize, HashMap<usize, usize>>,
                 clust1: usize, clust2: usize, count: usize) {
    let entry = shared_read_counts.entry(clust1).or_insert(HashMap::new());
    entry.insert(clust2, count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use estimation::contig_variants::*;
    use std::collections::HashSet;
    use model::variants;

    fn create_base(ref_sequence: &Vec<u8>, var_char: u8, pos: i64, sample_count: usize) -> Base {
        Base {
            tid: 0,
            pos,
            refr: ref_sequence[pos as usize..(pos as usize + 1)].to_vec(),
            variant: Variant::SNV(var_char),
            filters: vec!(),
            depth: vec![0, 5],
            truedepth: vec![0, 5],
            totaldepth: vec![5, 5],
            genotypes: HashSet::new(),
            quals: vec![0.; sample_count],
            referencedepth: vec![5, 0],
            physicalcov: vec![0; sample_count],
            baseq: vec![0; sample_count],
            mapq: vec![0; sample_count],
            conf: vec![0; sample_count],
            nucs: HashMap::new(),
            pernucs: HashMap::new(),
            ic: vec![0; sample_count],
            dc: vec![0; sample_count],
            xc: vec![0; sample_count],
            ac: vec![0; sample_count],
            af: vec![0.; sample_count],
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
        var_stats.add_contig(Some(&mut variant_abundances), 0, 0,
                             b"test".to_vec(), ref_sequence.len(), 0,
                             vec![10., 10., 0.], ups_and_downs);
        var_mat.add_contig(var_stats, 2, 0);

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
        var_stats.add_contig(Some(&mut variant_abundances), 0, 0,
                             b"test".to_vec(), ref_sequence.len(), 1,
                             vec![10., 10., 0.], ups_and_downs);

        var_mat.add_contig(var_stats, 2, 1);

        var_mat.generate_distances(0, "test");

        var_mat.run_fuzzy_scan(0.01, 0.05, 0.01, 0.01,
                               0., 0, 0., 0)

    }
}