use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use estimation::contig_variants::*;
use estimation::linkage::*;
use model::variants::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::fs::File;
use dbscan::fuzzy;
use itertools::{Itertools};
use rust_htslib::bam::HeaderView;


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
        contigs: HashMap<i32, Vec<u8>>,
        target_names: HashMap<i32, String>,
        target_lengths: HashMap<i32, f64>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        clusters: HashMap<i32, HashMap<i32, BTreeMap<String, (i32, usize)>>>,
        clusters_mean: HashMap<i32, f64>,
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
            contigs: HashMap::new(),
            target_names: HashMap::new(),
            target_lengths: HashMap::new(),
            sample_names: vec!["".to_string(); sample_count],
            kfrequencies: BTreeMap::new(),
            clusters: HashMap::new(),
            clusters_mean: HashMap::new(),
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

    /// Returns the variants at the current position
    /// as a mutable reference
    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>>;

    fn variants_of_contig(&mut self, tid: i32) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>>;

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(&mut self,
                  variant_stats: VariantStats,
                  sample_count: usize,
                  sample_idx: usize,
                  contig: Vec<u8>);

    /// Converts all variants into fuzzy::Var format
    fn generate_distances(&mut self, threads: usize, output_prefix: &str);

    /// Perform fuzzy DBSCAN clustering using proportionality
    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64, phi: f64);

    /// Takes clusters from DBSCAN and linkage method and writes variants to file as genotype
    fn generate_genotypes(&mut self,
                          output_prefix: &str);

    fn print_variant_stats(&self, output_prefix: &str);


}

impl VariantMatrixFunctions for VariantMatrix {
    fn setup(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut contigs,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ref mut kfrequencies,
                ..
            } => {
                *coverages = HashMap::new();
                *average_genotypes = HashMap::new();
                *all_variants = HashMap::new();
                *contigs = HashMap::new();
                *target_names = HashMap::new();
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
                ref mut depths,
                ..
            } => {
                sample_names[sample_idx] = sample_name;
                let target_count = header.target_count();

                for tid in (0..target_count).into_iter() {
                    // Initialize contig id in variant hashmap
                    let contig_variants = all_variants.entry(tid as i32)
                        .or_insert(HashMap::new());
                    let placeholder = HashMap::new();
                    let variants = match variant_records.get(&(tid as i32)) {
                        Some(map) => map,
                        _ => &placeholder,
                    };
                    let target_len = header.target_len(tid).unwrap();
                    // Apppend the sample index to each variant abundance
                    // Initialize the variant position index
                    // Also turns out to be the total number of variant positions
                    for pos in (0..target_len).into_iter() {
                        let position_variants = contig_variants.entry(pos as i64)
                            .or_insert(HashMap::new());
                        let mut variant_depth: f64 = 0.;

                        if variants.contains_key(&(pos as i64)) {
                            let abundance_map = variants.get(&(pos as i64)).unwrap();
                            for (variant, base_info) in abundance_map.iter() {
//                                info!("{:?}", base_info);

                                let sample_map = position_variants.entry(variant.clone())
                                    .or_insert(base_info.clone());
                                sample_map.combine_sample(base_info, sample_idx, 0);
                            }
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
                  sample_idx: usize,
                  contig: Vec<u8>) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut contigs,
                ref mut target_names,
                ref mut target_lengths,
                ref mut variances,
                ref mut depths,
                ..
            } => {
                match variant_stats {
                    VariantStats::VariantContigStats {
                        tid,
                        target_name,
                        target_len,
                        coverage,
                        variance,
                        variants,
                        depth,
//                        variations_per_n,
                        ..
                    } => {
                        let ag = average_genotypes.entry(tid).or_insert(
                            vec![0.0 as f64; sample_count]);
//                        ag[sample_idx] = mean_genotypes;
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

                        // copy across depths
                        let mut contig_variants = all_variants.entry(tid)
                            .or_insert(HashMap::new());
                        for (pos, d) in depth.iter().enumerate() {
                            let position_variants = contig_variants.entry(pos as i64)
                                .or_insert(HashMap::new());
                            for (variant, base_info) in position_variants.iter_mut() {
                                base_info.add_depth(sample_idx, *d);
                            }
                        }
                        depths.entry(tid).or_insert(depth);
                        contigs.entry(tid).or_insert(contig);
                        
                    }
                }
            }
        }
    }

    fn generate_distances(&mut self, _threads: usize, _output_prefix: &str) {
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
                    let contig_coverages = coverages.get(tid)
                        .expect("Unable to retrieve contig coverage");

                    let _max_coverage = contig_coverages.iter().cloned().fold1(f64::max)
                        .expect("Unable to retrieve max coverage");

                    variant_abundances.par_iter_mut().for_each(
                        |(position, hash)| {
                            // loop through each position that has variants ignoring positions that
                            // only contained the reference in all samples
                            if hash.keys().len() > 0 {

                                for (variant, base_info) in hash.iter_mut() {
                                    match variant {
                                        Variant::None => {},
                                        _ => {
                                            let _abundance: f64 = 0.;
                                            let _mean_var: f64 = 0.;
                                            let mut rel_abund = vec![0.0; sample_count as usize];

                                            // Get the mean abundance across samples
                                            let mut sample_idx: usize = 0;
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
                                                    debug!("No var depth {:?} {:?}", base_info.variant, base_info.truedepth);
                                                    var_depth
                                                        = base_info.truedepth[index] as f64;
                                                    base_info.depth[index] = base_info.truedepth[index];
                                                }
                                                let mut total_depth
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
                                                rel_abunds: rel_abund,
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

                info!("geoms {:?} {:?} {:?}", geom_mean_d, geom_mean_v, geom_mean_f);

                let geom_mean_v = geom_mean_v.lock().unwrap().clone();
                let geom_mean_v = geom_mean(&geom_mean_v);
                debug!("Geom Mean Var {:?}", geom_mean_v);

                let geom_mean_d = geom_mean_d.lock().unwrap().clone();
                let geom_mean_d = geom_mean(&geom_mean_d);
                debug!("Geom Mean Dep {:?}", geom_mean_d);

                let geom_mean_f = geom_mean_f.lock().unwrap().clone();
                let geom_mean_f = geom_mean(&geom_mean_f);
                debug!("Geom Mean Frq {:?}", geom_mean_f);
                info!("geoms {:?} {:?} {:?}", geom_mean_d, geom_mean_v, geom_mean_f);

                *variant_info = variant_info_all;
                *geom_mean_var = geom_mean_v;
                *geom_mean_dep = geom_mean_d;
                *geom_mean_frq = geom_mean_f;

            }
        }
    }

    fn run_fuzzy_scan(&mut self, e_min: f64, e_max: f64, pts_min: f64, pts_max: f64, phi: f64) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut geom_mean_frq,
                ref mut pred_variants,
                all_variants,
                target_lengths,
                ..
            } => {

                info!("Running fuzzyDBSCAN with {} Variants", variant_info.len());
                if variant_info.len() == 1 {
                    info!("Where did the variants go? {:?}", variant_info);
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
                    geom_var: geom_mean_var.clone(),
                    geom_dep: geom_mean_dep.clone(),
                    geom_frq: geom_mean_frq.clone(),
                };

                // run fuzzy DBSCAN
                let mut clusters = fuzzy_scanner.cluster(&variant_info[..]);

                // Sort the clusters by smallest to largest
                clusters.sort_by(
                    |a, b| a.len().cmp(&b.len())
                );

                // Perform read phasing clustering and return new clusters
                // and shared read info between clusters
//                linkage_clustering_of_variants(&clusters,
//                                               &variant_info,
//                                               all_variants);
                let (clusters, shared_read_count, mut condensed)
                    = linkage_clustering_of_clusters(&clusters,
                                           &variant_info,
                                           all_variants);

                // Collapse clusters with enough shared read info starting with smallest cluster
//                if clusters.len() > 1 {
//                    for (cluster_index, cluster) in clusters.iter().enumerate() {
////                    println!("{}", cluster.len());
//                        let cluster_map = shared_read_count.get(&cluster_index)
//                            .unwrap();
//                        for (clust2, jaccard) in cluster_map {
//                            debug!("Cluster {} Cluster {} Jaccard's {}",
//                                   cluster_index+1, *clust2+1, jaccard);
//                        };
//                    };
//                };

//                let dend = linkage(&mut condensed,
//                                   clusters.len(), Method::Ward);
//                info!("Dendogram {:?}", dend);

                let mut genome_length = 0.;
                for (tid, length) in target_lengths {
                    genome_length += *length;
                };
//                let required_variants = (1. - 0.9999) * genome_length;
//                debug!("Genome Length {} Required Variants {}", genome_length, required_variants);
//
//                let mut clusters_kept = Vec::new();
//                let mut noise = Vec::new();
//
//                for (index, cluster) in clusters.iter().enumerate() {
//                    if cluster.len() as f64 >= required_variants {
//                        debug!("Cluster {} Variants {}", index+1, cluster.len());
//                        clusters_kept.push(cluster);
//                    } else {
//                        noise.par_extend(cluster.par_iter().cloned());
//                    };
//                };
//                clusters_kept.push(&noise);

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

                        let mut prediction_count = prediction_count.lock().unwrap();
                        let count = prediction_count.entry(rank + 1)
                            .or_insert(HashSet::new());
                        count.insert(assignment.index);


//                        *category += 1;
                        let variant = &variant_info[assignment.index];
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

                        let mut prediction_features = prediction_features.lock().unwrap();
                        let feature = prediction_features.entry(rank+1).or_insert(HashSet::new());
                        feature.insert(variant.var.to_owned());

                    });
                });

                let mut prediction_count = prediction_count.lock().unwrap();
                let prediction_features = prediction_features.lock().unwrap();
                let mut prediction_variants = prediction_variants.lock().unwrap();


                for (cluster, pred_set) in prediction_features.iter() {
                    if pred_set.len() > 1 {
                        let count = prediction_count[cluster].len() as f64;
                        info!("Cluster {} Sites {}", cluster, prediction_count[cluster].len());
                    } else if !pred_set.contains(&Variant::None) {
                        info!("Cluster {} Sites {}", cluster, prediction_count[cluster].len());
                    } else {
//                        prediction_variants.remove_entry(cluster);
                    }
                }
//                debug!("Prediction count {:?}", prediction_count);
//                debug!("Prediction categories {:?}", prediction_features);
                *pred_variants = prediction_variants.clone();
            }
        }
    }

    fn generate_genotypes(&mut self, output_prefix: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                contigs,
                ref mut pred_variants,
                ..
            } => {

                pred_variants.par_iter().for_each(|(strain_index, genotype)|{
                    let file_name = format!("{}_strain_{}.fna", output_prefix.to_string(), strain_index);

                    let file_path = Path::new(&file_name);

                    // Open haplotype file or create one
                    let mut file_open = File::create(file_path)
                        .expect("No Read or Write Permission in current directory");

                    let mut genotype = genotype.clone();
                    let mut multivariant_sites = 0;
                    let mut tot_variations = 0;

                    // Generate the variant genome
                    for (tid, original_contig) in contigs.iter() {
                        let mut contig = String::new();

                        let mut skip_n = 0;
                        let mut skip_cnt = 0;
//                            let char_cnt = 0;
                        let mut variations = 0;

                        for (pos, base) in original_contig.iter().enumerate() {
                            if skip_cnt < skip_n {
                                skip_cnt += 1;
                            } else {
                                skip_n = 0;
                                skip_cnt = 0;
                                if genotype.contains_key(&tid) {
                                    let tid_genotype = genotype.get_mut(&tid).unwrap();

                                    if tid_genotype.contains_key(&(pos as i64)) {
                                        let categories = &genotype[tid][&(pos as i64)];
                                        let mut hash = HashSet::new();
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
                                            multivariant_sites += 1;
                                            debug!("Multi hash {:?} {:?}", hash, max_var)
                                        }
                                        match max_var {
                                            Variant::Deletion(size) => {
                                                // Skip the next n bases but rescue the reference prefix
//                                                debug!("DEL {:?}", size);

                                                skip_n = size - 1;
                                                skip_cnt = 0;
//                                                let first_byte = max_var.as_bytes()[0];
//                                                contig = contig + str::from_utf8(
//                                                    &[first_byte]).unwrap();
                                                variations += 1;
                                            },
                                            Variant::Insertion(alt) => {
//                                                debug!("INS {:?}", alt);

                                                // Insertions have a reference prefix that needs to be removed
                                                let removed_first_base = str::from_utf8(
                                                    &alt[1..]).unwrap();
                                                contig = contig + removed_first_base;
                                                variations += 1;
                                            },
                                            Variant::None => {
                                                debug!("None variant {:?}", max_var);
                                                contig = contig + str::from_utf8(&[*base]).unwrap();
                                            },
                                            Variant::MNV(alt) => {
//                                                debug!("MNV {:?}", alt);
                                                skip_n = alt.len() as u32 - 1;
                                                skip_cnt = 0;
                                                let removed_first_base = str::from_utf8(
                                                    &alt[1..]).unwrap();
                                                contig = contig + removed_first_base;
                                                variations += 1;
                                            },
                                            Variant::SNV(alt) => {
//                                                debug!("SNV {:?}", alt);

                                                contig = contig + str::from_utf8(&[alt]).unwrap();
                                                variations += 1;
                                            },
                                            _ => {
                                                debug!("Unknown variant {:?}", max_var);
                                                contig = contig + str::from_utf8(&[*base]).unwrap();
                                            }
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
                        tot_variations += variations;
                    }
                    debug!("{} Multivariant sites and single variant sites {} for Strain {}",
                          multivariant_sites, tot_variations, strain_index);
                });
            }
        }
    }

    fn print_variant_stats(&self, output_prefix: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                sample_names,
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

/// Add read count entry to cluster hashmap
pub fn add_entry(shared_read_counts: &mut HashMap<usize, HashMap<usize, usize>>,
                 clust1: usize, clust2: usize, count: usize) {
    let entry = shared_read_counts.entry(clust1).or_insert(HashMap::new());
    entry.insert(clust2, count);
}