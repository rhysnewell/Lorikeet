use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use estimation::contig_variants::*;
use model::variants::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::fs::File;
use dbscan::fuzzy;
use kodama::{Method, linkage};
use itertools::{Itertools};


#[derive(Debug)]
/// Container for all variants within a genome and associated clusters
pub enum VariantMatrix {
    VariantContigMatrix {
        coverages: HashMap<i32, Vec<f64>>,
        average_genotypes: HashMap<i32, Vec<f64>>,
        variances: HashMap<i32, Vec<f64>>,
        // TID, Position, Base, Var Depth, Total Depth
        all_variants: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
        variants_map: HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>,
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
            variants_map: HashMap::new(),
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
                  variant_records: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>);

    /// Returns the variants at the current position as mutable reference
    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>>;

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

    /// Connects fuzzy DBSCAN clusters based on shared read information
    fn linkage_clustering(clusters: &Vec<Vec<fuzzy::Assignment>>,
                          variant_info: &Vec<fuzzy::Var>,
                          variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>)
                          -> (Vec<Vec<fuzzy::Assignment>>, HashMap<usize, HashMap<usize, f64>>, Vec<f64>) ;

    /// Get all of the associated read ids for a given cluster
    fn get_read_set(variants: &fuzzy::Cluster,
                    variant_info: &Vec<fuzzy::Var>,
                    variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>) -> HashSet<i64>;

    /// Extract the read ids associated with a particular variant
    fn get_variant_set(variant: &fuzzy::Var,
                       variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>) -> HashSet<i64>;

    fn print_variant_stats(&self, output_prefix: &str);


}

impl VariantMatrixFunctions for VariantMatrix {
    fn setup(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut variants_map,
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
                *variants_map = HashMap::new();
                *contigs = HashMap::new();
                *target_names = HashMap::new();
                *target_lengths = HashMap::new();
                *sample_names = vec!();
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn add_sample(&mut self, sample_name: String, sample_idx: usize,
                  variant_records: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut depths,
                ..
            } => {
                sample_names[sample_idx] = sample_name;

                for (tid, depth) in depths.iter() {
                    // Initialize contig id in variant hashmap
                    let contig_variants = all_variants.entry(*tid)
                        .or_insert(HashMap::new());
                    let placeholder = HashMap::new();
                    let variants = match variant_records.get(tid) {
                        Some(map) => map,
                        _ => &placeholder,
                    };

                    // Apppend the sample index to each variant abundance
                    // Initialize the variant position index
                    // Also turns out to be the total number of variant positions
                    for (pos, total_depth) in depth.iter().enumerate() {
                        let position_variants = contig_variants.entry(pos as i64)
                            .or_insert(HashMap::new());
                        let mut variant_depth: f64 = 0.;

                        if variants.contains_key(&(pos as i64)) {
                            let abundance_map = variants.get(&(pos as i64)).unwrap();
                            for (variant, base_info) in abundance_map.iter() {
                                let sample_map = position_variants.entry(variant.clone())
                                    .or_insert(base_info.clone());
                                sample_map.combine_sample(base_info, sample_idx, *total_depth);
                            }
                        }
                    }
                }
                *depths = HashMap::new();
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
                        return contig_variants.get_mut(&pos)
                    },
                    _ => return None,
                };
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
                ref mut variants_map,
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
                all_variants.par_iter().for_each(|(tid, variant_abundances)| {
                    let contig_coverages = coverages.get(tid)
                        .expect("Unable to retrieve contig coverage");

                    let _max_coverage = contig_coverages.iter().cloned().fold1(f64::max)
                        .expect("Unable to retrieve max coverage");

                    variant_abundances.par_iter().for_each(
                        |(position, hash)| {
                            // loop through each position that has variants ignoring positions that
                            // only contained the reference in all samples
                            if hash.keys().len() > 0 {

                                for (variant, base_info) in hash.iter() {

                                    match variant {
                                        Variant::None => {},
                                        _ => {
                                            let _abundance: f64 = 0.;
                                            let _mean_var: f64 = 0.;

                                            let mut rel_abund = vec![0.0; sample_count as usize];

                                            // Get the mean abundance across samples
                                            let mut sample_idx: usize = 0;
                                            (0..sample_count).into_iter().for_each(|index| {
                                                let mut geom_mean_v =
                                                    geom_mean_v.lock().unwrap();
                                                let mut geom_mean_d =
                                                    geom_mean_d.lock().unwrap();
                                                let mut geom_mean_f =
                                                    geom_mean_f.lock().unwrap();


                                                let var_depth
                                                    = base_info.depth[index] as f64 + 1.;
                                                let total_depth
                                                    = base_info.totaldepth[index] as f64 + 1.;
                                                rel_abund[index] =
                                                    var_depth / total_depth;
                                                geom_mean_v[index] += (var_depth).ln();
                                                geom_mean_d[index] += (total_depth).ln();
                                                geom_mean_f[index] += (var_depth
                                                    / total_depth).ln();
                                            });

                                            let mut variant_info_all = variant_info_all
                                                .lock().unwrap();
                                            let point = fuzzy::Var {
                                                pos: *position,
                                                var: variant.clone(),
                                                deps: base_info.totaldepth.clone(),
                                                vars: base_info.depth.clone(),
                                                rel_abunds: rel_abund,
                                                tid: *tid,
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

                let geom_mean_v = geom_mean_v.lock().unwrap().clone();
                let geom_mean_v = geom_mean(&geom_mean_v);
                debug!("Geom Mean Var {:?}", geom_mean_v);

                let geom_mean_d = geom_mean_d.lock().unwrap().clone();
                let geom_mean_d = geom_mean(&geom_mean_d);
                debug!("Geom Mean Dep {:?}", geom_mean_d);

                let geom_mean_f = geom_mean_f.lock().unwrap().clone();
                let geom_mean_f = geom_mean(&geom_mean_f);
                debug!("Geom Mean Frq {:?}", geom_mean_f);

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
//                variant_map,
                target_lengths,
                ..
            } => {

                info!("Running fuzzyDBSCAN with {} Variants", variant_info.len());
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
//                let (clusters, shared_read_count, mut condensed)
//                    = Self::linkage_clustering(&clusters,
//                                               &variant_info,
//                                               &snps_map,
//                                               &indels_map);

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
                        prediction_variants.remove_entry(cluster);
                    }
                }
//                debug!("Prediction count {:?}", prediction_count);
                debug!("Prediction categories {:?}", prediction_features);
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
                                            match max_var {

                                                Variant::Deletion(size) => {
                                                    // Skip the next n bases but rescue the reference prefix
                                                    skip_n = size - 1;
                                                    skip_cnt = 0;
    //                                                let first_byte = max_var.as_bytes()[0];
    //                                                contig = contig + str::from_utf8(
    //                                                    &[first_byte]).unwrap();
                                                    variations += 1;
                                                },
                                                Variant::Insertion(alt) | Variant::Insertion(alt) => {
                                                    // Insertions have a reference prefix that needs to be removed
                                                    let removed_first_base = str::from_utf8(
                                                        &alt[1..]).unwrap();
                                                    contig = contig + removed_first_base;
                                                    variations += 1;
                                                },
                                                Variant::None => {
                                                    contig = contig + str::from_utf8(&[*base]).unwrap();
                                                },
                                                Variant::SNV(alt) => {
                                                    contig = contig + str::from_utf8(&[alt]).unwrap();
                                                    variations += 1;
                                                },
                                                _ => {
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
                        }
                    }
                });
            }
        }
    }

    /// Connects fuzzy DBSCAN clusters based on shared read information
    fn linkage_clustering(clusters: &Vec<Vec<fuzzy::Assignment>>,
                          variant_info: &Vec<fuzzy::Var>,
                          variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>)
                          -> (Vec<Vec<fuzzy::Assignment>>, HashMap<usize, HashMap<usize, f64>>, Vec<f64>) {

        if clusters.len() > 1 {
            let clusters_changed = Arc::new(Mutex::new(Vec::new()));
            let clusters_shared_reads = Arc::new(Mutex::new(HashMap::new()));
            let jaccard_distances = Arc::new(Mutex::new(
                vec![0.; (clusters.len().pow(2) - clusters.len()) / 2])
            );
            // Loop through each permutation of 2 clusters and observe shared variants in reads
            (0..clusters.len()).into_iter()
                .permutations(2)
                .collect::<Vec<Vec<usize>>>().into_par_iter().for_each(|(indices)| {

                // Get clusters by index
                let clust1 = &clusters[indices[0]];
                let clust2 = &clusters[indices[1]];

                // Total read set for each cluster
                let clust1_set: Arc<Mutex<HashSet<i64>>> = Arc::new(Mutex::new(HashSet::new()));
                let clust2_set: Arc<Mutex<HashSet<i64>>> = Arc::new(Mutex::new(HashSet::new()));

                // Cluster assignment index to avoid reappending to sets
                let clust1_index = Arc::new(Mutex::new(HashSet::new()));
                let clust2_index = Arc::new(Mutex::new(HashSet::new()));

                // Loop through first cluster
                clust1.par_iter().for_each(|assignment1| {
                    // Loop through second cluster
                    clust2.par_iter().for_each(|assignment2| {
                        if assignment1.index != assignment2.index {
                                let var1 = &variant_info[assignment1.index];
                                let var2 = &variant_info[assignment2.index];

                                // Read ids of first variant
                                let set1 = Self::get_variant_set(var1,
                                                                 variant_map);

                                // Read ids of second variant
                                let set2 = Self::get_variant_set(var2,
                                                                 variant_map);

                                // Extend each cluster read set
                                {
                                    let mut clust1_index =
                                        clust1_index.lock().unwrap();
                                    if !clust1_index.contains(&assignment1.index) {
                                        let mut clust1_set =
                                            clust1_set.lock().unwrap();
                                        clust1_set.par_extend(&set1);

                                        clust1_index.insert(assignment1.index);
                                    };
                                    let mut clust2_index =
                                        clust2_index.lock().unwrap();
                                    if !clust2_index.contains(&assignment2.index) {
                                        let mut clust2_set =
                                            clust2_set.lock().unwrap();
                                        clust2_set.par_extend(&set2);

                                        clust2_index.insert(assignment2.index);
                                    };
                                }
                                // Intersection of two read sets
//                            let intersection: BTreeSet<_> = set1.intersection(&set2).collect();
//
//                            if intersection.len() >= 1 {
//
//                                let mut clusters_changed =
//                                   clusters_changed.lock().unwrap();
//
//                                clusters_changed[indices[0]].push(assignment2.clone());
//                                clusters_changed[indices[1]].push(assignment1.clone());
//
//                                let mut clust1_set =
//                                    clust1_set.lock().unwrap();
//                                clust1_set.par_extend(&set2);
//
//                                let mut clust2_set =
//                                    clust2_set.lock().unwrap();
//                                clust2_set.par_extend(&set1);
//
//                           }
                        }
                    });
                });

                // Add the jaccard's similarity to the hashmap for the two clusters
                let clust1_set = clust1_set.lock().unwrap();
                let clust2_set = clust2_set.lock().unwrap();
                let intersection: HashSet<_> = clust1_set
                    .intersection(&clust2_set).collect();
//            let indices_set: HashSet<_> = [indices[0], indices[1]].iter().cloned().collect();
                let mut clusters_shared_reads = clusters_shared_reads
                    .lock().unwrap();
                let mut cluster_map = clusters_shared_reads.entry(indices[0])
                    .or_insert(HashMap::new());
//                let jaccard = intersection.len() as f64 /
//                    ((clust1_set.len() + clust2_set.len() - intersection.len()) as f64);
                let jaccard = intersection.len() as f64 /
                    std::cmp::min(clust1_set.len(), clust2_set.len()) as f64;
                let mut jaccard_distances = jaccard_distances.lock().unwrap();
                // Check to see if we have two or more samples
                if jaccard_distances.len() > 1 {
                    // Place each measure at appropriate index
                    jaccard_distances[
                        get_condensed_index(
                            indices[0], indices[1], clusters.len()).unwrap_or_else(|| 0)]
                        = 1. - jaccard;
                } else {
                    // If only 1 sample, then just put in first index
                    jaccard_distances[0] = 1. - jaccard;
                }
                // Add jaccard similarity to hashmap
                cluster_map.entry(indices[1]).or_insert(jaccard);
            });

            let mut clusters_shared_reads = clusters_shared_reads.lock().unwrap().clone();
            let mut jaccard_distances = jaccard_distances.lock().unwrap().clone();
            // Perform HAC using kodama
            let dend = linkage(&mut jaccard_distances,
                               clusters.len(), Method::Ward);

            info!("Dendogram {:?}", &dend);
            let changed = Arc::new(Mutex::new(0));
            // Step through each step in the dendrogram and combine clusters
            // that are under a certain dissimilarity
            dend.steps().into_par_iter().for_each(|step| {
                // Check to see that these are leaf clusters
                if step.cluster1 <= clusters.len() - 1 && step.cluster2 <= clusters.len() - 1 {
                    // combine clusters
                    if step.dissimilarity < 0.2 {
                        let mut new_cluster = Vec::new();
                        new_cluster.par_extend(clusters[step.cluster1].par_iter().cloned());
                        new_cluster.par_extend(clusters[step.cluster2].par_iter().cloned());

                        let mut clusters_changed
                            = clusters_changed.lock().unwrap();
//                        if clusters[step.cluster1].len() >= clusters[step.cluster2].len() {
//                            clusters_changed.push(clusters[step.cluster1].clone());
//                        } else {
//                            clusters_changed.push(clusters[step.cluster2].clone());
//                        }
                        let mut changed = changed.lock().unwrap();
                        *changed += 1;
                        clusters_changed.push(new_cluster);
                    } else { // cluster is by itself
                        let mut clusters_changed
                            = clusters_changed.lock().unwrap();
                        clusters_changed.push(clusters[step.cluster1].clone());
                        clusters_changed.push(clusters[step.cluster2].clone());
                    }
                // Check individually for leaf clusters
                } else if step.cluster1 <= clusters.len() - 1 {
                    let mut clusters_changed
                        = clusters_changed.lock().unwrap();
                    clusters_changed.push(clusters[step.cluster1].clone());
                } else if step.cluster2 <= clusters.len() - 1 {
                    let mut clusters_changed
                        = clusters_changed.lock().unwrap();
                    clusters_changed.push(clusters[step.cluster2].clone());
                }
            });

            let mut clusters_changed
                = clusters_changed.lock().unwrap().clone();
            let changed = changed.lock().unwrap();
            // If the number of clusters changed, then we rerun linkage clustering
            // First use of recursion properly, nice.
            if *changed > 0 {
                let (clusters_changed, clusters_shared_reads, jaccard_distances)
                    = Self::linkage_clustering(&clusters_changed,
                                               variant_info,
                                               variant_map);
                return (clusters_changed, clusters_shared_reads, jaccard_distances)
            } else {
                return (clusters.clone(), clusters_shared_reads, jaccard_distances)
            }
        } else {
            // create placeholder jaccard hashmap when there is only one cluster to prevent nothing
            // from being returned
            let mut jaccard_map = HashMap::new();
            let mut input = HashMap::new();
            input.insert(1, 1.0);
            jaccard_map.insert(0, input);

            return (clusters.clone(), jaccard_map, vec![1.0])
        }
    }

    /// Get all of the associated read ids for a given cluster
    fn get_read_set(variants: &fuzzy::Cluster,
                    variant_info: &Vec<fuzzy::Var>,
                    variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>) -> HashSet<i64> {

        let read_set = Arc::new(Mutex::new(HashSet::new()));

        variants.par_iter().for_each(|assignment|{
            let variant = &variant_info[assignment.index];
            let variant_set = &variant_map[&variant.tid][&variant.pos][&variant.var];
            variant_set.par_iter().for_each(|id|{
                let mut read_set = read_set.lock().unwrap();
                read_set.insert(*id);
            });
        });
        let read_set = read_set.lock().unwrap().clone();
        read_set
    }

    /// Extract the read ids associated with a particular variant
    fn get_variant_set(variant: &fuzzy::Var,
                       variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, HashSet<i64>>>>) -> HashSet<i64> {

        let mut variant_set = HashSet::new();

        variant_set = variant_map[&variant.tid][&variant.pos][&variant.var].clone();

        return variant_set
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

/// helper function to get the index of condensed matrix from it square form
fn get_condensed_index(i: usize, j: usize, n: usize) -> Option<usize>{
    if i == j {
        return None
    } else {
        return Some(n*i - i*(i+1)/2 + j - 1 - i)
    }
}