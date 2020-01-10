use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use pileup_structs::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use ndarray::{Array2, Array1, ArrayView};
use cogset::{Euclid, Dbscan, BruteScan};
use kodama::{Method, nnchain, Dendrogram};
use std::sync::{Arc, Mutex, MutexGuard};
use haplotypes_and_genotypes::*;
use ndarray_npy::write_npy;
use std::fs::File;
use std::process;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use tempfile;


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

    fn generate_distances(&mut self);

//    fn dbscan_cluster(&mut self, eps: f64, min_cluster_size: usize);

    fn generate_genotypes(&mut self, output_prefix: &str);

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
                            // Get the geometric means of the variant, depth, and reference counts
                            // at each variant position
//                            let var_geom: f32 = contig_sums[0].iter().product::<f32>()
//                                .powf((1 / variant_index) as f32);
//                            let dep_geom: f32 = contig_sums[1].iter().product::<f32>()
//                                .powf((1 / variant_index) as f32);
//                            let ref_geom: f32 = contig_sums[2].iter().product::<f32>()
//                                .powf((1 / variant_index) as f32);
//
//                            debug!("Ref CLR {:?}", contig_sums[2]);
//
//                            contig_sums[0] = contig_sums[0].iter()
//                                .map(|var| { (*var / var_geom).ln() }).collect();
//                            contig_sums[1] = contig_sums[1].iter()
//                                .map(|dep| { (*dep / dep_geom).ln() }).collect();
//                            contig_sums[2] = contig_sums[2].iter()
//                                .map(|refr| { (*refr / ref_geom).ln() }).collect();

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

    fn generate_distances(&mut self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                indels_map,
                snps_map,
                target_names,
                sample_names,
                coverages,
                ..
            } => {
                let sample_count = sample_names.len() as f32;


                let variant_info_all =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                // get basic variant info
                variants.par_iter().for_each(|(tid, variant_abundances)|{
                    variant_abundances.par_iter().for_each(
                        |(position, hash)| {
                        // loop through each position that has variants



                        for (var, abundances_vector) in hash.iter() {
                            let mut abundance: f32 = 0.;
                            let mut mean_var: f32 = 0.;
                            let mut total_d: f32 = 0.;
                            let mut freqs = Vec::new();
                            if !var.contains("R") {
                                // Get the mean abundance across samples
                                abundances_vector.iter().for_each(|(var, d)| {
                                    mean_var += *var;
                                    // Total depth of location
                                    total_d += *d;
                                    freqs.push(*var / *d);
                                    abundance += *var / *d;
                                });

                                mean_var = mean_var / sample_count;
                                abundance = abundance / sample_count;

                                let mut variant_info_all = variant_info_all
                                    .lock().unwrap();


                                variant_info_all.push(
                                    (position, var.to_string(),
                                     (abundance, freqs), tid));
                            }
                        }
                    });
                });

                let variant_info_all = variant_info_all.lock().unwrap();

                info!("Generating Variant Distances with {} Variants", variant_info_all.len());

                let mut variant_distances =
                    get_condensed_distances(&variant_info_all, indels_map, snps_map);

                let mut variant_distances = variant_distances
                    .lock()
                    .unwrap();

//                let strings: Vec<String> = variant_distances.iter().map(|n| n.to_string()).collect();

                let tmp_dir = TempDir::new("lorikeet_fifo")
                    .expect("Unable to create temporary directory");
                let fifo_path = tmp_dir.path().join("foo.pipe");

                // create new fifo and give read, write and execute rights to the owner.
                // This is required because we cannot open a Rust stream as a BAM file with
                // rust-htslib.
                unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
                    .expect(&format!("Error creating named pipe {:?}", fifo_path));

                let mut distances_file = tempfile::Builder::new()
                    .prefix("lorikeet-distances-vec")
                    .tempfile_in(tmp_dir.path())
                    .expect(&format!("Failed to create distances tempfile"));
//                writeln!(distances_file, "{:?}", variant_distances).expect("Unable to write to tempfile");
                let tmp_path = distances_file.path().to_str()
                    .expect("Failed to convert tempfile path to str").to_string();
                write_npy(&tmp_path, variant_distances.to_owned()).expect("Unable to write to tempfile");

//                println!("{:?}", variant_distances);
                let max_rank = 10;
                let mut ranks_rss = Arc::new(Mutex::new(vec![0.; max_rank]));

                (0..max_rank).into_par_iter().for_each(|rank| {
                    let cmd_string = format!(
                        "set -e -o pipefail; \
                     nmf.py {} True {} {}",
                        // NMF
                        rank+1,
                        100,
                        tmp_path);
                    info!("Queuing cmd_string: {}", cmd_string);
                    let mut python = std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&cmd_string)
                        .stderr(process::Stdio::piped())
                        .stdout(process::Stdio::piped())
                        .spawn()
                        .expect("Unable to execute bash");

                    let es = python.wait().expect("Unable to discern exit status");
                    if !es.success() {
                        error!("Error when running NMF: {:?}", cmd_string);
                        let mut err = String::new();
                        python.stderr.expect("Failed to grab stderr from NMF")
                            .read_to_string(&mut err).expect("Failed to read stderr into string");
                        error!("The overall STDERR was: {:?}", err);

                        process::exit(1);
                    } else {
                        let mut out = String::new();
                        python.stdout.expect("Failed to grab stdout from NMF").read_to_string(&mut out)
                            .expect("Failed to read stdout to string");
                        let mut ranks_rss = ranks_rss.lock().expect("Unable to lock RSS vec");
                        let rss: f32 = out.trim().parse().unwrap();
                        ranks_rss[rank as usize] = rss;
                    }
                });

                let ranks_rss = ranks_rss.lock().expect("unable to lock rss vec");
                let mut best_rank = 0;
                let mut best_rss = 0.;
                debug!("RSS Values {:?}", ranks_rss);

                for (rank, rss) in ranks_rss.iter().enumerate() {
                    if best_rank == 0 && best_rss == 0. && rank == 0 {
                        best_rank = rank + 1;
                        best_rss = *rss;
                    } else if &best_rss > rss {
                        best_rss = *rss;
                        best_rank = rank + 1;
                    } else if rss >= &best_rss {
                        break
                    }
                }

                let cmd_string = format!(
                    "set -e -o pipefail; \
                     nmf.py {} False {} {}",
                    // NMF
                    best_rank,
                    100,
                    tmp_path);
                info!("Queuing cmd_string: {}", cmd_string);
                let mut python = std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&cmd_string)
                    .stderr(process::Stdio::piped())
                    .stdout(process::Stdio::piped())
                    .spawn()
                    .expect("Unable to execute bash");

                let es = python.wait().expect("Unable to discern exit status");
                if !es.success() {
                    error!("Error when running NMF: {:?}", cmd_string);
                    let mut err = String::new();
                    python.stderr.expect("Failed to grab stderr from NMF")
                        .read_to_string(&mut err).expect("Failed to read stderr into string");
                    error!("The overall STDERR was: {:?}", err);

                    process::exit(1);
                } else {
                    let mut out = String::new();
                    python.stdout.expect("Failed to grab stdout from NMF").read_to_string(&mut out)
                        .expect("Failed to read stdout to string");
                    println!("{}", out);
                }


                tmp_dir.close().expect("Unable to close temp directory");


//                let py = gil.python();
//                let nmfpy = PyModule::import(py, "nmf.py").unwrap();
//                let variant_distances = variant_distances.to_owned().into_pyarray(py);
//                nmfpy.call1("perform_nmf", (variant_distances,)).unwrap();
            }
        }
    }

//    fn dbscan_cluster(&mut self, eps: f64, min_cluster_size: usize) {
//        match self {
//            PileupMatrix::PileupContigMatrix {
//                variants,
//                indels_map,
//                snps_map,
//                ref mut clusters,
//                ref mut clusters_mean,
//                ref mut dendrogram,
//                target_names,
//                sample_names,
//                coverages,
//                ..
//            } => {
//                let sample_count = sample_names.len() as f32;
//
//                let mut abundance_euclid =
//                    Arc::new(
//                        Mutex::new(
//                            Vec::new()));
//
//                let mut abundance_float =
//                    Arc::new(
//                        Mutex::new(
//                            Vec::new()));
//
//                let variant_info =
//                    Arc::new(
//                        Mutex::new(
//                            Vec::new()));
//
//                let variant_info_all =
//                    Arc::new(
//                        Mutex::new(
//                            Vec::new()));
//
//                let mut cluster_map =
//                    Arc::new(
//                        Mutex::new(
//                            HashMap::new()));
//
//                let mut cluster_hierarchies =
//                    Arc::new(
//                        Mutex::new(
//                            HashMap::new()));
//
//                let mut noise_set = Arc::new(
//                    Mutex::new(
//                        HashSet::new()));
//
//
//                for (tid, variant_abundances) in variants.iter() {
//                    variant_abundances.iter().for_each(
//                        |(position, hash)| {
//                            // loop through each position that has variants
//                            let mut abundance_euclid = abundance_euclid
//                                .lock().unwrap();
//                            let mut abundance_float = abundance_float
//                                .lock().unwrap();
//                            let mut variant_info = variant_info
//                                .lock().unwrap();
//                            let mut variant_info_all = variant_info_all
//                                .lock().unwrap();
//
//                            for (var, abundances_vector) in hash.iter() {
//                                let mut abundance: f32;
//                                let mut mean_var: f32 = 0.;
//                                let mut mean_d: f32 = 0.;
//                                if !var.contains("R") {
//                                    // Get the mean abundance across samples
//                                    // Get the mean abundance across samples
//                                    abundances_vector.iter().for_each(|(var, d)| {
//                                        mean_var += *var;
//                                        // Total depth of location
//                                        total_d += *d;
//                                        freq = *var / *d;
//                                        abundance += *var / *d;
//                                    });
//
//                                    mean_var = mean_var / sample_count;
//                                    abundance = abundance / sample_count;
//
//                                    let mut variant_info_all = variant_info_all
//                                        .lock().unwrap();
//
//                                    variant_info_all.push(
//                                        (position, var.to_string(),
//                                         (abundance, freqs), tid));
//
//
//                                    // This is a big hack but saves so much time
//                                    // Basically, any variant that has an abundance less than 0.05
//                                    // Automatically gets thrown into bin 0
//                                    // Downsides:
//                                    // Variants on border of eps look like they should cluster
//                                    // with bin 0, but instead create a new bin
//                                    // Upsides:
//                                    // Massive speed increase, massive decrease in memory usage
//                                    // Low abundant variants are all kind of in the same cluster any way
//                                    if mean_d >= 0. {
//                                        variant_info.push((position, var.to_string(), abundance, tid));
//                                        abundance_float.push(abundance);
//                                        abundance_euclid.push(Euclid([abundance.ln() as f64, mean_var.ln() as f64]));
//                                    } else {
//                                        let mut noise_set = noise_set.lock().unwrap();
//                                        noise_set.insert(variant_info_all.len());
//                                    }
//
//                                }
//                            }
//                        });
//                }
//
//                let abundance_euclid = abundance_euclid.lock().unwrap();
//                let variant_info = variant_info.lock().unwrap();
//                let variant_info_all = variant_info_all.lock().unwrap();
//                let abundance_float = abundance_float.lock().unwrap();
//                let scanner = BruteScan::new(&abundance_euclid);
//                if variant_info_all.len() > 0 {
//                    debug!("Beginning clustering of {} variants out of {}", abundance_euclid.len(),
//                           variant_info_all.len());
//                    let mut dbscan = Dbscan::new(scanner,
//                                                 eps as f64,
//                                                 min_cluster_size);
//
//
//                    let db_clusters = dbscan.by_ref().collect::<Vec<_>>();
//                    debug!("{} True clusters found for abundance > eps", db_clusters.len());
//
//                    let noise_points = dbscan.noise_points()
//                        .iter().cloned().collect::<Vec<_>>();
//
//                    // get the number of clusters so we can rescue orphan variants
//                    let mut number_of_clusters =
//                        Arc::new(
//                            Mutex::new(
//                                db_clusters.len() as i32));
//
//                    debug!("Sorting DBSCAN Clusters");
//
//                    db_clusters.par_iter().enumerate().for_each(|(cluster, index_vec)| {
//                        // Deal with clustered points
//                        index_vec.par_iter().for_each(|index| {
//                            let info = &variant_info[*index];
//                            let abundance = abundance_float[*index];
//
//                            let mut cluster_map = cluster_map.lock().unwrap();
//                            let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
//                            let mut number_of_clusters = number_of_clusters.lock().unwrap();
//
//                            let contig = cluster_map.entry(*info.3)
//                                .or_insert(HashMap::new());
//
//                            let position = contig
//                                .entry(*info.0).or_insert(BTreeMap::new());
//                            position.insert(info.1.to_string(), cluster as i32);
//
//                            let cluster_mean = cluster_hierarchies
//                                .entry(cluster as i32).or_insert(Vec::new());
//                            cluster_mean.push(abundance.to_owned());
//                        });
//                    });
//
//                    noise_points.par_iter().for_each(|noise| {
//                        let mut cluster_map = cluster_map.lock().unwrap();
//                        let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
//                        let mut number_of_clusters = number_of_clusters.lock().unwrap();
//
//
//                        // Deal with unclustered points by finding which cluster they are closest to
//                        let noise_info = &variant_info[*noise];
//                        let noise_abundance = noise_info.2;
//                        let mut curr_diff = 1.0;
//                        let mut curr_clus = 0;
//                        for (cluster, abundance) in cluster_hierarchies.iter() {
//                            let mean = abundance.iter().sum::<f32>() / abundance.len() as f32;
//
//                            let diff = (mean - noise_abundance).abs();
//                            if diff < curr_diff {
//                                curr_diff = diff;
//                                curr_clus = *cluster;
//                            }
//                        }
//
//                        let contig_map = cluster_map.entry(*noise_info.3)
//                            .or_insert(HashMap::new());
//
//                        let noise_position = contig_map.entry(*noise_info.0)
//                            .or_insert(BTreeMap::new());
//
//                        noise_position.insert(noise_info.1.to_string(), curr_clus);
//
//                        let noise_mean = cluster_hierarchies
//                            .entry(curr_clus).or_insert(Vec::new());
//                        noise_mean.push(noise_abundance.to_owned());
//                    });
//
//                    let noise_set = noise_set.lock().unwrap();
//
//                    noise_set.par_iter().for_each(|noise| {
//                        let mut cluster_map = cluster_map.lock().unwrap();
//                        let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
//                        let mut number_of_clusters = number_of_clusters.lock().unwrap();
//
//
//                        // Deal with unclustered points by finding which cluster they are closest to
//                        let noise_info = &variant_info_all[*noise];
//                        let noise_abundance = noise_info.2;
//                        let noise_frac = noise_abundance.0 / noise_abundance.1;
//                        let mut curr_diff = 1.0;
//                        let mut curr_clus = 0;
//
//                        for (cluster, abundance) in cluster_hierarchies.iter() {
//                            let mean = abundance.iter().sum::<f32>() / abundance.len() as f32;
//
//                            let diff = (mean - noise_frac).abs();
//                            if diff < curr_diff {
//                                curr_diff = diff;
//                                curr_clus = *cluster;
//                            }
//                        }
//
//                        let contig_map = cluster_map.entry(*noise_info.3)
//                            .or_insert(HashMap::new());
//
//                        let noise_position = contig_map.entry(*noise_info.0)
//                            .or_insert(BTreeMap::new());
//
//                        noise_position.insert(noise_info.1.to_string(), curr_clus);
//
//                        let noise_mean = cluster_hierarchies
//                            .entry(curr_clus).or_insert(Vec::new());
//                        noise_mean.push(noise_frac.to_owned());
//                    });
//
//                    let mut variant_distances =
//                        get_condensed_distances(&variant_info_all, indels_map, snps_map);
//
//                    let mut contig_coverage_means = HashMap::new();
//                    coverages.iter().map(|(tid, cov_vec)| {
//                        let contig = contig_coverage_means.entry(tid).or_insert(cov_vec.iter().sum::<f32>());
//                    });
//
//                    let mut variant_distances = variant_distances
//                        .lock()
//                        .unwrap();
//                    debug!("Performing HAC");
//
//
////                    let dend = nnchain(
////                        variant_distances.as_slice().expect("Unable to convert array to slice"),
////                        variant_info_all.len(),
////                        Method::Ward.into_method_chain()
////                            .expect("Incompatible linkage method"));
////                    debug!("Dendrogram {:?}", dend);
////
////                    let mut cluster_map = cluster_map.lock().unwrap();
////                    let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
////                    let mut means = HashMap::new();
////
////                    for (cluster, abundance) in cluster_hierarchies.iter() {
////                        let mean = abundance.iter().sum::<f32>() / abundance.len() as f32;
////
////                        means.insert(*cluster as i32, mean);
////                    }
////
////                    // Combine info of both clustering methods for easy access
////                    let mut full_cluster_map =
////                        Arc::new(
////                            Mutex::new(
////                                HashMap::new()));
////
////                    variant_info_all
////                        .par_iter().enumerate().for_each(|(index, (position, variant, abundance, tid))| {
////                        let db_cluster = cluster_map
////                            .get(tid).expect("No contig found when it should be")
////                            .get(position).expect("Position not found when it should be")
////                            .get(variant).expect("Variant not found when it should be");
////                        let mut full_cluster_map = full_cluster_map.lock().unwrap();
////
////                        let contig_map = full_cluster_map.entry(**tid)
////                            .or_insert(HashMap::new());
////
////                        let position_map = contig_map.entry(**position)
////                            .or_insert(BTreeMap::new());
////
////                        position_map.insert(variant.to_string(), (*db_cluster, index));
////                    });
////
////
////                    *dendrogram = dend;
////                    *clusters = full_cluster_map.lock().unwrap().clone();
////                    *clusters_mean = means;
////                    info!("{} Distinct variant frequency clusters found on {} contigs at eps {}",
////                          clusters_mean.len(), target_names.len(), eps);
//                }
//            }
//        }
//    }

    fn generate_genotypes(&mut self, output_prefix: &str) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut variants,
                ref mut clusters,
                ref mut clusters_mean,
                contigs,
                target_names,
                dendrogram,
                sample_names,
                ..
            } => {
                let sample_count = sample_names.len() as f32;
                // First we need to convert DBSCAN cluster ids into taxonomic rank ids
                // Clusters of higher abundance mutations will be closer to root
                // Ordered from root to leaf
                // We can then check how dissimilar two variants are on the hierarchical cluster
//                let mut ordered_clusters: Vec<_> = clusters_mean.iter().collect();
//                ordered_clusters
//                    .sort_by(|a, b|
//                        b.1.partial_cmp(a.1).unwrap());
//                debug!("{:?}", ordered_clusters);

                // Clusters is set up as HashMap<tid, HashMap<Position, BTreeMap<Variant, (Cluster_ID, Dendro_index)>>>
                // In this case we want to rearrange to HashMap<dendro_ID, HashMap<Position, (Variant, db_clust, tid)>>
                // This will allow us to disentangle positions where more than one variant is possible
                if dendrogram.len() > 0 {
                    let mut haplotypes_vec: Vec<Haplotype> = vec!();

                    debug!("Beginning haplotyping of dendrogram of length: {}", dendrogram.len());
                    let mut dendro_ids = Arc::new(Mutex::new(HashMap::new()));
                    for (tid, cluster_tid) in clusters.iter() {
                        cluster_tid.par_iter().for_each(|(position, variant_map)| {
                            for (variant, cluster) in variant_map.iter() {
                                let mut dendro_ids = dendro_ids.lock().unwrap();
                                let clust = dendro_ids.entry(cluster.1)
                                    .or_insert(BTreeMap::new());

                                clust.entry(*position)
                                    .or_insert((variant.to_string(), cluster.0, *tid));
                            }
                        });
                    }

                    // Numer of minimum clusters as inferred from DBSCAN
                    let k = clusters_mean.keys().len();

                    // Beginning roots (indices) of each cluster
                    // Since there are N - 1 steps in the dendrogram, to get k clusters we need the
                    // range of indices [N - 1 - 2k; N - 1 - k)
                    let n_1 = dendrogram.len();
                    // get the first k root labels
                    let mut cluster_root_labels = vec!();
                    let mut step_i = &dendrogram[n_1 - 1];
                    if k != 1 {
                        while cluster_root_labels.len() < k {
                            if cluster_root_labels.len() == 0 {
                                cluster_root_labels.push(step_i.cluster1);
                                cluster_root_labels.push(step_i.cluster2);
                            } else {
                                let mut cluster_to_check = cluster_root_labels
                                    .iter().max().unwrap().clone();

                                step_i = &dendrogram[cluster_to_check - n_1 - 1];
                                cluster_root_labels.push(step_i.cluster1);
                                cluster_root_labels.push(step_i.cluster2);

                                let cluster_to_check_i = cluster_root_labels.iter()
                                    .position(|x| x == &cluster_to_check).unwrap();
                                cluster_root_labels.remove(cluster_to_check_i);
                            }
                        };
                    } else {
                        cluster_root_labels.push(n_1+n_1);
                    }
//                let cluster_roots = (n_1 + 1 - 2 * (k)..n_1 + 1 - k);
                    let mut position_count: HashSet<usize> = HashSet::new();

                    for (index, cluster_root) in cluster_root_labels.into_iter().enumerate() {
                        if cluster_root > n_1 {
                            let cluster_root_id = cluster_root - n_1 - 1;
                            let hap_root = &dendrogram[cluster_root_id];
                            let mut new_haplotype = Haplotype::start(
                                hap_root.size, cluster_root_id, index);
                            let mut dendro_ids = dendro_ids.lock().unwrap();
                            new_haplotype.add_variants_per_genome(dendrogram, &dendro_ids, clusters);
                            debug!("{} {:?} {:?}",
                                   cluster_root_id,
                                   new_haplotype.node_size,
                                   new_haplotype.variants.len());

                            position_count.extend(&new_haplotype.variant_indices);
                            haplotypes_vec.push(new_haplotype);
                        } else {
                            let mut dendro_ids = dendro_ids.lock().unwrap();
                            let variant_pos = dendro_ids.get(&cluster_root).expect("Label not found");
                            let mut variant_map = HashMap::new();
                            for (pos, variant) in variant_pos.iter() {
                                let captured_tid = variant_map.entry(variant.2)
                                    .or_insert(HashMap::new());

                                let captured_var = captured_tid.entry(*pos)
                                    .or_insert(BTreeMap::new());

                                captured_var.entry(variant.0.clone())
                                    .or_insert((variant.1, cluster_root));

                                let clusters_tid = clusters.entry(variant.2)
                                    .or_insert(HashMap::new());

                                let cluster_pos = clusters_tid.entry(*pos)
                                    .or_insert(BTreeMap::new());

                                cluster_pos.insert(variant.0.clone(), (variant.1, cluster_root));
                            }

                            let mut new_haplotype = Haplotype {
                                root_cluster_id: cluster_root,
                                variant_indices: [cluster_root].into_iter().cloned().collect(),
                                variants: HashMap::new(),
                                variants_genome: variant_map,
                                node_size: 1,
                                haplotype_index: index,
                            };
                            position_count.extend(&new_haplotype.variant_indices);
                            haplotypes_vec.push(new_haplotype);
                        }
                    }
                    debug!("Variants found in tree {} {:?}", position_count.len(), position_count);

                    for (hap_index, haplotype) in haplotypes_vec.iter().enumerate() {
                        let file_name = format!("{}_strain_{}.fna", output_prefix.to_string(), hap_index);

                        let file_path = Path::new(&file_name);

                        // Open haplotype file or create one
                        let mut file_open = File::create(file_path)
                            .expect("No Read or Write Permission in current directory");

                        // Generate the consensus genome by checking each variant
                        // Variant has to be in more than 0.5 of population
                        for (tid, original_contig) in contigs.iter() {
                            let mut contig = String::new();

                            let mut skip_n = 0;
                            let mut skip_cnt = 0;
                            let mut char_cnt = 0;
                            let mut max_abund = 0.0;
                            let mut variations = 0;

                            for (pos, base) in original_contig.iter().enumerate() {
                                if skip_cnt < skip_n {
                                    skip_cnt += 1;
                                } else {
                                    let mut max_var = "";

                                    skip_n = 0;
                                    skip_cnt = 0;
                                    if haplotype.variants_genome.contains_key(&tid) {
                                        if haplotype.variants_genome[tid].contains_key(&(pos as i32)) {
                                            let hash = &haplotype.variants_genome[tid][&(pos as i32)];
                                            for (var, clusters) in hash.iter() {
                                                max_var = var;
                                                let abundance_map = &variants[tid][&(pos as i32)][var];
                                                let mut mean_var: f32 = 0.;
                                                let mut mean_d: f32 = 0.;
                                                abundance_map.iter().map(|(var, d)|{
                                                    mean_var += *var;
                                                    mean_d += *d;
                                                });
                                                mean_var = mean_var / sample_count;
                                                mean_d = mean_d / sample_count;
                                                let abundance= mean_var / mean_d;

                                                if abundance > max_abund {
                                                    max_abund = abundance;
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
                            writeln!(file_open, ">{}_strain_{}\t#max_variant_abundance_{}\t#variants_{}",
                                     target_names[tid],
                                     hap_index,
                                     max_abund,
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

fn get_condensed_distances(variant_info_all: &MutexGuard<Vec<(&i32, String, (f32, Vec<f32>), &i32)>>,
                           indels_map: &mut HashMap<i32, HashMap<i32, BTreeMap<String, BTreeSet<i64>>>>,
                           snps_map: &mut HashMap<i32, HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>) -> Arc<Mutex<Array1<f32>>> {
    let mut variant_distances: Arc<Mutex<Array1<f32>>>
        = Arc::new(
        Mutex::new(
            Array1::<f32>::zeros(((variant_info_all.len().pow(2) as usize
                - variant_info_all.len()) / 2 as usize))));
    debug!("Filling matrix of size {}",
           ((variant_info_all.len().pow(2) as usize - variant_info_all.len()) / 2 as usize));

    // produced condensed pairwise distances
    // described here: https://docs.rs/kodama/0.2.2/kodama/
    (0..variant_info_all.len()-1)
        .into_par_iter().for_each(|row_index| {
        let mut row_variant_set = &BTreeSet::new();
        let row_info = &variant_info_all[row_index];
        // lazily get the row variant read id set
        if indels_map[&row_info.3].contains_key(&row_info.0) {
            if indels_map[&row_info.3][&row_info.0].contains_key(&row_info.1) {
                row_variant_set = &indels_map[&row_info.3][&row_info.0][&row_info.1];
            } else if snps_map[&row_info.3].contains_key(&row_info.0) {
                let var_char = row_info.1.as_bytes()[0] as char;
                if snps_map[&row_info.3][&row_info.0].contains_key(&var_char) {
                    row_variant_set = &snps_map[&row_info.3][&row_info.0][&var_char];
                }
            }
        } else if snps_map[&row_info.3].contains_key(&row_info.0) {
            let var_char = row_info.1.as_bytes()[0] as char;
            if snps_map[&row_info.3][&row_info.0].contains_key(&var_char) {
                row_variant_set = &snps_map[&row_info.3][&row_info.0][&var_char];
            }
        }

        let row_start = *row_info.0 as usize;
        let row_end = row_start + row_info.1.len() - 1;

        (row_index+1..variant_info_all.len())
            .into_par_iter().for_each(|col_index| {
            if row_index == col_index {
//                let mut variant_distances = variant_distances.lock().unwrap();
//                variant_distances[[row_index, col_index]] = 0.;
            } else {
                let mut col_variant_set = &BTreeSet::new();
                let col_info = &variant_info_all[col_index];
                if indels_map[&col_info.3].contains_key(&col_info.0) {
                    if indels_map[&col_info.3][&col_info.0].contains_key(&col_info.1) {
                        col_variant_set = &indels_map[&col_info.3][&col_info.0][&col_info.1];
                    } else if snps_map[&col_info.3].contains_key(&col_info.0) {
                        let var_char = col_info.1.as_bytes()[0] as char;
                        if snps_map[&col_info.3][&col_info.0].contains_key(&var_char) {
                            col_variant_set = &snps_map[&col_info.3][&col_info.0][&var_char];
                        }
                    }
                } else if snps_map[&col_info.3].contains_key(&col_info.0) {
                    let var_char = col_info.1.as_bytes()[0] as char;
                    if snps_map[&col_info.3][&col_info.0].contains_key(&var_char) {
                        col_variant_set = &snps_map[&col_info.3][&col_info.0][&var_char];
                    }
                }

                let col_start = *col_info.0 as usize;
                let col_end = col_start + col_info.1.len() - 1;

                let mut distance: f32;

                // If the variants share positions, then instantly they can't be in the same
                // gentoype so max distance
                if row_start <= col_end && col_start <= row_end {
                    distance = 1.;
                } else {
//                    let intersection_len = (row_variant_set
//                        .intersection(&col_variant_set).collect::<HashSet<_>>().len()) as f32;
//
//                    let union_len = (row_variant_set
//                        .union(&col_variant_set).collect::<HashSet<_>>().len()) as f32;
//
//                    // Jaccard Similarity Modified for total read depth
//                    // Calculates the observed frequency of two variants together
//                    let jaccard = intersection_len /
//                        ((row_info.2).1 + (col_info.2).1 - intersection_len);
//
////                            let row_cov = contig_coverage_means[&row_info.3];
////                            let col_cov = contig_coverage_means[&col_info.3];
//
//                    let p: f32 = (row_info.2).0;
//                    let q: f32 = (col_info.2).0;
//                    // p * q = expected frequency
//                    let pq = p * q;
//
//                    // linkage disequilibrium
//                    let d = jaccard - pq;
//
//                    // correlation
//                    let rho = d.powf(2.) / (pq * (p - 1.) * (q - 1.));
//
//                    debug!("p {} q {} x {} D {} rho {} distance {}", p, q, jaccard, d, rho, 1. - rho);
                    // Distance will be defined as the mean between jaccard dist
                    // and dist_f
                    let mut corr = 0.;
                    let mut w = 0;
                    if (row_info.2).1.len() > 1 {
                        // Calculate pearson's correlation between two variants based on abundances
                        let x_bar = (row_info.2).0;
                        let y_bar = (col_info.2).0;
                        let mut s_x: f32 = 0.;
                        let mut s_y: f32 = 0.;
                        let mut top: f32 = 0.;

                        for (x, y) in (row_info.2).1.iter().zip((col_info.2).1.iter()) {
                            if (x > &0.) || (y > &0.) {
                                w += 1;
                            }
                            let x_diff = x - x_bar;
                            let y_diff = y - y_bar;
                            top += x_diff * y_diff;
                            s_x += x_diff.powf(2.);
                            s_y += y_diff.powf(2.)
                        }
                        corr = top / (s_x.powf(1. / 2.) * s_y.powf(1. / 2.));
                        distance = 1. - corr;

                    } else {
                        distance = ((row_info.2).0 - (col_info.2).0).abs();

                    }

                }

                match condensed_index(
                    row_index, col_index, variant_info_all.len()) {
                    Some(index) => {
                        let mut variant_distances = variant_distances.lock().unwrap();
                        variant_distances[[index]] = distance;
                    }
                    None => {
                        debug!("No corresponding index for row {} and col {}",
                               row_index, col_index);
                    }
                };
//                let mut variant_distances = variant_distances.lock().unwrap();
//                variant_distances[[row_index, col_index]] = distance;
//                variant_distances[[col_index, row_index]] = distance;
            }
        });
    });
    debug!("{:?}", variant_distances);
    return variant_distances
}

// helper function to get the index of condensed matrix from it square form
fn condensed_index(i: usize, j: usize, n: usize) -> Option<usize>{
    if i == j {
        return None
    } else {
        return Some(n*i - i*(i+1)/2 + j - 1 - i)
//        return Some(n*(n-1)/2 - (n - row_i)*(n - row_i - 1)/2 + col_j - row_i - 1)
    }
}