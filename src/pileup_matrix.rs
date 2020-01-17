use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use pileup_structs::*;
use matrix_handling::*;
use std::str;
use std::path::Path;
use std::io::prelude::*;
use rayon::prelude::*;
use ndarray::{Array2, Array1, ArrayView};
use cogset::{Euclid, Dbscan, BruteScan};
use kodama::{Method, nnchain, Dendrogram};
use std::sync::{Arc, Mutex, MutexGuard};
use haplotypes_and_genotypes::*;
use std::fs::File;
use std::process;
use std::cmp;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use tempfile;
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

                // A vector the length of the number of samples
                // Cumulatively calculates the product of abundaces from each sample
                let geom_mean_vec =
                    Arc::new(
                        Mutex::new(
                            vec![1.; sample_count as usize]));

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
                                let mut total_d: f32 = 0.;
                                let mut freqs = Vec::new();
                                if !var.contains("R") {
                                    // Get the mean abundance across samples
                                    let mut sample_idx: usize = 0;
                                    abundances_vector.iter().for_each(|(var, d)| {
                                        mean_var += *var;
                                        // Total depth of location
                                        total_d += *d;
                                        if var > &0. {
                                            let freq = (*var + 1.) / (*d + 1.);
                                            let mut geom_mean_vec = geom_mean_vec.lock().unwrap();
                                            geom_mean_vec[sample_idx] = geom_mean_vec[sample_idx] * freq;
                                            let sample_coverage = contig_coverages[sample_idx];

                                            freqs.push(freq * (sample_coverage / max_coverage));
                                            abundance += *var / *d;
                                        } else {
                                            freqs.push(0.);
                                        }
                                    });

                                    mean_var = mean_var / sample_count;
                                    abundance = abundance / sample_count;

                                    let mut variant_info_all = variant_info_all
                                        .lock().unwrap();


                                    variant_info_all.push(
                                        (position, var.to_string(),
                                         (total_d, freqs), tid));
                                }
                            }
                        });
                });
                let variant_info_all = variant_info_all.lock().unwrap();
                if variant_info_all.len() > 1 {

                    info!("Generating Variant Distances with {} Variants", variant_info_all.len());

                    let mut geom_mean_vec = geom_mean_vec.lock().unwrap();
                    let geom_means = geom_mean_vec.iter()
                        .map(|prod| {
                            prod / variant_info_all.len() as f32
                        }).collect::<Vec<f32>>();

                    let mut variant_distances =
                        get_condensed_distances(&variant_info_all[..],
                                                indels_map,
                                                snps_map,
                                                &geom_means,
                                                sample_count as i32);

                    let mut variant_distances =  variant_distances.lock().unwrap();

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

                    variant_distances.write_npy(&tmp_path);

//                println!("{:?}", variant_distances);

                    let cmd_string = format!(
                        "set -e -o pipefail; \
                        nmf.py {} {} {} {} {}",
                        // NMF
                        cmp::min(5, variant_info_all.len()),
                        cmp::min(15, variant_info_all.len()),
                        100,
                        tmp_path,
                        sample_count as i32);
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
                        println!("{}", sample_names[0]);
                        println!("{}", out);
                    }


                    tmp_dir.close().expect("Unable to close temp directory");


//                let py = gil.python();
//                let nmfpy = PyModule::import(py, "nmf.py").unwrap();
//                let variant_distances = variant_distances.to_owned().into_pyarray(py);
//                nmfpy.call1("perform_nmf", (variant_distances,)).unwrap();
                } else {
                    debug!("Not enough variants found {:?}, Non heterogeneous population", variant_info_all);
                }
            }
        }
    }

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
