use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use std::str;
use std::sync::{Arc, Mutex};
use std::io::prelude::*;
use rayon::prelude::*;
use permutation::*;
use codon_structs::*;
use haplotypes_and_genotypes::*;
use bio_types::strand;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
//use rusty_machine::learning::dbscan::DBSCAN;
//use rusty_machine::learning::UnSupModel;
//use rusty_machine::linalg::Matrix;
use cogset::{Euclid, Dbscan, BruteScan};
use kodama::{Method, nnchain, Dendrogram};
use nalgebra as na;
use itertools::{izip, Itertools};
use itertools::EitherOrBoth::{Both, Left, Right};
use std::path::Path;
use std::fs::File;

#[macro_use(array)]

pub enum PileupStats {
    PileupContigStats {
        nucfrequency: HashMap<i32, BTreeMap<char, BTreeSet<i32>>>,
        variants_in_reads: HashMap<i32, BTreeMap<i32, String>>,
        variant_abundances: HashMap<i32, BTreeMap<String, f64>>,
        variant_count: Vec<f64>,
        depth: Vec<f64>,
        indels: HashMap<i32, BTreeMap<String, BTreeSet<i32>>>,
        genotypes_per_position: HashMap<usize, usize>,
        mean_genotypes: f32,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: f32,
        variations_per_base: usize,
        coverage: f32,
        variance: f32,
        observed_contig_length: u32,
        num_covered_bases: i32,
        num_mapped_reads: u64,
        total_mismatches: u32,
        contig_end_exclusion: u32,
        min: f32,
        max: f32,
        method: String,
        read_error_rate: f64,
        // Clusters hashmap:
        // Key = Position
        // Value = K: Variant, V: (DBSCAN Cluster, HAC Index/initial cluster)
        clusters: HashMap<i32, BTreeMap<String, (i32, usize)>>,
        clusters_mean: HashMap<i32, f64>,
        dendrogram: Dendrogram<f64>,
        haplotypes: Vec<Haplotype>,
    }
}

impl PileupStats {
    pub fn new_contig_stats(min: f32, max: f32,
                            contig_end_exclusion: u32) -> PileupStats {
        PileupStats::PileupContigStats {
            nucfrequency: HashMap::new(),
            variants_in_reads: HashMap::new(),
            variant_abundances: HashMap::new(),
            variant_count: Vec::new(),
            depth: vec!(),
            indels: HashMap::new(),
            genotypes_per_position: HashMap::new(),
            mean_genotypes: 0.0,
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0.0,
            variations_per_base: 0,
            coverage: 0.00,
            variance: 0.00,
            observed_contig_length: 0,
            num_covered_bases: 0,
            contig_end_exclusion: contig_end_exclusion,
            num_mapped_reads: 0,
            total_mismatches: 0,
            min: min,
            max: max,
            method: "".to_string(),
            read_error_rate: 0.0,
            clusters: HashMap::new(),
            clusters_mean: HashMap::new(),
            dendrogram: Dendrogram::new(0),
            haplotypes: Vec::new(),
        }
    }
}

pub trait PileupFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i32>>>,
                  indels_positions: HashMap<i32, BTreeMap<String, BTreeSet<i32>>>,
                  tid: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  method: &str,
                  coverages: Vec<f32>,
                  ups_and_downs: Vec<i32>);

    fn calc_error(&mut self);

    fn calc_variants(&mut self,
                     min_variant_depth: usize,
                     coverage_fold: f32);

    fn filter_variants(&mut self, max_iter: usize);

    fn generate_variant_contig(&mut self,
                               original_contig: &Vec<u8>,
                               output_prefix: &str);

    fn generate_minimum_genotypes(&mut self) -> HashMap<usize, usize>;

    fn generate_svd(&mut self);

    fn calc_gene_mutations(&mut self,
                           gff_map: &HashMap<String, Vec<bio::io::gff::Record>>,
                           ref_sequence: &Vec<u8>,
                           codon_table: &CodonTable);

    fn cluster_variants(&mut self);

    fn print_variants(&mut self, ref_sequence: &Vec<u8>, sample_idx: i32);
}

impl PileupFunctions for PileupStats {
    fn setup(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variants_in_reads,
                ref mut variant_abundances,
                ref mut variant_count,
                ref mut depth,
                ref mut indels,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut variations_per_base,
                ref mut coverage,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ref mut clusters,
                ..
            } => {
                *nucfrequency = HashMap::new();
                *variants_in_reads = HashMap::new();
                *variant_abundances = HashMap::new();
                *variant_count = Vec::new();
                *depth = vec!();
                *indels = HashMap::new();
                *tid = 0;
                *total_indels = 0;
                *target_name = vec!();
                *target_len = 0.0;
                *variations_per_base = 0;
                *coverage = 0.00;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
                *clusters = HashMap::new();
            }
        }
    }

    fn add_contig(&mut self, nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i32>>>,
                  indel_positions: HashMap<i32, BTreeMap<String, BTreeSet<i32>>>,
                  target_id: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  method: &str,
                  coverages: Vec<f32>,
                  ups_and_downs: Vec<i32>) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variant_count,
                ref mut depth,
                ref mut indels,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut variance,
                ref mut coverage,
                ref mut method,
                ..
            } => {
                *nucfrequency = nuc_freq;
                *indels = indel_positions;
                *tid = target_id;
                *total_indels = total_indels_in_contig;
                *target_name = contig_name;
                *target_len = coverages[0];
                *coverage = coverages[1];
                *variance = coverages[2];
                *method = method.to_string();
                let variant_count_safe = Arc::new(Mutex::new(vec![0.; ups_and_downs.len()]));
                *depth = vec![0.; ups_and_downs.len()];
                let mut cumulative_sum = 0;
                for (pos, current) in ups_and_downs.iter().enumerate() {
                    cumulative_sum += *current;
                    depth[pos] = cumulative_sum as f64;
                }
//                let adjusted_depth = Arc::new(Mutex::new(depth.clone()));

                // Calculate how many reads have variant at each position
                // to go into linear regression predicting read error rate
                depth.par_iter().enumerate().for_each(
                    |(pos, _d)|{
                        let mut variant_count_safe = variant_count_safe.lock().unwrap();
//                        let mut adjusted_depth = adjusted_depth.lock().unwrap();
                        let snp_map = match nucfrequency.get(&(pos as i32)) {
                            Some(map) => map.to_owned(),
                            None => BTreeMap::new(),
                        };
                        let indel_map = match indels.get(&(pos as i32)){
                            Some(map) => map.to_owned(),
                            None => BTreeMap::new(),
                        };
                        if snp_map.len() > 1 {
                            for (key, reads) in snp_map.iter() {
                                if key != &"R".chars().collect::<Vec<char>>()[0] {
                                    variant_count_safe[pos] += reads.len() as f64;
                                }
                            }
                        }
                        if indel_map.len() > 0 {
                            for reads in indel_map.values() {
                                variant_count_safe[pos] += reads.len() as f64;
                            }
                        }
                });
                *variant_count = variant_count_safe.lock().unwrap().to_vec();

                debug!("new contig added {} with coverage {} and variance {}", tid, coverage, variance);
            }
        }
    }

    fn calc_error(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                variant_count,
                depth,
                read_error_rate,
                ..
            } => {
                let data = vec![("Y", variant_count.clone()), ("X", depth.clone())];
                let data = RegressionDataBuilder::new()
                    .build_from(data).expect("Unable to build regression from data");
                let formula = "Y ~ X";
                let model = FormulaRegressionBuilder::new()
                    .data(&data)
                    .formula(formula)
                    .fit()
                    .expect("Unable to fit data to formula");
                let parameters = model.parameters;
                let standard_errors = model.se;
                let pvalues = model.pvalues;
                debug!("Linear regression results: \n params {:?} \n se {:?} \n p-values {:?}",
                         parameters,
                         standard_errors,
                         pvalues.pairs());
                *read_error_rate = parameters.pairs()[0].1.clone() + 2.*standard_errors.pairs()[0].1.clone();

            }
        }
    }

    fn calc_variants(&mut self, min_variant_depth: usize, coverage_fold: f32){
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variants_in_reads,
                ref mut variant_abundances,
                depth,
                ref mut indels,
                target_len,
                ref mut variations_per_base,
                ref mut coverage,
                tid,
                read_error_rate,
                ..
            } => {
                let variants = Arc::new(Mutex::new(HashMap::new())); // The relative abundance of each variant
                let read_variants = Arc::new(Mutex::new(HashMap::new())); // The reads with variants and their positions
                let variant_count = Arc::new(Mutex::new(0));
                let indels = Arc::new(Mutex::new(indels));
                let nucfrequency = Arc::new(Mutex::new(nucfrequency));

                // for each location calculate if there is a variant based on read depth
                // Uses rayon multithreading
                depth.par_iter_mut().enumerate().for_each(|(i, d)| {
//                    let read_variants = Arc::clone(&read_variants);
//                    let variant_count = Arc::clone(&variant_count);
                    let mut rel_abundance = BTreeMap::new();
                    if (*coverage * (1.0 - coverage_fold) <= *d as f32
                        && *d as f32 <= *coverage * (1.0 + coverage_fold))
                        || (coverage_fold == 0.0) {
//                        if d >= &mut min_variant_depth.clone() {
                        // INDELS act differently to normal variants
                        // The reads containing this variant don't contribute to coverage values
                        // So we need to readjust the depth at these locations to show true
                        // read depth. i.e. variant depth + reference depth
                        let mut indels
                            = indels.lock().unwrap();
                        let indel_map = match indels.get(&(i as i32)) {
                            Some(map) => map.to_owned(),
                            None => BTreeMap::new(),
                        };
                        if indel_map.len() > 0 {
                            for (indel, read_ids) in indel_map.iter() {
                                let count = read_ids.len();
                                if indel.contains("N") {
                                    *d += count as f64;
                                }
                                if (count >= min_variant_depth) & (count as f64 / *d > *read_error_rate){
                                    rel_abundance.insert(indel.to_owned(), count as f64 / *d);
                                    for read in read_ids {
                                        let mut read_variants
                                            = read_variants.lock().unwrap();

                                        let read_vec = read_variants
                                            .entry(read.clone())
                                            .or_insert(BTreeMap::new());
                                        read_vec.insert(i as i32, indel.clone());
                                    }
                                    let mut variant_count = variant_count.lock().unwrap();
                                    *variant_count += 1;
                                } else {

                                    let indel_map_back = indels
                                        .entry(i as i32).or_insert(BTreeMap::new());
                                    indel_map_back.remove(indel);
                                }
                            }
                        }
                        let mut nucfrequency
                            = nucfrequency.lock().unwrap();

                        let nuc_map = match nucfrequency.get(&(i as i32)) {
                            Some(map) => map.to_owned(),
                            None => BTreeMap::new(),
                        };
                        if nuc_map.len() > 0 {
                            for (base, read_ids) in nuc_map.iter() {
                                let count = read_ids.len();

                                if (count >= min_variant_depth) & ((count as f64 / *d) > *read_error_rate) {
                                    rel_abundance.insert(base.to_string(), count as f64 / *d);

                                    for read in read_ids {
                                        let mut read_variants
                                            = read_variants.lock().unwrap();
                                        let read_vec = read_variants
                                            .entry(read.clone())
                                            .or_insert(BTreeMap::new());
                                        read_vec.insert(i as i32, base.to_string());
                                    }
                                    if base != &"R".chars().collect::<Vec<char>>()[0] {
                                        let mut variant_count = variant_count.lock().unwrap();
                                        *variant_count += 1;
                                    }
                                } else {

                                    let nuc_map_back = nucfrequency
                                        .entry(i as i32).or_insert(BTreeMap::new());
                                    nuc_map_back.remove(base);
                                }
                            }
                        };

//                        }
                    }

                    if rel_abundance.len() > 1 {
                        let mut variants = variants.lock().unwrap();
//                        debug!("Relative Abundances: {:?}", rel_abundance);
                        variants.insert(i as i32, rel_abundance);
                    } else if rel_abundance.len() == 1 {
                        if !(rel_abundance.contains_key(&"R".to_string())) {
                            let mut variants = variants.lock().unwrap();
//                            debug!("Relative Abundances: {:?}", rel_abundance);
                            variants.insert(i as i32, rel_abundance);
                        } else {
                            // Need to remove R from whitelist
                            let mut nucfrequency
                                = nucfrequency.lock().unwrap();
                            let nuc_map = nucfrequency
                                .entry(i as i32).or_insert(BTreeMap::new());
                            if nuc_map.len() > 0 {
                                for (base, read_ids) in nuc_map.iter() {
                                    let count = read_ids.len();
                                    if base == &("R".as_bytes()[0] as char) {
                                        for read in read_ids {
                                            let mut read_variants
                                                = read_variants.lock().unwrap();
                                            let read_vec = read_variants
                                                .entry(read.clone())
                                                .or_insert(BTreeMap::new());
                                            read_vec.remove(&(i as i32));
                                        }
                                    }
                                }
                            } else {
                                nucfrequency.remove(&(i as i32));
                            }
                        }
                    }

                });

                let read_variants = read_variants.lock().unwrap();
                *variants_in_reads = read_variants.to_owned();
                let variants = variants.lock().unwrap();
                *variant_abundances = variants.to_owned();
                let variant_count = variant_count.lock().unwrap();
                debug!("Total variants for {}: {:?}", tid, variant_count);
                *variations_per_base = *variant_count;
                let mut nucfrequency = nucfrequency.lock().unwrap();
                **nucfrequency = nucfrequency.to_owned();
                let mut indels = indels.lock().unwrap();
                **indels = indels.to_owned();
            }
        }
    }

    fn filter_variants(&mut self, max_iter: usize) {
//        match self{
//            PileupStats::PileupContigStats {
//                ref mut variant_abundances,
//                ..
//            } => {
//                let mut iteration = 0;
//                // Get most abundant and second most abundant variant values
//                let mut max_a = HashMap::new();
//                let mut max_b = HashMap::new();
//                variant_abundances.par_iter().enumerate().for_each(|(pos, var_map)|{
//                    if var_map.len() > 0 {
////                        let variant_values = var_map.values().collect::<Vec<f32>>();
////                        let variant_sum = variant_values.iter().sum();
////                        let ref_abundance = 1.0 - variant_sum;
////                        if ref_abundance >= 0.5 {
//
//                        }
//                    }
//                });
//            }
//        }
    }

    fn generate_variant_contig(&mut self,
                               original_contig: &Vec<u8>,
                               output_prefix: &str){
        match self {
            PileupStats::PileupContigStats {
                ref mut variant_abundances,
                ref mut clusters,
                clusters_mean,
                dendrogram,
                target_name,
                ref mut haplotypes,
                ..
            } => {
                // First we need to convert DBSCAN cluster ids into taxonomic rank ids
                // Clusters of higher abundance mutations will be closer to root
                // Ordered from root to leaf
                // We can then check how dissimilar two variants are on the hierarchical cluster
//                let mut ordered_clusters: Vec<_> = clusters_mean.iter().collect();
//                ordered_clusters
//                    .sort_by(|a, b|
//                        b.1.partial_cmp(a.1).unwrap());
//                debug!("{:?}", ordered_clusters);

                // Clusters is set up as HashMap<Position, BTreeMap<Variant, (Cluster_ID, Dendro_index)>>
                // In this case we want to rearrange to HashMap<dendro_ID, HashMap<Position, (Variant, db_clust)>>
                // This will allow us to disentangle positions where more than one variant is possible
                let mut dendro_ids = Arc::new(Mutex::new(HashMap::new()));
                clusters.par_iter().for_each(|(position, variant_map)| {
                    for (variant, cluster) in variant_map.iter() {
                        let mut dendro_ids = dendro_ids.lock().unwrap();
                        let clust = dendro_ids.entry(cluster.1)
                            .or_insert(BTreeMap::new());

                        clust.entry(*position)
                             .or_insert((variant.to_string(), cluster.0));
                    }
                });

                // Numer of minimum clusters as inferred from DBSCAN
                let k = clusters_mean.keys().len();

                // Beginning roots (indices) of each cluster
                // Since there are N - 1 steps in the dendrogram, to get k clusters we need the
                // range of indices [N - 1 - 2k; N - 1 - k)
                let n_1 = dendrogram.len();
                let cluster_roots = (n_1 + 1 - 2 * (k)..n_1 + 1 - k);
                let mut haplotypes_vec = vec![Haplotype::new(); k];
                let mut position_count: HashSet<usize> = HashSet::new();

                for (index, cluster_root_id) in cluster_roots.into_iter().enumerate() {
                    let hap_root = &dendrogram[cluster_root_id];
                    let mut new_haplotype = Haplotype::start(
                        hap_root.size, cluster_root_id, index);
                    let mut dendro_ids = dendro_ids.lock().unwrap();
                    new_haplotype.add_variants(dendrogram, &dendro_ids);
                    debug!("{:?}", new_haplotype);

                    for (pos, variants) in new_haplotype.variants.iter(){
                        let cluster_pos = clusters.entry(*pos)
                            .or_insert(variants.clone());
                        for (variant, clust) in variants.iter(){
                            cluster_pos.insert(variant.to_string(), *clust);
                        }

                    }

                    position_count.extend(&new_haplotype.variant_indices);
                    haplotypes_vec[index] = new_haplotype;
                }
                debug!("Variants found in tree {}", position_count.len());

//                print!("[");
//                for step in dendrogram.steps(){
//                    println!("[{}, {}, {}, {}],", step.cluster1, step.cluster2, step.dissimilarity, step.size);
//                }
//                println!("]");
//                print!("{{");
//                for (pos, cluster) in clusters.iter() {
//                    println!("{:?}: {:?},", pos, cluster);
//                }
//                println!("}}");
//                let file_name = output_prefix.to_string() + &"_".to_owned()
//                    + &"genotypes".to_owned()
//                    + &".fna".to_owned();
//
//                let file_path = Path::new(&file_name);
//
//                // Open haplotype file or create one
//                let mut file_open = match File::open(file_path) {
//                    Ok(fasta) => fasta,
//                    Err(_e) => {
//                        match File::create(file_path) {
//                            Ok(fasta) => fasta,
//                            Err(e) => {
//                                println!("Cannot create file {:?}", e);
//                                std::process::exit(1)
//                            },
//                        }
//                    },
//                };
//
//                for (hap_index, haplotype) in haplotypes_vec.iter().enumerate() {
//                    writeln!(file_open, ">{}_haplotype_{}",
//                             str::from_utf8(target_name).expect("UTF-8 error"),
//                             hap_index);
//                    let mut contig = original_contig.clone();
//
//                    // Generate the new potential genotype
//                    for (pos, variant_set) in haplotype.variants.iter() {
//
//                        if variant_abundances[&(pos as i32)].len() > 0 {
//                            let hash = &variant_abundances[&(pos as i32)];
//                            for (var, abundance) in hash.iter() {
//                                if abundance > &max_abund {
//                                    max_var = var;
//                                    max_abund = *abundance;
//                                }
//                            }
//                            if max_abund >= 0.5 {
//                                if max_var.contains("N") {
//                                    skip_n = max_var.len() - 1;
//                                    skip_cnt = 0;
//                                } else {
//                                    contig = contig + max_var;
//                                }
//                            } else {
//                                contig = contig + str::from_utf8(&[*base]).unwrap();
//                            }
//                        } else {
//                            contig = contig + str::from_utf8(&[*base]).unwrap();
//                        }
//                    };
//                    contig = contig + "\n \n";
//                    match file_open.write_all(contig.as_bytes()) {
//                        Ok(consensus_genome) => consensus_genome,
//                        Err(e) => {
//                            println!("Cannot write to file {:?}", e);
//                            std::process::exit(1)
//                        }
//                    };
                *haplotypes = haplotypes_vec;
            }
        }
    }

    fn generate_minimum_genotypes(&mut self) -> HashMap<usize, usize>{
        match self {
            PileupStats::PileupContigStats {
                ref mut variant_abundances,
                ref mut nucfrequency,
                ref mut indels,
                ref mut variants_in_reads,
                ref mut genotypes_per_position,
                ref mut mean_genotypes,
                tid,
                target_len,
                variations_per_base,
                coverage,
                ..
            } => {
                let genotypes =
                    Arc::new(Mutex::new(HashMap::new()));
                let variant_count =
                    Arc::new(Mutex::new(0));
                let total_genotype_count =
                    Arc::new(Mutex::new(0));

                debug!("starting genotyping of tid {}, of length {}, and var per b {} at {} times coverage",
                        tid, target_len, variations_per_base, coverage);

                variant_abundances.par_iter().for_each(|(position, variants)| {
                    // For each variant we calculate the minimum number of genotypes possible
                    // based on variants in reads mapping to this variant location
                    if variants.len() > 0 {

                        let genotypes = Arc::clone(&genotypes);
                        let variant_count = Arc::clone(&variant_count);
                        let total_genotype_count = Arc::clone(&total_genotype_count);

                        let mut genotype_pos = HashMap::new();

                        for (var, _abundance) in variants.iter() {
                            let genotype_count = genotype_pos.entry(var.to_string())
                                .or_insert(0);

                            let genotype_vec = Arc::new(
                                Mutex::new(
                                    Vec::new()));

                            let mut read_ids = BTreeSet::new();

                            let indel_map = match indels
                                .get(position){
                                Some(map) => map.to_owned(),
                                None => BTreeMap::new()
                            };

                            let nuc_map = match nucfrequency
                                .get(position){
                                Some(map) => map.to_owned(),
                                None => BTreeMap::new()
                            };

                            if indel_map.contains_key(var) {
                                read_ids =
                                    match indel_map.get(var) {
                                        Some(ids) => ids.clone(),
                                        None => {
                                            println!("Variant not in indel hash");
                                            std::process::exit(1)
                                        },
                                    };
                            } else if nuc_map.contains_key(
                                &(var.chars().collect::<Vec<char>>()[0])) {
                                read_ids =
                                    match nuc_map
                                        .get(&(var.chars().collect::<Vec<char>>()[0])) {
                                        Some(ids) => ids.clone(),
                                        None => {
                                            println!("Variant not in frequency Hash");
                                            std::process::exit(1)
                                        },
                                    };
                            }
                            let mut left_most_variants: Vec<i32> = Vec::new();
                            let read_vec = read_ids.into_iter().collect::<Vec<i32>>();


                            for read_id in read_vec.iter() {
                                let position_map = match variants_in_reads.get(read_id) {
                                    Some(positions) => {
                                        positions},
                                    None => {
                                        debug!("read id not recorded in variant map {}, {}", var, read_id);
                                        break
                                    },
                                };
                                left_most_variants.push(
                                    *position_map.keys().cloned().collect::<Vec<i32>>().iter().min().unwrap());
                            }
                            // Generate the permutation of read id indices that create this list ordering
                            let permuted = permutation::sort(&left_most_variants[..]);
                            // Order the read vec by this permutation
                            let read_vec_sorted = permuted.apply_slice(&read_vec[..]);

                            // Loop through reads based on their left most variant position
                            read_vec_sorted.par_iter().for_each(|read_id|{
                                let mut genotype_record= Genotype::start(*position as usize);;
                                let position_map = variants_in_reads
                                    .get(read_id).expect("Read ID not found");
                                let mut genotype_vec = genotype_vec.lock().unwrap();

                                if genotype_vec.len() == 0 {
                                    // No genotype observed yet, so create one
                                    genotype_record.new(*read_id, position_map, &variant_abundances);

                                    genotype_vec.push(genotype_record.clone());
                                } else {
//                                    let position_map_variants: Vec<String> = position_map.values().cloned().collect();
//                                    let position_map_positions: Vec<i32> = position_map.keys().cloned().collect();
                                    let position_set = position_map.keys().cloned().collect::<HashSet<i32>>();

                                    let mut new_genotype;

                                    'genotypes: for genotype in genotype_vec.iter_mut() {

                                        // Create HashSets of variant positions to use intersection
                                        let genotype_position_set =
                                            genotype.ordered_variants.keys().cloned().collect::<HashSet<i32>>();

                                        let diff: Vec<i32> = position_set
                                            .symmetric_difference(&genotype_position_set).cloned().collect();

                                        // If there is a difference between positions of current read
                                        // and previous genotype
                                        if diff.len() > 0  {
                                            // Unique positions in both current genotype and
                                            // current read
                                            // check for internal differences except edges
                                            if diff
                                                .iter()
                                                .any(|read_pos|
                                                    (Some(read_pos) > genotype_position_set.iter().min())
                                                    && (Some(read_pos) < genotype_position_set.iter().max())){
                                                // difference against current

                                                genotype_record.new(*read_id, position_map, &variant_abundances);
                                                continue

                                            } else if diff
                                                .iter()
                                                .any(|read_pos|
                                                    (Some(read_pos) > position_set.iter().min())
                                                        && (Some(read_pos) < position_set.iter().max())){
                                                // difference against read
                                                // possible new genotype detected
                                                genotype_record.new(*read_id, position_map, &variant_abundances);
                                                continue

                                                // Check if variant location sits outside current genotype bounds
                                            // But contains no differences within
                                            } else if (genotype.ordered_variants.keys().max() <= diff.iter().max())
                                                || (diff.iter().min() <= genotype.ordered_variants.keys().min()) {
                                                // check variants against stored variants for a genotype
                                                let shared_positions: Vec<i32> = position_set
                                                    .intersection(&genotype_position_set).cloned().collect();
                                                new_genotype = genotype.check_variants(position_map, shared_positions);
                                                if new_genotype {
                                                    // possible new genotype detected
                                                    genotype_record.new(*read_id, position_map, &variant_abundances);
                                                    continue

                                                } else {
                                                    // Extension of previous genotype
                                                    genotype.read_ids.insert(*read_id);

                                                    for new_position in diff.iter() {
                                                        if !(genotype.ordered_variants.contains_key(&new_position)) {
                                                            let variant_map = &variant_abundances.get(new_position);
                                                            match variant_map {
                                                                Some(_variant_map) => {
                                                                    genotype.ordered_variants
                                                                        .insert(new_position.clone(), position_map[new_position].clone());
                                                                },
                                                                None => continue,
                                                            }
                                                        }
                                                    };

                                                    // reset genotype_record and move to next read
                                                    genotype_record = Genotype::start(*position as usize);
                                                    break 'genotypes;
                                                }
                                            }
                                        } else {
                                            // check variants against stored variants for a genotype
                                            let shared_positions: Vec<i32> = position_set
                                                .intersection(&genotype_position_set).cloned().collect();
                                            new_genotype = genotype.check_variants(position_map, shared_positions);

                                            if new_genotype {
                                                // possible new genotype detected
                                                genotype_record.new(*read_id, position_map, &variant_abundances);
                                                continue


                                            } else {
                                                // No difference with a previous genotype, reset current
                                                genotype_record = Genotype::start(*position as usize);
                                                genotype.read_ids.insert(*read_id);
                                                break
                                            }
                                        }
                                    }
                                    if genotype_record.ordered_variants.len() > 0 {
                                        // New genotype detected
                                        genotype_vec.push(genotype_record);
//                                        debug!("genotypes {:?}", genotype_vec);

//                                        genotype_record = Genotype::start(*position as usize);
                                    }
                                }
                            });

                            let genotype_vec = genotype_vec.lock().unwrap();

                            *genotype_count += genotype_vec.len();

                            let mut total_genotype_count = total_genotype_count.lock().unwrap();
                            *total_genotype_count += genotype_vec.len();
                            let mut variant_count = variant_count.lock().unwrap();
                            *variant_count += 1;
                        }
                        let mut genotypes = genotypes.lock().unwrap();
                        let mut genotypes_at_position: usize = 0;
                        for value in genotype_pos.values(){
                            genotypes_at_position += *value;
                        };
                        genotypes.insert(*position as usize, genotypes_at_position);
                    }
                });
                //Calc the mean number of genotypes per variant
                let variant_count = variant_count.lock().unwrap();
                let total_genotype_count = total_genotype_count.lock().unwrap();
                if *variant_count > 0 {
                    *mean_genotypes = *total_genotype_count as f32 / *variant_count as f32;
                } else {
                    *mean_genotypes = 0.0 as f32;
                }
                *genotypes_per_position = genotypes.lock().unwrap().clone();
                let to_return = genotypes.lock().unwrap().clone();
                return to_return
            }
        }
    }

    fn generate_svd(&mut self) {

        match self {
            PileupStats::PileupContigStats {
                indels,
                nucfrequency,
                variant_abundances,
                variants_in_reads,
                variations_per_base,
                ..
            } => {
                // Here we are gonna calculate how 'distant' one variant is from another
                // Based on which read ids they share
                // So we need a v by v matrix where each cell is the count of how many
                // read ids they share

                // total reads
                let n = variants_in_reads.keys().len();

                let mut variant_distances: na::base::DMatrix<f64>
                    = na::base::DMatrix::zeros(*variations_per_base, *variations_per_base);

                // Set up the distance matrix of size n*n
                let mut variant_indices = HashMap::new();
                let mut variant_index = 0usize;
                let mut read_indices = HashMap::new();
                let mut read_index = 0usize;

                for (read_id, read_map) in variants_in_reads.iter() {
                    let row_index = read_indices.entry(read_id).or_insert(read_index);
                    for (genomic_pos, variant) in read_map.iter() {
                        if !(variant.contains("R")){
                            let mut column_index;
                            if variant_indices.contains_key(&(genomic_pos, variant)) {
                                column_index = variant_indices
                                    .entry((genomic_pos, variant))
                                    .or_insert(variant_index);
                            } else {
                                column_index = variant_indices
                                    .entry((genomic_pos, variant))
                                    .or_insert(variant_index);
                                variant_index += 1;
                            }

                            variant_distances[(*row_index, *column_index)] = 1.;
                        }
                    }
                    read_index += 1;
                }
                // Use SVD from ndarray_linalg
                let svd_array = variant_distances.svd(false, false);
                // Get the v_t singular vector matrix from SVD
                let v_t = svd_array.singular_values.into_iter().cloned().collect::<Vec<_>>();

            }
        }

    }

    fn calc_gene_mutations(&mut self,
                           gff_map: &HashMap<String, Vec<bio::io::gff::Record>>,
                           ref_sequence: &Vec<u8>,
                           codon_table: &CodonTable) {
        match self {
            PileupStats::PileupContigStats {
                indels,
                variant_abundances,
                tid: _,
                target_name,
                depth,
                ..
            } => {
                let contig_name = String::from_utf8(target_name.clone())
                    .expect("Cannot create string from target_name");
                let placeholder = Vec::new();
                let gff_records = match gff_map.get(&contig_name){
                    Some(records) => records,
                    None => &placeholder,
                };
                debug!("Calculating population dN/dS from reads for {} genes", gff_records.len());
                let mut print_stream = &mut Mutex::new(std::io::stdout());
                gff_records.par_iter().enumerate().for_each(|(_id, gene)| {
                    let dnds = codon_table.find_mutations(gene, variant_abundances, indels, ref_sequence, depth);
                    let strand = gene.strand().expect("No strandedness found");
                    let frame: usize = gene.frame().parse().unwrap();
                    let start = gene.start().clone() as usize - 1;
                    let end = gene.end().clone() as usize;

                    let strand_symbol = match strand {
                        strand::Strand::Forward | strand::Strand::Unknown => {
                            '+'.to_string()
                        },
                        strand::Strand::Reverse => {
                            '-'.to_string()
                        }

                    };
                    let gene_id = &gene.attributes()["ID"].split("_").collect::<Vec<&str>>()[1];
                    let mut contig = gene.seqname().to_owned();
                    contig = contig + "_";
                    let mut indel_map = BTreeMap::new();

                    for cursor in start..end+1 {
                        let variant_map = match variant_abundances.get(&(cursor as i32)){
                            Some(map) => map,
                            None => continue,
                        };
                        if indels.contains_key(&(cursor as i32)){
                            indel_map = indels.get(&(cursor as i32))
                                .expect("No INDEL at this location").to_owned();
                        };
                        if variant_map.len() > 0 {
                            let mut print_stream = print_stream.lock().unwrap();


                            for (variant, abundance) in variant_map {
                                if variant.to_owned().contains("R"){
                                    continue
                                }
                                write!(print_stream, "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t",
                                         contig.clone()+gene_id, gene.start(),
                                         gene.end(), frame, strand_symbol, dnds, cursor).expect("Unable to write to stream");


                                if variant.to_owned().contains("N") {
                                    writeln!(print_stream, "{}\t{}\t{:.3}\t{}\t{}",
                                           variant,
                                           str::from_utf8(
                                               &ref_sequence[cursor..cursor
                                                   + variant.len() as usize]).unwrap(),
                                           abundance, depth[cursor], "D").expect("Unable to write to stream");

                                } else if indel_map.contains_key(variant) {
                                     writeln!(print_stream,"{}\t{}\t{:.3}\t{}\t{}",
                                           variant,
                                           str::from_utf8(
                                               &[ref_sequence[cursor]]).unwrap(),
                                           abundance, depth[cursor], "I").expect("Unable to write to stream");
                                } else {
                                    writeln!(print_stream, "{}\t{}\t{:.3}\t{}\t{}",
                                             variant,
                                             ref_sequence[cursor] as char,
                                             abundance,
                                             depth[cursor], "S").expect("Unable to write to stream");
                                }
                            }
                        }
                    }
                });

            }
        }

    }

    fn cluster_variants(&mut self) {
        match self{
            PileupStats::PileupContigStats {
                variant_abundances,
                ref mut clusters,
                ref mut clusters_mean,
                indels,
                nucfrequency,
                tid,
                ref mut dendrogram,
                ..} => {

                let mut abundance_euclid =
                    Arc::new(
                    Mutex::new(
                        Vec::new()));

                let mut abundance_float =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                let variant_info =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                let variant_info_all =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                let mut cluster_map =
                    Arc::new(
                    Mutex::new(
                        HashMap::new()));

                let mut cluster_hierarchies =
                    Arc::new(
                        Mutex::new(
                            HashMap::new()));


                let eps = 0.05;
                let min_cluster_size = 2;

                variant_abundances.iter().for_each(
                    |(position, hash)|{
                        // loop through each position that has variants
                        let mut abundance_euclid = abundance_euclid
                            .lock().unwrap();
                        let mut abundance_float = abundance_float
                            .lock().unwrap();
                        let mut variant_info = variant_info
                            .lock().unwrap();
                        let mut variant_info_all = variant_info_all
                            .lock().unwrap();

                        for (var, abundance) in hash.iter() {
                            if !var.contains("R") {
                                // This is a big hack but saves so much time
                                // Basically, any variant that has an abundance less than 0.05
                                // Automatically gets thrown into bin 0
                                // Downsides:
                                // Variants on border of eps look like they should cluster
                                // with bin 0, but instead create a new bin
                                // Upsides:
                                // Massive speed increase, massive decrease in memory usage
                                // Low abundant variants are all kind of in the same cluster any way
                                if abundance >= &eps {
                                    variant_info.push((position, var.to_string(), abundance));
                                    abundance_float.push(*abundance);
                                    abundance_euclid.push(Euclid([*abundance]));
                                } else {
                                    let mut cluster_map =
                                        cluster_map.lock().unwrap();
                                    let mut cluster_hierarchies =
                                        cluster_hierarchies.lock().unwrap();

                                    let position = cluster_map
                                        .entry(*position).or_insert(BTreeMap::new());
                                    position.insert(var.to_string(), 0);

                                    let cluster_mean = cluster_hierarchies
                                        .entry(0).or_insert(Vec::new());
                                    cluster_mean.push(abundance.to_owned());
                                }
                                variant_info_all.push((position, var.to_string(), abundance));
                            }
                        }
                });

                let abundance_euclid = abundance_euclid.lock().unwrap();
                let variant_info = variant_info.lock().unwrap();
                let variant_info_all = variant_info_all.lock().unwrap();
                let abundance_float = abundance_float.lock().unwrap();
                let scanner = BruteScan::new(&abundance_euclid);
                debug!("Beginning clustering of {} variants out of {}", abundance_euclid.len(),
                       variant_info_all.len());
                let mut dbscan = Dbscan::new(scanner,
                                             eps,
                                             min_cluster_size);


                let db_clusters = dbscan.by_ref().collect::<Vec<_>>();
                debug!("{} True clusters found for abundance > eps", db_clusters.len());

                let noise_points = dbscan.noise_points()
                    .iter().cloned().collect::<Vec<_>>();

                // get the number of clusters so we can rescue orphan variants
                let mut number_of_clusters =
                    Arc::new(
                        Mutex::new(
                            clusters.len() as i32 + 1));


                // All cluster ids are + 1, because we have the variants with abundances < eps
                // outside of the clustering algorithm already as cluster id 0
                db_clusters.par_iter().enumerate().for_each(|(cluster, index_vec)|{
                    // Deal with clustered points
                    index_vec.par_iter().for_each(|index|{
                        let info = &variant_info[*index];
                        let abundance = abundance_float[*index];

                        let mut cluster_map = cluster_map.lock().unwrap();
                        let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
                        let mut number_of_clusters = number_of_clusters.lock().unwrap();

                        let position = cluster_map
                            .entry(*info.0).or_insert(BTreeMap::new());
                        position.insert(info.1.to_string(), cluster as i32 + 1);

                        let cluster_mean = cluster_hierarchies
                            .entry(cluster as i32 + 1).or_insert(Vec::new());
                        cluster_mean.push(abundance.to_owned());
                    });
                });

                noise_points.par_iter().for_each(|noise|{

                    let mut cluster_map = cluster_map.lock().unwrap();
                    let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
                    let mut number_of_clusters = number_of_clusters.lock().unwrap();

                    // Deal with unclustered points
                    let noise_info = &variant_info[*noise];
                    let noise_abundance = abundance_float[*noise];

                    let noise_position = cluster_map.entry(*noise_info.0)
                        .or_insert(BTreeMap::new());
                    noise_position.insert(noise_info.1.to_string(), *number_of_clusters);

                    let noise_mean = cluster_hierarchies
                        .entry(*number_of_clusters).or_insert(Vec::new());
                    noise_mean.push(noise_abundance.to_owned());
                    *number_of_clusters += 1;
                });

                let mut variant_distances: Arc<Mutex<Vec<f64>>>
                    = Arc::new(
                    Mutex::new(
                        vec![0.; (variant_info_all.len().pow(2) as usize - variant_info_all.len()) / 2 as usize]));

                // produced condensed pairwise distances
                // described here: https://docs.rs/kodama/0.2.2/kodama/
                (0..variant_info_all.len())
                    .into_par_iter().enumerate().for_each(|(row_index, row_info)|{
                    let mut row_variant_set = &BTreeSet::new();
                    let row_info = &variant_info_all[row_index];
                    // lazily get the row variant read id set
                    if indels.contains_key(&row_info.0) {
                        if indels[row_info.0].contains_key(&row_info.1){
                            row_variant_set = &indels[&row_info.0][&row_info.1];
                        } else if nucfrequency.contains_key(&row_info.0) {
                            let var_char = row_info.1.as_bytes()[0] as char;
                            if nucfrequency[&row_info.0].contains_key(&var_char){
                                row_variant_set = &nucfrequency[&row_info.0][&var_char];
                            }
                        }
                    } else if nucfrequency.contains_key(&row_info.0) {
                        let var_char = row_info.1.as_bytes()[0] as char;
                        if nucfrequency[&row_info.0].contains_key(&var_char){
                            row_variant_set = &nucfrequency[&row_info.0][&var_char];
                        }
                    }
                    (row_index..variant_info_all.len())
                        .into_par_iter().enumerate().for_each(|(col_index, col_info)|{
                        let mut col_variant_set= &BTreeSet::new();
                        let col_info = &variant_info_all[col_index];
                        if indels.contains_key(&col_info.0) {
                            if indels[&col_info.0].contains_key(&col_info.1){
                                col_variant_set = &indels[&col_info.0][&col_info.1];
                            } else if nucfrequency.contains_key(&col_info.0) {
                                let var_char = col_info.1.as_bytes()[0] as char;
                                if nucfrequency[&col_info.0].contains_key(&var_char){
                                    col_variant_set = &nucfrequency[&col_info.0][&var_char];
                                }
                            }
                        } else if nucfrequency.contains_key(&col_info.0) {
                            let var_char = col_info.1.as_bytes()[0] as char;
                            if nucfrequency[&col_info.0].contains_key(&var_char){
                                col_variant_set = &nucfrequency[&col_info.0][&var_char];
                            }
                        }

                        let intersection_len = (row_variant_set
                            .intersection(&col_variant_set).collect::<HashSet<_>>().len()) as f64;

                        let union_len = (row_variant_set
                            .union(&col_variant_set).collect::<HashSet<_>>().len()) as f64;

                        // Jaccard Similarity
                        let jaccard = intersection_len / union_len;

                        // Distance between abundance values
                        let dist_f = (row_info.2 - col_info.2).abs();

                        // Distance will be defined as the mean between jaccard dist
                        // and dist_f
                        let distance = ((1. - jaccard) + dist_f) / 2.;

                        match condensed_index(
                            row_index, col_index, variant_info.len()) {
                            Some(index) => {
                                let mut variant_distances = variant_distances.lock().unwrap();
                                variant_distances[index] = distance;
                            }
                            None => {
                                debug!("No corresponding index for row {} and col {}",
                                       row_index, col_index);
                            }
                        };
                    });
                });


                let mut variant_distances = variant_distances
                    .lock()
                    .unwrap();
                debug!("Variant Distance Vector {:?}", variant_distances);


                let dend = nnchain(
                    &mut variant_distances,
              variant_info_all.len(),
                  Method::Ward.into_method_chain()
                        .expect("Incompatible linkage method"));
                debug!("Dendrogram {:?}", dend);

                let mut cluster_map = cluster_map.lock().unwrap();
                let mut cluster_hierarchies = cluster_hierarchies.lock().unwrap();
                let mut means = HashMap::new();

                for (cluster, abundances) in cluster_hierarchies.iter() {
                    let mean = abundances.iter().sum::<f64>()/abundances.len() as f64;
                    means.insert(*cluster as i32, mean);
                }

                // Combine info of both clustering methods for easy access
                let mut full_cluster_map =
                    Arc::new(
                        Mutex::new(
                            HashMap::new()));

                variant_info_all
                    .par_iter().enumerate().for_each(|(index, (position, variant, abundance))|{
                    let db_cluster = cluster_map
                        .get(position).expect("Position not found when it should be")
                        .get(variant).expect("Variant not found when it should be");
                    let mut full_cluster_map = full_cluster_map.lock().unwrap();

                    let position_map = full_cluster_map.entry(**position)
                        .or_insert(BTreeMap::new());

                    position_map.insert(variant.to_string(), (*db_cluster, index));
                });


                *dendrogram = dend;
                *clusters = full_cluster_map.lock().unwrap().clone();
                *clusters_mean = means;
                info!("{} Distinct variant frequency clusters found on contig {} at eps {}",
                      clusters_mean.len(), tid, eps);
            }
        }
    }

    fn print_variants(&mut self, ref_sequence: &Vec<u8>, sample_idx: i32){
        match self {
            PileupStats::PileupContigStats {
                indels,
                variant_abundances,
                depth,
                tid,
                genotypes_per_position,
                clusters,
                haplotypes,
                ..

            } => {
                info!("Outputting {} variant locations", variant_abundances.keys().len());
                for (position, d) in depth.iter().enumerate() {

                    // loop through each position that has variants
                    let hash = match variant_abundances.get(&(position as i32)) {
                        Some(hash) => hash,
                        None => continue,
                    };

                    let cluster_map = match clusters.get(&(position as i32)) {
                        Some(map) => map,
                        None => continue,
                    };

                    for (var, abundance) in hash.iter() {
                        // for each variant at a location
                        let indel_map = match indels.get(&(position as i32)) {
                            Some(map) => map.to_owned(),
                            None => BTreeMap::new(),
                        };
                        let cluster_val = match cluster_map.get(var) {
                            Some(i) => *i,
                            None => (-1, 0),
                        };

                        if indel_map.contains_key(var) {
                            // How does this print N for insertions?
                            if var.to_owned().contains("N"){
                                print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position,
                                       var,
                                       str::from_utf8(
                                           &ref_sequence[position..position
                                               + var.len() as usize]).unwrap(),
                                       abundance, d);

                            } else {
                                print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position,
                                       var,
                                       str::from_utf8(
                                           &[ref_sequence[position]]).unwrap(),
                                       abundance, d);
                            }

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_count) => {
                                    print!("{}\t", gtype_count);
                                },
                                None => {
                                    print!("0\t");
                                },
                            };
                            println!("{}\t{}", cluster_val.0, cluster_val.1);

                        } else if var.len() == 1 && var != &"R".to_string(){
                            print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position,
                                   var,
                                   ref_sequence[position] as char,
                                   abundance, d);

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_count) => {
                                    print!("{}\t", gtype_count);
                                },
                                None => {
                                    print!("0\t");
                                },
                            };
                            println!("{}\t{}", cluster_val.0, cluster_val.1);
                        }
                    }
                };
            }
        }
    }
}

// helper function to get the index of condensed matrix from it square form
fn condensed_index(row_i: usize, col_j: usize, n: usize) -> Option<usize>{
    if row_i == col_j {
        return None
    } else {
        if row_i < col_j {
            return Some(n*row_i - row_i*(row_i+1)/2 + col_j - 1 - row_i)
        } else {
            return Some(n*col_j - col_j*(col_j+1)/2 + row_i - 1 - col_j)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotypes_two_positions() {
        let mut contig = PileupStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i32>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());

        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [0, 1, 2].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [3, 4, 5].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                   [0, 1].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [2, 3].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [4, 5].iter().cloned().collect::<BTreeSet<i32>>());

        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          "mean",
                          vec![2., 6., 0.],
                          vec![6,0]);

        match contig {
            PileupStats::PileupContigStats {
                ref mut read_error_rate,
                ..
            } => {
                *read_error_rate = 0.0;
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);
        let genotypes = contig.generate_minimum_genotypes();
        println!("{:?}", genotypes);

        assert_eq!(genotypes[&0], 4);
        assert_eq!(genotypes[&1], 4);
    }

    #[test]
    fn test_genotypes_three_positions() {
        let mut contig = PileupStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i32>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());
        nuc_freq.insert(2, BTreeMap::new());

        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [0, 1, 2].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [3, 4, 5].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                   [0, 1].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [2, 3].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [4, 5].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(2).or_insert(BTreeMap::new());

        pos.insert("G".chars().collect::<Vec<char>>()[0],
                   [1, 2, 4, 5].iter().cloned().collect::<BTreeSet<i32>>());


        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          "mean",
                          vec![2., 6., 0.],
                          vec![6,4,0]);

        match contig {
            PileupStats::PileupContigStats {
                ref mut read_error_rate,
                ..
            } => {
                *read_error_rate = 0.0;
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);
        let genotypes = contig.generate_minimum_genotypes();
        println!("{:?}", genotypes);

        assert_eq!(genotypes[&0], 4);
        assert_eq!(genotypes[&1], 4);
        assert_eq!(genotypes[&2], 3);
    }

    #[test]
    fn test_genotypes_four_positions() {
        let mut contig = PileupStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i32>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());
        nuc_freq.insert(2, BTreeMap::new());
        nuc_freq.insert(3, BTreeMap::new());

        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [0].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                           [5].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                           [0, 2].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [1].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                           [3, 5].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("C".chars().collect::<Vec<char>>()[0],
                           [4].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(2).or_insert(BTreeMap::new());

        pos.insert("G".chars().collect::<Vec<char>>()[0],
                           [0, 1, 2, 3].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [4, 5].iter().cloned().collect::<BTreeSet<i32>>());

        let mut pos = nuc_freq.entry(3).or_insert(BTreeMap::new());

        pos.insert("C".chars().collect::<Vec<char>>()[0],
                           [2].iter().cloned().collect::<BTreeSet<i32>>());
        pos.insert("T".chars().collect::<Vec<char>>()[0],
                           [3, 4, 5].iter().cloned().collect::<BTreeSet<i32>>());


        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          "mean",
                          vec![2., 6., 0.],
                          vec![3,3,0,-6]);

        match contig {
            PileupStats::PileupContigStats {
                ref mut read_error_rate,
                ..
            } => {
                *read_error_rate = 0.0;
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);
        let genotypes = contig.generate_minimum_genotypes();
        println!("{:?}", genotypes);

        assert_eq!(genotypes[&0], 2);
        assert_eq!(genotypes[&1], 5);
        assert_eq!(genotypes[&2], 5);
        assert_eq!(genotypes[&3], 4);
    }
}