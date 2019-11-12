use std::collections::{HashMap, HashSet};
use std::collections::BTreeMap;
use std::str;
use std::sync::{Arc, Mutex};
use std::io::prelude::*;
use rayon::prelude::*;
use permutation::*;
use codon_structs::*;
use bio_types::strand;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use rusty_machine::learning::dbscan::DBSCAN;
use rusty_machine::learning::UnSupModel;
use rusty_machine::linalg::Matrix;



#[derive(Debug, Clone)]
pub struct Genotype {
    read_ids: HashSet<i32>,
    start_var_pos: usize,
    ordered_variants: HashMap<i32, String>,
    frequencies: Vec<f64>,
}

impl Genotype {
    fn start(position: usize) -> Genotype {
        Genotype {
            read_ids: HashSet::new(),
            start_var_pos: position,
            ordered_variants: HashMap::new(),
            frequencies: Vec::new(),
        }
    }

    fn new(&mut self, read_id: i32, position_map: &BTreeMap<i32, String>, variant_abundances: &Vec<HashMap<String, f64>>) {
        self.read_ids = HashSet::new();
        self.read_ids.insert(read_id);
        for (pos, variant) in position_map.iter() {
            let variant_map = &variant_abundances[*pos as usize];
            if variant_map.len() > 0 {
                debug!("Fetching position: {:?} {:?} {:?}", pos, variant_map, variant);
                self.frequencies.push(variant_map[variant].clone());
                self.ordered_variants.insert(pos.clone(), variant.to_string());
            }
        }
    }

    fn check_variants(&self, position_map: &BTreeMap<i32, String>) -> bool {
        // check variants against stored variants for a genotype
        let mut new_var= false;
        for (check_pos, check_var) in position_map.iter() {
            if self.ordered_variants.contains_key(&check_pos) {
                let current_var = match self
                    .ordered_variants
                    .get(&check_pos) {
                    Some(var) => var,
                    None => {
                        println!("Position not recorded in variant map");
                        std::process::exit(1)
                    }
                };
                if current_var != check_var {
                    //Then this is a new genotype
                    new_var = true;
                    break
                } else {
                    new_var = false;
                }
            } else {
                // current variant does not contain this position, must be new genotype
                new_var = true;
                break
            }
        }
        return new_var
    }
}

pub enum PileupStats {
    PileupContigStats {
        nucfrequency: Vec<HashMap<char, HashSet<i32>>>,
        variants_in_reads: HashMap<i32, BTreeMap<i32, String>>,
        variant_abundances: Vec<HashMap<String, f64>>,
        variant_count: Vec<f64>,
        depth: Vec<f64>,
        indels: Vec<HashMap<String, HashSet<i32>>>,
        genotypes_per_position: HashMap<usize, HashMap<String, usize>>,
        mean_genotypes: f32,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: f32,
        variations_per_base: f32,
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
    }
}

impl PileupStats {
    pub fn new_contig_stats(min: f32, max: f32,
                            contig_end_exclusion: u32) -> PileupStats {
        PileupStats::PileupContigStats {
            nucfrequency: vec!(),
            variants_in_reads: HashMap::new(),
            variant_abundances: Vec::new(),
            variant_count: Vec::new(),
            depth: vec!(),
            indels: vec!(),
            genotypes_per_position: HashMap::new(),
            mean_genotypes: 0.0,
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0.0,
            variations_per_base: 0.00,
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
        }
    }
}

pub trait PileupFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
                  indels_positions: Vec<HashMap<String, HashSet<i32>>>,
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
                               original_contig: Vec<u8>,
                               consensus_genome: std::fs::File);

    fn generate_minimum_genotypes(&mut self);

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
                ..
            } => {
                *nucfrequency = vec!();
                *variants_in_reads = HashMap::new();
                *variant_abundances = Vec::new();
                *variant_count = Vec::new();
                *depth = vec!();
                *indels = vec!();
                *tid = 0;
                *total_indels = 0;
                *target_name = vec!();
                *target_len = 0.0;
                *variations_per_base = 0.00;
                *coverage = 0.00;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
        }
    }

    fn add_contig(&mut self, nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
                  indel_positions: Vec<HashMap<String, HashSet<i32>>>,
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
                let mut variant_count_safe = Arc::new(Mutex::new(vec![0.; ups_and_downs.len()]));
                *depth = vec![0.; ups_and_downs.len()];
                let mut cumulative_sum = 0;
                for (pos, current) in ups_and_downs.iter().enumerate() {
                    cumulative_sum += *current;
                    depth[pos] = cumulative_sum as f64;
                }
                let mut adjusted_depth = Arc::new(Mutex::new(depth.clone()));

                // Calculate how many reads have variant at each position
                // to go into linear regression predicting read error rate
                nucfrequency.par_iter().zip(indels.par_iter()).enumerate().for_each(
                    |(pos, (snp_map, indel_map))|{
                        let mut variant_count_safe = variant_count_safe.lock().unwrap();
                        let mut adjusted_depth = adjusted_depth.lock().unwrap();
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
                let mut data = vec![("Y", variant_count.clone()), ("X", depth.clone())];
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
                let variants = Arc::new(Mutex::new(vec![HashMap::new(); depth.len()])); // The relative abundance of each variant
                let read_variants = Arc::new(Mutex::new(HashMap::new())); // The reads with variants and their positions
                let variant_count = Arc::new(Mutex::new(0));
                let indels_backup = Arc::new(Mutex::new(indels.clone()));
                let nucfrequency_backup = Arc::new(Mutex::new(nucfrequency.clone()));

                // for each location calculate if there is a variant based on read depth
                // Uses rayon multithreading
                depth.par_iter_mut().enumerate().for_each(|(i, d)| {
                    let read_variants = Arc::clone(&read_variants);
                    let variant_count = Arc::clone(&variant_count);
                    let mut rel_abundance = HashMap::new();
                    if (*coverage * (1.0 - coverage_fold) <= *d as f32
                        && *d as f32 <= *coverage * (1.0 + coverage_fold))
                        || (coverage_fold == 0.0) {
//                        if d >= &mut min_variant_depth.clone() {
                            // INDELS act differently to normal variants
                            // The reads containing this variant don't contribute to coverage values
                            // So we need to readjust the depth at these locations to show true
                            // read depth. i.e. variant depth + reference depth
                            if indels[i].len() > 0 {
                                for (indel, read_ids) in indels[i].iter() {
                                    let count = read_ids.len();
                                    *d += count as f64;
                                    if (count >= min_variant_depth) & (count as f64 / *d > *read_error_rate){
                                        rel_abundance.insert(indel.clone(), count as f64 / *d);
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
                                        let mut indels_backup
                                            = indels_backup.lock().unwrap();
                                        indels_backup[i].remove(indel);
                                    }
                                }
                            }
                            if nucfrequency[i].len() > 0 {
                                for (base, read_ids) in nucfrequency[i].iter() {
                                    let count = read_ids.len();
                                    if (count >= min_variant_depth) & (count as f64 / *d > *read_error_rate) {
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
                                        let mut nucfrequency_backup
                                            = nucfrequency_backup.lock().unwrap();
                                        nucfrequency_backup[i].remove(base);
                                    }
                                }
                            };

//                        }
                    }

                    if rel_abundance.len() > 1 {
                        let mut variants = variants.lock().unwrap();
                        debug!("Relative Abundances: {:?}", rel_abundance);
                        variants[i] = rel_abundance;
                    } else {

                    }

                });

                let read_variants = read_variants.lock().unwrap();
                *variants_in_reads = read_variants.clone();
                let variants = variants.lock().unwrap();
                *variant_abundances = variants.clone();
                let variant_count = variant_count.lock().unwrap();
                debug!("Total variants for {}: {:?}", tid, variant_count);
                *variations_per_base = (*variant_count) as f32/target_len.clone() as f32;
                let nucfrequency_backup = nucfrequency_backup.lock().unwrap();
                *nucfrequency = nucfrequency_backup.clone();
                let indels_backup = indels_backup.lock().unwrap();
                *indels = indels_backup.clone();
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
                               original_contig: Vec<u8>,
                               mut consensus_genome: std::fs::File){
        match self {
            PileupStats::PileupContigStats {
                ref mut variant_abundances,
                ..
            } => {
                let mut contig = String::new();

                let mut skip_n = 0;
                let mut skip_cnt = 0;
                // Generate the consensus genome by checking each variant
                // Variant has to be in more than 0.5 of population
                for (pos, base) in original_contig.iter().enumerate() {
                    if skip_cnt < skip_n {
                        skip_cnt += 1;
                    } else {
                        let mut max_var = "";
                        let mut max_abund = 0.0;
                        skip_n = 0;
                        skip_cnt = 0;
                        if variant_abundances[pos].len() > 0 {
                            let hash = &variant_abundances[pos];
                            for (var, abundance) in hash.iter() {
                                if abundance > &max_abund {
                                    max_var = var;
                                    max_abund = *abundance;
                                }
                            }
                            if max_abund >= 0.5{
                                if max_var.contains("N") {
                                    skip_n = max_var.len() - 1;
                                    skip_cnt = 0;
                                } else {
                                    contig = contig + max_var;
                                }
                            } else {
                                contig = contig + str::from_utf8(&[*base]).unwrap();
                            }
                        } else {
                            contig = contig + str::from_utf8(&[*base]).unwrap();
                        }
                    }
                };
                contig = contig + "\n";
                match consensus_genome.write_all(contig.as_bytes()) {
                    Ok(consensus_genome) => consensus_genome,
                    Err(e) => {
                        println!("Cannot write to file {:?}", e);
                        std::process::exit(1)}
                };
            }
        }
    }

    fn generate_minimum_genotypes(&mut self) {
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

                variant_abundances.par_iter().enumerate().for_each(|(position, variants)| {
                    // For each variant we calculate the minimum number of genotypes possible
                    // based on variants in reads mapping to this variant location
                    if variants.len() > 0 {
                        let mut genotype_record;

                        let genotypes = Arc::clone(&genotypes);
                        let variant_count = Arc::clone(&variant_count);
                        let total_genotype_count = Arc::clone(&total_genotype_count);

                        let mut genotype_pos = HashMap::new();

                        for (var, _abundance) in variants.iter() {
                            let genotype_count = genotype_pos.entry(var.to_string())
                                .or_insert(0);

                            let mut genotype_vec = Vec::new();

                            genotype_record = Genotype::start(position);
                            let mut read_ids = HashSet::new();
                            if indels[position].contains_key(var) {
                                read_ids =
                                    match indels[position].get(var) {
                                        Some(ids) => ids.clone(),
                                        None => {
                                            println!("Variant not in indel hash");
                                            std::process::exit(1)
                                        },
                                    };
                            } else if nucfrequency[position].contains_key(
                                &(var.clone().into_bytes()[0] as char)) {
                                read_ids =
                                    match nucfrequency[position]
                                        .get(&(var.clone().into_bytes()[0] as char)) {
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
                            for read_id in read_vec_sorted.iter() {
                                let position_map = match variants_in_reads.get(read_id) {
                                    Some(positions) => positions,
                                    None => {
                                        debug!("read id not recorded in variant map {}, {}", var, read_id);
                                        break
                                    },
                                };

                                if genotype_vec.len() == 0 {
                                    // No genotype observed yet, so create one
                                    debug!("Read ID : {} Position Map: {:?}", read_id, position_map);
                                    genotype_record.new(*read_id, position_map, &variant_abundances);

                                    genotype_vec.push(genotype_record.clone());
                                } else {
//                                    let position_map_variants: Vec<String> = position_map.values().cloned().collect();
                                    let position_map_positions: Vec<i32> = position_map.keys().cloned().collect();
                                    let position_set = position_map.keys().cloned().collect::<HashSet<i32>>();

                                    let mut new_genotype = false;

                                    'outer: for genotype in genotype_vec.iter_mut() {

                                        // Create HashSets of variant positions to use intersection
                                        let genotype_position_set =
                                            genotype.ordered_variants.keys().cloned().collect::<HashSet<i32>>();

                                        let diff: Vec<i32> = genotype_position_set
                                            .symmetric_difference(&position_set).cloned().collect();

                                        if diff.len() > 0 {
                                            // Positional difference found
                                            // Check if new genotype
                                            for pos in diff.iter() {
                                                if (genotype.ordered_variants.keys().min() < Some(pos))
                                                    && (Some(pos) < genotype.ordered_variants.keys().max()) {
                                                    // possible new genotype detected
                                                    genotype_record.new(*read_id, position_map, &variant_abundances);
                                                    break 'outer;
                                                }
                                            }
                                            if (genotype.ordered_variants.keys().min() == diff.iter().min())
                                                || (diff.iter().max() > genotype.ordered_variants.keys().max()) {
                                                // check variants against stored variants for a genotype
                                                new_genotype = genotype.check_variants(position_map);
                                                if new_genotype {
                                                    // possible new genotype detected
                                                    genotype_record.new(*read_id, position_map, &variant_abundances);
                                                } else {
                                                    // Extension of previous genotype
                                                    genotype.read_ids.insert(*read_id);

                                                    for (new_position, new_variant) in position_map.iter() {
                                                        if !(genotype.ordered_variants.contains_key(&new_position)) {
                                                            let variant_map = &variant_abundances[*new_position as usize];
                                                            genotype.ordered_variants
                                                                .insert(new_position.clone(), new_variant.clone());
                                                            genotype_record.frequencies.push(
                                                                variant_map[new_variant])
                                                        }
                                                    };

                                                    // reset genotype_record
                                                    genotype_record = Genotype::start(position);
                                                    genotype.read_ids.insert(*read_id);
                                                    break
                                                }
                                            }
                                        } else {
                                            // check variants against stored variants for a genotype
                                            new_genotype = genotype.check_variants(position_map);

                                            if new_genotype {
                                                // possible new genotype detected
                                                genotype_record.new(*read_id, position_map, &variant_abundances);

                                            } else {
                                                // No difference with a previous genotype, reset current
                                                genotype_record = Genotype::start(position);
                                                genotype.read_ids.insert(*read_id);
                                                break
                                            }
                                        }
                                    }
                                    if genotype_record.ordered_variants.len() > 0 {
                                        // New genotype detected
                                        genotype_vec.push(genotype_record);
//                                        debug!("genotypes {:?}", genotype_vec);

                                        genotype_record = Genotype::start(position);
                                    }
                                }
                            }

                            *genotype_count += genotype_vec.len();
                            let mut total_genotype_count = total_genotype_count.lock().unwrap();
                            *total_genotype_count += genotype_vec.len();
                            let mut variant_count = variant_count.lock().unwrap();
                            *variant_count += 1;
                        }
                        let mut genotypes = genotypes.lock().unwrap();
                        genotypes.insert(position, genotype_pos);
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
                tid,
                target_name,
                depth,
                ..
            } => {
                let contig_name = String::from_utf8(target_name.clone())
                    .expect("Cannot create string from target_name");
                let placeholder = Vec::new();
                let mut gff_records = match gff_map.get(&contig_name){
                    Some(records) => records,
                    None => &placeholder,
                };
                debug!("Calculating population dN/dS from reads for {} genes", gff_records.len());
                let mut print_stream = &mut Mutex::new(std::io::stdout());
                gff_records.par_iter().enumerate().for_each(|(id, gene)| {
                    let dnds = codon_table.find_mutations(gene, variant_abundances, ref_sequence, depth);
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

                    for (cursor, variant_map) in
                        variant_abundances[start..end].to_vec().iter().enumerate() {
                        let cursor = cursor + start;
                        if variant_map.len() > 0 {
                            let mut print_stream = print_stream.lock().unwrap();
                            for (variant, abundance) in variant_map {
                                write!(print_stream, "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t",
                                         contig.clone()+gene_id, gene.start(),
                                         gene.end(), frame, strand_symbol, dnds, cursor);
                                if variant.to_owned().contains("N") {
                                    writeln!(print_stream, "{}\t{}\t{}\t{}",
                                           variant,
                                           str::from_utf8(
                                               &ref_sequence[cursor..cursor
                                                   + variant.len() as usize]).unwrap(),
                                           abundance, depth[cursor]);

                                } else if variant.len() > 1 {
                                     writeln!(print_stream,"{}\t{}\t{:.3}\t{}",
                                           str::from_utf8(
                                               &[ref_sequence[cursor-1]]).unwrap().to_owned() + variant,
                                           str::from_utf8(
                                               &[ref_sequence[cursor-1]]).unwrap(),
                                           abundance, depth[cursor]);
                                } else {
                                    writeln!(print_stream, "{}\t{}\t{:.3}\t{}",
                                             variant,
                                             ref_sequence[cursor] as char, abundance, depth[cursor]);
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
                genotypes_per_position,
                depth,
                ..} => {

                let mut abundance_lnddivg =
                    Arc::new(
                    Mutex::new(
                        Vec::new()));

                let mut variant_info =
                    Arc::new(
                        Mutex::new(
                            Vec::new()));

                variant_abundances.iter().enumerate().for_each(
                    |(position, hash)|{
                        // loop through each position that has variants
                        let mut abundance_lnddivg = abundance_lnddivg.lock().unwrap();
                        let mut variant_info = variant_info.lock().unwrap();

                        let d = depth[position];

                        for (var, abundance) in hash.iter() {
                            variant_info.push((position, var.clone()));
                            // for each variant at a location

                            // genotypes associated with that position and variant
                            let mut g_count = 0.;
                            match genotypes_per_position.get(&position) {
                                Some(gtype_hash) => {
                                    match gtype_hash.get(var) {
                                        Some(gtype_count) => {
                                            g_count = gtype_count.clone() as f64;
                                        },
                                        None => {
                                            g_count = 0.;
                                        }
                                    }
                                },
                                None => {
                                    g_count = 0.;
                                },
                            };
                            abundance_lnddivg.push(abundance.clone());
                            abundance_lnddivg.push((d/g_count).ln());
                        }
                });
                let abundance_lnddivg = abundance_lnddivg.lock().unwrap();
                let variant_info = variant_info.lock().unwrap();
                // Put inputs into matrix form
                let inputs = Matrix::new(
                    abundance_lnddivg.len()/2, 2, abundance_lnddivg.to_vec());

                let mut model = DBSCAN::new(0.1, 2);
                model.train(&inputs).unwrap();
                let clustering = model.clusters().unwrap();
                for (cluster, info) in clustering.iter().zip(variant_info.iter()){
                    let cluster = match cluster {
                        Some(c) => *c as i32,
                        None => -1,
                    };
//                    println!("Cluster: {} Pos: {} Var: {}", cluster, info.0, info.1);
                }
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
                ..

            } => {
                variant_abundances.iter().enumerate().for_each(|(position, hash)|{
                    // loop through each position that has variants
//                    let position = *position as usize;
                    let d = depth[position];

                    for (var, abundance) in hash.iter() {
                        // for each variant at a location
                        if indels[position].contains_key(var) {
                            // How does this print N for insertions?
                            if var.to_owned().contains("N"){
                                print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position,
                                       var,
                                       str::from_utf8(
                                           &ref_sequence[position..position
                                               + var.len() as usize]).unwrap(),
                                       abundance, d);

                            } else {
                                print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position-1,
                                       str::from_utf8(
                                           &[ref_sequence[position-1]]).unwrap().to_owned() + var,
                                       str::from_utf8(
                                           &[ref_sequence[position-1]]).unwrap(),
                                       abundance, d);
                            }

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_hash) => {
                                    match gtype_hash.get(&var.to_string()) {
                                        Some(gtype_count) => {
                                            print!("{}\t", gtype_count);
                                        },
                                        None => {
                                            print!("0\t");
                                        }
                                    }
                                },
                                None => {
                                    print!("0\t");
                                },
                            };
                            println!{"{}", sample_idx};

                        } else if var.len() == 1{
                            print!("{}\t{}\t{}\t{}\t{:.3}\t{}\t", tid, position,
                                   var,
                                   ref_sequence[position] as char,
                                   abundance, d);

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_hash) => {
                                    match gtype_hash.get(&var.to_string()) {
                                        Some(gtype_count) => {
                                            print!("{}\t", gtype_count);
                                        },
                                        None => {
                                            print!("0\t");
                                        }
                                    }
                                },
                                None => {
                                    print!("0\t");
                                },
                            };
                            println!{"{}", sample_idx};
                        }
                    }
                });
            }
        }
    }
}
