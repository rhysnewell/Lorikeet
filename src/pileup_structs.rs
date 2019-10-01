use std::collections::{HashMap, HashSet};
use std::collections::BTreeMap;
//use ndarray::Array2;
//use simhash::*;
//use mash::*;
//use distance::*;
//use itertools::Itertools;
use std::str;
//use std::fs::File;
use std::sync::{Arc, Mutex};
use std::io::prelude::*;
use rayon::prelude::*;
use permutation::*;
//use std::iter::FromIterator;
//use rust_htslib::bam::record::{Cigar, CigarStringView};

#[derive(Debug, Clone)]
pub struct Genotype {
    read_ids: HashSet<i32>,
    base_positions: Vec<i32>,
    start_var_pos: usize,
    ordered_variants: HashMap<i32, String>,
}

pub enum PileupStats {
    PileupContigStats {
        nucfrequency: Vec<HashMap<char, HashSet<i32>>>,
        variants_in_reads: HashMap<i32, BTreeMap<i32, String>>,
        variant_abundances: BTreeMap<i32, HashMap<String, f32>>,
        depth: Vec<usize>,
        indels: Vec<HashMap<String, HashSet<i32>>>,
        genotypes_per_position: HashMap<usize, HashMap<String, usize>>,
        mean_genotypes: f32,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: usize,
        variations_per_base: f32,
        coverage: f32,
        variance: f32,
        observed_contig_length: u32,
        num_covered_bases: i32,
        contig_end_exclusion: u32,
        min_fraction_covered_bases: f32,
        min: f32,
        max: f32,
    }
}

impl PileupStats {
    pub fn new_contig_stats(min: f32, max: f32, min_fraction_covered_bases: f32,
                            contig_end_exclusion: u32) -> PileupStats {
        PileupStats::PileupContigStats {
            nucfrequency: vec!(),
            variants_in_reads: HashMap::new(),
            variant_abundances: BTreeMap::new(),
            depth: vec!(),
            indels: vec!(),
            genotypes_per_position: HashMap::new(),
            mean_genotypes: 0.0,
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0,
            variations_per_base: 0.00,
            coverage: 0.00,
            variance: 0.00,
            observed_contig_length: 0,
            num_covered_bases: 0,
            contig_end_exclusion: contig_end_exclusion,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max,
        }
    }
}

pub trait PileupFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
                  read_depth: Vec<usize>,
                  indels_positions: Vec<HashMap<String, HashSet<i32>>>,
                  tid: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize);

    fn calc_variants(&mut self,
                     depth_thresh: usize,
                     variant_fraction: f64);

    fn generate_variant_matrix(&mut self);

    fn generate_variant_contig(&mut self,
                               original_contig: Vec<u8>,
                               consensus_genome: std::fs::File);

    fn generate_genotypes(&mut self);

    fn calc_coverage(&mut self) -> f32;

    fn print_variants(&mut self, ref_sequence: Vec<u8>, depth_thresh: usize);
}

impl PileupFunctions for PileupStats {
    fn setup(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variants_in_reads,
                ref mut variant_abundances,
                ref mut depth,
                ref mut indels,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut variations_per_base,
                ref mut coverage,
                ref mut num_covered_bases,
                ..
            } => {
                *nucfrequency = vec!();
                *variants_in_reads = HashMap::new();
                *variant_abundances = BTreeMap::new();
                *depth = vec!();
                *indels = vec!();
                *tid = 0;
                *total_indels = 0;
                *target_name = vec!();
                *target_len = 0;
                *variations_per_base = 0.00;
                *coverage = 0.00;
                *num_covered_bases = 0;
            }
        }
    }

    fn add_contig(&mut self, nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
                  read_depth: Vec<usize>,
                  indel_positions: Vec<HashMap<String, HashSet<i32>>>,
                  target_id: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut depth,
                ref mut indels,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ..
            } => {
                *nucfrequency = nuc_freq;
                *depth = read_depth;
                *indels = indel_positions;
                *tid = target_id;
                *total_indels = total_indels_in_contig;
                *target_name = contig_name;
                *target_len = contig_len;

            }
        }
    }

    fn calc_variants(&mut self, depth_thresh: usize, variant_fraction: f64){
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variants_in_reads,
                ref mut variant_abundances,
                depth,
                ref mut indels,
                total_indels,
                target_len,
                ref mut variations_per_base,
                ref mut coverage,
                ..
            } => {
                let mut variants = Arc::new(Mutex::new(BTreeMap::new())); // The relative abundance of each variant
                let mut read_variants = Arc::new(Mutex::new(HashMap::new())); // The reads with variants and their positions
                let mut variant_count = Arc::new(Mutex::new(0));

                // for each location calculate if there is a variant based on read depth
                // Uses rayon multithreading
                depth.into_par_iter().enumerate().for_each(|(i, d)| {
                    let mut variant = Arc::clone(&variants);
                    let mut read_variants = Arc::clone(&read_variants);
                    let mut variant_count = Arc::clone(&variant_count);
                    let mut rel_abundance = HashMap::new();
                    if *coverage * 0.75 < *d as f32 && *d as f32 <= *coverage * 1.25 {
                        if d >= &mut depth_thresh.clone() {
                            if nucfrequency[i].len() > 0 {
                                for (base, read_ids) in nucfrequency[i].iter() {
                                    let count = read_ids.len();
                                    if count as f64 / depth_thresh as f64 >= variant_fraction {
                                        rel_abundance.insert(base.to_string(), count as f32 / *d as f32);
                                        for read in read_ids {
                                            let mut read_variants = read_variants.lock().unwrap();
                                            let read_vec = read_variants
                                                .entry(read.clone())
                                                .or_insert(BTreeMap::new());
                                            read_vec.insert(i as i32, base.to_string());
                                        }
                                    }
                                    let mut variant_count = variant_count.lock().unwrap();
                                    *variant_count += 1;

                                }
                            };
                            if indels[i].len() > 0 {
                                for (indel, read_ids) in indels[i].iter() {
                                    let count = read_ids.len();
                                    if count as f64 / depth_thresh as f64 >= variant_fraction {
                                        rel_abundance.insert(indel.clone(), count as f32 / *d as f32);
                                        for read in read_ids {
                                            let mut read_variants = read_variants.lock().unwrap();

                                            let read_vec = read_variants
                                                .entry(read.clone())
                                                .or_insert(BTreeMap::new());
                                            read_vec.insert(i as i32, indel.clone());
                                        }
                                    }
                                    let mut variant_count = variant_count.lock().unwrap();
                                    *variant_count += 1;
                                }
                            }
                        }
                    }

                    if rel_abundance.len() > 0 {
                        let mut variants = variants.lock().unwrap();
                        variants.insert(i as i32, rel_abundance);
                    }

                });

                let mut read_variants = read_variants.lock().unwrap();
                debug!("read variants {:?}", read_variants);
                *variants_in_reads = read_variants.clone();
                let mut variants = variants.lock().unwrap();
                debug!("variants abundances {:?}", variants);
                *variant_abundances = variants.clone();
                let mut variant_count = variant_count.lock().unwrap();
                *variations_per_base = (*variant_count+*total_indels as i32) as f32/target_len.clone() as f32;
            }
        }
    }

    fn generate_variant_matrix(&mut self){
        match self {
            PileupStats::PileupContigStats {
                variants_in_reads: _,
                target_len: _,
                ..
            } => {
//                let mut position_covariance = BTreeMap::new();
//                let mut covar_array = Array2::<usize>::zeros((*target_len, *target_len));
//
//                for (_read_id, positions) in variants_in_reads{
//                    for position in positions.keys().collect().iter().combinations(2){
//
//                        covar_array[[*position[0] as usize, *position[1] as usize]] += 1 as usize;
//                        covar_array[[*position[1] as usize, *position[0] as usize]] += 1 as usize;
//
//                    }
//                }
////                println!("{:?}", covar_array);
//                for row in covar_array.genrows(){
////                    for col in row{
////                        print!("{},")
////                    }
//                    println!("{}", row.iter().fold(String::new(), |acc, &num| acc + &num.to_string() + "\t"));
////                    row.iter().fold(String::new(),|mut s,&n| {print!(format!((s,"{},",n).ok()); s});
//                }
            }
        }
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
                        if variant_abundances.contains_key(&(pos as i32)) {
                            let hash = variant_abundances.get(&(pos as i32)).unwrap();
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

    fn generate_genotypes(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut variant_abundances,
                ref mut nucfrequency,
                ref mut indels,
                ref mut variants_in_reads,
                ref mut genotypes_per_position,
                ref mut mean_genotypes,
                ..
            } => {
                let mut genotypes =
                    Arc::new(Mutex::new(HashMap::new()));
                let mut variant_count =
                    Arc::new(Mutex::new(0));
                let mut total_genotype_count =
                    Arc::new(Mutex::new(0));

                variant_abundances.par_iter().for_each(|(position, variants)| {
                    let position = *position as usize;
                    let mut genotype_record;

                    let mut genotypes = Arc::clone(&genotypes);
                    let mut variant_count = Arc::clone(&variant_count);
                    let mut total_genotype_count = Arc::clone(&total_genotype_count);

                    let mut genotype_pos = HashMap::new();

                    for (var, _abundance) in variants.iter() {
                        let genotype_count = genotype_pos.entry(var.to_string())
                            .or_insert(0);

                        let mut genotype_vec = Vec::new();

                        genotype_record = Genotype {
                            read_ids: HashSet::new(),
                            base_positions: Vec::new(),
                            start_var_pos: position,
                            ordered_variants: HashMap::new(),
                        };
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

                        } else if var.len() == 1 {
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
                        let mut read_vec = read_ids.into_iter().collect::<Vec<i32>>();

                        for read_id in read_vec.iter() {
                            let mut position_map = match variants_in_reads.get(read_id) {
                                Some(positions) => positions,
                                None => {
                                    debug!("read id not recorded in variant map {}, {}", var, read_id);
                                    break
                                },
                            };
                            left_most_variants.push(
                                *position_map.keys().cloned().collect::<Vec<i32>>().iter().min().unwrap());
                        }
                        let permuted = permutation::sort(&left_most_variants[..]);
                        let read_vec_sorted = permuted.apply_slice(&read_vec[..]);

                        debug!("read vec {:?}", read_vec_sorted);
                        for read_id in read_vec_sorted.iter() {
                            let mut position_map = match variants_in_reads.get(read_id) {
                                Some(positions) => positions,
                                None => {
                                    debug!("read id not recorded in variant map {}, {}", var, read_id);
                                    break
                                },
                            };
                            debug!("read id {} positions {:?}", read_id, position_map);
                            if genotype_vec.len() == 0 {
                                genotype_record.read_ids.insert(*read_id);
                                for (pos, variant) in position_map.iter(){
                                    genotype_record.base_positions.push(pos.clone());
                                    genotype_record.ordered_variants.insert(pos.clone(), variant.to_string());
                                }
                                genotype_record.base_positions.sort();
                                genotype_vec.push(genotype_record.clone());

                            } else {
//                                    let position_map_variants: Vec<String> = position_map.values().cloned().collect();
                                let position_map_positions: Vec<i32> = position_map.keys().cloned().collect();
                                let position_set = position_map.keys().cloned().collect::<HashSet<i32>>();

                                let mut new_genotype = false;

                                for genotype in genotype_vec.iter_mut(){
                                    let genotype_position_set =
                                        genotype.base_positions.iter().cloned().collect::<HashSet<i32>>();

                                    let diff: Vec<i32> = genotype_position_set
                                        .symmetric_difference(&position_set).cloned().collect();
                                    if diff.len() > 0 {
                                        // Positional difference found
                                        // Check if new genotype
                                        for pos in diff.iter() {

                                            if (genotype.base_positions.iter().min() < Some(pos))
                                                && (Some(pos) < genotype.base_positions.iter().max()) {
                                                // possible new genotype detected
                                                genotype_record.read_ids = HashSet::new();
                                                genotype_record.read_ids.insert(*read_id);
                                                for (pos, variant) in position_map.iter(){
                                                    genotype_record.base_positions.push(pos.clone());
                                                    genotype_record.ordered_variants.insert(pos.clone(), variant.to_string());
                                                }
                                                genotype_record.base_positions.sort();
//                                                    new_genotype = true;
                                                break
                                            }
                                        }
                                        if (genotype.base_positions.iter().min() == diff.iter().min())
                                            || (diff.iter().max() > genotype.base_positions.iter().max()) {
                                            // check variants against stored variants for a genotype
                                            for (check_pos, check_var) in position_map.iter() {
                                                if genotype.ordered_variants.contains_key(&check_pos) {
                                                    let current_var = match genotype
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
                                                        new_genotype = true;
                                                    }
                                                }
                                            }
                                            if new_genotype{
                                                // possible new genotype detected
                                                genotype_record.read_ids = HashSet::new();
                                                genotype_record.read_ids.insert(*read_id);
                                                for (pos, variant) in position_map.iter(){
                                                    genotype_record.base_positions.push(pos.clone());
                                                    genotype_record.ordered_variants.insert(pos.clone(), variant.to_string());
                                                }
                                                genotype_record.base_positions.sort();
                                            } else {
                                                // Extension of previous genotype
                                                genotype.read_ids.insert(*read_id);
                                                for base_position in position_map_positions.iter(){
                                                    if !(genotype.base_positions.contains(&base_position)){
                                                        genotype.base_positions.push(*base_position);
                                                    }
                                                };
                                                genotype.base_positions.sort();
                                                for (new_position, new_variant) in position_map.iter() {
                                                    if !(genotype.ordered_variants.contains_key(&new_position)) {
                                                        genotype.ordered_variants
                                                                .insert(new_position.clone(), new_variant.clone());
                                                    }
                                                };

                                                // reset genotype_record
                                                genotype_record = Genotype {
                                                    read_ids: HashSet::new(),
                                                    base_positions: Vec::new(),
                                                    start_var_pos: position,
                                                    ordered_variants: HashMap::new(),
                                                };
                                                genotype.read_ids.insert(*read_id);
                                                break
                                            }
                                        }

                                    } else {
                                        // check variants against stored variants for a genotype
                                        for (check_pos, check_var) in position_map.iter() {
                                            if genotype.ordered_variants.contains_key(&check_pos) {
                                                let current_var = match genotype
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
                                                    new_genotype = true;
                                                }
                                            }
                                        }
                                        if new_genotype{
                                            // possible new genotype detected
                                            genotype_record.read_ids = HashSet::new();
                                            genotype_record.read_ids.insert(*read_id);
                                            for (pos, variant) in position_map.iter(){
                                                genotype_record.base_positions.push(pos.clone());
                                                genotype_record.ordered_variants.insert(pos.clone(), variant.to_string());
                                            }
                                            genotype_record.base_positions.sort();
                                        } else {
                                            // No difference with a previous genotype, reset current
                                            genotype_record = Genotype {
                                                read_ids: HashSet::new(),
                                                base_positions: Vec::new(),
                                                start_var_pos: position,
                                                ordered_variants: HashMap::new(),
                                            };
                                            genotype.read_ids.insert(*read_id);
                                            break
                                        }
                                    }
                                }
                                if genotype_record.base_positions.len() > 0 {
                                    genotype_vec.push(genotype_record);

                                    genotype_record = Genotype {
                                        read_ids: HashSet::new(),
                                        base_positions: Vec::new(),
                                        start_var_pos: position,
                                        ordered_variants: HashMap::new(),
                                    };
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
                });
                //Calc the mean number of genotypes per variant
                let mut variant_count = variant_count.lock().unwrap();
                let mut total_genotype_count = total_genotype_count.lock().unwrap();
                if *variant_count > 0 {
                    *mean_genotypes = *total_genotype_count as f32 / *variant_count as f32;
                } else {
                    *mean_genotypes = 0.0 as f32;
                }
                *genotypes_per_position = genotypes.lock().unwrap().clone();
            }
        }
    }

    fn calc_coverage(&mut self) -> f32 {
        match self {
            PileupStats::PileupContigStats {
                ref mut depth,
                target_len,
                ref mut coverage,
                ref mut variance,
                observed_contig_length,
                num_covered_bases,
                contig_end_exclusion,
                min_fraction_covered_bases,
                min,
                max,
                ..

            } => {
                let len1 = target_len;
                match *contig_end_exclusion * 2 < *len1 as u32 {
                    true => {
                        debug!("Adding len1 {}", len1);
                        *observed_contig_length += *len1 as u32 - 2 * *contig_end_exclusion
                    },
                    false => {
                        debug!("Contig too short - less than twice the contig-end-exclusion");
                    }
                }

                debug!("Total observed length now {}", *observed_contig_length);
                let mut counts: Vec<usize> = vec!();
                let start_from = *contig_end_exclusion as usize;
                let end_at = *len1 - *contig_end_exclusion as usize - 1;
                let mut cumulative_sum;
                for (i, current) in depth.iter().enumerate() {
                    cumulative_sum = *current;
                    if i >= start_from && i <= end_at {
                        if cumulative_sum > 0 {
                            *num_covered_bases += 1
                        }
                        if counts.len() <= cumulative_sum {
                            (counts).resize(cumulative_sum + 1, 0);
                        }
                        (counts)[cumulative_sum] += 1
                    }
                }

                let total_bases = *observed_contig_length;
                debug!("Calculating coverage with num_covered_bases {}, observed_length {} and counts {:?}",
                       num_covered_bases, observed_contig_length, counts);
                let answer = match total_bases {
                    0 => {
                        *variance = 0.0;
                        0.0},
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
                            *variance = 0.0;
                            0.0
                        } else {
                            let min_index: usize = (*min * total_bases as f32).floor() as usize;
                            let max_index: usize = (*max * total_bases as f32).ceil() as usize;
                            if *num_covered_bases == 0 { return 0.0; }
//                            counts[0] += 0;

                            let mut num_accounted_for: usize = 0;
                            let mut total: usize = 0;
                            let mut started = false;
                            let mut i = 0;

                            let mut k = 0;
                            // Ensure K is within the range of coverages - take the
                            // lowest coverage.
                            while counts[k] == 0 {
                                k += 1;
                            }
                            let mut ex = 0;
                            let mut ex2 = 0;

                            for (x, num_covered) in counts.iter().enumerate() {
                                num_accounted_for += num_covered.clone() as usize;

                                if num_covered > &0 {
                                    let nc = *num_covered as usize;
                                    ex += (x-k) * nc;
                                    ex2 += (x-k)*(x-k) * nc;
                                }
                                debug!("start: i {}, num_accounted_for {}, total {}, min {}, max {}",
                                       i, num_accounted_for, total, min_index, max_index);
                                if num_accounted_for >= min_index {
                                    debug!("inside");
                                    if started {
                                        if num_accounted_for > max_index {
                                            debug!("num_accounted_for {}, *num_covered {}",
                                                   num_accounted_for, *num_covered);
                                            let num_excess = num_accounted_for - *num_covered as usize;
                                            let num_wanted = match max_index >= num_excess {
                                                true => max_index - num_excess + 1,
                                                false => 0
                                            };
                                            debug!("num wanted1: {}", num_wanted);
                                            total += num_wanted * i;
                                            break;
                                        } else {
                                            total += *num_covered as usize * i;
                                        }
                                    } else {
                                        if num_accounted_for > max_index {
                                            // all coverages are the same in the trimmed set
                                            total = (max_index - min_index + 1) * i;
                                            started = true
                                        } else if num_accounted_for < min_index {
                                            debug!("too few on first")
                                        } else {
                                            let num_wanted = num_accounted_for - min_index + 1;
                                            debug!("num wanted2: {}", num_wanted);
                                            total = num_wanted * i;
                                            started = true;
                                        }
                                    }
                                }
                                debug!("end i {}, num_accounted_for {}, total {}", i, num_accounted_for, total);

                                i += 1;
                            }
                            // Return sample variance not population variance since
                            // almost all MAGs are incomplete.
                            *variance =
                                (ex2 as f32 - (ex*ex) as f32/total_bases as f32) / (total_bases - 1) as f32;
                            total as f32 / (max_index - min_index) as f32
                        }
                    }
                };
                *coverage = answer.clone();
                return answer
            }
        }
    }

    fn print_variants(&mut self, ref_sequence: Vec<u8>, _depth_thresh: usize){
        match self {
            PileupStats::PileupContigStats {
                nucfrequency,
                indels,
                variant_abundances,
                variants_in_reads,
                depth,
                tid,
                genotypes_per_position,
                ..

            } => {
                println!("tid\tpos\tvariant\treference\tabundance\tdepth\tgenotypes");
                for (position, hash) in variant_abundances.iter() {
                    // loop through each position that has variants
                    let position = *position as usize;
                    let d = depth[position];

                    for (var, abundance) in hash.iter() {
                        // for each variant at a location
                        if indels[position].contains_key(var) {
                            // How does this print N for insertions?
                            print!("{}\t{}\t{}\t{}\t{}\t{}\t", tid, position,
                                   var,
                                   str::from_utf8(
                                       &ref_sequence[position..position
                                           + var.len() as usize]).unwrap(),
                                   abundance, d);

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_hash) => {
                                    match gtype_hash.get(&var.to_string()) {
                                        Some(gtype_count) => {
                                            print!("{}\n", gtype_count);
                                        },
                                        None => {
                                            print!("0\n");
                                        }
                                    }
                                },
                                None => {
                                    print!("0\n");
                                },
                            };


//                            let read_ids =
//                                match indels[position].get(var) {
//                                    Some(ids) => ids,
//                                    None => {
//                                        println!("Variant not in indel hash");
//                                        std::process::exit(1)
//                                    },
//                                };
//                            let mut connected_bases = HashSet::new();
//                            for read_id in read_ids {
//                                if variants_in_reads.contains_key(&read_id) {
//                                    connected_bases = connected_bases.union(
//                                        &variants_in_reads[&read_id].keys().cloned().collect::<HashSet<i32>>())
//                                        .cloned().collect::<HashSet<i32>>();
//                                }
//                            }
//                            print!("{:?}\n", connected_bases)
                        } else if var.len() == 1{
                            print!("{}\t{}\t{}\t{}\t{}\t{}\t", tid, position,
                                   var,
                                   ref_sequence[position] as char,
                                   abundance, d);

                            // Print number of genotypes associated with that position and variant
                            match genotypes_per_position.get(&position) {
                                Some(gtype_hash) => {
                                    match gtype_hash.get(&var.to_string()) {
                                        Some(gtype_count) => {
                                            print!("{}\n", gtype_count);
                                        },
                                        None => {
                                            print!("0\n");
                                        }
                                    }
                                },
                                None => {
                                    print!("0\n");
                                },
                            };

//                            let read_ids =
//                                match nucfrequency[position]
//                                    .get(&(var.clone().into_bytes()[0] as char)){
//                                Some(ids) => ids,
//                                None => {
//                                    println!("Variant not in frequency Hash");
//                                    std::process::exit(1)},
//                            };

//                            let mut connected_bases = HashSet::new();
//                            for read_id in read_ids{
//                                if variants_in_reads.contains_key(&read_id) {
//                                    connected_bases = connected_bases.union(
//                                        &variants_in_reads[&read_id].keys().cloned().collect::<HashSet<i32>>())
//                                        .cloned().collect::<HashSet<i32>>();
//                                }
//                            }
//                            print!("{:?}\n", connected_bases)
                        }
                    }
                };
            }
        }
    }
}