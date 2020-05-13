use std::collections::{HashMap, HashSet, BTreeMap};
use std::str;
use std::sync::{Arc, Mutex};
use std::io::prelude::*;
use rayon::prelude::*;
use estimation::codon_structs::*;
use bio_types::strand;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use rust_htslib::bcf::record;

use std::path::Path;
use std::fs::OpenOptions;

use model::variants::*;

pub enum VariantStats {
    VariantContigStats {
        variants: HashMap<i64, HashMap<Variant, Base>>,
//        nucfrequency: HashMap<i32, BTreeMap<char, BTreeSet<i64>>>,
//        variants_in_reads: HashMap<i64, BTreeMap<i32, String>>,
//        variant_abundances: HashMap<i32, BTreeMap<String, (f64, f64)>>,
        variant_count: Vec<f64>,
        depth: Vec<i32>,
//        indels: HashMap<i32, BTreeMap<String, BTreeSet<i64>>>,
//        genotypes_per_position: HashMap<usize, usize>,
//        mean_genotypes: f64,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: f64,
//        variations_per_n: usize,
//        total_variants: usize,
        coverage: f64,
        variance: f64,
        observed_contig_length: u32,
        num_covered_bases: i32,
        num_mapped_reads: u64,
        total_mismatches: u64,
        contig_end_exclusion: u64,
        min: f64,
        max: f64,
        method: String,
        regression: (f64, f64, f64),
        // Clusters hashmap:
        // Key = Position
        // Value = K: Variant, V: (DBSCAN Cluster, HAC Index/initial cluster)
//        clusters: HashMap<i32, BTreeMap<String, (i32, usize)>>,
//        clusters_mean: HashMap<i32, f64>,
    }
}

impl VariantStats {
    pub fn new_contig_stats(min: f64, max: f64,
                            contig_end_exclusion: u64) -> VariantStats {
        VariantStats::VariantContigStats {
            variants: HashMap::new(),
            tid: 0,
            variant_count: vec!(),
            depth: vec!(),
            total_indels: 0,
            target_name: vec!(),
            target_len: 0.0,
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
            regression: (0., 0., 0.),
        }
    }
}

pub trait VariantFunctions {
    fn setup(&mut self);

    fn len(&mut self) -> usize;

    fn add_contig(&mut self,
                  variant_map: Option<&mut HashMap<i64, HashMap<Variant, Base>>>,
                  tid: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize,
                  sample_idx: usize,
                  coverages: Vec<f64>,
                  ups_and_downs: Vec<i32>);

    /// Perform linear regression between total mismacthes and read depth
    fn calc_error(&mut self, ani: f32) -> usize;

//    /// Filter out variants from potential sequencing or mapping errors
//    fn filter_variants(&mut self,
//                       min_variant_depth: usize,
//                       coverage_fold: f64);

    /// Replace reference variants with most dominant variant observed within the reads
    fn polish_contig(&mut self,
                     original_contig: &Vec<u8>,
                     output_prefix: &str);

    /// Perform dN/dS calculations based on read mapping using modified Jukes-Cantor method
    fn calc_gene_mutations(&mut self,
                           gff_map: &HashMap<String, Vec<bio::io::gff::Record>>,
                           ref_sequence: &Vec<u8>,
                           codon_table: &CodonTable);

//    /// Prints out variant info for current contig
//    fn print_variants(&mut self, ref_sequence: &Vec<u8>, stoit_name: &str);
}

impl VariantFunctions for VariantStats {
    fn setup(&mut self) {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                ref mut variant_count,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut coverage,
                ref mut num_covered_bases,
                ref mut num_mapped_reads, ..
            } => {
                *variants = HashMap::new();
                *variant_count = Vec::new();
                *depth = vec!();
                *tid = 0;
                *total_indels = 0;
                *target_name = vec!();
                *target_len = 0.0;
                *coverage = 0.00;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
        }
    }

    fn len(&mut self) -> usize {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                ..
            } => {
                variants.len()
            }
        }
    }

    fn add_contig(&mut self,
                  variant_map: Option<&mut HashMap<i64, HashMap<Variant, Base>>>,
                  target_id: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize,
                  sample_idx: usize,
                  coverages: Vec<f64>,
                  ups_and_downs: Vec<i32>) {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                ref mut variant_count,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut variance,
                ref mut coverage,
                ref mut method,
                ..
            } => {

                *tid = target_id;
                *total_indels = total_indels_in_contig;
                *target_name = contig_name;
                *target_len = contig_len as f64;
                *coverage = coverages[1] as f64;
                *variance = coverages[2] as f64;
                *method = method.to_string();
                *variants = match variant_map {
                    Some(map) => map.clone(),
                    _ => HashMap::new(),
                };

                let variant_count_safe = Arc::new(Mutex::new(vec![0.; ups_and_downs.len()]));
                *depth = vec![0; ups_and_downs.len()];
                let mut cumulative_sum = 0;
                for (pos, current) in ups_and_downs.iter().enumerate() {
                    cumulative_sum += *current;
                    depth[pos] = cumulative_sum;
                }

                // Calculate how many reads have variant at each position
                // to go into linear regression predicting read error rate
                depth.par_iter_mut().enumerate().for_each(|(pos, d)| {
                    let mut variant_count_safe = variant_count_safe.lock().unwrap();
                    let mut var_set = match variants.get(&(pos as i64)) {
                        Some(set) => set.to_owned(),
                        None => HashMap::new(),
                    };
                    // Add depth of variant to count file if variant is present
                    for (var, base) in var_set.iter_mut() {
                        let var_depth = match base.variant {
                            Variant::None => 0,
                            _ => base.depth.par_iter().sum(),
                        };
                        if var_depth != 0 {
                            variant_count_safe[pos] += var_depth as f64;
                        }
                    }
                });

                *variant_count = variant_count_safe.lock().unwrap().to_vec();

                debug!("new contig added {} with coverage {} and variance {}", tid, coverage, variance);
            }
        }
    }

    fn calc_error(&mut self, ani: f32) -> usize {
        match self {
            VariantStats::VariantContigStats {
                variant_count,
                depth,
                regression,
                target_len,
                ..
            } => {
                // convert depth to f64 for this to work
                let depth_64: Vec<f64> = depth.par_iter().map(|x| *x as f64).collect();
                let data = vec![("Y", variant_count.clone()), ("X", depth_64.clone())];

                let data = RegressionDataBuilder::new()
                    .build_from(data).expect("Unable to build regression from data");
                let formula = "Y ~ X";
                let model = FormulaRegressionBuilder::new()
                    .data(&data)
                    .formula(formula)
                    .fit()
                    .expect("Unable to fit data to formula");
                let parameters = model.parameters;
                let standard_errors = model.se.pairs();
                let pvalues = model.pvalues;
                debug!("Linear regression results: \n params {:?} \n se {:?} \n p-values {:?}",
                         parameters,
                         standard_errors,
                         pvalues.pairs());
                // return results as intercept, effect size, se
                *regression = (parameters.intercept_value as f64,
                               parameters.regressor_values[0] as f64,
                               standard_errors[0].1 as f64);

                if ani > 0. {
                    // calculate the minimum number of variant sites needed for specified ANI
                    // ANI is not divided by 100 here as it is already betwen 0 and 1
                    let variants_for_ani = (*target_len as f32 * (1. - ani)) as usize;
                    info!("Variants needed for ANI {}", variants_for_ani);

                    // Sort variant count vec descending
                    variant_count.sort_by(|a, b| b.partial_cmp(a).unwrap());

                    variant_count[variants_for_ani] as usize
                } else {
                    0
                }

            }
        }
    }

//    fn calc_variants(&mut self, min_variant_depth: usize, coverage_fold: f64){
//        match self {
//            VariantStats::VariantContigStats {
//                ref mut variants,
//                depth,
//                ref mut coverage,
//                tid,
//                regression,
//                ..
//            } => {
//                let variants = Arc::new(Mutex::new(HashMap::new())); // The relative abundance of each variant
//                let read_variants = Arc::new(Mutex::new(HashMap::new())); // The reads with variants and their positions
//                let variant_count = Arc::new(Mutex::new(0));
//                let indels = Arc::new(Mutex::new(indels));
////                let mut outside_coverage = Arc::new(Mutex::new(HashMap::new()));
//                let nucfrequency = Arc::new(Mutex::new(nucfrequency));
////                let min_variant_fraction = min_variant_depth as f64 / 100.;
//                // for each location calculate if there is a variant based on read depth
//                // Uses rayon multithreading
//                depth.par_iter_mut().enumerate().for_each(|(i, d)| {
////                    let read_variants = Arc::clone(&read_variants);
////                    let variant_count = Arc::clone(&variant_count);
//                    let mut rel_abundance = BTreeMap::new();
////                        if d >= &mut min_variant_depth.clone() {
//                    // INDELS act differently to normal variants
//                    // The reads containing this variant don't contribute to coverage values
//                    // So we need to readjust the depth at these locations to show true
//                    // read depth. i.e. variant depth + reference depth
//                    let mut indels
//                        = indels.lock().unwrap();
//                    let indel_map = match indels.get(&(i as i32)) {
//                        Some(map) => map.to_owned(),
//                        None => BTreeMap::new(),
//                    };
//                    if indel_map.len() > 0 {
//                        for (indel, read_ids) in indel_map.iter() {
//                            let count = read_ids.len();
////                                if indel.contains("N") {
////                                    *d += count as f64;
////                                }
//                            *d += count as f64;
//                            if (*coverage * (1.0 - coverage_fold) <= *d as f64
//                                && *d as f64 <= *coverage * (1.0 + coverage_fold))
//                                || (coverage_fold == 0.0) {
//                                if (count >= min_variant_depth)
////                                    && (count as f64 / *d >= min_variant_fraction)
//                                    && ((count as f64) > *d as f64 * 6. * (regression.1 + regression.2)) {
//                                    rel_abundance.insert(indel.to_owned(), (count as f64, *d as f64));
//                                    for read in read_ids {
//                                        let mut read_variants
//                                            = read_variants.lock().unwrap();
//
//                                        let read_vec = read_variants
//                                            .entry(read.clone())
//                                            .or_insert(BTreeMap::new());
//                                        read_vec.insert(i as i32, indel.clone());
//                                    }
//                                    let mut variant_count = variant_count.lock().unwrap();
//                                    *variant_count += 1;
//                                } else {
//                                    let indel_map_back = indels
//                                        .entry(i as i32).or_insert(BTreeMap::new());
//                                    indel_map_back.remove(indel);
//                                }
//                            }
//                        }
//                    }
//                    if (*coverage * (1.0 - coverage_fold) <= *d as f64
//                        && *d as f64 <= *coverage * (1.0 + coverage_fold))
//                        || (coverage_fold == 0.0) {
//                        let mut nucfrequency
//                            = nucfrequency.lock().unwrap();
//
//                        let nuc_map = match nucfrequency.get(&(i as i32)) {
//                            Some(map) => map.to_owned(),
//                            None => BTreeMap::new(),
//                        };
//                        if nuc_map.len() > 1 {
//                            for (base, read_ids) in nuc_map.iter() {
//                                let count = read_ids.len();
//
//                                if base != &"R".chars().collect::<Vec<char>>()[0] {
//                                    if (count >= min_variant_depth)
//    //                                    && (count as f64 / *d >= min_variant_fraction)
//                                        && ((count as f64) > *d as f64 * 6. * (regression.1 + regression.2)) {
//                                        rel_abundance.insert(base.to_string(), (count as f64, *d as f64));
//
//                                        for read in read_ids {
//                                            let mut read_variants
//                                                = read_variants.lock().unwrap();
//                                            let read_vec = read_variants
//                                                .entry(read.clone())
//                                                .or_insert(BTreeMap::new());
//                                            read_vec.insert(i as i32, base.to_string());
//                                        }
//                                        if base != &"R".chars().collect::<Vec<char>>()[0] {
//                                            let mut variant_count = variant_count.lock().unwrap();
//                                            *variant_count += 1;
//                                        }
//                                    } else {
//                                        let nuc_map_back = nucfrequency
//                                            .entry(i as i32).or_insert(BTreeMap::new());
//                                        nuc_map_back.remove(base);
//                                    }
//                                }
//                            };
//                        } else {
//                            let nuc_map_back = nucfrequency
//                                .entry(i as i32).or_insert(BTreeMap::new());
//                            nuc_map_back.remove(&"R".chars().collect::<Vec<char>>()[0]);
//                        };
//
////                        }
//                    } else {
//
//                    }
//
//                    if rel_abundance.len() > 1 {
//                        let mut variants = variants.lock().unwrap();
////                        debug!("Relative Abundances: {:?}", rel_abundance);
//                        variants.insert(i as i32, rel_abundance);
//                    } else if rel_abundance.len() == 1 {
//                        if !(rel_abundance.contains_key(&"R".to_string())) {
//                            let mut variants = variants.lock().unwrap();
////                            debug!("Relative Abundances: {:?}", rel_abundance);
//                            variants.insert(i as i32, rel_abundance);
//                        } else {
//                            // Need to remove R from whitelist
//                            let mut nucfrequency
//                                = nucfrequency.lock().unwrap();
//                            let nuc_map = nucfrequency
//                                .entry(i as i32).or_insert(BTreeMap::new());
//                            if nuc_map.len() > 1 {
//                                for (base, read_ids) in nuc_map.iter() {
//                                    let count = read_ids.len();
//                                    if base == &("R".as_bytes()[0] as char) {
//                                        for read in read_ids {
//                                            let mut read_variants
//                                                = read_variants.lock().unwrap();
//                                            let read_vec = read_variants
//                                                .entry(read.clone())
//                                                .or_insert(BTreeMap::new());
//                                            read_vec.remove(&(i as i32));
//                                        }
//                                    }
//                                }
//                            } else {
//                                nucfrequency.remove(&(i as i32));
//                            }
//                        }
//                    }
//
//                });
//
//                let read_variants = read_variants.lock().unwrap();
//                *variants_in_reads = read_variants.to_owned();
//                let variants = variants.lock().unwrap();
//                *variant_abundances = variants.to_owned();
//                let variant_count = variant_count.lock().unwrap();
//                debug!("Total variants for {}: {:?}", tid, variant_count);
//                // Total variants is the actual amount of variants that passed all thresholds
//                *total_variants = *variant_count;
//                let mut nucfrequency = nucfrequency.lock().unwrap();
//                **nucfrequency = nucfrequency.to_owned();
//                let mut indels = indels.lock().unwrap();
//                **indels = indels.to_owned();
//            }
//        }
//    }


    fn polish_contig(&mut self,
                     original_contig: &Vec<u8>,
                     output_prefix: &str) {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                target_name,
                ..
            } => {
                let file_name = output_prefix.to_string()
                    + &"_polished.fna".to_owned();

                let file_path = Path::new(&file_name);

                // Open haplotype file or create one

                let mut file_open = OpenOptions::new().append(true).create(true)
                    .open(file_path).expect("No Read or Write Permission in current directory");
                file_open.write(b">").unwrap();
                file_open.write(target_name).unwrap();
                file_open.write(b"\n").unwrap();


                let mut contig = String::new();

                let mut variations = 0;
                let mut skip_n = 0;
                let mut skip_cnt = 0;
                // Generate the consensus genome by checking each variant
                // Variant has to be in more than 0.5 of population
                for (pos, base) in original_contig.iter().enumerate() {
                    if skip_cnt < skip_n {
                        skip_cnt += 1;
                    } else {
                        let mut max_var = Variant::None;
                        let mut max_abund = 0.0;
                        skip_n = 0;
                        skip_cnt = 0;
                        if variants.contains_key(&(pos as i64)){
                            let alleles = &variants[&(pos as i64)];
                            for (var, base) in alleles.iter() {
                                if base.freq[0] > max_abund {
                                    max_var = var.clone();
                                    max_abund = base.freq[0];
                                }
                            }
                            if max_abund >= 0.5 {
                                match max_var {
                                    Variant::Deletion(size) => {
                                        // Skip the next n bases but rescue the reference prefix
                                        skip_n = size - 1;
                                        skip_cnt = 0;
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
                            contig = contig + str::from_utf8(&[*base]).unwrap();
                        }
                    }
                };
//                contig = contig + "\n";
                for line in contig.as_bytes().to_vec()[..].chunks(60).into_iter(){
                    file_open.write(line).unwrap();
                    file_open.write(b"\n").unwrap();
                };
//                writeln!(file_open, "{:?}", contig);
            }
        }
    }

    fn calc_gene_mutations(&mut self,
                           gff_map: &HashMap<String, Vec<bio::io::gff::Record>>,
                           ref_sequence: &Vec<u8>,
                           codon_table: &CodonTable) {
        match self {
            VariantStats::VariantContigStats {
                variants,
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
                let print_stream = Arc::new(Mutex::new(std::io::stdout()));
                gff_records.par_iter().enumerate().for_each(|(_id, gene)| {
                    let dnds = codon_table.find_mutations(gene, variants, ref_sequence, depth);
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

//                    for cursor in start..end+1 {
//                        let variant_map = match variants.get(&(cursor as i64)){
//                            Some(map) => map,
//                            None => continue,
//                        };
//                        if variant_map.len() > 0 {
//                            let mut print_stream = print_stream.lock().unwrap();
//
//
//                            for variant in variant_map {
//                                if variant.to_owned().contains("R"){
//                                    continue
//                                }
//                                write!(print_stream, "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t",
//                                         contig.clone()+gene_id, gene.start(),
//                                         gene.end(), frame, strand_symbol, dnds, cursor).expect("Unable to write to stream");
//
//
//                                if variant.to_owned().contains("N") {
//                                    writeln!(print_stream, "{}\t{}\t{}\t{}\t{}",
//                                           variant,
//                                           str::from_utf8(
//                                               &ref_sequence[cursor..cursor
//                                                   + variant.len() as usize]).unwrap(),
//                                           abundance.0, abundance.1, "D").expect("Unable to write to stream");
//
//                                } else if indel_map.contains_key(variant) {
//                                     writeln!(print_stream,"{}\t{}\t{:.3}\t{}\t{}",
//                                           variant,
//                                           str::from_utf8(
//                                               &[ref_sequence[cursor]]).unwrap(),
//                                           abundance.0, abundance.1, "I").expect("Unable to write to stream");
//                                } else {
//                                    writeln!(print_stream, "{}\t{}\t{}\t{}\t{}",
//                                             variant,
//                                             ref_sequence[cursor] as char,
//                                             abundance.0, abundance.1, "S").expect("Unable to write to stream");
//                                }
//                            }
//                        }
//                    }
                });

            }
        }

    }

//    fn print_variants(&mut self, ref_sequence: &Vec<u8>, stoit_name: &str){
//        match self {
//            VariantStats::VariantContigStats {
//                variants,
//                depth,
//                tid,
//                ..
//
//            } => {
//                info!("Outputting {} variant locations", variant_abundances.keys().len());
//                for (position, _d) in depth.iter().enumerate() {
//
//                    // loop through each position that has variants
//                    let hash = match variant_abundances.get(&(position as i32)) {
//                        Some(hash) => hash,
//                        None => continue,
//                    };
//
//                    let cluster_map = match clusters.get(&(position as i32)) {
//                        Some(map) => map.clone(),
//                        None => BTreeMap::new(),
//                    };
//
//                    for (var, abundance) in hash.iter() {
//                        // for each variant at a location
//                        let indel_map = match indels.get(&(position as i32)) {
//                            Some(map) => map.to_owned(),
//                            None => BTreeMap::new(),
//                        };
//
//                        let cluster_val = match cluster_map.get(var) {
//                            Some(i) => *i,
//                            None => (-1, 0),
//                        };
//
//                        if indel_map.contains_key(var) {
//                            // How does this print N for insertions?
//                            if var.to_owned().contains("N"){
//                                print!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
//                                       stoit_name,
//                                       tid,
//                                       position,
//                                       var,
//                                       str::from_utf8(
//                                           &ref_sequence[position..position
//                                               + var.len() as usize]).unwrap(),
//                                       abundance.0, abundance.1);
//
//                            } else {
//                                print!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
//                                       stoit_name,
//                                       tid, position,
//                                       var,
//                                       str::from_utf8(
//                                           &[ref_sequence[position]]).unwrap(),
//                                       abundance.0, abundance.1);
//                            }
//
//                            // Print number of genotypes associated with that position and variant
//                            match genotypes_per_position.get(&position) {
//                                Some(gtype_count) => {
//                                    print!("{}\t", gtype_count);
//                                },
//                                None => {
//                                    print!("0\t");
//                                },
//                            };
//                            println!("{}\t{}", cluster_val.0, cluster_val.1);
//
//                        } else if var.len() == 1 && var != &"R".to_string(){
//                            print!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
//                                   stoit_name,
//                                   tid, position,
//                                   var,
//                                   ref_sequence[position] as char,
//                                   abundance.0, abundance.1);
//
//                            // Print number of genotypes associated with that position and variant
//                            match genotypes_per_position.get(&position) {
//                                Some(gtype_count) => {
//                                    print!("{}\t", gtype_count);
//                                },
//                                None => {
//                                    print!("0\t");
//                                },
//                            };
//                            println!("{}\t{}", cluster_val.0, cluster_val.1);
//                        }
//                    }
//                };
//            }
//        }
//    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileups::*;
    use std::str;
    use std::fs::File;
    use coverm::mapping_parameters::*;
    use coverm::shard_bam_reader::*;
    use coverm::genome_exclusion::*;
    use coverm::bam_generator::*;

//    fn test_with_stream<R: NamedBamReader + Send,
//        G: NamedBamReaderGenerator<R> + Send>(
//        expected: &str,
//        bam_readers: Vec<G>,
//        mut reference: bio::io::fasta::IndexedReader<File>,
//        proper_pairs_only: bool,
//        n_threads: usize,
//        coverage_fold: f32,
//        min_var_depth: usize,
//        min: f32,
//        max: f32,
//        mode: &str,
//        include_indels: bool,
//        include_soft_clipping: bool) {
////        let mut stream = Cursor::new(Vec::new());
//        {
//            reads_mapped_vec = pileup_variants(
//                bam_readers,
//                &mut coverage_taker,
//                coverage_estimators,
//                print_zero_coverage_contigs,
//                flag_filters,
//                false,
//                );
//        }
////        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
//    }
    #[test]
    fn test_genotypes_two_positions() {
        let mut contig = VariantStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i64>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());

        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [0, 1, 2].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [3, 4, 5].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                   [0, 1].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [2, 3].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [4, 5].iter().cloned().collect::<BTreeSet<i64>>());

        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          2,
                          "mean",
                          vec![2., 6., 0.],
                          vec![6,0]);

        match contig {
            VariantStats::VariantContigStats {
                ref mut regression,
                ..
            } => {
                *regression = (0., 0., 0.);
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);
    }

    #[test]
    fn test_genotypes_three_positions() {
        let mut contig = VariantStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i64>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());
        nuc_freq.insert(2, BTreeMap::new());

        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [0, 1, 2].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [3, 4, 5].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                   [0, 1].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                   [2, 3].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                   [4, 5].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(2).or_insert(BTreeMap::new());

        pos.insert("G".chars().collect::<Vec<char>>()[0],
                   [1, 2, 4, 5].iter().cloned().collect::<BTreeSet<i64>>());


        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          3,
                          "mean",
                          vec![2., 6., 0.],
                          vec![6,4,0]);

        match contig {
            VariantStats::VariantContigStats {
                ref mut regression,
                ..
            } => {
                *regression = (0., 0., 0.);
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);


    }

    #[test]
    fn test_genotypes_four_positions() {
        let mut contig = VariantStats::new_contig_stats(0.,
                                                       1.,
                                                       0);
        let mut nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i64>>> = HashMap::new();
        nuc_freq.insert(0, BTreeMap::new());
        nuc_freq.insert(1, BTreeMap::new());
        nuc_freq.insert(2, BTreeMap::new());
        nuc_freq.insert(3, BTreeMap::new());
        let reference = "TGTG".as_bytes();
        let mut pos = nuc_freq.entry(0).or_insert(BTreeMap::new());

        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [0].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                           [5].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(1).or_insert(BTreeMap::new());

        pos.insert("T".chars().collect::<Vec<char>>()[0],
                           [0, 2].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [1].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("R".chars().collect::<Vec<char>>()[0],
                           [3, 5].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("C".chars().collect::<Vec<char>>()[0],
                           [4].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(2).or_insert(BTreeMap::new());

        pos.insert("G".chars().collect::<Vec<char>>()[0],
                           [0, 1, 2, 3].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("A".chars().collect::<Vec<char>>()[0],
                           [4, 5].iter().cloned().collect::<BTreeSet<i64>>());

        let mut pos = nuc_freq.entry(3).or_insert(BTreeMap::new());

        pos.insert("C".chars().collect::<Vec<char>>()[0],
                           [2].iter().cloned().collect::<BTreeSet<i64>>());
        pos.insert("T".chars().collect::<Vec<char>>()[0],
                           [3, 4, 5].iter().cloned().collect::<BTreeSet<i64>>());


        contig.add_contig(nuc_freq,
                          HashMap::new(),
                          0,
                          0,
                          "contig_name".as_bytes().to_vec(),
                          4,
                          "mean",
                          vec![2., 6., 0.],
                          vec![3,3,0,-6]);

        match contig {
            VariantStats::VariantContigStats {
                ref mut regression,
                ..
            } => {
                *regression = (0., 0., 0.);
            }
        }

        // filters variants across contig
        contig.calc_variants(
            1,
            0.);

        let snps = match &contig {
            VariantStats::VariantContigStats {
                variant_abundances,
                ..
            } => {
                variant_abundances
            }
        };
        assert_eq!(snps.len(), 4);

        assert_eq!(contig.print_variants(&reference.to_vec(), "test"),
        println!("test    0       0       A       T       1       3       0       -1      0
                    test    0       1       A       G       1       6       0       -1      0
                    test    0       1       C       G       1       6       0       -1      0
                    test    0       1       T       G       2       6       0       -1      0
                    test    0       2       A       T       2       6       0       -1      0
                    test    0       2       G       T       4       6       0       -1      0
                    test    0       3       C       G       1       0       0       -1      0
                    test    0       3       T       G       3       0       0       -1      0"));
    }
}