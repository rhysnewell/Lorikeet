use std::collections::{HashMap, HashSet};
use std::collections::{BTreeMap, BTreeSet};
use std::cmp::min;
use std::str;
use rm::linalg::Matrix;
use simhash::*;
use mash::*;
use distance::*;
use rust_htslib::bam::record::{Cigar, CigarStringView};

#[derive(Debug)]
pub struct VariantBase {
    read_ids: HashSet<i32>,
    connected_bases: HashSet<i32>,
    indel: bool,
    seen: bool,
    pos: usize,
    seq: String,
    abundance: f32,
}

pub enum PileupStats {
    PileupContigStats {
        nucfrequency: Vec<HashMap<char, HashSet<i32>>>,
        kfrequency: BTreeMap<Vec<u8>, usize>,
        variants_in_reads: HashMap<i32, HashSet<i32>>,
        variant_abundances: BTreeMap<i32, HashMap<String, f32>>,
        depth: Vec<usize>,
        indels: Vec<HashMap<String, HashSet<i32>>>,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: usize,
        variations_per_base: f32,
        coverage: f32,
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
            kfrequency: BTreeMap::new(),
            variants_in_reads: HashMap::new(),
            variant_abundances: BTreeMap::new(),
            depth: vec!(),
            indels: vec!(),
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0,
            variations_per_base: 0.00,
            coverage: 0.00,
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
                  k_freq: BTreeMap<Vec<u8>, usize>,
                  read_depth: Vec<usize>,
                  indels_positions: Vec<HashMap<String, HashSet<i32>>>,
                  tid: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize);

    fn calc_variants(&mut self,
                     depth_thresh: usize,
                     variant_fraction: f64);

    fn generate_variant_contig(&mut self,
                               original_contig: Vec<u8>,
                               depth_thresh: usize,
                               variant_fraction: f64);

    fn calc_coverage(&mut self) -> f32;
}

impl PileupFunctions for PileupStats {
    fn setup(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut kfrequency,
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
                *kfrequency = BTreeMap::new();
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
                  k_freq: BTreeMap<Vec<u8>, usize>,
                  read_depth: Vec<usize>,
                  indel_positions: Vec<HashMap<String, HashSet<i32>>>,
                  target_id: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize) {
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut kfrequency,
                ref mut depth,
                ref mut indels,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ..
            } => {
                *nucfrequency = nuc_freq;
                *kfrequency = k_freq;
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
                kfrequency,
                ref mut variants_in_reads,
                ref mut variant_abundances,
                ref mut depth,
                ref mut indels,
                tid,
                total_indels,
                target_name,
                target_len,
                ref mut variations_per_base,
                ref mut coverage,
                ..
            } => {
                let mut variants = BTreeMap::new(); // The relative abundance of each variant
                let mut read_variants = HashMap::new(); // The reads with variants and their positions
                let mut variant_count = 0;
                let mut cursor = 0;
                let mut depth_sum = 0;


                for (i, d) in depth.iter().enumerate(){
                    let mut rel_abundance = HashMap::new();
                    if d >= &depth_thresh {
                        if nucfrequency[i].len() > 0 {
                            for (base, read_ids) in nucfrequency[i].iter() {
                                variant_count += 1;
                                let count = read_ids.len();
                                if count as f64 / *d as f64 >= variant_fraction {
                                    rel_abundance.insert(base.to_string(), count as f32 / *d as f32);
                                    for read in read_ids {
                                        let mut read_vec = read_variants
                                            .entry(read.clone())
                                            .or_insert(HashSet::new());
                                        read_vec.insert(cursor as i32);
                                    }
                                }
                            }
                        };
                        if indels[i].len() > 0 {
                            for (indel, read_ids) in indels[i].iter(){
                                variant_count += 1;
                                let count = read_ids.len();
                                if count as f64 / *d as f64 >= variant_fraction {
                                    rel_abundance.insert(indel.clone(), count as f32 / *d as f32);
                                    for read in read_ids {
                                        let mut read_vec = read_variants
                                            .entry(read.clone())
                                            .or_insert(HashSet::new());
                                        read_vec.insert(cursor as i32);
                                    }
                                }
                            }
                        }
                    }

                    if rel_abundance.len() > 0 {
                        variants.insert(cursor, rel_abundance);
                    }

                    cursor += 1;
                    depth_sum += d;
                }

                debug!("read variants {:?}", read_variants);
                *variants_in_reads = read_variants;
                debug!("variants abundances {:?}", variants);
                *variant_abundances = variants;
                *variations_per_base = (variant_count+*total_indels as i32) as f32/target_len.clone() as f32;
            }
        }
    }

    fn generate_variant_contig(&mut self,
                               original_contig: Vec<u8>,
                               depth_thresh: usize,
                               variant_fraction: f64){
        match self {
            PileupStats::PileupContigStats {
                ref mut nucfrequency,
                ref mut variants_in_reads,
                ref mut depth,
                ref mut indels,
                ref mut variations_per_base,
                ref mut coverage,
                ref mut variant_abundances,
                ..
            } => {
                let mut variant_cooccurrences = HashMap::new();
                let mut variant_record;
                for (i, d) in depth.iter().enumerate(){
                    if d >= &depth_thresh {
                        if nucfrequency[i].len() > 0 {
                            for (k , v) in nucfrequency[i].iter() {
                                let count = v.len();
                                if count as f64 / *d as f64 >= variant_fraction {

                                    let mut occ = variant_cooccurrences.entry(i)
                                        .or_insert(vec!());

                                    variant_record = VariantBase {
                                        read_ids: HashSet::new(),
                                        connected_bases: HashSet::new(),
                                        indel: false,
                                        seen: false,
                                        pos: i,
                                        seq: k.clone().to_string(),
                                        abundance: variant_abundances[&(i as i32)][&k.to_string()],
                                    };
                                    variant_record.read_ids = variant_record.read_ids
                                                                .union(&v)
                                                                .cloned()
                                                                .collect::<HashSet<i32>>();
                                    for read_id in v {
                                        variant_record.connected_bases = variant_record
                                                                .connected_bases
                                                                .union(&variants_in_reads[&read_id])
                                                                .cloned()
                                                                .collect::<HashSet<i32>>();
                                    };
                                    occ.push(variant_record);
                                }
                            }
                        };
                        if indels[i].len() > 0 {
                            for (k , v) in indels[i].iter() {
                                let count = v.len();
                                if count as f64 / *d as f64 >= variant_fraction {
                                    let mut occ = variant_cooccurrences.entry(i)
                                        .or_insert(vec!());

                                    variant_record = VariantBase {
                                        read_ids: HashSet::new(),
                                        connected_bases: HashSet::new(),
                                        indel: true,
                                        seen: false,
                                        pos: i,
                                        seq: k.clone().to_string(),
                                        abundance: variant_abundances[&(i as i32)][k],
                                    };

                                    variant_record.read_ids = variant_record.read_ids.union(&v).cloned().collect::<HashSet<i32>>();
                                    for read_id in v {
                                        variant_record.connected_bases = variant_record
                                            .connected_bases
                                            .union(&variants_in_reads[&read_id]).cloned()
                                            .collect::<HashSet<i32>>();
                                    };
                                    occ.push(variant_record);
                                }
                            }
                        }
                    }
                };
                debug!("{:?}", variant_cooccurrences);

                let mut contig_variants = HashMap::new();
                for (base, variants) in variant_cooccurrences.iter(){
                    let mut variant_track = vec![HashSet::new(); original_contig.len()];
                    let mut base_char;
                    for variant in variants{
                        if variant.indel {
                            if variant.seq.contains("N"){
                                let mut pos = variant.pos + 1;
                                for del in 0..variant.seq.len(){
                                    variant_track[pos].insert('N'.to_string());
                                    pos += 1;
                                }
                            } else {
                                base_char = (original_contig[variant.pos] as char).to_string();
                                base_char.push_str(&variant.seq.clone());
                                variant_track[variant.pos].insert(base_char);
                            }
                        } else{
                            variant_track[variant.pos].insert(variant.seq.clone());
                        };
                        for connected_base in variant.connected_bases.clone(){
                            let mut connected_variants = &variant_cooccurrences[&(connected_base as usize)];
                            for covariant in connected_variants{
                                if !variant.read_ids.is_disjoint(&covariant.read_ids) {
                                    if covariant.indel {
                                        if covariant.seq.contains("N"){
                                            let mut pos = covariant.pos + 1;
                                            for del in 0..covariant.seq.len(){
                                                variant_track[pos].insert('N'.to_string());
                                                pos += 1;
                                            }
                                        } else {
                                            base_char = (original_contig[covariant.pos] as char).to_string();
                                            base_char.push_str(&covariant.seq.clone());
                                            variant_track[covariant.pos].insert(base_char);
                                        }
                                    } else{
                                        variant_track[covariant.pos].insert(covariant.seq.clone());
                                    };
                                }
                            }
                        }
                    }
//                    debug!("Variant Track: {:?}", variant_track);
                    let mut modified_contig = String::new();
                    for (var, reference) in variant_track.iter().zip(original_contig.iter()){
                        if var.len() > 0{
                            for v in var.iter(){
                                if v != &'N'.to_string() {
                                    modified_contig.push_str(v);
                                } else {
                                    break
                                }
                            }
                        } else {
                            modified_contig.push(*reference as char);
                        }
                    }
//                    debug!("Modified Contig: {:?}", modified_contig);
                    contig_variants.insert(modified_contig.clone(), mash_sequence(&modified_contig.into_bytes(),
                            1000,
                            21,
                            0)); // parse to mash
                };
                debug!("Contig Variants: {:?}", contig_variants);
//                let mut hashes =  vec!();
//                hashes.extend(contig_variants.keys().map(|x| x+1));
                for el1 in contig_variants.values() {
                    for el2 in contig_variants.values(){
//                        println!("{}", hash_similarity(*el1, *el2));
                        if el1 != el2 {
                            let dist = distance(&el1,
                                                &el2,
                                                false);
                            println!("{:?}", dist)
                        }

                    }

                }
            }
        }
    }

    fn calc_coverage(&mut self) -> f32 {
        match self {
            PileupStats::PileupContigStats {
                nucfrequency,
                kfrequency,
                ref mut depth,
                tid,
                total_indels,
                target_name,
                target_len,
                variations_per_base,
                ref mut coverage,
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
                let mut cumulative_sum = 0;
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
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
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
                            for num_covered in counts.iter() {
                                num_accounted_for += *num_covered as usize;
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
                            total as f32 / (max_index - min_index) as f32
                        }
                    }
                };
                *coverage = answer.clone();
                return answer
            }
        }
    }
}

pub enum PileupMatrix {
    PileupContigMatrix {
        variants: Vec<f32>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        coverages: Vec<f32>,
        tids: Vec<i32>,
        total_indels_across_contigs: Vec<usize>,
        target_names: Vec<Vec<u8>>,
    },
    PileupVariantMatrix {
        variants: Vec<f32>,
        read_variants: Vec<HashMap<i32, HashSet<i32>>>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        coverages: Vec<f32>,
        tids: Vec<i32>,
        total_indels_across_contigs: Vec<usize>,
        target_names: Vec<Vec<u8>>,
    }
}

impl PileupMatrix {
    pub fn new_contig_stats() -> PileupMatrix {
        PileupMatrix::PileupContigMatrix {
            variants: vec!(),
            kfrequencies: BTreeMap::new(),
            coverages: vec!(),
            tids: vec!(),
            total_indels_across_contigs: vec!(),
            target_names: vec!(),
        }
    }
    pub fn new_contig_variants() -> PileupMatrix {
        PileupMatrix::PileupVariantMatrix {
            variants: vec!(),
            read_variants: vec!(),
            kfrequencies: BTreeMap::new(),
            coverages: vec!(),
            tids: vec!(),
            total_indels_across_contigs: vec!(),
            target_names: vec!(),
        }
    }
}

pub trait PileupMatrixFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  pileup_stats: PileupStats,
                  number_of_contigs: usize);

    fn print_matrix(&self);

}

impl PileupMatrixFunctions for PileupMatrix{
    fn setup(&mut self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut variants,
                ref mut kfrequencies,
                ref mut coverages,
                ref mut tids,
                ref mut total_indels_across_contigs,
                ref mut target_names,
            } => {
                *variants = vec!();
                *kfrequencies = BTreeMap::new();
                *coverages = vec!();
                *tids = vec!();
                *total_indels_across_contigs = vec!();
                *target_names = vec!();
            },
            PileupMatrix::PileupVariantMatrix {
                ref mut variants,
                ref mut read_variants,
                ref mut kfrequencies,
                ref mut coverages,
                ref mut tids,
                ref mut total_indels_across_contigs,
                ref mut target_names,
            } => {
                *variants = vec!();
                *read_variants = vec!();
                *kfrequencies = BTreeMap::new();
                *coverages = vec!();
                *tids = vec!();
                *total_indels_across_contigs = vec!();
                *target_names = vec!();
            }
        }
    }

    fn add_contig(&mut self, mut pileup_stats: PileupStats, number_of_contigs: usize) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut variants,
                ref mut kfrequencies,
                ref mut coverages,
                ref mut tids,
                ref mut total_indels_across_contigs,
                ref mut target_names,
            } => {
                match pileup_stats {
                    PileupStats::PileupContigStats {
                        ref mut kfrequency,
                        ref mut variants_in_reads,
                        ref mut tid,
                        ref mut total_indels,
                        ref mut target_name,
                        ref mut variations_per_base,
                        ref mut coverage,
                        ref mut num_covered_bases,
                        ..
                    } => {
                        variants.push(*variations_per_base);
                        let contig_order_id = variants.len() as usize - 1;
                        for (tet, count) in kfrequency.iter() {
                            let mut count_vec = kfrequencies.entry(tet.to_vec())
                                .or_insert(vec![0; number_of_contigs]);
                            count_vec[contig_order_id] = *count;
                        }
                        coverages.push(*coverage);
                        tids.push(*tid);
                        total_indels_across_contigs.push(*total_indels);
                        target_names.push(target_name.clone());
                    }
                }
            },
            PileupMatrix::PileupVariantMatrix {
                ref mut variants,
                ref mut read_variants,
                ref mut kfrequencies,
                ref mut coverages,
                ref mut tids,
                ref mut total_indels_across_contigs,
                ref mut target_names,
            } => {
                match pileup_stats {
                    PileupStats::PileupContigStats {
                        ref mut kfrequency,
                        ref mut variants_in_reads,
                        ref mut tid,
                        ref mut total_indels,
                        ref mut target_name,
                        ref mut variations_per_base,
                        ref mut coverage,
                        ref mut num_covered_bases,
                        ..
                    } => {
                        variants.push(*variations_per_base);
                        read_variants.push(variants_in_reads.clone());
                        let contig_order_id = variants.len() as usize - 1;
                        for (tet, count) in kfrequency.iter() {
                            let mut count_vec = kfrequencies.entry(tet.to_vec())
                                                            .or_insert(vec![0; number_of_contigs]);
                            count_vec[contig_order_id] = *count;
                        }
                        coverages.push(*coverage);
                        tids.push(*tid);
                        total_indels_across_contigs.push(*total_indels);
                        target_names.push(target_name.clone());
                    }
                }
            }
        }
    }

    fn print_matrix(&self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variants,
                kfrequencies,
                coverages,
                tids,
                total_indels_across_contigs,
                target_names,
            } => {
                for i in 0..variants.len() {
                    print!("{}\t{}\t{}\t{}\t",
                             std::str::from_utf8(&target_names[i][..]).unwrap(),
                             variants[i],
                             coverages[i],
                             total_indels_across_contigs[i]);
                    for (kmer, counts) in kfrequencies.iter(){
                        print!("{}\t", counts[i]);
                    }
                    print!("\n");
                }

            },
            PileupMatrix::PileupVariantMatrix {
                variants,
                read_variants,
                kfrequencies,
                coverages,
                tids,
                total_indels_across_contigs,
                target_names,
            } => {
                for i in 0..variants.len() {
                    print!("{}\t{}\t{}\t{}\t",
                           std::str::from_utf8(&target_names[i][..]).unwrap(),
                           variants[i],
                           coverages[i],
                           total_indels_across_contigs[i]);
                    for (kmer, counts) in kfrequencies.iter(){
                        print!("{}\t", counts[i]);
                    }
                    print!("\n");
                }

            }
        }
    }
}