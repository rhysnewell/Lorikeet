use std::collections::{HashMap, HashSet};
use std::collections::BTreeMap;
use pileup_structs::*;
use std::str;
use std::fs::File;
use std::io::prelude::*;

#[derive(Debug, Clone)]
pub enum PileupMatrix {
    PileupContigMatrix {
        coverages: HashMap<i32, Vec<f32>>,
        average_genotypes: HashMap<i32, Vec<f32>>,
        variances: HashMap<i32, Vec<f32>>,
        target_names: HashMap<i32, String>,
        target_lengths: HashMap<i32, usize>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
    }
}

impl PileupMatrix {
    pub fn new_matrix() -> PileupMatrix {
        PileupMatrix::PileupContigMatrix {
            coverages: HashMap::new(),
            average_genotypes: HashMap::new(),
            variances: HashMap::new(),
            target_names: HashMap::new(),
            target_lengths: HashMap::new(),
            sample_names: vec!(),
            kfrequencies: BTreeMap::new(),
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
                  pileup_stats: PileupStats);

    fn print_stats(&self);

    fn print_kmers(&self);

}

impl PileupMatrixFunctions for PileupMatrix{
    fn setup(&mut self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut variances,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ref mut kfrequencies,
            } => {
                *coverages = HashMap::new();
                *average_genotypes = HashMap::new();
                *variances = HashMap::new();
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

    fn add_contig(&mut self, mut pileup_stats: PileupStats) {
        match self {
            PileupMatrix::PileupContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut variances,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                match pileup_stats {
                    PileupStats::PileupContigStats {
                        ref mut tid,
                        ref mut target_name,
                        ref mut target_len,
                        ref mut coverage,
                        ref mut variations_per_base,
                        ref mut mean_genotypes,
                        ref mut variance,
                        ..
                    } => {
                        let ag = average_genotypes.entry(*tid).or_insert(vec!());
                        ag.push(*mean_genotypes);
                        let var = variances.entry(*tid).or_insert(vec!());
                        var.push(*variance);
                        let cov = coverages.entry(*tid).or_insert(vec!());
                        cov.push(*coverage);
                        target_names.insert(*tid,
                                            str::from_utf8(target_name).unwrap().to_string());
                        target_lengths.insert(*tid, target_len.clone());
                    }
                }
            }
        }
    }

    fn print_stats(&self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                variances,
                average_genotypes,
                coverages,
                target_names,
                target_lengths,
                sample_names,
                ..
            } => {
//                for i in 0..variants.len() {
//                    print!("{}\t",
//                           std::str::from_utf8(&target_names[i][..]).unwrap(),
////                             variants[i],
////                             coverages[i],
////                             total_indels_across_contigs[i]
//                    );
//                    print!("\n");
//                }
            }
        }
    }

    fn print_kmers(&self) {
        match self {
            PileupMatrix::PileupContigMatrix {
                kfrequencies,
                target_names,
                ..
            } => {
                for (tid, name) in target_names.iter() {
                    print!("{}\t",
                           name);
                    for (_kmer, counts) in kfrequencies.iter(){
                        print!("{}\t", counts[*tid as usize]);
                    }
                    print!("\n");
                }
            }
        }
    }
}