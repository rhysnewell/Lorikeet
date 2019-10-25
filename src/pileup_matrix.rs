use std::collections::{HashMap, BTreeMap};
use pileup_structs::*;
use std::str;
use std::path::Path;
use std::fs::File;
use std::io::prelude::*;
use rayon::prelude::*;
use itertools::izip;

#[derive(Debug, Clone)]
pub enum PileupMatrix {
    PileupContigMatrix {
        coverages: HashMap<i32, Vec<f32>>,
        average_genotypes: HashMap<i32, Vec<f32>>,
        variances: HashMap<i32, Vec<f32>>,
        target_names: HashMap<i32, String>,
        target_lengths: HashMap<i32, f32>,
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
                  pileup_stats: PileupStats,
                  sample_count: usize,
                  sample_idx: usize);

    fn print_stats(&self, output_prefix: &str);

    fn print_kmers(&self, output_prefix: &str, kmer_size: &usize);

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

    fn add_contig(&mut self, mut pileup_stats: PileupStats, sample_count: usize, sample_idx: usize) {
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
                        ref mut mean_genotypes,
                        ref mut variance,
                        ..
                    } => {
                        let ag = average_genotypes.entry(*tid).or_insert(
                            vec![0.0 as f32; sample_count]);
                        ag[sample_idx] = *mean_genotypes;
                        let var = variances.entry(*tid).or_insert(
                            vec![0.0 as f32; sample_count]
                        );
                        var[sample_idx] = *variance;
                        let cov = coverages.entry(*tid).or_insert(
                            vec![0.0 as f32; sample_count]
                        );
                        cov[sample_idx] = *coverage;
                        target_names.insert(*tid,
                                            str::from_utf8(target_name).unwrap().to_string());
                        target_lengths.insert(*tid, target_len.clone());
                    }
                }
            }
        }
    }

    fn print_stats(&self, output_prefix: &str) {
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
                let file_name = output_prefix.to_string() + &"_".to_owned()
                    + &"contig_stats".to_owned()
                    + &".tsv".to_owned();
                let file_path = Path::new(&file_name);
                let mut file_open = match File::create(file_path) {
                    Ok(fasta) => fasta,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    },
                };
                write!(file_open, "contigName\tcontigLen\ttotalAvgDepth\ttotalAvgGeno").unwrap();
                for sample_name in sample_names.iter(){
                    write!(file_open, "\t{}.bam\t{}.bam-var\t{}.bam-gen",
                           &sample_name, &sample_name, &sample_name).unwrap();
                }
                write!(file_open, "\n").unwrap();
                for (tid, contig_name) in target_names.iter() {
                    write!(file_open, "{}\t{}", contig_name, target_lengths[tid]).unwrap();
                    let placeholder = vec![0.0 as f32; sample_names.len() as usize];
                    let coverage_vec = match coverages.get(tid) {
                        Some(vector) => vector,
                        None => &placeholder,
                    };
                    let coverage_sum: f32 = coverage_vec.par_iter().sum();
                    write!(file_open, "\t{}", coverage_sum/coverage_vec.len() as f32).unwrap();
                    let variance_vec = match variances.get(tid) {
                        Some(vector) => vector,
                        None => &placeholder,
                    };
                    let genotype_vec = match average_genotypes.get(tid) {
                        Some(vector) => vector,
                        None => &placeholder,
                    };
                    let genotype_sum: f32 = genotype_vec.par_iter().sum();
                    write!(file_open, "\t{}", genotype_sum/genotype_vec.len() as f32).unwrap();

                    for (coverage, variance, genotypes) in izip!(coverage_vec,
                                                                 variance_vec,
                                                                 genotype_vec){
                        write!(file_open, "\t{}\t{}\t{}", coverage, variance, genotypes).unwrap();
                    }
                    write!(file_open, "\n").unwrap();
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