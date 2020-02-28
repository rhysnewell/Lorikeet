pub mod pileups;
pub mod pileup_structs;
pub mod pileup_matrix;
pub mod cli;
pub mod external_command_checker;
pub mod codon_structs;
pub mod matrix_handling;
pub mod variants;
pub mod estimation;
pub mod model;
pub mod utils;
pub mod factorization;
pub mod dbscan;

extern crate bio;
extern crate clap;
extern crate seq_io;
extern crate bio_types;
extern crate linregress;
extern crate coverm;
extern crate bird_tool_utils;
extern crate csv;
extern crate statrs;
extern crate ordered_float;
extern crate rust_htslib;
extern crate env_logger;
extern crate nix;
extern crate tempdir;
extern crate tempfile;
extern crate rand;
extern crate itertools;
extern crate rayon;
#[macro_use]
extern crate ndarray;
extern crate ndarray_linalg;
extern crate strum;

#[macro_use]
extern crate approx;
#[macro_use]
extern crate log;
#[macro_use]
extern crate strum_macros;
#[macro_use]
extern crate derive_builder;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate derive_new;
#[macro_use]
extern crate pest_derive;

use std::path::Path;
use coverm::genomes_and_contigs::GenomesAndContigs;
use std::io::Read;


pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";


/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_to_genome(){
        let mut contig_to_genome = GenomesAndContigs::new();
        let genome = String::from("genome0");
        let index = contig_to_genome.establish_genome(genome);
        contig_to_genome.insert(String::from("contig1"), index);
        assert_eq!(
            String::from("genome0"),
            *(contig_to_genome.genome_of_contig(&String::from("contig1")).unwrap()));
    }

    #[test]
    fn test_read_genome_fasta_files_one_genome(){
        let contig_to_genome = read_genome_fasta_files(&vec!["tests/data/genome1.fna"]);
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq1")).unwrap());
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq2")).unwrap());
    }
}



