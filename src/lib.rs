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
pub mod dbscan;

extern crate bio;
extern crate clap;
extern crate seq_io;
extern crate bio_types;
extern crate linregress;
extern crate coverm;
extern crate bird_tool_utils;
extern crate nymph;
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
extern crate ndarray;
extern crate ndarray_linalg;
extern crate openblas_src;
extern crate strum;

extern crate approx;
#[macro_use]
extern crate log;
extern crate strum_macros;
extern crate derive_builder;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate derive_new;
extern crate pest_derive;

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