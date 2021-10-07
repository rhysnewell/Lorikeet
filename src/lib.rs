#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case,
    unused
)]
pub mod activity_profile;
pub mod annotator;
pub mod assembly;
pub mod cli;
pub mod estimation;
pub mod external_command_checker;
pub mod genotype;
pub mod graphs;
pub mod haplotype;
pub mod model;
pub mod pair_hmm;
pub mod read_error_corrector;
pub mod read_orientation;
pub mod read_threading;
pub mod reads;
pub mod reference;
pub mod smith_waterman;
pub mod test_utils;
pub mod utils;

// HTS and bio files
extern crate bio;
extern crate bio_types;
extern crate rust_htslib;
extern crate seq_io;

// Birds and CoverM
extern crate bird_tool_utils;
extern crate coverm;
extern crate galah;

// Stats
extern crate kodama;
extern crate linregress;
#[macro_use]
extern crate ndarray;
extern crate mathru;
extern crate ndarray_npy;
extern crate statrs;

// Utilities
#[macro_use]
extern crate approx;
extern crate clap;
extern crate compare;
extern crate csv;
extern crate env_logger;
extern crate glob;
extern crate hashlink;
extern crate indexmap;
extern crate itertools;
extern crate libm;
extern crate multimap;
extern crate nix;
extern crate num;
extern crate ordered_float;
extern crate petgraph;
extern crate rand;
extern crate rayon;
extern crate scoped_threadpool;
extern crate strum;
extern crate tempdir;
extern crate tempfile;

#[macro_use]
extern crate log;
extern crate derive_builder;
extern crate indicatif;
extern crate strum_macros;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate lazy_static;
extern crate derive_new;
extern crate pest_derive;
#[macro_use]
extern crate enum_ordinalize;
extern crate term;

use clap::*;
use std::process;

pub fn parse_percentage(m: &ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage: f32 = m.value_of(parameter).unwrap().parse().unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                process::exit(1);
            }
            info!("Using {} {}%", parameter, percentage * 100.0);
            percentage
        }
        false => 0.0,
    }
}

// Enum for exclusion out here so long read can find it
pub enum GenomeExclusionTypes {
    SeparatorType,
    NoneType,
    GenomesAndContigsType,
}
