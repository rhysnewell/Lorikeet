pub mod cli;
pub mod dbscan;
pub mod estimation;
pub mod external_command_checker;
pub mod model;
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
extern crate ndarray;
extern crate ndarray_npy;
extern crate statrs;
extern crate libm;

// Utilities
extern crate clap;
extern crate csv;
extern crate env_logger;
extern crate glob;
extern crate itertools;
extern crate nix;
extern crate ordered_float;
extern crate rand;
extern crate rayon;
extern crate scoped_threadpool;
extern crate tempdir;
extern crate tempfile;

//extern crate plotly;
extern crate strum;

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

use clap::*;
use std::process;

pub fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(m.value_of(parameter), f32).unwrap();
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
