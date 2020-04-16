pub mod pileups;
pub mod pileup_structs;
pub mod pileup_matrix;
pub mod cli;
pub mod external_command_checker;
pub mod codon_structs;
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
extern crate csv;
extern crate itertools;
extern crate kodama;
extern crate statrs;
extern crate ordered_float;
extern crate env_logger;
extern crate nix;
extern crate tempdir;
extern crate tempfile;
extern crate rand;
extern crate rayon;
extern crate rust_htslib;

//extern crate plotly;
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

use clap::*;
use std::process;

pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";

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