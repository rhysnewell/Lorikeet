#![allow(
    non_upper_case_globals,
    non_snake_case,
    non_camel_case_types
)]
// extern crate openssl;
// extern crate openssl_sys;

pub mod abundance;
pub mod activity_profile;
pub mod ani_calculator;
pub mod annotator;
pub mod assembly;
pub mod bam_parsing;
pub mod cli;
pub mod evolve;
pub mod external_command_checker;
pub mod genotype;
pub mod graphs;
pub mod haplotype;
pub mod linkage;
pub mod model;
pub mod pair_hmm;
pub mod processing;
pub mod read_error_corrector;
pub mod read_orientation;
pub mod read_threading;
pub mod reads;
pub mod reference;
pub mod smith_waterman;
pub mod test_utils;
pub mod utils;

// Stats
#[macro_use]
extern crate ndarray;
// Utilities
#[macro_use]
extern crate enum_ordinalize;
#[macro_use]
extern crate log;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate approx;

use std::process;

pub const AUTHOR: &str =
    "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology";
pub const AUTHOR_AND_EMAIL: &str =
    "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <rhys.newell94 near gmail.com>";

pub fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.contains_id(parameter) {
        true => {
            let mut percentage: f32 = *m.get_one(parameter).unwrap();
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
