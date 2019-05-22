extern crate strainm;
use strainm::bam_generator::*;
use strainm::bam_generator::*;
use strainm::filter;
use strainm::external_command_checker;
use strainm::mapping_parameters::*;
use strainm::FlagFilter;
use strainm::CONCATENATED_FASTA_FILE_SEPARATOR;
use strainm::genomes_and_contigs::GenomesAndContigs;


extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::env;
use std::str;
use std::process;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate tempfile;
use tempfile::NamedTempFile;
#[macro_use]
extern crate lazy_static;


fn main() {

}
