#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate coverm;
extern crate lorikeet_genome;
extern crate rayon;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;
extern crate bio;
extern crate itertools;
extern crate petgraph;
extern crate rand;
extern crate term;
#[macro_use]
extern crate ntest;

use coverm::bam_generator::generate_indexed_named_bam_readers_from_bam_files;
use coverm::bam_generator::IndexedNamedBamReader;
use rust_htslib::bam::Record;

#[test]
fn test_fetch_coordinates() {
    let mut bam_generated = generate_indexed_named_bam_readers_from_bam_files(
        vec!["tests/data/ben/bams_short/lorikeet-genome.random10000.20_differences.sim_reads.1.fq.bam"],
        1,
    )
        .into_iter()
        .next()
        .unwrap();

    let start = 1000;
    let stop = 2000;

    bam_generated.fetch((
        0, start, // fetch is possibly 1-based
        stop,
    ));

    let mut records = Vec::new();
    let mut record = Record::new();
    let mut first = true;
    while bam_generated.read(&mut record) == true {
        assert!(record.pos() <= stop);
        records.push(record.clone());
    }
}
