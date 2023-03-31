#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::bam_parsing::bam_generator::generate_indexed_named_bam_readers_from_bam_files;
use lorikeet_genome::bam_parsing::bam_generator::IndexedNamedBamReader;
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

    bam_generated
        .fetch((
            0, start, // fetch is possibly 1-based
            stop,
        ))
        .expect("Unable to fetch coordinates");

    let mut records = Vec::new();
    let mut record = Record::new();
    let mut first = true;
    while bam_generated.read(&mut record) == true {
        assert!(record.pos() <= stop);
        records.push(record.clone());
    }
}
