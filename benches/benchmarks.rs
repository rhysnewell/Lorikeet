extern crate assert_cli;
extern crate criterion;
extern crate rayon;
extern crate tempdir;

use assert_cli::Assert;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rayon::prelude::*;

fn test_small_variant_call_no_threads() {
    let td = tempdir::TempDir::new("tmp/").unwrap();

    Assert::main_binary()
        .with_args(&[
            "call",
            "-r",
            "tests/data/ben/random10000.fna",
            "-c",
            "tests/data/ben/four_strains/random10000.10_differences.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.10_differences.sim_reads.2.fq",
            "tests/data/ben/four_strains/random10000.20_differences.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.20_differences.sim_reads.2.fq",
            "tests/data/ben/four_strains/random10000.combined_differences.1.fq",
            "tests/data/ben/four_strains/random10000.combined_differences.2.fq",
            "tests/data/ben/four_strains/random10000.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.sim_reads.2.fq",
            "tests/data/ben/random10000.5_differences.shared_locations1.fq",
            "tests/data/ben/random10000.5_differences.shared_locations2.fq",
            "-o",
            td.path().to_str().unwrap(),
        ])
        .succeeds()
        .unwrap()
}

fn bench_small_variant_call_no_threads(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("small variant call: no threads", |b| {
        b.iter(|| test_small_variant_call_no_threads())
    })
}

fn test_small_variant_call_eight_threads() {
    let td = tempdir::TempDir::new("tmp/").unwrap();

    Assert::main_binary()
        .with_args(&[
            "call",
            "-r",
            "tests/data/ben/random10000.fna",
            "-c",
            "tests/data/ben/four_strains/random10000.10_differences.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.10_differences.sim_reads.2.fq",
            "tests/data/ben/four_strains/random10000.20_differences.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.20_differences.sim_reads.2.fq",
            "tests/data/ben/four_strains/random10000.combined_differences.1.fq",
            "tests/data/ben/four_strains/random10000.combined_differences.2.fq",
            "tests/data/ben/four_strains/random10000.sim_reads.1.fq",
            "tests/data/ben/four_strains/random10000.sim_reads.2.fq",
            "tests/data/ben/random10000.5_differences.shared_locations1.fq",
            "tests/data/ben/random10000.5_differences.shared_locations2.fq",
            "-t",
            "8",
            "-o",
            td.path().to_str().unwrap(),
        ])
        .succeeds()
        .unwrap()
}

fn bench_small_variant_call_eight_threads(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("small variant call: 8 threads", |b| {
        b.iter(|| test_small_variant_call_eight_threads())
    })
}

criterion_group!(
    benches,
    bench_small_variant_call_no_threads,
    bench_small_variant_call_eight_threads
);
criterion_main!(benches);
