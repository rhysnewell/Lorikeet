extern crate criterion;
extern crate rayon;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rayon::prelude::*;

fn bench_single_iter(c: &mut Criterion) -> &mut Criterion {
    let test_vec = (0..1000).into_iter().collect::<Vec<usize>>();
    c.bench_function("single_iter", |b| {
        b.iter(|| test_single_iter(black_box(&test_vec)))
    })
}

fn test_single_iter(test_vec: &Vec<usize>) {
    let mut sum = 0;
    let mut mult = 1;

    for a in test_vec.iter() {
        sum += a;
        mult += a;
    }
}

fn bench_par_iters(c: &mut Criterion) -> &mut Criterion {
    let test_vec = (0..1000).into_iter().collect::<Vec<usize>>();
    c.bench_function("par_iters", |b| {
        b.iter(|| test_par_iters(black_box(&test_vec)))
    })
}

fn test_par_iters(test_vec: &Vec<usize>) {
    let _sum = test_vec.par_iter().sum::<usize>();
    let _sum = test_vec.par_iter().sum::<usize>();
}

fn bench_collect(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("collect", |b| b.iter(|| test_into_iter()))
}
fn test_into_iter() {
    let _test_vec = (0..100000).into_iter().collect::<Vec<usize>>();
}

fn bench_par_collect(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("par_collect", |b| b.iter(|| test_into_par_iter()))
}
fn test_into_par_iter() {
    let _test_vec = (0..100000).into_par_iter().collect::<Vec<usize>>();
}

fn bench_map_collect(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("map_collect", |b| b.iter(|| test_map_collect()))
}
fn test_map_collect() {
    let _test_vec = (0..100000)
        .into_iter()
        .map(|i| i * 1000)
        .collect::<Vec<usize>>();
}

fn bench_par_map_collect(c: &mut Criterion) -> &mut Criterion {
    c.bench_function("par_map_collect", |b| {
        b.iter(|| test_into_par_map_collect())
    })
}
fn test_into_par_map_collect() {
    let _test_vec = (0..100000)
        .into_par_iter()
        .map(|i| i * 1000)
        .collect::<Vec<usize>>();
}

criterion_group!(
    benches,
    bench_single_iter,
    bench_par_iters,
    bench_collect,
    bench_par_collect,
    bench_map_collect,
    bench_par_map_collect,
);
criterion_main!(benches);
