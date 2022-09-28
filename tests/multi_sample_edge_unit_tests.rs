#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

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

use lorikeet_genome::graphs::base_edge::BaseEdge;
use lorikeet_genome::graphs::multi_sample_edge::MultiSampleEdge;
use lorikeet_genome::utils::utils::make_permutations;
use std::cmp::max;

struct MultiplicityTestProvider {
    counts_per_sample: Vec<usize>,
    num_samples_pruning: usize,
}

impl MultiplicityTestProvider {
    pub fn new(counts_per_sample: Vec<usize>, num_samples_pruning: usize) -> Self {
        Self {
            counts_per_sample,
            num_samples_pruning,
        }
    }
}

fn test_multiplicity(cfg: MultiplicityTestProvider) {
    let mut edge = MultiSampleEdge::new(false, 0, cfg.num_samples_pruning);
    assert_eq!(edge.get_multiplicity(), 0);
    assert_eq!(edge.get_pruning_multiplicity(), 0);

    let mut copy = edge.clone();
    assert_eq!(
        edge.get_current_single_sample_multiplicity(),
        copy.get_current_single_sample_multiplicity()
    );
    assert_eq!(edge.get_dot_label(), copy.get_dot_label());
    assert_eq!(
        edge.get_pruning_multiplicity(),
        copy.get_pruning_multiplicity()
    );
    assert_eq!(edge.get_multiplicity(), copy.get_multiplicity());

    let mut total = 0;
    for i in 0..cfg.counts_per_sample.len() {
        let mut count_for_sample = 0;
        for count in 0..cfg.counts_per_sample[i] {
            edge.inc_multiplicity(1);
            total += 1;
            count_for_sample += 1;
            assert_eq!(edge.get_multiplicity(), total);
            assert_eq!(
                edge.get_current_single_sample_multiplicity(),
                count_for_sample
            );
        }
        edge.flush_single_sample_multiplicity();
    }

    let mut counts = cfg.counts_per_sample.clone();
    counts.push(0);
    counts.sort();
    let prune = counts[max(counts.len() as i32 - cfg.num_samples_pruning as i32, 0) as usize];
    assert_eq!(edge.get_multiplicity(), total);
    assert_eq!(edge.get_pruning_multiplicity(), prune);
}

#[test]
fn make_multiplicity_data() {
    let counts_per_sample = vec![0, 1, 2, 3, 4, 5];
    for num_samples_pruning in vec![1, 2, 3] {
        for n_samples in vec![1, 2, 3, 4, 5] {
            for perm in make_permutations(&counts_per_sample, n_samples, false) {
                test_multiplicity(MultiplicityTestProvider::new(perm, num_samples_pruning))
            }
        }
    }
}
