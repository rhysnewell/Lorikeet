#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;
extern crate rand;

use lorikeet_genome::genotype::genotype_allele_counts::GenotypeAlleleCounts;
use lorikeet_genome::genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use lorikeet_genome::model::allele_list::AlleleList;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reads::read_clipper::ReadClipper;
use lorikeet_genome::test_utils::read_clipper_test_utils::ReadClipperTestUtils;
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::math_utils::{MathUtils, RunningAverage};
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use rand::distributions::{Distribution, Normal};
use rand::rngs::ThreadRng;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::ops::Deref;
use std::sync::Mutex;

#[test]
fn test_running_average() {
    let numbers = vec![1, 2, 4, 5, 3, 128, 25678, -24];
    let mut r = RunningAverage::new();

    for b in numbers.iter() {
        r.add(*b as f64)
    }

    assert_eq!(numbers.len(), r.obs_count());
    assert!(r.mean() - 3224.625 < 2e-10);
    assert!(
        r.stddev() - 9072.6515881128 < 2e-10,
        "value {}",
        r.stddev() - 9072.6515881128
    );
}

#[test]
fn test_approximate_log_sum_log() {
    let required_precision = 1e-4;
    assert!(
        relative_eq!(
            MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 3),
            3.0_f64.log10(),
            epsilon = required_precision
        ),
        "Value {}",
        MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 3)
    );
    assert!(relative_eq!(
        MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 0),
        std::f64::NEG_INFINITY,
        epsilon = required_precision
    ));
    assert!(relative_eq!(
        MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 3),
        3.0_f64.log10(),
        epsilon = required_precision
    ));
    assert!(relative_eq!(
        MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 2),
        2.0_f64.log10(),
        epsilon = required_precision
    ));
    assert!(relative_eq!(
        MathUtils::approximate_log10_sum_log10_vec(vec![0.0, 0.0, 0.0].as_slice(), 0, 1),
        0.0_f64,
        epsilon = required_precision
    ));
    assert!(relative_eq!(
        MathUtils::approximate_log10_sum_log10_vec(
            vec![
                std::f64::NEG_INFINITY,
                std::f64::NEG_INFINITY,
                std::f64::NEG_INFINITY,
            ]
            .as_slice(),
            0,
            3
        ),
        std::f64::NEG_INFINITY,
        epsilon = required_precision
    ));

    let mut rnd = ThreadRng::default();
    let normal = Normal::new(0.0, 1.0);

    for j in 0..5 {
        for i in 0..5 {
            let a: f64 = (1 + 3 * j) as f64 * normal.sample(&mut rnd);
            let b: f64 = (1 + 3 * j) as f64 * normal.sample(&mut rnd);
            let c: f64 = (1 + 3 * j) as f64 * normal.sample(&mut rnd);

            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10_vec(vec![a].as_slice(), 0, 1),
                a,
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10_vec(
                    vec![a, std::f64::NEG_INFINITY].as_slice(),
                    0,
                    1
                ),
                a,
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10_vec(vec![a, b].as_slice(), 0, 2),
                (10.0_f64.powf(a) + 10.0_f64.powf(b)).log10(),
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10(a, b),
                (10.0_f64.powf(a) + 10.0_f64.powf(b)).log10(),
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10(b, a),
                (10.0_f64.powf(a) + 10.0_f64.powf(b)).log10(),
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::approximate_log10_sum_log10_vec(vec![a, b, c].as_slice(), 0, 3),
                (10.0_f64.powf(a) + 10.0_f64.powf(b) + 10.0_f64.powf(c)).log10(),
                epsilon = required_precision
            ));
        }
    }
}

#[test]
fn test_approximate_log_sum_log_on_slice() {
    let required_precision = 1e-4;
    assert!(
        relative_eq!(
            MathUtils::approximate_log10_sum_log10_vec(vec![-32.0, -39.0, -46.0].as_slice(), 0, 3),
            -31.9999,
            epsilon = required_precision
        ),
        "Value {}",
        MathUtils::approximate_log10_sum_log10_vec(vec![-32.0, -39.0, -46.0].as_slice(), 0, 3)
    );
    assert!(
        relative_eq!(
            MathUtils::approximate_log10_sum_log10_vec(
                vec![-35.0, -32.0, -39.0, -46.0, -48.0].as_slice(),
                1,
                4
            ),
            -31.9999,
            epsilon = required_precision
        ),
        "Value {}",
        MathUtils::approximate_log10_sum_log10_vec(
            vec![-35.0, -32.0, -39.0, -46.0, -48.0].as_slice(),
            1,
            4
        )
    );
}

#[test]
fn test_log_10_sum_log10() {
    let required_precision = 1e-14;

    let log3 = 0.477121254719662;
    assert!(relative_eq!(
        MathUtils::log10_sum_log10(&vec![0.0, 0.0, 0.0], 0, 3),
        log3,
        epsilon = required_precision
    ));

    assert!(relative_eq!(
        MathUtils::log10_sum_log10(&vec![-5.15], 0, 1),
        -5.15,
        epsilon = required_precision
    ));
    assert!(relative_eq!(
        MathUtils::log10_sum_log10(&vec![130.0], 0, 1),
        130.0,
        epsilon = required_precision
    ));

    assert!(relative_eq!(
        MathUtils::log10_sum_log10(&vec![0.0, 0.0], 0, 2),
        (10.0_f64.powf(0.0) + 10.0_f64.powf(0.0)).log10(),
        epsilon = required_precision
    ));

    let mult_partition_factor = vec![0.999, 0.98, 0.95, 0.90, 0.8, 0.5, 0.3, 0.1, 0.05, 0.001];
    let n_partitions = vec![2, 4, 8, 16, 32, 64, 128, 256, 512, 1028];

    for alpha in mult_partition_factor {
        let log_alpha = (alpha as f64).log10();
        let log_one_minus_alpha = (1.0 - alpha as f64).log10();
        for n_part in n_partitions.clone() {
            let mut multiplicative = vec![0.0; n_part];
            let mut equal = vec![0.0; n_part];

            let mut remaining_log = 0.0; // realspace = 1
            for i in 0..n_part - 1 {
                equal[i] = -((n_part as f64).log10());
                let piece = remaining_log + log_alpha; // take a*remaining, leaving remaining-a*remaining = (1-a)*remaining
                multiplicative[i] = piece;
                remaining_log = remaining_log + log_one_minus_alpha;
            }
            equal[n_part - 1] = -((n_part as f64).log10());
            multiplicative[n_part - 1] = remaining_log;
            assert!(relative_eq!(
                MathUtils::log10_sum_log10(&equal, 0, equal.len()),
                0.0,
                epsilon = required_precision
            ));
            assert!(relative_eq!(
                MathUtils::log10_sum_log10(&multiplicative, 0, multiplicative.len()),
                0.0,
                epsilon = required_precision
            ));
        }
    }
}

#[test]
fn test_normalize() {
    let error = 1e-6;

    let normalized_log10 = MathUtils::normalize_log10(
        vec![3.0_f64.log10(), 2.0_f64.log10(), 1.0_f64.log10()],
        true,
    );

    let normalized_log10_expected = vec![
        (3.0_f64 / 6.0_f64).log10(),
        (2.0_f64 / 6.0_f64).log10(),
        (1.0_f64 / 6.0_f64).log10(),
    ];

    assert_equals_double_array(&normalized_log10, &normalized_log10_expected, error);
}

#[test]
fn test_log10_factorial() {
    assert!(relative_eq!(
        MathUtils::log10_factorial(4.0),
        1.3802112,
        epsilon = 1e-6
    ));
    assert!(relative_eq!(
        MathUtils::log10_factorial(10.0),
        6.559763,
        epsilon = 1e-6
    ));
    assert!(relative_eq!(
        MathUtils::log10_factorial(200.0),
        374.896888,
        epsilon = 1e-6
    ));
    assert!(relative_eq!(
        MathUtils::log10_factorial(12342.0),
        45138.2626503,
        epsilon = 1e-6
    ));
}

fn assert_equals_double_array(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "Arrays of differing lengths actual {} -> expected {}",
        actual.len(),
        expected.len()
    );
    for (i, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            relative_eq!(a, e, epsilon = tolerance),
            "Actual {}, Expected {} at index {}",
            a,
            e,
            i
        );
    }
}
