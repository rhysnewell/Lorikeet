#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use lorikeet_genome::utils::quality_utils::QualityUtils;
use mathru::statistics::distrib::Beta;
use statrs::function::beta::ln_beta;

fn test_leading_order_in_error_rate(num_ref: usize, num_alt: usize, error_rate: f64) {
    let qual = QualityUtils::error_prob_to_qual(error_rate);
    let alt_quals = vec![qual; num_alt];

    let calculated =
        HaplotypeCallerEngine::log_likelihood_ratio_constant_error(num_ref, num_alt, error_rate);
    let expected =
        ln_beta((num_ref + 1) as f64, (num_alt + 1) as f64) - (num_alt as f64) * error_rate.ln();
    println!("calc {} expected {}", calculated, expected);
    assert!(
        relative_eq!(calculated, expected, epsilon = 0.07),
        "expected {} calculated {}",
        expected,
        calculated
    );
}

#[test]
fn get_leading_order_data() {
    test_leading_order_in_error_rate(100, 5, 0.0001);
    test_leading_order_in_error_rate(100, 20, 0.0001);
    test_leading_order_in_error_rate(100, 200, 0.0001);
    test_leading_order_in_error_rate(10, 2, 0.0001);
    test_leading_order_in_error_rate(1000, 2, 0.00001);
    test_leading_order_in_error_rate(10000, 2, 0.000001);
}
