extern crate gkl;
extern crate lorikeet_genome;
extern crate rust_htslib;

use gkl::pairhmm::{forward};
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::allele_likelihoods::AlleleLikelihoods;
use lorikeet_genome::pair_hmm::pair_hmm::PairHMM;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::{
    AVXMode, PairHMMInputScoreImputator,
};
use lorikeet_genome::reads::read_utils::ReadUtils;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use rust_htslib::bam::record::{Cigar, CigarString};
use std::cmp::max;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

#[test]
fn test_likelihoods_avx() {
    let file = File::open("tests/resources/pairhmm-testdata.txt").unwrap();
    let lines = BufReader::new(file).lines();
    for line in lines {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let hap_bases = tokens.next().unwrap().as_bytes();
        let hap = Haplotype::new(hap_bases, true);
        let bases = tokens.next().unwrap().as_bytes();
        let parse_qual = |tokens: &mut std::str::SplitWhitespace, min| -> Vec<u8> {
            tokens
                .next()
                .unwrap()
                .as_bytes()
                .iter()
                .copied()
                .map(|b| max(min, b - 33))
                .collect()
        };
        let base_quals = &parse_qual(&mut tokens, 6);
        let insertion_quals = &parse_qual(&mut tokens, 0);
        let deletion_quals = &parse_qual(&mut tokens, 0);
        let gcp = &parse_qual(&mut tokens, 0);
        let expected_result = tokens.next().unwrap().parse::<f64>().unwrap();
        let read_length = bases.len();
        // println!("")
        let avx_function = forward().unwrap();
        let result = avx_function(
            hap_bases,
            bases,
            base_quals,
            insertion_quals,
            deletion_quals,
            gcp,
        );

        println!("Direct result {}", result);
        // This should work but doesn't??
        assert!((result - expected_result).abs() < 1e-5);

        println!(
            "quals {:?} {:?} {:?} {:?}",
            &base_quals, &insertion_quals, &deletion_quals, &gcp
        );
        let mut read = ArtificialReadUtils::create_artificial_read(
            bases,
            base_quals,
            CigarString::from(vec![Cigar::Match(read_length as u32)]),
        );
        ReadUtils::set_insertion_base_qualities(&mut read, insertion_quals);
        ReadUtils::set_deletion_base_qualities(&mut read, deletion_quals);

        let mut read_map = HashMap::new();
        read_map.insert(0, vec![read.clone()]);

        let hap_vec = vec![hap.clone()];
        let mut hmm = PairHMM::initialize(&hap_vec, &read_map, AVXMode::AVX);
        let mut likelihoods =
            AlleleLikelihoods::new(hap_vec.clone(), vec!["sample".to_string()], read_map);

        let score_imputator = PairHMMInputScoreImputator::new(gcp[0]);

        hmm.compute_log10_likelihoods(0, &mut likelihoods, vec![read.clone()], &score_imputator);
        let la = hmm.get_log_likelihood_array();

        assert!((la[0] - expected_result).abs() < 1e-5, "Likelihood not in expected range for PairHMM implementation: got {} expected {}, diff {}", la[0], expected_result, la[0] - expected_result);
    }
}
