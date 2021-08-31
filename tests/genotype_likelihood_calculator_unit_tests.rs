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

use lorikeet_genome::genotype::genotype_allele_counts::GenotypeAlleleCounts;
use lorikeet_genome::genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reads::read_clipper::ReadClipper;
use lorikeet_genome::test_utils::read_clipper_test_utils::ReadClipperTestUtils;
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::ops::Deref;
use std::sync::Mutex;

lazy_static! {
    static ref MAXIMUM_ALLELE: Vec<usize> = vec![1, 2, 5, 6];
    static ref PLOIDY: Vec<usize> = vec![1, 2, 3, 20];
    static ref READ_COUNTS: Vec<Vec<usize>> = vec![
        vec![10, 100, 50],
        vec![0, 100, 10, 1, 50],
        vec![1, 2, 3, 4, 20],
        vec![10, 0],
    ];
}

fn test_ploidy_and_maximum_allele(ploidy: usize, allele_count: usize) {
    let mut calculator = GenotypeLikelihoodCalculators::get_instance(ploidy, allele_count);
    assert_eq!(calculator.ploidy, ploidy);
    assert_eq!(calculator.allele_count, allele_count);
    assert_eq!(
        calculator.genotype_count as usize,
        calculate_genotype_count(ploidy, allele_count),
        "Ploidy {} Allele count {}",
        ploidy,
        allele_count
    );
    let genotype_count = calculator.genotype_count;
    let test_genotype_count = min(30000, genotype_count as usize);
    for i in 0..test_genotype_count {
        let allele_counts = calculator.genotype_allele_counts_at(i);
        if i > 0 {
            // assert!(calculator.genotype_allele_counts_at(i - 1) < allele_counts);
        };

        // if i % 100 == 0 {
        //     println!("i {}, allele_counts {}", i, allele_counts.distinct_allele_count());
        // }
        let mut allele_array = vec![0; ploidy];
        let mut index = 0;
        for j in 0..allele_counts.distinct_allele_count() {
            for allele_count in
                allele_array[index..(index + allele_counts.allele_count_at(j))].iter_mut()
            {
                *allele_count = allele_counts.allele_index_at(j);
            }
            index += allele_counts.allele_count_at(j);
        }
        let mut allele_count_array = vec![0; allele_counts.distinct_allele_count() << 1];
        allele_counts.copy_allele_counts(&mut allele_count_array, 0);
        assert_eq!(index, ploidy);
        assert_eq!(calculator.alleles_to_index(allele_array.as_slice()), i);
        assert_eq!(
            calculator.allele_counts_to_index(allele_count_array.as_slice()),
            i
        );
    }
}

// fn test_likelihood_calculation(ploidy: usize, allele_count: usize, read_count: &[usize]) {
//     let mut read_likelihoods = ReadLikelihoodsUnitTester::read_likelihoods(allele_count, read_count);
//
//     let mut calculator = GenotypeLikelihoodCalculators::get_instance(ploidy, allele_count);
//     let genotype_count = calculator.genotype_count;
//     let test_genotype_count = min(30000, genotype_count);
//     let sample_count = read_count.len();
//
//     let permutation = AlleleLikelihoodMatrixMapper::new(read_likelihoods.get_allele_list().permutation(read_likelihoods.get_allele_list()));
//     // println!("Read Counts {:?} read likelihoods {:?}", read_count, &read_likelihoods);
//     for s in 0..sample_count {
//         let mut sample_likelihoods = read_likelihoods.sample_matrix(sample_count);
//         let mut genotype_likelihoods = calculator.genotype_likelihoods(&sample_likelihoods, &permutation);
//         let genotype_likelihoods_doubles = genotype_likelihoods.get_likelihoods();
//         assert_eq!(genotype_likelihoods_doubles.len(), genotype_count as usize);
//         for i in 0..test_genotype_count {
//             let mut genotype_allele_counts = calculator.genotype_allele_counts_at(i as usize);
//             let mut read_genotype_likelihoods = vec![0.0; sample_likelihoods.ncols()];
//             for r in 0..sample_likelihoods.ncols() {
//                 let mut components = vec![0.0; genotype_allele_counts.distinct_allele_count()];
//                 for ar in 0..genotype_allele_counts.distinct_allele_count() {
//                     let a = genotype_allele_counts.allele_index_at(ar);
//                     let a_count = genotype_allele_counts.allele_count_at(ar);
//                     let read_lk = sample_likelihoods[[a, r]];
//                     components[ar] = read_lk + (a_count as f64).log10();
//                 }
//                 read_genotype_likelihoods[r] = MathUtils::approximate_log10_sum_log10_vec(components.as_slice(), 0, components.len()) - (ploidy as f64).log10();
//             }
//             let genotype_likelihood = read_genotype_likelihoods.iter().sum::<f64>();
//             assert!(relative_eq!(genotype_likelihoods_doubles[i as usize], genotype_likelihood, epsilon=1e-4));
//         }
//
//     }
// }

fn calculate_genotype_count(ploidy: usize, allele_count: usize) -> usize {
    if ploidy == 0 {
        return 0;
    } else if ploidy == 1 {
        return allele_count;
    } else if ploidy == 2 {
        return (allele_count * (allele_count + 1)) >> 1;
    } else if allele_count == 0 {
        return 0;
    } else {
        return calculate_genotype_count(ploidy - 1, allele_count)
            + calculate_genotype_count(ploidy, allele_count - 1);
    }
}

#[test]
fn ploidy_and_maximum_allele_data() {
    for i in PLOIDY.iter() {
        for j in MAXIMUM_ALLELE.iter() {
            test_ploidy_and_maximum_allele(*i, *j)
        }
    }
}

// #[test]
// fn ploidy_and_maximum_allele_and_read_counts_data() {
//     for i in PLOIDY.iter() {
//         for j in MAXIMUM_ALLELE.iter() {
//             for k in READ_COUNTS.iter() {
//                 test_likelihood_calculation(*i, *j, k)
//             }
//         }
//     }
// }
