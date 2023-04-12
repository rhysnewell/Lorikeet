#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::math_utils::MathUtils;
use std::cmp::min;

lazy_static! {
    static ref MAXIMUM_ALLELE: Vec<usize> = vec![1, 2, 5, 6];
    // static ref MAXIMUM_ALLELE: Vec<usize> = vec![5];
    static ref PLOIDY: Vec<usize> = vec![1, 2, 3, 20];
    // static ref PLOIDY: Vec<usize> = vec![20];
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

fn test_likelihood_calculation(ploidy: usize, allele_count: usize, read_count: &[usize]) {
    let mut read_likelihoods =
        ReadLikelihoodsUnitTester::read_likelihoods(allele_count, read_count);

    let mut calculator = GenotypeLikelihoodCalculators::get_instance(ploidy, allele_count);
    let genotype_count = calculator.genotype_count;
    let test_genotype_count = min(30000, genotype_count);
    let sample_count = read_count.len();
    // println!("ploidy {} allele count {} read count {:?}", ploidy, allele_count, &read_count);
    let permutation = AlleleLikelihoodMatrixMapper::new(
        read_likelihoods
            .get_allele_list()
            .permutation(read_likelihoods.get_allele_list()),
    );
    for s in 0..sample_count {
        let number_of_evidences = read_likelihoods.sample_evidence_count(s);
        let sample_likelihoods = read_likelihoods.sample_matrix(s);

        let genotype_likelihoods =
            calculator.genotype_likelihoods(&sample_likelihoods, &permutation, number_of_evidences);
        let genotype_likelihoods_doubles = genotype_likelihoods.get_likelihoods();
        assert_eq!(genotype_likelihoods_doubles.len(), genotype_count as usize);
        // println!("sample {}", s);
        for i in 0..test_genotype_count {
            let genotype_allele_counts = calculator.genotype_allele_counts_at(i as usize);
            
            if i == 999 || i == 998 || i == 10 {
                println!("i {} distinct {} genotype count {} index {}", i, genotype_allele_counts.distinct_allele_count(), genotype_count, genotype_allele_counts.index());
            }
            // println!("i {} distinct {} genotype count {}", i, genotype_allele_counts.distinct_allele_count(), genotype_count);
            let mut read_genotype_likelihoods = vec![0.0; number_of_evidences];
            for r in 0..number_of_evidences {
                let mut components = vec![0.0; genotype_allele_counts.distinct_allele_count()];
                for ar in 0..genotype_allele_counts.distinct_allele_count() {
                    let a = genotype_allele_counts.allele_index_at(ar);
                    let a_count = genotype_allele_counts.allele_count_at(ar);
                    
                    let read_lk = sample_likelihoods[[a, r]];
                    if i == 999  {
                        println!("ar {} a {} a_count {} read_lk {}", ar, a, a_count, read_lk);
                    }
                    components[ar] = read_lk + (a_count as f64).log10();
                }
                read_genotype_likelihoods[r] = MathUtils::approximate_log10_sum_log10_vec(
                    components.as_slice(),
                    0,
                    components.len(),
                ) - (ploidy as f64).log10();
            }
            let genotype_likelihood = read_genotype_likelihoods.iter().sum::<f64>();
            if i == 999  {
                println!("read_genotype_likelihoods {:?}", read_genotype_likelihoods);
                println!("Before test expect {} actual {} index {}", genotype_likelihood, genotype_likelihoods_doubles[i as usize], i);
                println!("Ploidy {} allele count {} read count {:?}", ploidy, allele_count, &read_count);
            }
            assert!(
                relative_eq!(
                    genotype_likelihood,
                    genotype_likelihoods_doubles[i as usize],
                    epsilon = 1e-4
                ),
                "Expected {} Actual {} index {}",
                genotype_likelihood,
                genotype_likelihoods_doubles[i as usize],
                i
            );
        }
    }
}

// Don't use genotypeIndexMap in code? Its only used in ReferenceConfidenceVariantContextMerger
// fn test_genotype_index_map(ploidy: usize, old_allele_count: usize, new_allele_count: usize) {
//     let mut rnd = ThreadRng::default();
//     let max_allele_count = max(old_allele_count, new_allele_count);
//     let mut allele_map = Vec::new();
//     let mut reverse_map = HashMap::new();
//
//     for i in 0..new_allele_count {
//         allele_map[i] = rnd.gen_range(0, old_allele_count);
//         let reverse_map_entry = reverse_map.entry(allele_map[i]).or_insert(HashSet::new());
//         reverse_map_entry.insert(i);
//     };
//
//     let mut calculators = GenotypeLikelihoodCalculators::default();
//     let mut calculator = GenotypeLikelihoodCalculators::get_instance(ploidy, max_allele_count);
//
//     let mut
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

#[test]
fn ploidy_and_maximum_allele_and_read_counts_data() {
    for i in PLOIDY.iter() {
        for j in MAXIMUM_ALLELE.iter() {
            for k in READ_COUNTS.iter() {
                test_likelihood_calculation(*i, *j, k)
            }
        }
    }
}
