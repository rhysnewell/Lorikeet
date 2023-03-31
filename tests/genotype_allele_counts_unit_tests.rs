#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::genotype::genotype_allele_counts::GenotypeAlleleCounts;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;




use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::simple_interval::{Locatable};



use std::cmp::{max, Ordering};
use std::collections::{HashSet};


use std::sync::Mutex;

lazy_static! {
    static ref PLOIDY: Vec<usize> = vec![1, 2, 3, 7, 10];
    static ref PLOIDY_WITH_ZERO: Vec<usize> = vec![0, 1, 2, 3, 7, 10];
}

const MAXIMUM_ALLELE_INDEX: usize = 10;

struct GenotypeAlleleCountsUnitTest {
    test_alleles: Vec<ByteArrayAllele>,
}

impl GenotypeAlleleCountsUnitTest {
    fn new() -> Self {
        let mut test_alleles = Vec::new();
        let mut string_builder = String::new();

        string_builder.push('A');
        for i in 0..=50 {
            test_alleles.push(ByteArrayAllele::new(string_builder.as_bytes(), i == 0));
            string_builder.push('A');
        }

        Self { test_alleles }
    }
}

fn test_first(ploidy: usize) {
    let mut subject = GenotypeAlleleCounts::first(ploidy);
    let alleles = GenotypeAlleleCountsUnitTest::new();
    let log10_combination_count = subject.log10_combination_count();
    assert_eq!(subject.ploidy(), ploidy);
    assert_eq!(subject.distinct_allele_count(), 1);
    assert!(
        relative_eq!(log10_combination_count, 0.0, epsilon = 1e-3),
        "Log10 combo count {} right {} +/- {}",
        log10_combination_count,
        0.0,
        1e-3
    );
    assert_eq!(subject.allele_count_at(0), ploidy);
    assert_eq!(subject.allele_count_for(0), ploidy);
    assert_eq!(subject.allele_rank_for(0), 0);
    assert_eq!(subject.allele_rank_for(1), -2);
    assert!(subject.contains_allele(0));
    assert!(!subject.contains_allele(1));
    assert_eq!(subject.allele_index_at(0), 0);
    assert_eq!(subject.maximum_allele_index(), 0);
    assert_eq!(subject.minimum_allele_index(), 0);
    assert_eq!(subject.cmp(&subject), Ordering::Equal);
    assert_eq!(&subject, &subject);
    assert_eq!(subject.index(), 0);
    assert_eq!(
        subject.as_allele_list(&alleles.test_alleles),
        vec![alleles.test_alleles[0].clone(); ploidy]
    );

    for maximum_allele_index in 0..=MAXIMUM_ALLELE_INDEX {
        let mut expected = vec![0; maximum_allele_index + 1];
        expected[0] = ploidy as i32;
        assert_eq!(
            subject.allele_counts_by_index(maximum_allele_index),
            expected
        );
    }
}

fn test_next(ploidy: usize) {
    if ploidy == 0 {
        test_next_zero_ploidy();
    } else if ploidy == 1 {
        test_next_one_ploidy()
    } else {
        test_ploidy_two_or_more(ploidy)
    }
}

fn test_next_zero_ploidy() {
    let alleles = GenotypeAlleleCountsUnitTest::new();
    let mut first = GenotypeAlleleCounts::first(0);
    let next = first.next();
    assert_eq!(&first, &next);
    assert_eq!(next.distinct_allele_count(), 0);
    assert_eq!(next.ploidy(), 0);
    assert_eq!(next.index(), 0);
    assert_eq!(next.as_allele_list(&alleles.test_alleles), Vec::new());

    for maximum_allele_index in 0..=10 {
        let expected = vec![0; maximum_allele_index + 1];
        assert_eq!(next.allele_counts_by_index(maximum_allele_index), expected);
    }
}

fn test_next_one_ploidy() {
    let alleles = GenotypeAlleleCountsUnitTest::new();
    let ploidy2 = GenotypeAlleleCounts::first(2);
    let mut current = GenotypeAlleleCounts::first(1);

    while !current.contains_allele(MAXIMUM_ALLELE_INDEX + 1) {
        let mut next = current.next();
        assert!(
            relative_eq!(next.log10_combination_count(), 0.0, epsilon = 1e-3),
            "Log10 combo count {} right {} +/- {}",
            next.log10_combination_count(),
            0.0,
            1e-3
        );
        assert_eq!(next.minimum_allele_index(), next.maximum_allele_index());
        assert_eq!(
            next.minimum_allele_index(),
            current.minimum_allele_index() + 1
        );
        assert_eq!(next.allele_count_at(0), 1);
        assert_eq!(
            next.allele_index_at(0),
            next.minimum_allele_index() as usize
        );
        assert_eq!(
            next.allele_rank_for(next.minimum_allele_index() as usize),
            0
        );
        assert_eq!(
            next.allele_rank_for(next.minimum_allele_index() as usize + 1),
            -2
        );
        assert_eq!(
            next.allele_count_for(next.minimum_allele_index() as usize),
            1
        );
        assert_eq!(
            next.allele_count_for(next.minimum_allele_index() as usize + 1),
            0
        );
        assert_eq!(next.ploidy(), 1);

        let mut dest = vec![0; next.distinct_allele_count() * 2];
        next.copy_allele_counts(&mut dest, 0);

        assert_eq!(dest, vec![next.index(), 1]);

        assert!(&next > &current);
        assert!(&current < &next);
        assert!(&next == &next);
        assert_ne!(&next, &ploidy2);
        assert_eq!(next.index(), current.index() + 1);
        assert_eq!(next.ploidy(), current.ploidy());
        assert_eq!(
            next.as_allele_list(&alleles.test_alleles),
            vec![alleles.test_alleles[next.maximum_allele_index() as usize].clone()]
        );

        for maximum_allele_index in 0..=MAXIMUM_ALLELE_INDEX {
            let mut expected = vec![0; maximum_allele_index + 1];
            if maximum_allele_index >= current.minimum_allele_index() as usize + 1 {
                expected[current.minimum_allele_index() as usize + 1] = 1;
            };
            assert_eq!(next.allele_counts_by_index(maximum_allele_index), expected);
        }

        current = next;
    }
}

fn test_ploidy_two_or_more(ploidy: usize) {
    if ploidy < 2 {
        panic!("Ploidy has to be two or more >:(");
    };
    let alleles = GenotypeAlleleCountsUnitTest::new();

    let mut current = GenotypeAlleleCounts::first(ploidy);
    loop {
        let mut next = current.next();
        if next.contains_allele(MAXIMUM_ALLELE_INDEX + 1) {
            break;
        };

        if ploidy == 2 {
            assert!(
                relative_eq!(
                    next.log10_combination_count(),
                    if next.distinct_allele_count() == 2 {
                        2.0_f64.log10()
                    } else {
                        0.0
                    },
                    epsilon = 1e-3
                ),
                "log10_combination_count {} is not expected {}",
                next.log10_combination_count(),
                if next.distinct_allele_count() == 2 {
                    2.0_f64.log10()
                } else {
                    0.0
                }
            );
        } else if ploidy == 3 {
            assert!(
                relative_eq!(
                    next.log10_combination_count(),
                    if next.distinct_allele_count() == 3 {
                        6.0_f64.log10()
                    } else if next.distinct_allele_count() == 2 {
                        6.0_f64.log10() - 2.0_f64.log10()
                    } else {
                        0.0
                    },
                    epsilon = 1e-3
                ),
                "log10_combination_count {} is not expected {}",
                next.log10_combination_count(),
                if next.distinct_allele_count() == 3 {
                    6.0_f64.log10()
                } else if next.distinct_allele_count() == 2 {
                    6.0_f64.log10() - 2.0_f64.log10()
                } else {
                    0.0
                }
            )
        } else {
            if next.distinct_allele_count() == 1 {
                assert!(relative_eq!(
                    next.log10_combination_count(),
                    0.0,
                    epsilon = 1e-3
                ));
            } else if next.distinct_allele_count() == ploidy {
                assert!(relative_eq!(
                    next.log10_combination_count(),
                    MathUtils::log10_factorial(ploidy as f64),
                    epsilon = 1e-3
                ))
            }
        }

        let allele_counts_as_list = next.sorted_allele_counts.clone();
        let absent_alleles = Mutex::new(HashSet::new());

        next.for_each_absent_allele_index(
            |allele_index: usize| {
                let mut absent_alleles = absent_alleles.lock().unwrap();
                absent_alleles.insert(allele_index);
            },
            MAXIMUM_ALLELE_INDEX + 1,
        );

        let mut actual_allele_counts = vec![0; next.distinct_allele_count() * 2];
        next.copy_allele_counts(&mut actual_allele_counts, 0);

        assert_eq!(&allele_counts_as_list, &actual_allele_counts);
        let absent_alleles = absent_alleles.lock().unwrap();
        assert_eq!(
            absent_alleles.len(),
            MAXIMUM_ALLELE_INDEX + 1 - next.distinct_allele_count()
        );
        next.for_each_allele_index_and_count(|index, _count| {
            assert!(!absent_alleles.contains(&index))
        });

        if current.distinct_allele_count() == 1 {
            assert_eq!(
                next.maximum_allele_index(),
                current.maximum_allele_index() + 1
            );
            assert_eq!(next.distinct_allele_count(), 2);
            assert_eq!(next.minimum_allele_index(), 0);
        } else {
            assert_eq!(next.maximum_allele_index(), current.maximum_allele_index());
            assert_eq!(
                next.minimum_allele_index(),
                if current.allele_count_at(0) > 1 {
                    0
                } else if current.allele_count_at(0) == 1 {
                    current.minimum_allele_index() + 1
                } else {
                    current.minimum_allele_index()
                }
            );
        }

        // Checking on 0's new count and current.minAllele + 1 alleles.
        assert_eq!(
            next.allele_count_for(0),
            current.allele_count_for(current.minimum_allele_index() as usize) - 1
        );
        assert_eq!(
            next.allele_count_for(current.minimum_allele_index() as usize + 1),
            current.allele_count_for(current.minimum_allele_index() as usize + 1) + 1
        );

        // Checks current.minAllele count
        assert_eq!(
            next.allele_count_for(current.minimum_allele_index() as usize),
            if current.minimum_allele_index() == 0 {
                current.allele_count_at(0) - 1
            } else {
                0
            }
        );

        let mut total_count_sum = 0;
        let mut expected_allele_counts_by_index =
            vec![0; max(MAXIMUM_ALLELE_INDEX as i64, next.maximum_allele_index()) as usize + 1];
        for i in 0..next.distinct_allele_count() {
            let count = next.allele_count_at(i);
            let index = next.allele_index_at(i);
            expected_allele_counts_by_index[index] = count as i32;
            // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
            assert_eq!(next.allele_count_for(index), count);
            total_count_sum += count;
            // Check on counts of, in theory, unaffected allele counts.
            if index as i64 > current.minimum_allele_index() + 1 {
                assert_eq!(
                    next.allele_count_for(index),
                    current.allele_count_for(index)
                )
            };
        }

        assert_eq!(
            next.allele_counts_by_index(max(
                MAXIMUM_ALLELE_INDEX as i64,
                next.maximum_allele_index()
            ) as usize),
            expected_allele_counts_by_index
        );
        assert_eq!(total_count_sum, ploidy);

        assert!(next > current);
        assert!(current < next);
        assert_eq!(&next, &next);
        assert_ne!(&current, &next);
        assert_eq!(next.index(), current.index() + 1);
        assert_eq!(next.ploidy(), ploidy);

        //check as_allele_list
        let mut expected_list = Vec::new();
        for i in 0..next.distinct_allele_count() {
            for _j in 0..next.allele_count_at(i) {
                expected_list.push(alleles.test_alleles[next.allele_index_at(i)].clone());
            }
        }
        assert_eq!(next.as_allele_list(&alleles.test_alleles), expected_list);

        current = next;
    }
}

#[test]
fn ploidy_data() {
    for ploidy in &*PLOIDY {
        test_first(*ploidy)
    }
}

#[test]
fn ploidy_data_with_zero() {
    for i in 0..PLOIDY_WITH_ZERO.len() {
        test_next(PLOIDY_WITH_ZERO[i])
    }
}

#[test]
fn ploidy_data_with_zero_increase() {
    for i in 0..PLOIDY_WITH_ZERO.len() {
        test_increase(PLOIDY_WITH_ZERO[i])
    }
}

fn test_increase(ploidy: usize) {
    if ploidy == 0 {
        test_next_zero_ploidy_increase();
    } else if ploidy == 1 {
        test_next_one_ploidy_increase()
    } else {
        test_ploidy_two_or_more_increase(ploidy)
    }
}

fn test_next_zero_ploidy_increase() {
    let mut first = GenotypeAlleleCounts::first(0);
    let next = first.next();
    let alleles = GenotypeAlleleCountsUnitTest::new();

    assert_eq!(&first, &next);
    assert_eq!(next.distinct_allele_count(), 0);
    assert_eq!(next.ploidy(), 0);
    assert_eq!(next.index(), 0);
    assert_eq!(next.as_allele_list(&alleles.test_alleles), Vec::new());
    for maximum_allele_index in 0..=10 {
        let expected = vec![0; maximum_allele_index + 1];
        assert_eq!(next.allele_counts_by_index(maximum_allele_index), expected);
    }

    first.increase(1);
    assert_eq!(&first, &next);
}

fn test_next_one_ploidy_increase() {
    let mut next = GenotypeAlleleCounts::first(1);
    let _alleles = GenotypeAlleleCountsUnitTest::new();

    while !next.contains_allele(MAXIMUM_ALLELE_INDEX + 1) {
        let current = next.clone();
        next.increase(1);
        assert_eq!(next.minimum_allele_index(), next.maximum_allele_index());
        assert_eq!(
            next.minimum_allele_index(),
            current.minimum_allele_index() + 1
        );
        assert_eq!(next.allele_count_at(0), 1);
        assert_eq!(
            next.allele_index_at(0),
            next.minimum_allele_index() as usize
        );
        assert_eq!(
            next.allele_rank_for(next.minimum_allele_index() as usize),
            0
        );
        assert_eq!(
            next.allele_rank_for(next.minimum_allele_index() as usize + 1),
            -2
        );
        assert_eq!(
            next.allele_count_for(next.minimum_allele_index() as usize),
            1
        );
        assert_eq!(
            next.allele_count_for(next.minimum_allele_index() as usize + 1),
            0
        );
        assert_eq!(next.ploidy(), 1);

        assert!(next > current);
        assert!(current < next);

        assert_eq!(next.index(), current.index() + 1);
        assert_eq!(next.ploidy(), current.ploidy());

        for maximum_allele_index in 0..=MAXIMUM_ALLELE_INDEX {
            let mut expected = vec![0; maximum_allele_index + 1];
            if maximum_allele_index >= current.minimum_allele_index() as usize + 1 {
                expected[current.minimum_allele_index() as usize + 1] = 1;
            }
            assert_eq!(next.allele_counts_by_index(maximum_allele_index), expected);
        }
    }
}

fn test_ploidy_two_or_more_increase(ploidy: usize) {
    assert!(ploidy >= 2);

    let mut next = GenotypeAlleleCounts::first(ploidy);

    while !next.contains_allele(MAXIMUM_ALLELE_INDEX + 1) {
        let current = next.clone();

        next.increase(1);
        if current.distinct_allele_count() == 1 {
            assert_eq!(
                next.maximum_allele_index(),
                current.maximum_allele_index() + 1
            );
            assert_eq!(next.distinct_allele_count(), 2);
            assert_eq!(next.minimum_allele_index(), 0);
        } else {
            assert_eq!(next.maximum_allele_index(), current.maximum_allele_index());
            assert_eq!(
                next.minimum_allele_index(),
                if current.allele_count_at(0) > 1 {
                    0
                } else if current.allele_count_at(0) == 1 {
                    current.minimum_allele_index() + 1
                } else {
                    current.minimum_allele_index()
                }
            );
        };

        // Checking on 0's new count and current.minAllele + 1 alleles.
        assert_eq!(
            next.allele_count_for(0),
            current.allele_count_for(current.minimum_allele_index() as usize) - 1
        );
        assert_eq!(
            next.allele_count_for(current.minimum_allele_index() as usize + 1),
            current.allele_count_for(current.minimum_allele_index() as usize + 1) + 1
        );

        // check current.min_allele count
        assert_eq!(
            next.allele_count_for(current.minimum_allele_index() as usize),
            if current.minimum_allele_index() == 0 {
                current.allele_count_at(0) - 1
            } else {
                0
            }
        );

        let mut total_count_sum = 0;
        let mut expected_allele_counts_by_index =
            vec![0; max(MAXIMUM_ALLELE_INDEX, next.maximum_allele_index() as usize) + 1];
        for i in 0..next.distinct_allele_count() {
            let count = next.allele_count_at(i);
            let index = next.allele_index_at(i);
            expected_allele_counts_by_index[index] = count as i32;
            // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
            assert_eq!(next.allele_count_for(index), count);
            total_count_sum += count;
            // Check on counts of, in theory, unaffected allele counts.
            if index > current.minimum_allele_index() as usize + 1 {
                assert_eq!(
                    next.allele_count_for(index),
                    current.allele_count_for(index)
                );
            }
        }

        assert_eq!(
            next.allele_counts_by_index(max(
                MAXIMUM_ALLELE_INDEX,
                next.maximum_allele_index() as usize
            )),
            expected_allele_counts_by_index
        );
        assert_eq!(total_count_sum, ploidy);
        assert!(next > current);
        assert!(current < next);
        assert!(current != next);
        assert!(next != current);
        assert_eq!(next.index(), current.index() + 1);
        assert_eq!(next.ploidy(), ploidy);
    }
}
