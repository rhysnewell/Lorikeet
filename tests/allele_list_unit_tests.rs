#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;

use lorikeet_genome::model::allele_list::{AlleleList, AlleleListPermutation};
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::test_utils::allele_list_unit_tester::AlleleListUnitTester;
use rand::prelude::ThreadRng;
use rand::Rng;

lazy_static! {
    static ref ALLELE_COUNT: Vec<usize> = vec![0, 1, 5, 10, 20];
    static ref MAX_ALLELE_LENGTH: Vec<usize> = vec![1, 2, 3, 10];
}

#[test]
fn test_empty_list() {
    let al: AlleleList<ByteArrayAllele> = AlleleList::empty_allele_list();
    assert_eq!(al.number_of_alleles(), 0);
}

#[test]
fn single_allele_list_data() {
    let mut allele_lists = vec![Vec::new(); ALLELE_COUNT.len() * MAX_ALLELE_LENGTH.len()];
    let mut allele_unit_tester = AlleleListUnitTester::new();
    let mut next_index = 0;
    for i in 0..ALLELE_COUNT.len() {
        for j in 0..MAX_ALLELE_LENGTH.len() {
            allele_lists[next_index] =
                allele_unit_tester.generate_random_alleles(ALLELE_COUNT[i], MAX_ALLELE_LENGTH[j]);
            next_index += 1;
        }
    }

    for list in allele_lists.into_iter() {
        test_self_permutation(&list);
        test_subset_permutation(&list);
        test_self_permutation(&list);
    }
}

// #[test]
// fn two_allele_list_data() {
//     let mut allele_lists = vec![Vec::new(); ALLELE_COUNT.len() * MAX_ALLELE_LENGTH.len()];
//     let mut allele_unit_tester = AlleleListUnitTester::new();
//     let mut next_index = 0;
//     for i in 0..ALLELE_COUNT.len() {
//         for j in 0..MAX_ALLELE_LENGTH.len() {
//             allele_lists[next_index] =
//                 allele_unit_tester.generate_random_alleles(ALLELE_COUNT[i], MAX_ALLELE_LENGTH[j]);
//             next_index += 1;
//         }
//     }

//     for i in 0..allele_lists.len() {
//         for j in 0..allele_lists.len() {
//             test_equals(&allele_lists[i], &allele_lists[j])
//         }
//     }
// }

// fn test_index_of_reference(alleles1: &Vec<ByteArrayAllele>) {
//
// }

fn test_self_permutation(alleles1: &Vec<ByteArrayAllele>) {
    let original_allele_list = AlleleList::new(alleles1);
    let self_permutation = original_allele_list
        .clone()
        .permutation(original_allele_list.clone());
    assert_eq!(
        self_permutation.from_size(),
        original_allele_list.number_of_alleles()
    );
    assert_eq!(
        self_permutation.to_size(),
        original_allele_list.number_of_alleles()
    );
    assert!(self_permutation.is_non_permuted());
    assert!(!self_permutation.is_partial());
    for i in 0..original_allele_list.number_of_alleles() {
        assert_eq!(
            self_permutation.get_allele(i),
            original_allele_list.get_allele(i).unwrap()
        );
        assert_eq!(self_permutation.from_index(i), i);
        assert_eq!(self_permutation.to_index(i).unwrap(), i);
        assert_eq!(self_permutation.from_list(), self_permutation.to_list());
        AlleleListUnitTester::assert_allele_list(
            &original_allele_list,
            self_permutation.from_list(),
        );
    }
    // assert!(AlleleList::equals(self_permutation.))
}

fn test_subset_permutation(alleles1: &Vec<ByteArrayAllele>) {
    let mut rnd = ThreadRng::default();
    let mut random_subset_alleles = Vec::new();
    for allele in alleles1.iter() {
        if rnd.gen_bool(0.5) {
            random_subset_alleles.push(allele.clone());
        };
    }

    let original_allele_list = AlleleList::new(alleles1);
    let target_random_allele_list = AlleleList::new_from_vec(random_subset_alleles);
    let subset = original_allele_list
        .clone()
        .permutation(target_random_allele_list.clone());
    if original_allele_list.number_of_alleles() == target_random_allele_list.number_of_alleles() {
        // stop because this is invalid for this test
    } else {
        assert!(subset.is_partial());
        assert!(!subset.is_non_permuted());
        assert_eq!(subset.from_size(), original_allele_list.number_of_alleles());
        assert_eq!(
            subset.to_size(),
            target_random_allele_list.number_of_alleles()
        );
        AlleleListUnitTester::assert_allele_list(&original_allele_list, subset.from_list());
        AlleleListUnitTester::assert_allele_list(&target_random_allele_list, subset.to_list());

        for i in 0..target_random_allele_list.number_of_alleles() {
            assert_eq!(
                Some(subset.from_index(i)),
                original_allele_list
                    .index_of_allele(target_random_allele_list.get_allele(i).unwrap())
            );
        }

        for j in 0..original_allele_list.number_of_alleles() {
            let allele = original_allele_list.get_allele(j);
            assert_eq!(
                subset.to_index(j),
                target_random_allele_list.index_of_allele(allele.unwrap())
            )
        }
    }
}

// fn test_shuffle_permutation(alleles1: &Vec<ByteArrayAllele>) {
//     let mut rnd = ThreadRng::default();
//     let original_allele_list = AlleleList::new(alleles1);
//     if original_allele_list.number_of_alleles() <= 1 {
//         // can't do anything so ignore
//     } else {
//         let mut target_allele_array = original_allele_list
//             .as_list_of_alleles()
//             .into_iter()
//             .cloned()
//             .collect::<Vec<ByteArrayAllele>>();
//         let mut from_index = vec![0; target_allele_array.len()];
//         for i in 0..from_index.len() {
//             from_index[i] = i;
//         }

//         for i in 0..target_allele_array.len() - 1 {
//             let swap_index = rnd.gen_range(0, target_allele_array.len() - i - 1);
//             let other_index = from_index[swap_index + i + 1];
//             let other = target_allele_array[swap_index + i + 1].clone();
//             from_index[swap_index + i + 1] = from_index[i];
//             from_index[i] = other_index;
//             target_allele_array[swap_index + i + 1] = target_allele_array[i].clone();
//             target_allele_array[i] = other;
//         }

//         let target_allele_list = AlleleList::new(&target_allele_array);

//         let permutation = original_allele_list
//             .clone()
//             .permutation(target_allele_list.clone());
//         assert!(!permutation.is_non_permuted());
//         AlleleListUnitTester::assert_allele_list(&original_allele_list, permutation.from_list());
//         AlleleListUnitTester::assert_allele_list(&target_allele_list, permutation.to_list());
//         assert!(!permutation.is_partial());
//         assert_eq!(
//             permutation.from_size(),
//             original_allele_list.number_of_alleles()
//         );
//         assert_eq!(
//             permutation.to_size(),
//             target_allele_list.number_of_alleles()
//         );
//         for i in 0..permutation.from_size() {
//             assert_eq!(
//                 permutation.to_index(i),
//                 target_allele_list.index_of_allele(original_allele_list.get_allele(i).unwrap())
//             );
//             assert_eq!(
//                 Some(permutation.from_index(i)),
//                 original_allele_list.index_of_allele(target_allele_list.get_allele(i).unwrap())
//             );
//             assert_eq!(permutation.from_index(i), from_index[i]);
//         }
//     }
// }

// fn test_equals(alleles1: &Vec<ByteArrayAllele>, alleles2: &Vec<ByteArrayAllele>) {
//     println!("test_equals: {:?} {:?}", alleles1, alleles2);
//     let allele_list1 = AlleleList::new(alleles1);
//     let allele_list2 = AlleleList::new(alleles2);

//     assert!(allele_list1.eq(&allele_list2));
//     assert!(allele_list2.eq(&allele_list1));

//     assert!(allele_list1.as_list_of_alleles() == allele_list2.as_list_of_alleles())
// }
