
use std::collections::HashMap;


use crate::model::allele_likelihoods::AlleleLikelihoods;
use crate::model::byte_array_allele::{ByteArrayAllele};
use crate::processing::lorikeet_engine::ReadType;
use crate::reads::bird_tool_reads::BirdToolRead;
use crate::test_utils::allele_list_unit_tester::AlleleListUnitTester;
use crate::utils::artificial_read_utils::ArtificialReadUtils;



/**
 * Constains utilities for tests that need to create read-likelihoods.
 */
pub struct ReadLikelihoodsUnitTester {}

impl ReadLikelihoodsUnitTester {
    pub fn read_likelihoods(
        allele_count: usize,
        read_count: &[usize],
    ) -> AlleleLikelihoods<ByteArrayAllele> {
        let sample_count = read_count.len();
        let allele_list = AlleleListUnitTester::allele_list(allele_count, 100, true).unwrap();
        let sample_list = (0..sample_count).map(|i| format!("SAMPLE_{}", i)).collect();
        let mut sample_to_reads = HashMap::new();
        for i in 0..sample_count {
            sample_to_reads.insert(i, Self::read_list(i, read_count[i]));
        }

        let mut likelihoods =
            AlleleLikelihoods::new_from_allele_list(allele_list, sample_list, sample_to_reads);

        for s in 0..sample_count {
            let sample_likelihoods = &mut likelihoods.values_by_sample_index[s];
            for a in 0..allele_count {
                for r in 0..read_count[s] {
                    sample_likelihoods[[a, r]] = Self::test_likelihood(s, a, r);
                }
            }
        }

        return likelihoods;
    }

    fn test_likelihood(sample_index: usize, allele_index: usize, read_index: usize) -> f64 {
        return -(((3 * (sample_index + 1) + 7 * (allele_index + 1) + 11 * (read_index + 1))
            as f64)
            .abs());
    }

    fn read_list(sample_index: usize, read_count: usize) -> Vec<BirdToolRead> {
        let reads = (0..read_count)
            .into_iter()
            .map(|read_index| {
                BirdToolRead::new(
                    ArtificialReadUtils::create_artificial_read_default(
                        &format!("READ_{}_{}", sample_index, read_index),
                        1,
                        1,
                        100,
                        false,
                    ),
                    sample_index,
                    ReadType::Short,
                )
            })
            .collect();

        return reads;
    }
}
