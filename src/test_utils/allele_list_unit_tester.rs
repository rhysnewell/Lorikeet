use haplotype::haplotype::Haplotype;
use hashlink::LinkedHashSet;
use model::allele_list::AlleleList;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use rand::{rngs::ThreadRng, Rng};
use std::collections::HashSet;
use std::hash::Hash;
use test_utils::random_dna::RandomDNA;
use utils::errors::BirdToolError;

/**
 * Helper class for those unit-test classes that test on implementations of SampleList.
 *
 * @author Rhys Newell rhys.newell at hdr.qut.edu.au;
 */

pub struct AlleleListUnitTester {
    random: ThreadRng,
    random_dna: RandomDNA,
}

impl AlleleListUnitTester {
    pub fn new() -> Self {
        Self {
            random: ThreadRng::default(),
            random_dna: RandomDNA::new(ThreadRng::default()),
        }
    }

    /**
     * Test that the contents of an allele-list are the ones expected.
     * <p/>
     * <p>
     * This method perform various consistency check involving all the {@link org.broadinstitute.hellbender.utils.genotyper.AlleleList} interface methods.
     * Therefore calling this method is equivalent to a thorough check of the {@link org.broadinstitute.hellbender.utils.genotyper.AlleleList} aspect of
     * the {@code actual} argument.
     * </p>
     *
     * @param actual   the sample-list to assess.
     * @param expected the expected sample-list.
     * @throws IllegalArgumentException if {@code expected} is {@code null} or contains
     *                                  {@code null}s which is an indication of an bug in the testing code.
     * @throws RuntimeException         if there is some testing assertion exception which
     *                                  is an indication of an actual bug the code that is been tested.
     */
    pub fn assert_allele_list<A: Allele + Send + Sync + Hash>(
        actual: &AlleleList<A>,
        expected: Vec<&A>,
    ) {
        let mut expected_allele_set = LinkedHashSet::new();
        assert_eq!(actual.number_of_alleles(), expected.len());
        for i in 0..expected.len() {
            let expected_allele = expected[i];
            let actual_allele = match actual.get_allele(i) {
                Some(actual_allele) => actual_allele,
                None => panic!("Can't retrieve actual allele {}", i)
            };
            assert!(!expected_allele_set.contains(actual_allele));
            assert_eq!(actual_allele, expected_allele);
            assert_eq!(
                actual.index_of_allele(actual_allele).unwrap(),
                i,
                "allele index mismatch"
            );
            expected_allele_set.insert(actual_allele.clone());
        }
    }

    /**
     * Generate testing alleles. Guarantees that alleles are unique.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    pub fn generate_random_unique_alleles(
        &mut self,
        allele_count: usize,
        max_allele_length: usize,
    ) -> Vec<ByteArrayAllele> {
        let mut set = HashSet::new();
        while set.len() < allele_count {
            let allele_length = self.random.gen_range(0, max_allele_length) + 1;
            let allele = self.random_dna.next_bases(allele_length);
            let result = ByteArrayAllele::new(allele.as_bytes(), false);
            set.insert(result);
        }

        return set.into_iter().collect::<Vec<ByteArrayAllele>>();
    }

    /**
     * Generate testing alleles.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    pub fn generate_random_alleles(
        &mut self,
        allele_count: usize,
        max_allele_length: usize,
    ) -> Vec<ByteArrayAllele> {
        let mut set = Vec::new();
        while set.len() < allele_count {
            let allele_length = self.random.gen_range(0, max_allele_length) + 1;
            let allele = self.random_dna.next_bases(allele_length);
            let result = ByteArrayAllele::new(allele.as_bytes(), false);
            set.push(result);
        }

        return set;
    }

    /**
     * Generate testing alleles.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     * @param skipIfRepeats throw an test-skip exception {@link SkipException} if the resulting allele-list
     *                     has repeats, thus is size is less than {@code alleleCount}
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    pub fn allele_list(
        allele_count: usize,
        max_allele_length: usize,
        skip_if_repeats: bool,
    ) -> Result<AlleleList<ByteArrayAllele>, BirdToolError> {
        let mut alleles =
            Self::new().generate_random_unique_alleles(allele_count, max_allele_length);
        if allele_count > 0 {
            alleles[0] = ByteArrayAllele::new(&alleles[0].bases[..], true);
        }

        let allele_list = AlleleList::new_from_vec(alleles);

        if skip_if_repeats && allele_list.number_of_alleles() != allele_list.len() {
            return Err(BirdToolError::SkipException(format!(
                "Repeated alleles, should be infrequent"
            )));
        };

        return Ok(allele_list);
    }
}
