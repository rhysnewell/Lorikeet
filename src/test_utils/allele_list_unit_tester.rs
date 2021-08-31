use haplotype::haplotype::Haplotype;
use model::allele_list::AlleleList;
use model::byte_array_allele::ByteArrayAllele;
use rand::{rngs::ThreadRng, Rng};
use std::collections::HashSet;
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
