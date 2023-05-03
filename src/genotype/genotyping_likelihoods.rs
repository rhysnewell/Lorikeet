use crate::genotype::genotype_likelihoods::GenotypeLikelihoods;
use crate::haplotype::homogenous_ploidy_model::PloidyModel;
use crate::model::allele_list::AlleleList;
use crate::model::byte_array_allele::Allele;

/**
 * Genotyping Likelihoods collection.
 */
pub struct GenotypingLikelihoods<P: PloidyModel, A: Allele> {
    likelihoods: Vec<GenotypeLikelihoods>,
    ploidy_model: P,
    alleles: AlleleList<A>,
}

impl<P: PloidyModel, A: Allele> GenotypingLikelihoods<P, A> {
    /**
     * Creates a new genotyping-likelihoods collection given the genotype alleles, the sample ploidy model and the
     *   likelihoods themselves.
     * <p>
     * Notice that this constructor does not check whether the likelihood array lengths corresponds to the sample plodies and
     * number of alleles.
     * </p>
     *
     * @param alleles the genotyping alleles.
     * @param ploidyModel the ploidy model.
     * @param likelihoods the actual genotype likelihoods, one element per sample.
     *
     * @throws IllegalArgumentException if any argument is {@code null}, or the number of samples in {@code ploidyModel}
     *  does not correspond with the number of likelihoods arrays in {@code likelihoods}
     */
    pub fn new(
        alleles: AlleleList<A>,
        ploidy_model: P,
        likelihoods: Vec<GenotypeLikelihoods>,
    ) -> GenotypingLikelihoods<P, A> {
        Self {
            alleles,
            ploidy_model,
            likelihoods,
        }
    }

    pub fn number_of_sample(&self) -> usize {
        self.ploidy_model.number_of_samples()
    }

    // pub fn index_of_sample(&self, sample: &String) -> usize {
    //
    // }

    /**
     * Returns the ploidy of the sample given its index in the collection.
     *
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid index for this collection:
     *   [0,{@link #numberOfSamples()).
     *
     * @return 0 or greater.
     */
    pub fn sample_ploidy(&self, sample_index: usize) -> usize {
        self.ploidy_model.sample_ploidy(sample_index)
    }

    /**
     * Returns the genotype-likelihoods of the sample given its index in the collection.
     *
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid index for this collection:
     *   [0,{@link #numberOfSamples()).
     *
     * @return never {@code null}.
     */
    pub fn sample_likelihoods(&self, sample_index: usize) -> &GenotypeLikelihoods {
        &self.likelihoods[sample_index]
    }

    pub fn number_of_alleles(&self) -> usize {
        self.alleles.len()
    }

    pub fn index_of_allele(&self, allele: &A) -> Option<usize> {
        self.alleles.index_of_allele(allele)
    }

    pub fn get_allele(&self, index: usize) -> Option<&A> {
        self.alleles.get_allele(index)
    }
}
