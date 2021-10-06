use itertools::Itertools;
use rayon::prelude::*;

/**
 * Information about the number of chromosome per sample at a given location.
 *
 */
pub trait PloidyModel {
    /**
     * Return the assumed ploidy for a sample given its index.
     *
     * @param sampleIndex target sample index.
     * @return 0 or greater.
     */
    fn sample_ploidy(&self, sample_index: usize) -> usize;

    /**
     * Checks whether the ploidy is homogeneous across all samples.
     *
     * @return {@code true} if all samples has the same ploidy.
     */
    fn is_homogenous(&self) -> bool;

    /**
     * Sum of all ploidy across all samples.
     * <p>
     *     It must match the sum of all ploidies across samples.
     * </p>
     *
     * @return 0 or greater.
     */
    fn total_ploidy(&self) -> usize;

    fn number_of_samples(&self) -> usize;
}

/**
* {@link PloidyModel} implementation tailored to work with a homogeneous constant ploidy
* across samples and positions.
*/
pub struct HomogeneousPloidyModel {
    pub(crate) sample_list: Vec<String>,
    pub(crate) ploidy: usize,
}

impl HomogeneousPloidyModel {
    pub fn new(sample_list: Vec<String>, ploidy: usize) -> Self {
        Self {
            sample_list,
            ploidy,
        }
    }
}

impl PloidyModel for HomogeneousPloidyModel {
    fn sample_ploidy(&self, sample_index: usize) -> usize {
        self.ploidy
    }

    fn is_homogenous(&self) -> bool {
        true
    }

    fn total_ploidy(&self) -> usize {
        self.ploidy * self.sample_list.len()
    }

    fn number_of_samples(&self) -> usize {
        self.sample_list.len()
    }
}

/**
 * General heterogeneous ploidy model.
 */
pub struct HeterogeneousPloidyModel {
    pub(crate) sample_list: Vec<String>,
    pub(crate) ploidies: Vec<usize>,
    pub(crate) ploidy_sum: usize,
    pub(crate) is_homogenous: bool,
}

impl HeterogeneousPloidyModel {
    pub fn new(sample_list: Vec<String>, ploidies: Vec<usize>) -> Self {
        assert!(
            sample_list.len() == ploidies.len(),
            "Sample list and ploidy array length must match"
        );
        let ploidy_sum = ploidies.par_iter().sum();
        let is_homogenous = ploidies.len() == 0 || ploidies.iter().all_equal();
        Self {
            sample_list,
            ploidies,
            ploidy_sum,
            is_homogenous,
        }
    }
}

impl PloidyModel for HeterogeneousPloidyModel {
    fn sample_ploidy(&self, sample_index: usize) -> usize {
        assert!(
            sample_index <= self.sample_list.len() - 1,
            "Index out of bounds for ploidy list"
        );
        self.ploidies[sample_index]
    }

    fn is_homogenous(&self) -> bool {
        self.is_homogenous
    }

    fn total_ploidy(&self) -> usize {
        self.ploidy_sum
    }

    fn number_of_samples(&self) -> usize {
        self.sample_list.len()
    }
}
