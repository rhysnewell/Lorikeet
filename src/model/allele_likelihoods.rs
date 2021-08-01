use model::byte_array_allele::ByteArrayAllele;
use model::variants::Allele;
use ndarray::{Array2, Array3};
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use std::collections::HashMap;
use std::hash::Hash;
use utils::math_utils::MathUtils;
use utils::simple_interval::{Locatable, SimpleInterval};

lazy_static! {
    pub static ref LOG_10_INFORMATIVE_THRESHOLD: f64 = 0.2;
    pub static ref NATURAL_LOG_INFORMATIVE_THRESHOLD: f64 =
        MathUtils::log10_to_log(*LOG_10_INFORMATIVE_THRESHOLD);
}

/**
 * Evidence-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
pub struct AlleleLikelihoods {
    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    pub(crate) evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    /**
     * Evidence disqualified by .
     */
    pub(crate) filtered_evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    /**
     * Indexed per sample, allele and finally evidence (within sample).
     * <p>
     *     valuesBySampleIndex[s][a][r] == lnLk(R_r | A_a) where R_r comes from Sample s.
     * </p>
     */
    pub(crate) values_by_sample_index: Vec<Array2<f64>>,
    /**
     * Keeps track of the maximum number of evidences and likelihood values that can be stored
     * stored across all alleles.
     */
    pub(crate) likelihoods_matrix_evidence_capacity_by_sample_index: Vec<usize>,
    /**
     * Holds the number of evidence per sample.
     */
    pub(crate) number_of_evidences: Vec<usize>,
    /**
     * Sample list.
     */
    pub(crate) samples: Vec<String>,
    /**
     * Allele list.
     */
    pub(crate) alleles: Vec<ByteArrayAllele>,
    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_INDEX}.
     */
    pub(crate) reference_allele_index: Option<usize>,
    // /**
    //  * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
    //  */
    // pub(crate) sample_matrices: Vec<Array2<f64>>,
    pub(crate) is_natural_log: bool,
}

impl AlleleLikelihoods {
    pub fn new(
        alleles: Vec<ByteArrayAllele>,
        samples: Vec<String>,
        evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    ) -> AlleleLikelihoods {
        let sample_count = samples.len();
        let allele_count = alleles.len();

        let allele_count = alleles.len();
        let mut values_by_sample_index = vec![Array2::<f64>::zeros((0, 0)); sample_count];
        let mut likelihoods_matrix_evidence_capacity_by_sample_index = vec![0; sample_count];
        let reference_allele_index = Self::find_reference_allele_index(&alleles);
        let mut number_of_evidences = vec![0; sample_count];
        Self::setup_indexes(
            &evidence_by_sample_index,
            sample_count,
            allele_count,
            &mut number_of_evidences,
            &mut likelihoods_matrix_evidence_capacity_by_sample_index,
            &mut values_by_sample_index,
        );

        // let sample_matrices = vec![Array2::<f64>::zeros((0, 0)); sample_count];

        Self {
            evidence_by_sample_index,
            filtered_evidence_by_sample_index: HashMap::new(),
            values_by_sample_index,
            likelihoods_matrix_evidence_capacity_by_sample_index,
            number_of_evidences,
            samples,
            alleles,
            reference_allele_index,
            is_natural_log: false,
            // sample_matrices
        }
    }

    fn new_private(
        alleles: Vec<ByteArrayAllele>,
        samples: Vec<String>,
        evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
        filtered_evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
        values: Vec<Vec<Vec<f64>>>,
    ) {
        let sample_count = samples.len();
    }

    fn get_informative_threshold(&self) -> f64 {
        if self.is_natural_log {
            return *NATURAL_LOG_INFORMATIVE_THRESHOLD;
        } else {
            return *LOG_10_INFORMATIVE_THRESHOLD;
        }
    }

    fn find_reference_allele_index(alleles: &Vec<ByteArrayAllele>) -> Option<usize> {
        alleles.par_iter().position_first(|a| a.is_ref)
    }

    fn setup_indexes(
        evidence_by_sample: &HashMap<usize, Vec<BirdToolRead>>,
        sample_count: usize,
        allele_count: usize,
        number_of_evidences: &mut Vec<usize>,
        likelihoods_matrix_evidence_by_sample_index: &mut Vec<usize>,
        values_by_sample_index: &mut Vec<Array2<f64>>,
    ) {
        for s in 0..sample_count {
            let sample_evidences = evidence_by_sample.get(&s);
            let sample_evidence_count = match sample_evidences {
                None => 0,
                Some(sample_evidences) => sample_evidences.len(),
            };
            number_of_evidences[s] = sample_evidence_count;

            let sample_values = vec![vec![0.0; sample_evidence_count]; allele_count];
            let sample_values = Array2::<f64>::zeros((allele_count, sample_evidence_count));
            likelihoods_matrix_evidence_by_sample_index[s] = sample_evidence_count;
            values_by_sample_index[s] = sample_values;
        }
    }

    /**
     * Adjusts likelihoods so that for each unit of evidence, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each unit of evidence based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    pub fn normalize_likelihoods(
        &mut self,
        maximum_likelihood_difference_cap: f64,
        symmetrically_normalize_alleles_to_reference: bool,
    ) {
        if maximum_likelihood_difference_cap != f64::NEG_INFINITY {
            let allele_count = self.alleles.len();
            if !(allele_count == 0 || allele_count == 1) {
                for s in 0..self.values_by_sample_index.len() {
                    let mut sample_values = &mut self.values_by_sample_index[s];
                    let evidence_count = self.evidence_by_sample_index.get(&s).unwrap().len();
                    for r in 0..evidence_count {}
                }
            }
        }
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    fn normalize_likelihoods_per_evidence(
        &mut self,
        maximum_best_alt_likelihood_difference: f64,
        sample_values: &mut Array2<f64>,
        sample_index: usize,
        evidence_index: usize,
        symmetrically_normalize_alleles_to_reference: bool,
    ) {
        //allow the best allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        let best_allele = self.search_best_allele(
            sample_values,
            sample_index,
            evidence_index,
            symmetrically_normalize_alleles_to_reference,
            &None,
        );

        let worst_likelihood_cap = best_allele.likelihood + maximum_best_alt_likelihood_difference;

        let allele_count = self.alleles.len();

        sample_values
            .column_mut(evidence_index)
            .par_mapv_inplace(|value| {
                if value < worst_likelihood_cap {
                    worst_likelihood_cap
                } else {
                    value
                }
            });
    }

    /**
     * Search the best allele for a unit of evidence.
     *
     * @param sampleIndex including sample index.
     * @param evidenceIndex  target evidence index.
     *
     * @param priorities An array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     *                   uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score.
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    fn search_best_allele(
        &self,
        sample_values: &mut Array2<f64>,
        sample_index: usize,
        evidence_index: usize,
        can_be_reference: bool,
        priorities: &Option<Vec<f64>>,
    ) -> BestAllele {
        let allele_count = self.alleles.len();
        if allele_count == 0
            || (allele_count == 1 && self.reference_allele_index.unwrap() == 0 && !can_be_reference)
        {
            return BestAllele::new(
                sample_index,
                evidence_index,
                None,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            );
        };

        let mut best_allele_index = if can_be_reference || self.reference_allele_index.unwrap() != 0
        {
            0
        } else {
            1
        };

        let mut second_best_index = 0;
        let mut best_likelihood = sample_values[[best_allele_index, evidence_index]];
        let mut second_best_likelihood = f64::NEG_INFINITY;

        for a in (best_allele_index + 1)..allele_count {
            if !can_be_reference && self.reference_allele_index.unwrap() == a {
                continue;
            };

            let candidate_likelihood = sample_values[[a, evidence_index]];
            if candidate_likelihood > best_likelihood {
                second_best_index = best_allele_index;
                best_allele_index = a;
                second_best_likelihood = best_likelihood;
                best_likelihood = candidate_likelihood;
            } else if candidate_likelihood > second_best_likelihood {
                second_best_index = a;
                second_best_likelihood = candidate_likelihood;
            }
        }

        match priorities {
            None => {
                // pass
            }
            Some(priorities) => {
                if (best_likelihood - second_best_likelihood) < self.get_informative_threshold() {
                    let mut best_priority = priorities[best_allele_index];
                    let mut second_best_priority = priorities[second_best_index];
                    for a in 0..allele_count {
                        let candidate_likelihood = sample_values[[a, evidence_index]];
                        if a == best_allele_index
                            || (!can_be_reference && a == self.reference_allele_index.unwrap())
                            || (best_likelihood - candidate_likelihood)
                                > self.get_informative_threshold()
                        {
                            continue;
                        };

                        let candidate_priority = priorities[a];
                        if candidate_priority > best_priority {
                            second_best_index = best_allele_index;
                            best_allele_index = a;
                            second_best_priority = best_priority;
                            best_priority = candidate_priority;
                        } else if candidate_priority > second_best_priority {
                            second_best_index = a;
                            second_best_priority = candidate_priority;
                        }
                    }
                }
            }
        };

        best_likelihood = sample_values[[best_allele_index, evidence_index]];
        second_best_likelihood = if second_best_index != best_allele_index {
            sample_values[[second_best_index, evidence_index]]
        } else {
            f64::NEG_INFINITY
        };

        return BestAllele::new(
            sample_index,
            evidence_index,
            Some(best_allele_index),
            best_likelihood,
            second_best_likelihood,
        );
    }

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param log10MinTrueLikelihood Function that returns the minimum likelihood that the best allele for a unit of evidence must have
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    pub fn filter_poorly_modeled_evidence(
        &mut self,
        log10_min_true_likelihood: Box<dyn Fn(&BirdToolRead) -> f64>,
    ) {
        let number_of_samples = self.samples.len();
        for sample_index in 0..number_of_samples {
            let mut indexes_to_remove;
            {
                let sample_evidence = self.evidence_by_sample_index.get(&sample_index).unwrap();
                let number_of_evidence = sample_evidence.len();
                indexes_to_remove = (0..number_of_evidence)
                    .into_iter()
                    .filter(|i| {
                        self.maximum_likelihood_over_all_alleles(sample_index, *i)
                            < (log10_min_true_likelihood)(&sample_evidence[*i])
                    })
                    .collect::<Vec<usize>>();
            }
            self.remove_evidence_by_index(sample_index, indexes_to_remove)
        }
    }

    // remove evidence and unset the {@code evidenceIndexBySampleIndex} cache for this sample
    // assumes that evidencesToRemove is sorted and without duplicates.
    fn remove_evidence_by_index(&mut self, sample_index: usize, evidences_to_remove: Vec<usize>) {
        // Retain the filtered evidence for later genotyping purposes
        let filtered = self
            .filtered_evidence_by_sample_index
            .entry(sample_index)
            .or_insert(Vec::new());
        let num_to_remove = evidences_to_remove.len();
        if num_to_remove > 0 {
            let old_evidence_count = self.number_of_evidences[sample_index];
            let new_evidence_count = old_evidence_count - num_to_remove;

            // update the list of evidence and evidence count
            let old_evidence = self.evidence_by_sample_index.entry(0).or_insert(Vec::new());
            let mut num_removed = 0;
            for n in 0..old_evidence_count {
                if num_removed < num_to_remove && n == evidences_to_remove[num_removed] {
                    num_removed += 1;
                } else {
                    filtered.push(old_evidence.remove(n));

                    // update the likelihoods arrays in place
                    let mut sample_values = &mut self.values_by_sample_index[sample_index];
                    for mut allele_values in sample_values.rows_mut() {
                        allele_values[n - num_removed] = allele_values[n]
                    }
                }
            }

            // set to NaN lks of the deleted positions in lk value arrays.
            let mut sample_values = &mut self.values_by_sample_index[sample_index];
            for mut allele_values in sample_values.rows_mut() {
                allele_values
                    .slice_mut(s![new_evidence_count..])
                    .par_mapv_inplace(|v| f64::NAN);
            }

            self.number_of_evidences[sample_index] = new_evidence_count;
            // invalidate evidence-index_by_sample_index occurs here. We don't have this struct field
            // Unsure if we need yet
            // TODO: Make sure this functions
        }
    }

    // The evidenceToIndex map becomes invalid when the evidence list is modified, for example by deleting evidence
    // When adding evidence it is simple enough to add new entries to the map, but we must be careful to do so.
    // fn invalidate_evidence_to_index_cache(&mut self, sample_index) {
    //     self.evide
    // }

    fn maximum_likelihood_over_all_alleles(
        &self,
        sample_index: usize,
        evidence_index: usize,
    ) -> f64 {
        let mut result = f64::NEG_INFINITY;
        let allele_count = self.alleles.len();
        let sample_values = &self.values_by_sample_index[sample_index];
        for a in 0..allele_count {
            if sample_values[[a, evidence_index]] > result {
                result = sample_values[[a, evidence_index]];
            };
        }

        return result;
    }
}

// pub struct LikelihoodMatrix<L: Locatable> {
//     pub(crate) evidence: Vec<L>,
//     pub(crate) alleles: Vec<Allele>,
//     pub(crate) values: Array2<f64>,
// }

pub struct BestAllele {
    allele_index: Option<usize>,
    /**
     * The containing sample.
     */
    sample_index: usize,
    /**
     * The query evidence.
     */
    evidence_index: usize,
    /**
     * If allele != null, the indicates the likelihood of the evidence.
     */
    likelihood: f64,
    /**
     * Confidence that the evidence actually was generated under that likelihood.
     * This is equal to the difference between this and the second best allele match.
     */
    confidence: f64,
}

impl BestAllele {
    pub fn new(
        sample_index: usize,
        evidence_index: usize,
        best_allele_index: Option<usize>,
        likelihood: f64,
        second_best_likelihood: f64,
    ) -> Self {
        let confidence = if likelihood == second_best_likelihood {
            0.0
        } else {
            likelihood - second_best_likelihood
        };
        Self {
            allele_index: best_allele_index,
            sample_index,
            evidence_index,
            likelihood,
            confidence,
        }
    }

    pub fn is_informative(&self) -> bool {
        return self.confidence > *LOG_10_INFORMATIVE_THRESHOLD;
    }
}
