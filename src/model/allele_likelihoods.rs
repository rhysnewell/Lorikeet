use assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use haplotype::haplotype::Haplotype;
use hashlink::linked_hash_map::LinkedHashMap;
use model::allele_list::AlleleList;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variants::NON_REF_ALLELE;
use ndarray::{Array2, ArrayBase, ArrayView, Axis};
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use statrs::statistics::Median;
use std::cmp::max;
use std::collections::HashMap;
use std::f64::NAN;
use std::hash::Hash;
use test_utils::allele_list_unit_tester::AlleleListUnitTester;
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
#[derive(Debug, Clone)]
pub struct AlleleLikelihoods<A: Allele> {
    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    pub(crate) evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    /**
     * Evidence disqualified by .
     */
    pub filtered_evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
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
    pub(crate) alleles: AlleleList<A>,
    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_INDEX}.
     */
    pub(crate) reference_allele_index: Option<usize>,
    // /**
    //  * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
    //  */
    // pub(crate) sample_matrices: Vec<Array2<f64>>,
    pub(crate) is_natural_log: bool,

    pub(crate) subsetted_genomic_loc: Option<SimpleInterval>,
}

impl<A: Allele> AlleleLikelihoods<A> {
    pub fn new(
        alleles: Vec<A>,
        samples: Vec<String>,
        evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    ) -> AlleleLikelihoods<A> {
        let allele_list = AlleleList::new_from_vec(alleles);
        Self::new_from_allele_list(allele_list, samples, evidence_by_sample_index)
    }

    pub fn new_from_allele_list(
        alleles: AlleleList<A>,
        samples: Vec<String>,
        evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    ) -> AlleleLikelihoods<A> {
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
            subsetted_genomic_loc: None,
            // sample_matrices
        }
    }

    pub fn new_from_likelihoods(
        alleles: AlleleList<A>,
        samples: Vec<String>,
        evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
        filtered_evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
        values_by_sample_index: Vec<Array2<f64>>,
    ) -> Self {
        let sample_count = samples.len();
        let number_of_evidences = (0..sample_count)
            .into_iter()
            .map(|i| evidence_by_sample_index.get(&i).unwrap().len())
            .collect::<Vec<usize>>();
        let likelihoods_matrix_evidence_capacity_by_sample_index = values_by_sample_index
            .iter()
            .map(|sample_values| {
                sample_values
                    .axis_iter(Axis(1))
                    .map(|row| row.into_par_iter().filter(|x| !x.is_nan()).count())
                    .min()
                    .unwrap_or(0)
            })
            .collect::<Vec<usize>>();

        let reference_allele_index = Self::find_reference_allele_index(&alleles);

        Self {
            evidence_by_sample_index,
            filtered_evidence_by_sample_index,
            values_by_sample_index,
            likelihoods_matrix_evidence_capacity_by_sample_index,
            number_of_evidences,
            samples,
            alleles,
            reference_allele_index,
            is_natural_log: false,
            subsetted_genomic_loc: None,
        }
    }

    pub fn samples(&self) -> &Vec<String> {
        &self.samples
    }

    pub fn alleles(&self) -> &AlleleList<A> {
        &self.alleles
    }

    /**
     * Returns sample name given its index.
     *
     * @param sampleIndex query index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is negative.
     *
     * @return never {@code null}.
     */
    pub fn index_of_sample(&self, sample: &String) -> usize {
        self.samples.iter().position(|s| s == sample).unwrap()
    }

    /**
     * Returns the index of an allele within the likelihood collection.
     *
     * @param allele the query allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return usize
     */
    pub fn index_of_allele(&self, allele: &A) -> Option<usize> {
        self.alleles.index_of_allele(allele)
    }

    pub fn number_of_samples(&self) -> usize {
        self.samples.len()
    }

    pub fn number_of_alleles(&self) -> usize {
        self.alleles.len()
    }

    pub fn get_allele_list(&self) -> AlleleList<A> {
        self.alleles.clone()
    }

    pub fn get_allele_list_byte_array(&self) -> AlleleList<ByteArrayAllele> {
        AlleleList::new_from_vec(
            self.alleles
                .list
                .iter()
                .map(|a| ByteArrayAllele::new(a.get_bases(), a.is_reference()))
                .collect::<Vec<ByteArrayAllele>>(),
        )
    }

    /**
     * Returns the quantity of evidence that belongs to a sample in the evidence-likelihood collection.
     * @param sampleIndex the query sample index.
     *
     *
     * @return 0 or greater.
     */
    pub fn sample_evidence_count(&self, sample_index: usize) -> usize {
        match self.evidence_by_sample_index.get(&sample_index) {
            Some(reads) => return reads.len(),
            None => 0,
        }
    }

    /**
     * Returns the units of evidence that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return None or perhaps a zero-length array if there is no evidence in sample. No element in
     *   the array will be null.
     */
    pub fn sample_evidence(&self, sample_index: usize) -> Option<&Vec<BirdToolRead>> {
        self.evidence_by_sample_index.get(&sample_index)
    }

    /**
     * Returns the index of evidence within a sample evidence-likelihood sub collection.
     * @param sampleIndex the sample index.
     * @param evidence the query evidence.
     * @return None if there is no such evidence in that sample, 0 or greater otherwise.
     */
    pub fn evidence_index(&self, sample_index: usize, evidence: &BirdToolRead) -> Option<usize> {
        let index = self.evidence_by_sample_index.get(&sample_index).unwrap();
        return index.iter().position(|read| read == evidence);
    }

    pub fn set(
        &mut self,
        sample_index: usize,
        allele_index: usize,
        evidence_index: usize,
        value: f64,
    ) {
        self.values_by_sample_index[sample_index][[allele_index, evidence_index]] = value
    }
    // fn new_private(
    //     alleles: Vec<ByteArrayAllele>,
    //     samples: Vec<String>,
    //     evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    //     filtered_evidence_by_sample_index: HashMap<usize, Vec<BirdToolRead>>,
    //     values: Vec<Vec<Vec<f64>>>,
    // ) {
    //     let sample_count = samples.len();
    // }

    fn get_informative_threshold(&self) -> f64 {
        if self.is_natural_log {
            return *NATURAL_LOG_INFORMATIVE_THRESHOLD;
        } else {
            return *LOG_10_INFORMATIVE_THRESHOLD;
        }
    }

    fn find_reference_allele_index(alleles: &AlleleList<A>) -> Option<usize> {
        alleles
            .as_list_of_alleles()
            .par_iter()
            .position_first(|a| a.is_reference())
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

    pub fn set_variant_calling_subset_used(&mut self, loc: &SimpleInterval) {
        self.subsetted_genomic_loc = Some(loc.clone())
    }

    /**
     * Returns the location used for subsetting. May be null.
     */
    pub fn get_variant_calling_subset_applied(&self) -> &Option<SimpleInterval> {
        &self.subsetted_genomic_loc
    }

    pub fn evidence_count(&self) -> usize {
        self.evidence_by_sample_index
            .values()
            .flat_map(|reads| reads.iter())
            .count()
    }

    pub fn sample_matrix(&mut self, sample_index: usize) -> &mut Array2<f64> {
        &mut self.values_by_sample_index[sample_index]
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
        debug!(
            "about to normalize with max cap {}",
            maximum_likelihood_difference_cap
        );
        if maximum_likelihood_difference_cap != f64::NEG_INFINITY {
            let allele_count = self.alleles.len();
            debug!("Allele count {}", allele_count);
            if !(allele_count == 0 || allele_count == 1) {
                for s in 0..self.values_by_sample_index.len() {
                    let evidence_count = self.evidence_by_sample_index.get(&s).unwrap().len();
                    for r in 0..evidence_count {
                        self.normalize_likelihoods_per_evidence(
                            maximum_likelihood_difference_cap,
                            s,
                            r,
                            symmetrically_normalize_alleles_to_reference,
                        );
                    }
                }
            };
        };
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    fn normalize_likelihoods_per_evidence(
        &mut self,
        maximum_best_alt_likelihood_difference: f64,
        sample_index: usize,
        evidence_index: usize,
        symmetrically_normalize_alleles_to_reference: bool,
    ) {
        //allow the best allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        let best_allele = self.search_best_allele(
            sample_index,
            evidence_index,
            symmetrically_normalize_alleles_to_reference,
            &None,
        );

        let worst_likelihood_cap = best_allele.likelihood + maximum_best_alt_likelihood_difference;
        debug!(
            "Worst likelihood cap {}, best likelihood {}, max best diff {}",
            worst_likelihood_cap, best_allele.likelihood, maximum_best_alt_likelihood_difference
        );

        let allele_count = self.alleles.len();

        let mut sample_values = &mut self.values_by_sample_index[sample_index];

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
        sample_index: usize,
        evidence_index: usize,
        can_be_reference: bool,
        priorities: &Option<Vec<i32>>,
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

        let sample_values = &self.values_by_sample_index[sample_index];
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
     * Updates the likelihood of the NonRef allele (if present) based on the likelihoods of a set of non-symbolic
     */
    pub fn update_non_ref_allele_likelihoods(
        &mut self,
        alleles_to_consider: AlleleList<A>,
        non_ref_allele_index: Option<usize>,
    ) {
        match non_ref_allele_index {
            None => return,
            Some(non_ref_allele_index) => {
                let allele_count = self.alleles.number_of_alleles();
                let non_symbolic_allele_count = allele_count - 1;
                // likelihood buffer reused across evidence:
                let mut qualified_allele_likelihoods = vec![0.0; non_symbolic_allele_count];
                for s in 0..self.samples.len() {
                    let evidence_count = self.evidence_by_sample_index.get(&s).unwrap().len();
                    for r in 0..evidence_count {
                        let sample_values = &self.values_by_sample_index[s];

                        let best_allele = self.search_best_allele(s, r, true, &None);
                        let mut number_of_qualified_allele_likelihoods = 0;
                        for i in 0..allele_count {
                            let allele_likelihood = sample_values[[i, r]];
                            if !allele_likelihood.is_nan() {
                                if i != non_ref_allele_index
                                    && allele_likelihood < best_allele.likelihood
                                    && alleles_to_consider
                                        .index_of_allele(self.alleles().get_allele(i))
                                        .is_some()
                                {
                                    qualified_allele_likelihoods
                                        [number_of_qualified_allele_likelihoods] =
                                        allele_likelihood;
                                    number_of_qualified_allele_likelihoods += 1;
                                }
                            }
                        }
                        let non_ref_likelihood = qualified_allele_likelihoods.median();
                        // when the median is NaN that means that all applicable likekihoods are the same as the best
                        // so the evidence is not informative at all given the existing alleles. Unless there is only one (or zero) concrete
                        // alleles with give the same (the best) likelihood to the NON-REF. When there is only one (or zero) concrete
                        // alleles we set the NON-REF likelihood to NaN.
                        let mut sample_values = &mut self.values_by_sample_index[s];

                        sample_values[[non_ref_allele_index, r]] = if !non_ref_likelihood.is_nan() {
                            non_ref_likelihood
                        } else if non_symbolic_allele_count <= 1 {
                            NAN
                        } else {
                            best_allele.likelihood
                        }
                    }
                }
            }
        }
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each evidence in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and evidence as the original.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this evidence-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    pub fn marginalize<'b>(
        &self,
        new_to_old_allele_map: &'b LinkedHashMap<usize, Vec<&'b A>>,
    ) -> Self {
        let new_alleles = new_to_old_allele_map.keys().collect::<Vec<&usize>>();
        let old_allele_count = self.alleles.len();
        let new_allele_count = new_alleles.len();

        // we get the index correspondence between new old -> new allele, None entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        let old_to_new_allele_index_map =
            self.old_to_new_allele_index_map(new_to_old_allele_map, old_allele_count, new_alleles);
        debug!("old to new allele map {:?}", &old_to_new_allele_index_map);
        // We calculate the marginal likelihoods.
        let new_likelihood_values = self.marginal_likelihoods(
            old_allele_count,
            new_allele_count,
            old_to_new_allele_index_map,
        );

        let sample_count = self.number_of_samples();
        debug!("new liklelihood values {:?}", &new_likelihood_values);

        let new_allele_list = new_to_old_allele_map
            .keys()
            .map(|i| self.alleles.get_allele(*i).clone())
            .collect::<Vec<A>>();

        let mut result = Self::new_from_likelihoods(
            AlleleList::new(&new_allele_list),
            self.samples.clone(),
            self.evidence_by_sample_index.clone(),
            self.filtered_evidence_by_sample_index.clone(),
            new_likelihood_values,
        );
        result.is_natural_log = self.is_natural_log;

        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    fn marginal_likelihoods(
        &self,
        old_allele_count: usize,
        new_allele_count: usize,
        old_to_new_allele_index_map: Vec<Option<usize>>,
    ) -> Vec<Array2<f64>> {
        let sample_count = self.samples.len();
        let mut result = vec![Array2::zeros((0, 0)); sample_count];
        debug!("new allele count {}", new_allele_count);

        for s in 0..sample_count {
            let sample_evidence_count = self.sample_evidence_count(s);
            let old_sample_values = &self.values_by_sample_index[s];

            let mut new_sample_values = Array2::zeros((new_allele_count, sample_evidence_count));
            // We initiate all likelihoods to -Inf.
            new_sample_values.fill(f64::NEG_INFINITY);

            // For each old allele and read we update the new table keeping the maximum likelihood.
            for r in 0..sample_evidence_count {
                for a in 0..old_allele_count {
                    let new_allele_index = old_to_new_allele_index_map[a];
                    match new_allele_index {
                        None => continue,
                        Some(new_allele_index) => {
                            let likelihood = old_sample_values[[a, r]];
                            if likelihood > new_sample_values[[new_allele_index, r]] {
                                new_sample_values[[new_allele_index, r]] = likelihood;
                            };
                        }
                    }
                }
            }
            result[s] = new_sample_values;
        }

        return result;
    }

    /**
     * Remove those reads that do not comply with a requirement.
     *
     * @param predicate the predicate representing the requirement.
     *
     * <p>
     *     This method modifies the current read-likelihoods collection.
     * </p>
     * <p>
     *     Any exception thrown by the predicate will be propagated to the calling code.
     * </p>
     *
     * @throws IllegalArgumentException if {@code predicate} is {@code null}.
     */
    pub fn retain_evidence(
        &mut self,
        predicate: &Box<dyn Fn(&BirdToolRead, &SimpleInterval) -> bool>,
        interval: &SimpleInterval,
    ) {
        let sample_count = self.samples.len();

        for s in 0..sample_count {
            // Remove evidence from the primary data
            let remove_indices = self
                .evidence_by_sample_index
                .get(&s)
                .unwrap()
                .into_iter()
                .enumerate()
                .filter(|(_, read)| !predicate(read, interval))
                .map(|(idx, _)| idx)
                .collect::<Vec<usize>>();
            self.remove_evidence_by_index(s, remove_indices);

            // If applicable also apply the predicate to the filters
            self.filtered_evidence_by_sample_index
                .get_mut(&s)
                .unwrap()
                .retain(|r| predicate(r, interval));
        }
    }

    /**
     * Add more evidence to the collection.
     *
     * @param evidenceBySample evidence to add.
     * @param initialLikelihood the likelihood for the new entries.
     *
     * @throws
     */
    pub fn add_evidence(
        &mut self,
        evidence_by_sample: HashMap<usize, Vec<BirdToolRead>>,
        initial_likelihood: f64,
    ) {
        for (sample_index, new_sample_evidence) in evidence_by_sample {
            if new_sample_evidence.is_empty() {
                continue;
            };

            let old_evidence_count = self
                .evidence_by_sample_index
                .get(&sample_index)
                .unwrap()
                .len();
            self.append_evidence(new_sample_evidence, sample_index);
            let new_evidence_count = self
                .evidence_by_sample_index
                .get(&sample_index)
                .unwrap()
                .len();
            self.extends_likelihood_arrays(
                initial_likelihood,
                sample_index,
                old_evidence_count,
                new_evidence_count,
            );
        }
    }

    /// Extends the likelihood array
    fn extends_likelihood_arrays(
        &mut self,
        initial_likelihood: f64,
        sample_index: usize,
        old_evidence_count: usize,
        new_evidence_count: usize,
    ) {
        let number_of_alleles = self.alleles.number_of_alleles();
        self.ensure_likelihoods_matrix_evidence_capacity(
            sample_index,
            new_evidence_count,
            number_of_alleles,
            initial_likelihood,
        )
    }

    /// Resizes the lk value holding arrays to be able to handle at least "X" amount of evidence
    fn ensure_likelihoods_matrix_evidence_capacity(
        &mut self,
        sample_index: usize,
        x: usize,
        number_of_alleles: usize,
        initial_likelihood: f64,
    ) {
        let current_capacity =
            self.likelihoods_matrix_evidence_capacity_by_sample_index[sample_index];
        if current_capacity < x {
            let new_capacity = x << 1; // we double it to avoid repetitive 1-element extensions resizing.
            for i in 0..(new_capacity - current_capacity) {
                self.values_by_sample_index[sample_index].push_column(ArrayView::from(
                    &vec![initial_likelihood; number_of_alleles],
                ));
            }
            self.likelihoods_matrix_evidence_capacity_by_sample_index[sample_index] = new_capacity;
        }
    }

    /// Append the new evidence reference into the structure per-sample, returning the count of evidence actually added (duplicates are not added)
    /// NOTE: the evidence-to-index cache is updated in place and not invalidated via {@link #invalidateEvidenceToIndexCache(int)} because adding new evidence
    /// to the cache, as opposed to removing evidence, is just a matter of appending entries
    fn append_evidence(&mut self, new_sample_evidence: Vec<BirdToolRead>, sample_index: usize) {
        let mut sample_evidence = self
            .evidence_by_sample_index
            .entry(sample_index)
            .or_insert(Vec::new());
        sample_evidence.par_extend(new_sample_evidence);
    }

    // calculates an old to new allele index map array.
    fn old_to_new_allele_index_map<'b>(
        &self,
        new_to_old_allele_map: &'b LinkedHashMap<usize, Vec<&'b A>>,
        old_allele_count: usize,
        new_alleles: Vec<&'b usize>,
    ) -> Vec<Option<usize>> {
        let mut old_to_new_allele_index_map = vec![None; old_allele_count];
        debug!("New alleles {:?}", &new_alleles);
        for new_index in 0..new_alleles.len() {
            let new_allele = new_alleles[new_index];
            for old_allele in new_to_old_allele_map.get(new_allele).unwrap() {
                let old_allele_index = self.index_of_allele(*old_allele);
                match old_allele_index {
                    None => {
                        panic!(
                            "Missing old {:?} allele in likelihood collection: New_allele {}",
                            old_allele, new_allele
                        );
                    }
                    Some(old_allele_index) => {
                        if old_to_new_allele_index_map[old_allele_index] != None {
                            panic!("Collision detected: Two new alleles refer to same old allele");
                        };
                        old_to_new_allele_index_map[old_allele_index] = Some(*new_allele);
                    }
                };
            }
        }

        return old_to_new_allele_index_map;
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
                        debug!(
                            "read value {} and thresh {} ",
                            self.maximum_likelihood_over_all_alleles(sample_index, *i),
                            (log10_min_true_likelihood)(&sample_evidence[*i])
                        );
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
            let mut old_evidence = self
                .evidence_by_sample_index
                .remove(&sample_index)
                .unwrap_or(Vec::new());
            let mut new_evidence = Vec::new();
            let mut num_removed = 0;
            for (n, read) in old_evidence.into_iter().enumerate() {
                if num_removed < num_to_remove && n == evidences_to_remove[num_removed] {
                    num_removed += 1;
                    filtered.push(read);
                } else {
                    new_evidence.push(read);

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
            self.evidence_by_sample_index
                .insert(sample_index, new_evidence);
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

    pub fn best_alleles_breaking_ties_main(
        &self,
        tie_breaking_priority: Box<dyn Fn(&A) -> i32>,
    ) -> Vec<BestAllele> {
        return (0..self.number_of_samples())
            .into_iter()
            .flat_map(|n| self.best_alleles_tie_breaking(n, &tie_breaking_priority))
            .collect::<Vec<BestAllele>>();
    }

    pub fn best_alleles_breaking_ties_for_sample(&self, sample_index: usize) -> Vec<BestAllele> {
        self.best_alleles_tie_breaking(
            sample_index,
            &AssemblyBasedCallerUtils::reference_tiebreaking_priority(),
        )
    }

    /**
     * Returns the collection of best allele estimates for one sample's evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit of evidence in the evidence-likelihoods collection.
     */
    fn best_alleles_tie_breaking(
        &self,
        sample_index: usize,
        tie_breaking_priority: &Box<dyn Fn(&A) -> i32>,
    ) -> Vec<BestAllele> {
        //TODO: this currently just does ref vs alt.  Really we want CIGAR complexity.
        let priorities = Some(
            self.alleles
                .list
                .iter()
                .map(|a| (tie_breaking_priority)(a))
                .collect::<Vec<i32>>(),
        );
        let evidence_count = self
            .evidence_by_sample_index
            .get(&sample_index)
            .unwrap()
            .len();

        return (0..evidence_count)
            .into_par_iter()
            .map(|r| self.search_best_allele(sample_index, r, true, &priorities))
            .collect::<Vec<BestAllele>>();
    }

    pub fn change_evidence(&mut self, evidence_replacements: HashMap<ReadIndexer, BirdToolRead>) {
        for (read_indexer, replacement_read) in evidence_replacements {
            let mut sample_evidence = self
                .evidence_by_sample_index
                .get_mut(&read_indexer.sample_index);

            match sample_evidence {
                Some(sample_evidence) => {
                    sample_evidence[read_indexer.evidence_index] = replacement_read;
                }
                None => continue,
            };
        }
    }
}

// pub struct LikelihoodMatrix<L: Locatable> {
//     pub(crate) evidence: Vec<L>,
//     pub(crate) alleles: Vec<Allele>,
//     pub(crate) values: Array2<f64>,
// }

#[derive(Debug)]
pub struct BestAllele {
    pub allele_index: Option<usize>,
    /**
     * The containing sample.
     */
    pub sample_index: usize,
    /**
     * The query evidence.
     */
    pub evidence_index: usize,
    /**
     * If allele != null, the indicates the likelihood of the evidence.
     */
    pub likelihood: f64,
    /**
     * Confidence that the evidence actually was generated under that likelihood.
     * This is equal to the difference between this and the second best allele match.
     */
    pub confidence: f64,
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

#[derive(Debug, Clone, Copy, Hash, Ord, PartialOrd, PartialEq, Eq)]
pub struct ReadIndexer {
    sample_index: usize,
    evidence_index: usize,
}

impl ReadIndexer {
    pub fn new(sample_index: usize, evidence_index: usize) -> ReadIndexer {
        Self {
            sample_index,
            evidence_index,
        }
    }
}
