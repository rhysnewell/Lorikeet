use genotype::genotype_allele_counts::GenotypeAlleleCounts;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use model::allele_list::AlleleListPermutation;
use model::byte_array_allele::Allele;
use ndarray::Array2;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::BinaryHeap;
use utils::math_utils::MathUtils;

#[derive(Clone, Debug)]
pub struct GenotypeLikelihoodCalculator {
    /**
     * Number of genotypes given this calculator {@link #ploidy} and {@link #alleleCount}.
     */
    pub genotype_count: i32,
    /**
     * Ploidy for this calculator.
     */
    pub ploidy: usize,
    /**
     * Number of genotyping alleles for this calculator.
     */
    pub allele_count: usize,
    /**
     * Offset table for this calculator.
     *
     * <p>
     *     This is a shallow copy of {@link GenotypeLikelihoodCalculators#alleleFirstGenotypeOffsetByPloidy} when the calculator was created
     *     thus it follows the same format as that array. Please refer to its documentation.
     * </p>
     *
     * <p>You can assume that this offset table contain at least (probably more) the numbers corresponding to the allele count and ploidy for this calculator.
     * However since it might have more than that and so you must use {@link #alleleCount} and {@link #ploidy} when
     * iterating through this array rather that its length or the length of its components.</p>.
     */
    pub allele_first_genotype_offset_by_ploidy: Array2<i32>,
    /**
     * Genotype table for this calculator.
     *
     * <p>It is ensure that it contains all the genotypes for this calculator ploidy and allele count, maybe more. For
     * that reason you must use {@link #genotypeCount} when iterating through this array and not relay on its length.</p>
     */
    pub genotype_allele_counts: Vec<GenotypeAlleleCounts>,
    /**
     * Buffer used as a temporary container for likelihood components for genotypes stratified by reads.
     *
     * <p>
     *     It is indexed by genotype index and then by read index. The read capacity is increased as needed by calling
     *     {@link #ensureReadCapacity(int) ensureReadCapacity}.
     * </p>
     */
    pub read_likelihoods_by_genotype_index: Vec<Vec<f64>>,
    /**
     * Buffer field use as a temporal container for sorted allele counts when calculating the likelihood of a
     * read in a genotype.
     * <p>
     *      This array follows the same format as {@link GenotypeAlleleCounts#sortedAlleleCounts}. Each component in the
     *      genotype takes up two positions in the array where the first indicate the allele index and the second its frequency in the
     *      genotype. Only non-zero frequency alleles are represented, taking up the first positions of the array.
     * </p>
     *
     * <p>
     *     This array is sized so that it can accommodate the maximum possible number of distinct alleles in any
     *     genotype supported by the calculator, value stored in {@link #maximumDistinctAllelesInGenotype}.
     * </p>
     */
    pub genotype_alleles_and_counts: Vec<usize>,
    /**
     * Maximum number of components (or distinct alleles) for any genotype with this calculator ploidy and allele count.
     */
    pub maximum_distinct_alleles_in_genotype: usize,
    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    pub last_overhead_counts: GenotypeAlleleCounts,
    /**
     * Buffer used as a temporary container for likelihood components for genotypes stratified by alleles, allele frequency and reads.
     *
     * <p>To improve performance we use a 1-dimensional array to implement a 3-dimensional one as some of those dimension
     * have typically very low depths (allele and allele frequency)</p>
     *
     * <p>
     *     The value contained in position <code>[a][f][r] == log10Lk(read[r] | allele[a]) + log10(f) </code>. Exception is
     *     for f == 0 whose value is undefined (in practice 0.0) and never used.
     * </p>
     *
     * <p>
     *     It is indexed by read, then by allele and then by the number of copies of the allele. For the latter
     *     there are as many entries as the ploidy of the calculator + 1 (to accommodate zero copies although is
     *     never used in practice).
     * </p>
     */
    pub(crate) read_allele_likelihood_by_allele_count: Vec<f64>,
    /**
     * Indicates how many reads the calculator supports.
     *
     * <p>This figure is increased dynamically as per the
     * calculation request calling {@link #ensureReadCapacity(int) ensureReadCapacity}.<p/>
     */
    pub read_capacity: i32,
    /**
     * Buffer field use as a temporal container for component likelihoods when calculating the likelihood of a
     * read in a genotype. It is stratified by read and the allele component of the genotype likelihood... that is
     * the part of the likelihood sum that correspond to a particular allele in the genotype.
     *
     * <p>
     *     It is implemented in a 1-dimensional array since typically one of the dimensions is rather small. Its size
     *     is equal to {@link #readCapacity} times {@link #maximumDistinctAllelesInGenotype}.
     * </p>
     *
     * <p>
     *     More concretely [r][i] == log10Lk(read[r] | allele[i]) + log(freq[i]) where allele[i] is the ith allele
     *     in the genotype of interest and freq[i] is the number of times it occurs in that genotype.
     * </p>
     */
    read_genotype_likelihood_components: Vec<f64>,

    /**
     * Max-heap for integers used for this calculator internally.
     */
    allele_heap: BinaryHeap<usize>,
}

impl GenotypeLikelihoodCalculator {
    pub fn new(
        ploidy: usize,
        allele_count: usize,
        allele_first_genotype_offset_by_ploidy: Array2<i32>,
        mut genotype_table_by_ploidy: Vec<Vec<GenotypeAlleleCounts>>,
    ) -> GenotypeLikelihoodCalculator {
        let genotype_count = allele_first_genotype_offset_by_ploidy[[ploidy, allele_count]];
        let maximum_distinct_alleles_in_genotype = std::cmp::min(ploidy, allele_count);
        // debug!("Genotype count {} Ploidy {} Allele Count {}", genotype_count, ploidy, allele_count);

        GenotypeLikelihoodCalculator {
            genotype_allele_counts: genotype_table_by_ploidy.remove(ploidy),
            read_likelihoods_by_genotype_index: vec![Vec::new(); genotype_count as usize],
            genotype_alleles_and_counts: vec![0; maximum_distinct_alleles_in_genotype as usize * 2],
            maximum_distinct_alleles_in_genotype,
            last_overhead_counts: GenotypeAlleleCounts::build_empty(),
            read_capacity: -1,
            genotype_count,
            allele_count,
            ploidy,
            allele_first_genotype_offset_by_ploidy,
            read_genotype_likelihood_components: vec![],
            allele_heap: BinaryHeap::with_capacity(ploidy),
            read_allele_likelihood_by_allele_count: Vec::new(),
        }
    }

    // /**
    //  * Makes sure that temporal arrays and matrices are prepared for a number of reads to process.
    //  * @param requestedCapacity number of read that need to be processed.
    //  */
    /**
     * Returns the genotype associated to a particular likelihood index.
     *
     * <p>If {@code index} is larger than {@link GenotypeLikelihoodCalculators#MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY},
     *  this method will reconstruct that genotype-allele-count iteratively from the largest strongly referenced count available.
     *  or the last requested index genotype.
     *  </p>
     *
     * <p> Therefore if you are iterating through all genotype-allele-counts you should do sequentially and incrementally, to
     * avoid a large efficiency drop </p>.
     *
     * @param index query likelihood-index.
     * @return never {@code null}.
     */
    pub fn genotype_allele_counts_at(&mut self, index: usize) -> &mut GenotypeAlleleCounts {
        if index >= self.genotype_count as usize {
            panic!(
                "Invalid likelihood index {} >= {} (Genotype count for n-alleles = {} and {}",
                index, self.genotype_count, self.allele_count, self.ploidy
            );
        } else if index < self.genotype_allele_counts.len() {
            return &mut self.genotype_allele_counts[index];
        } else if self.last_overhead_counts.is_null() || self.last_overhead_counts.index() > index {
            let mut result = self.genotype_allele_counts
                [GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY - 1]
                .clone();

            // let mut result = &mut self.genotype_allele_counts[index];
            //
            result.increase(
                index as i32
                    - GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY as i32
                    + 1,
            );

            // self.genotype_allele_counts.push(result);
            // result.increase(1);
            //
            self.last_overhead_counts = result;
            // return result;
            return &mut self.last_overhead_counts;
        } else {
            self.last_overhead_counts
                .increase(index as i32 - self.last_overhead_counts.index() as i32);
            return &mut self.last_overhead_counts;
            // return &mut self.genotype_allele_counts[index]
        }
    }

    /**
     * Give a list of alleles, returns the likelihood array index.
     * @param alleleIndices the indices of the alleles in the genotype, there should be as many repetition of an
     *                      index as copies of that allele in the genotype. Allele indices do not need to be sorted in
     *                      any particular way.
     *
     * @return never {@code null}.
     */
    pub fn alleles_to_index(&mut self, allele_indices: &[usize]) -> usize {
        if self.ploidy == 0 {
            return 0;
        }

        self.allele_heap.clear();
        for i in 0..allele_indices.len() {
            self.allele_heap.push(allele_indices[i]);
        }

        return self.allele_heap_to_index();
    }

    /**
     * Returns the likelihood index given the allele counts.
     *
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts} with 0 offset.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is not a valid {@code allele count array}:
     *  <ul>
     *      <li>is {@code null},</li>
     *      <li>or its length is not even,</li>
     *      <li>or it contains any negatives,
     *      <li>or the count sum does not match the calculator ploidy,</li>
     *      <li>or any of the alleles therein is negative or greater than the maximum allele index.</li>
     *  </ul>
     *
     * @return 0 or greater but less than {@link #genotypeCount}.
     */
    pub fn allele_counts_to_index(&mut self, allele_count_array: &[usize]) -> usize {
        if allele_count_array.len() % 2 != 0 {
            panic!("The allele counts array cannot have odd length")
        }
        self.allele_heap.clear();
        for i in (0..allele_count_array.len()).step_by(2) {
            let index = allele_count_array[i];
            let count = allele_count_array[i + 1];
            for _ in (0..count).into_iter() {
                self.allele_heap.push(index)
            }
        }

        return self.allele_heap_to_index();
    }

    /**
     * Transforms the content of the heap into an index.
     *
     * <p>
     *     The heap contents are flushed as a result, so is left ready for another use.
     * </p>
     *
     * @return a valid likelihood index.
     */
    fn allele_heap_to_index(&mut self) -> usize {
        if self.allele_heap.len() != self.ploidy {
            panic!("The sum of allele counts must be equal to the ploidy of the calculator");
        }

        if self.allele_heap.peek().unwrap() >= &self.allele_count {
            panic!(
                "Invalid allele {:?} more than the maximum {}",
                self.allele_heap.peek(),
                self.allele_count - 1
            )
        }

        let mut result = 0;
        for p in (0..=self.ploidy).into_iter().rev() {
            if p == 0 {
                continue;
            };
            let allele = self.allele_heap.pop().unwrap();
            result += self.allele_first_genotype_offset_by_ploidy[[p, allele]]
        }
        return result as usize;
    }

    /**
     * Calculate the likelihoods given the list of alleles and the likelihood map.
     *
     * @param likelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code likelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code likelihoods}.
     *
     * @return never {@code null}.
     */
    pub fn genotype_likelihoods<A: Allele>(
        &mut self,
        likelihoods: &Array2<f64>,
        permutation: &AlleleLikelihoodMatrixMapper<A>,
        number_of_evidences: usize,
    ) -> GenotypeLikelihoods {
        debug!("Single likelihoods array {:?}", likelihoods);
        debug!("Number of evidences {}", number_of_evidences);
        debug!("permutation {:?}", permutation);
        let read_likelihoods_by_genotype_index = self
            .get_read_raw_read_likelihoods_by_genotype_index(
                likelihoods,
                permutation,
                number_of_evidences,
            );
        return GenotypeLikelihoods::from_log10_likelihoods(read_likelihoods_by_genotype_index);
    }

    /**
     * A helper method that actually does the matrix operations but returns the raw values.
     *
     * @return the raw array (in log10 likelihoods space) of the GL for each genotype
     */
    fn get_read_raw_read_likelihoods_by_genotype_index<A: Allele>(
        &mut self,
        likelihoods: &Array2<f64>,
        permutation: &AlleleLikelihoodMatrixMapper<A>,
        number_of_evidences: usize,
    ) -> Vec<f64> {
        assert!(
            permutation.permutation.number_of_alleles() == self.allele_count,
            "Mismatch between likelihood matrix and allele_count {} -> {}",
            permutation.permutation.number_of_alleles(),
            self.allele_count
        );

        self.ensure_read_capcity(likelihoods.ncols());

        // [x][y][z] = z * LnLk(Read_x | Allele_y)
        self.read_likelihood_components_by_allele_count(
            likelihoods,
            permutation,
            number_of_evidences,
        );
        self.genotype_likelihood_by_read(number_of_evidences);

        return self.genotype_likelihoods_private(number_of_evidences);
    }

    /**
     * Calculates the final genotype likelihood array out of the likelihoods for each genotype per read.
     *
     * @param readLikelihoodsByGenotypeIndex <i>[g][r]</i> likelihoods for each genotype <i>g</i> and <i>r</i>.
     * @param readCount number of reads in the input likelihood arrays in {@code genotypeLikelihoodByRead}.
     * @return never {@code null}, one position per genotype where the <i>i</i> entry is the likelihood of the ith
     *   genotype (0-based).
     */
    fn genotype_likelihoods_private(&self, read_count: usize) -> Vec<f64> {
        let mut result = vec![0.0; self.genotype_count as usize];
        let denominator = (read_count as f64) * (self.ploidy as f64).log10();

        // instead of dividing each read likelihood by ploidy ( so subtract log10(ploidy) )
        // we multiply them all and the divide by ploidy^readCount (so substract readCount * log10(ploidy) )
        for g in 0..self.genotype_count as usize {
            result[g] = self.read_likelihoods_by_genotype_index[g][0..read_count]
                .iter()
                .sum::<f64>()
                - denominator;
        }

        return result;
    }

    /**
     * Calculates the likelihood component of each read on each genotype.
     *
     * NOTE: this is not actually the read likelihood component for each genotype, it is the sum of the log read likelihoods components
     *       for each genotype without having been normalized by the the denominator of the ploidy, that happens in the final step
     *
     * @param readLikelihoodComponentsByAlleleCount [a][f][r] likelihood stratified by allele <i>a</i>, frequency in genotype <i>f</i> and
     *                                              read <i>r</i>.
     * @param readCount number of reads in {@code readLikelihoodComponentsByAlleleCount}.
     * @return never {@code null}.
     */
    fn genotype_likelihood_by_read(
        &mut self,
        // equivalent vec of self.allele_count
        read_count: usize,
    ) {
        // Here we don't use the convenience of {@link #genotypeAlleleCountsAt(int)} within the loop to spare instantiations of
        // GenotypeAlleleCounts class when we are dealing with many genotypes.

        let mut allele_counts_index = 0;
        // let mut allele_counts = &mut self.genotype_allele_counts[allele_counts_index];
        for genotype_index in 0..(self.genotype_count as usize) {
            // let mut read_likelihoods = &mut self.read_likelihoods_by_genotype_index[genotype_index];
            let component_count =
                self.genotype_allele_counts[allele_counts_index].distinct_allele_count();

            match component_count {
                1 => {
                    self.single_component_genotype_likelihood_by_read(
                        genotype_index,
                        allele_counts_index,
                        read_count,
                    );
                }
                2 => {
                    self.two_component_genotype_likelihood_by_read(
                        genotype_index,
                        allele_counts_index,
                        read_count,
                    );
                }
                _ => {
                    self.many_component_genotype_likelihood_by_read(
                        genotype_index,
                        allele_counts_index,
                        read_count,
                    );
                }
            }
            if genotype_index < (self.genotype_count - 1) as usize {
                allele_counts_index = self.next_genotype_allele_counts(allele_counts_index);
                // allele_counts = &mut self.genotype_allele_counts[allele_counts_index];
            }
        }
    }

    fn next_genotype_allele_counts(
        &mut self,
        // allele_counts: &mut GenotypeAlleleCounts,
        index: usize,
    ) -> usize {
        // let index = allele_counts.index();
        let mut result;
        let cmp = (index + 1)
            .checked_sub(GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY);
        match cmp {
            None => {
                // result = &mut self.genotype_allele_counts[index + 1];
                return index + 1;
            }
            Some(cmp) => {
                if cmp == 0 {
                    self.genotype_allele_counts[index].increase(1);
                    return index;
                } else {
                    // allele_counts.increase(1);
                    // result = allele_counts;
                    // return result
                    result = self.genotype_allele_counts[index].clone();
                    result.increase(1);
                    self.genotype_allele_counts.push(result);
                    return index + 1;
                }
                // result = self.genotype_allele_counts[index].clone();
                // result.increase(1);
                // self.genotype_allele_counts.push(result);
                // return index + 1;
            }
        }
    }

    /**
     * Calculates the likelihood component by read for a given genotype allele count assuming that there are
     * exactly one allele present in the genotype.
     */
    fn single_component_genotype_likelihood_by_read(
        &mut self,
        genotype_index: usize,
        genotype_allele_counts_index: usize,
        read_count: usize,
    ) {
        let mut likelihood_by_read = &mut self.read_likelihoods_by_genotype_index[genotype_index];
        let mut genotype_allele_counts =
            &mut self.genotype_allele_counts[genotype_allele_counts_index];
        let allele = genotype_allele_counts.allele_index_at(0);
        // the count of the only component must be = ploidy.
        // let mut offset = (allele * (self.ploidy + 1) + self.ploidy) * read_count;
        likelihood_by_read[..read_count].clone_from_slice(
            &self.read_allele_likelihood_by_allele_count[((allele * (self.ploidy + 1)
                + self.ploidy)
                * read_count)
                ..(read_count + (allele * (self.ploidy + 1) + self.ploidy) * read_count)],
        );
    }

    /**
     * Calculates the likelihood component by read for a given genotype allele count assuming that there are
     * exactly two alleles present in the genotype (with arbitrary non-zero counts each).
     */
    fn two_component_genotype_likelihood_by_read(
        &mut self,
        genotype_index: usize,
        genotype_allele_counts_index: usize,
        read_count: usize,
    ) {
        let mut genotype_allele_counts =
            &mut self.genotype_allele_counts[genotype_allele_counts_index];
        let allele_0 = genotype_allele_counts.allele_index_at(0);
        let freq_0 = genotype_allele_counts.allele_count_at(0);
        let allele_1 = genotype_allele_counts.allele_index_at(1);
        let freq_1 = self.ploidy - freq_0; // no need to get it from genotypeAlleleCounts.

        let mut allele_0_LnL_offset = read_count * ((self.ploidy + 1) * allele_0 + freq_0);
        let mut allele_1_LnL_offset = read_count * ((self.ploidy + 1) * allele_1 + freq_1);
        let mut likelihood_by_read = &mut self.read_likelihoods_by_genotype_index[genotype_index];

        for r in 0..read_count {
            let ln_lk_0 = self.read_allele_likelihood_by_allele_count[allele_0_LnL_offset];
            allele_0_LnL_offset += 1;
            let ln_lk_1 = self.read_allele_likelihood_by_allele_count[allele_1_LnL_offset];
            allele_1_LnL_offset += 1;
            likelihood_by_read[r] = MathUtils::approximate_log10_sum_log10(ln_lk_0, ln_lk_1);
        }
    }

    /**
     * General genotype likelihood component by read calculator. It does not make any assumption in the exact
     * number of alleles present in the genotype.
     */
    fn many_component_genotype_likelihood_by_read(
        &mut self,
        genotype_index: usize,
        genotype_allele_counts_index: usize,
        read_count: usize,
    ) {
        // First we collect the allele likelihood component for all reads and place it
        // in readGenotypeLikelihoodComponents for the final calculation per read.
        let mut genotype_allele_counts =
            &mut self.genotype_allele_counts[genotype_allele_counts_index];

        genotype_allele_counts.copy_allele_counts(&mut self.genotype_alleles_and_counts, 0);

        let component_count = genotype_allele_counts.distinct_allele_count();
        let allele_data_size = (self.ploidy + 1) * read_count;
        let mut cc = 0;
        for c in 0..component_count {
            let allele_index = self.genotype_alleles_and_counts[cc];
            cc += 1;
            let allele_count = self.genotype_alleles_and_counts[cc];
            cc += 1;
            // alleleDataOffset will point to the index of the first read likelihood for that allele and allele count.
            let mut allele_data_offset =
                allele_data_size * allele_index + allele_count * read_count;

            let mut read_data_offset = c;
            for r in 0..read_count {
                self.read_genotype_likelihood_components[read_data_offset] =
                    self.read_allele_likelihood_by_allele_count[allele_data_offset];
                allele_data_offset += 1;
                read_data_offset += self.maximum_distinct_alleles_in_genotype;
            }
        }

        let mut likelihood_by_read = &mut self.read_likelihoods_by_genotype_index[genotype_index];
        // Calculate the likelihood per read.
        let mut read_data_offset = 0;
        for r in 0..read_count {
            likelihood_by_read[r] = MathUtils::approximate_log10_sum_log10_vec(
                &self.read_genotype_likelihood_components,
                read_data_offset,
                read_data_offset + component_count,
            );
            read_data_offset += self.maximum_distinct_alleles_in_genotype;
        }
    }

    /**
     * Makes sure that temporal arrays and matrices are prepared for a number of reads to process.
     * @param requestedCapacity number of read that need to be processed.
     */
    pub fn ensure_read_capcity(&mut self, requested_capacity: usize) {
        if self.read_capacity == -1 {
            let minimum_capacity = max(requested_capacity, 10); // Never go too small, 10 is the minimum.
            self.read_allele_likelihood_by_allele_count =
                vec![0.0; minimum_capacity * self.allele_count * (self.ploidy + 1)];
            for i in 0..self.genotype_count as usize {
                self.read_likelihoods_by_genotype_index[i] = vec![0.0; minimum_capacity];
            }
            self.read_genotype_likelihood_components = vec![0.0; self.ploidy * minimum_capacity];
            self.read_capacity = minimum_capacity as i32;
        } else if (self.read_capacity as usize) < requested_capacity {
            let double_capacity = (requested_capacity as usize) << 1;
            self.read_allele_likelihood_by_allele_count =
                vec![0.0; double_capacity * self.allele_count * (self.ploidy + 1)];
            for i in 0..self.genotype_count as usize {
                self.read_likelihoods_by_genotype_index[i] = vec![0.0; double_capacity];
            }
            self.read_genotype_likelihood_components =
                vec![0.0; self.maximum_distinct_alleles_in_genotype * double_capacity];
            self.read_capacity = double_capacity as i32;
        }
    }

    /**
     * Returns a 3rd matrix with the likelihood components.
     *
     * <pre>
     *     result[y][z][x] :=  z * lnLk ( read_x | allele_y ).
     * </pre>
     *
     * @return never {@code null}.
     */
    fn read_likelihood_components_by_allele_count<A: Allele>(
        &mut self,
        likelihoods: &Array2<f64>,
        permutation: &AlleleLikelihoodMatrixMapper<A>,
        number_of_evidences: usize,
    ) {
        let read_count = number_of_evidences;
        let allele_data_size = read_count * (self.ploidy + 1);

        // frequency1Offset = readCount to skip the useless frequency == 0. So now we are at the start frequency == 1
        // frequency1Offset += alleleDataSize to skip to the next allele index data location (+ readCount) at each iteration.
        let mut frequency_1_offset = read_count;
        for a in 0..self.allele_count {
            debug!(
                "no slice cloning from slice offset {} evidences {} dest size {} source size {}",
                frequency_1_offset,
                number_of_evidences,
                self.read_allele_likelihood_by_allele_count.len(),
                &likelihoods.row(permutation.permutation.from_index(a)).len()
            );
            debug!(
                "cloning from slice offset {} evidences {} dest size {} source size {}",
                frequency_1_offset,
                number_of_evidences,
                self.read_allele_likelihood_by_allele_count
                    [frequency_1_offset..frequency_1_offset + number_of_evidences]
                    .len(),
                &likelihoods.row(permutation.permutation.from_index(a)).len()
            );
            self.read_allele_likelihood_by_allele_count
                [frequency_1_offset..frequency_1_offset + number_of_evidences]
                .clone_from_slice(
                    &likelihoods
                        .row(permutation.permutation.from_index(a))
                        .as_slice()
                        .unwrap()[0..number_of_evidences],
                );

            // p = 2 because the frequency == 1 we already have it.
            let mut destination_offset = frequency_1_offset + read_count;
            for frequency in 2..=self.ploidy {
                let log10_frequency = (frequency as f64).log10();
                let mut source_offset = frequency_1_offset;
                for r in 0..read_count {
                    self.read_allele_likelihood_by_allele_count[destination_offset] = self
                        .read_allele_likelihood_by_allele_count[source_offset]
                        + log10_frequency;
                    destination_offset += 1;
                    source_offset += 1;
                }
            }
            frequency_1_offset += allele_data_size;
        }
    }

    /**
     * Composes a genotype index map given a allele index recoding.
     *
     * @param oldToNewAlleleIndexMap allele recoding. The ith entry indicates the index of the allele in original encoding
     *                               that corresponds to the ith allele index in the final encoding.
     *
     * @throws IllegalArgumentException if this calculator cannot handle the recoding provided. This is
     * the case when either {@code oldToNewAlleleIndexMap}'s length or any of its element (+ 1 as they are 0-based) is larger
     * this calculator's {@link #alleleCount()}. Also if any {@code oldToNewAllelesIndexMap} element is negative.
     *
     * @return never {@code null}.
     */
    pub fn genotype_index_map(&mut self, old_to_new_allele_index_map: &[usize]) -> Vec<usize> {
        let result_allele_count = old_to_new_allele_index_map.len();
        assert!(result_allele_count <= self.allele_count, "This calculator does not have enough capacity for handling {} alleles", result_allele_count);

        let result_length = if result_allele_count == self.allele_count {
            self.genotype_count as usize
        } else {
            GenotypeLikelihoodCalculators::genotype_count(self.ploidy, result_allele_count) as usize
        };

        let mut result = vec![0; result_length];
        let mut sorted_allele_counts = vec![0; max(self.ploidy, self.allele_count) << 1];
        self.allele_heap.clear();
        let mut allele_counts_index = 0;
        for i in 0..result_length {
            self.genotype_index_map_per_genotype_index(i, allele_counts_index, old_to_new_allele_index_map, &mut result, &mut sorted_allele_counts);
            if i < result_length - 1 {
                allele_counts_index = self.next_genotype_allele_counts(allele_counts_index);
            }
        }

        result
    }

    /**
     * Performs the genotype mapping per new genotype index.
     *
     * @param newGenotypeIndex the target new genotype index.
     * @param alleleCounts tha correspond to {@code newGenotypeIndex}.
     * @param oldToNewAlleleIndexMap the allele mapping.
     * @param destination where to store the new genotype index mapping to old.
     * @param sortedAlleleCountsBuffer a buffer to re-use to get the genotype-allele-count's sorted allele counts.
     */
    fn genotype_index_map_per_genotype_index(
        &mut self,
        new_genotype_index: usize,
        // old_allele_counts_index: Option<usize>,
        allele_counts_index: usize,
        old_to_new_allele_index_map: &[usize],
        destination: &mut Vec<usize>,
        sorted_allele_counts_buffer: &mut Vec<usize>,
    ) {
        // let distinct_allele_count = match old_allele_counts_index {
        //     None => {
        //         let distinct_allele_count = self.genotype_allele_counts[allele_counts_index].distinct_allele_count();
        //         self.genotype_allele_counts[allele_counts_index].copy_allele_counts(destination, 0);
        //         distinct_allele_count
        //     },
        //     Some(old_index) => {
        //         if old_index == allele_counts_index
        //     }
        // }
        let distinct_allele_count = self.genotype_allele_counts[allele_counts_index].distinct_allele_count();
        self.genotype_allele_counts[allele_counts_index].copy_allele_counts(destination, 0);

        let mut jj = 0;
        for _ in 0..distinct_allele_count {
            let old_index = sorted_allele_counts_buffer[jj];
            jj += 1;
            let repeats = sorted_allele_counts_buffer[jj];
            jj += 1;
            let new_index = old_to_new_allele_index_map[jj];
            jj += 1;

            if new_index >= self.allele_count {
                panic!("Found invalid new allele index ({}) for old index ({})", new_index, old_index);
            };

            for k in 0..repeats {
                self.allele_heap.push(new_index);
            }
        }

        let genotype_index = self.allele_heap_to_index();
        destination[new_genotype_index] = genotype_index;
    }
}
