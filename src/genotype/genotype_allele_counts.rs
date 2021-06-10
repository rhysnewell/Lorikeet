/**
The following code is adapted from the broadinstitute GATK HaplotypeCaller program
*/

use ordered_float::{NotNan, OrderedFloat};
use ndarray::{Array, Array2, ArrayBase, OwnedRepr};
use utils::math_utils::MathUtils;
use rayon::prelude::*;


/**
 * Collection of allele counts for a genotype. It encompasses what alleles are present in the genotype and in what number.</p>
 *
 * <p>Alleles are represented herein by their indices running from <b>0</b> to <b>N-1</b> where <i>N</i> is the number of alleles.</p>
 *
 * <p>Genotypes are represented as a single array of alternating alleles and counts, where only alleles with non-zero counts are included:
 * [allele 1, count1, allele 2, count2. . .]</p>
 *
 * <p>Each allele present in a genotype (count != 0) has a <i>rank</i>, that is the 0-based ordinal of
 * that allele amongst the ones present in the genotype as sorted by their index.</p>
 *
 * <p>For example:</p>
 *
 * <p><b>[0,1,2,1]</b> has two alleles with indices <b>0</b> and <b>2</b>, both with count 1, corresponding to diploid genotype 0/2.
 * The rank of <b>0</b> is <i>0</i> whereas the rank of <b>2</b> is <i>1</i>.</p>
 *
 * <p><b>[2,1,4,2,7,1]</b> has three alleles with indices <b>2</b>, <b>4</b> and <b>7</b>. <b>2</b> and <b>7</b> have count 1 whereas <b>4</b> has count 2.
 * It corresponds to tetraploid genotype 2/4/4/7
 * The rank of <b>2</b> is <i>0</i>, the rank of <b>4</b> is <i>1</i>. and the rank of <b>7</b> is <i>2</i>.</p>
 *
 * <p>In contrast, in both examples above both <b>3</b> and <b>10</b> (and many others) are absent thus they have no rank (represented by <i>-1</i> whenever applies).</p>
 *
 * <p><b>[0,0,1,2]</b> is not valid because allele 0 has a count of 0 and should be absent from the array.</p>
 *
 * <p><b>[1,1,0,1]</b> is not valid because allele 1 comes before allele 0.</p>
 *
 * <p>{@link GenotypeAlleleCounts} instances have themselves their own index (returned by {@link #index() index()}, that indicate their 0-based ordinal within the possible genotype combinations with the same ploidy.</p>
 *
 * <p>For example, for ploidy 3:</p>
 *
 * <table>
 *     <th>Index</th><th>Genotype</th>
 *     <tr><td>0</td><td><b>0/0/0</b></td></tr>
 *     <tr><td>1</td><td><b>0/0/1</b></td></tr>
 *     <tr><td>2</td><td><b>0/1/1</b></td></tr>
 *     <tr><td>3</td><td><b>1/1/1</b></td></tr>
 *     <tr><td>4</td><td><b>0/0/2</b></td></tr>
 *     <tr><td>6</td><td><b>0/1/2</b></td></tr>
 *     <tr><td>7</td><td><b>1/1/2</b></td></tr>
 *     <tr><td>8</td><td><b>0/2/2</b></td></tr>
 *     <tr><td>9</td><td><b>1/2/2</b></td></tr>
 *     <tr><td>10</td><td><b>2/2/2</b></td></tr>
 *     <tr><td>11</td><td><b>0/0/3</b></td></tr>
 *     <tr><td>12</td><td><b>0/1/3</b></td></tr>
 *     <tr><td>13</td><td><b>1/1/3</b></td></tr>
 *     <tr><td>14</td><td><b>0/2/3</b></td></tr>
 *     <tr><td>15</td><td><b>1/2/3</b></td></tr>
 *     <tr><td>16</td><td><b>2/2/3</b></td></tr>
 *     <tr><td>17</td><td><b>0/3/3</b></td></tr>
 *     <tr><td>...</td><td>...</td></tr>
 * </table>
 *
 * The total number of possible genotypes is only bounded by the maximum allele index.
 */

#[derive(Clone, Debug)]
pub struct GenotypeAlleleCounts {
    ploidy: usize,
    distinct_allele_count: usize,
    sorted_allele_counts: Vec<usize>,
    index: usize,
    log10_combination_count: f64,
}

impl GenotypeAlleleCounts {
    const UNCOMPUTED_LOG_10_COMBINATION_COUNT: f64 = -1.;

    pub fn build_empty() -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy: 0,
            index: 0,
            sorted_allele_counts: Vec::new(),
            distinct_allele_count: 0,
            log10_combination_count: -1.,
        }
    }

    pub fn is_null(&self) -> bool {
        if self.ploidy == 0
            && self.index == 0
            && self.sorted_allele_counts.len() == 0
            && self.distinct_allele_count == 0
            && self.log10_combination_count == -1 {
            true
        }
        false
    }

    /**
     * Creates a new unphased genotype.
     *
     * <p>This method assumes that the invoker is passing a well formatted and sorted allele frequency array.
     * Not checks are done for the sake of performance.</p>
     *
     * <p>
     *     The input argument {@code sortedAlleleCounts} list the index of alleles included in the unphased genotype
     *     and their frequency in the genotype in a single array using consecutive pairs:<br/>
     *
     *     <pre> [allele_1,freq_1,allele_2,freq_2, ... , allele_i, freq_i, ... , allele_n, freq_n]</pre>
     *
     *     <br/>
     *     No entry can have frequency == 0 (these must be omitted) and entries are sorted by allele index without
     *     any repetitions so that if <i>i < j</i> then <i>allele_i < allele_j</i>.
     *
     * </p>
     *
     * <p>
     *     The {@code ploidy} provided must be equal to the sum of all frequencies in {@code sortedAlleleCounts}
     * </p>
     * @param ploidy the genotype ploidy.
     * @param sortedAlleleCounts sorted allele counts following the restrictions above.
     * @param index the genotype index.
     */
    pub fn build(ploidy: usize, index: usize, sorted_allele_counts: &[usize]) -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy,
            index,
            sorted_allele_counts: Vec::from(sorted_allele_counts),
            distinct_allele_count: sorted_allele_counts.len() >> 1,
            log10_combination_count: -1.,
        }
    }

    pub fn build_with_count(
        ploidy: usize,
        index: usize,
        sorted_allele_counts: Vec<usize>,
        distinct_allele_count: usize
    ) -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy,
            index,
            sorted_allele_counts,
            distinct_allele_count,
            log10_combination_count: -1.,
        }
    }

    pub fn ploidy(&self) -> usize { self.ploidy }

    pub fn index(&self) -> usize { self.index }
    /**
     * Gets the log10 combination count, computing it if uninitialized.  Note that the invoked MathUtils method uses fast cached
     * log10 values of integers for any reasonable ploidy.
     *
     * This method should be invoked on instances of {@link GenotypeAlleleCounts} cached in {@link GenotypeLikelihoodCalculators::genotypeTableByPloidy}.
     * Such usage allows the result of this computation to be cached once for an entire run of HaplotypeCaller.
     * @return
     */
    pub fn log10_combination_count(&mut self) -> f64 {
        if self.log10_combination_count == GenotypeAlleleCounts::UNCOMPUTED_LOG_10_COMBINATION_COUNT {
            self.log10_combination_count =
                MathUtils::log10_factorial(self.ploidy as f64)
                    - (0..self.distinct_allele_count).into_par_iter().map(|i| {
                    &self.sorted_allele_counts[i]
                }).sum()
        }
        return self.log10_combination_count
    }

    pub fn sum_over_allele_indices_and_counts<F>(&self, f: F) -> f64
        where F: Fn(usize, f64, &Vec<f64>) -> f64 {
        // let result = (0..self.distinct_allele_count)
        0.
    }

    /**
     * Increases the allele counts a number of times.
     *
     * <p>
     *     This method must not be invoked on cached genotype-allele-counts that are meant to remain constant,
     *     such as the ones contained in {@link GenotypeLikelihoodCalculators#genotypeTableByPloidy}.
     * </p>
     *
     * @param times the number of times to increase.
     *
     * @throws IllegalArgumentException if {@code times} is negative.
     */
    pub fn increase(&mut self, times: i32) {
        if times >= 0 {
            for _ in 0..times {
                // if the ploidy is zero there is only one possible genotype.
                if self.distinct_allele_count == 0 {
                    break
                }

                if self.distinct_allele_count == 1 {
                    if self.ploidy == 1 {
                        self.sorted_allele_counts[0] += 1;
                    } else {
                        if self.sorted_allele_counts.len() < 4 {
                            self.sorted_allele_counts.resize(4, 0);
                        }

                        self.sorted_allele_counts[2] = self.sorted_allele_counts[0] + 1;
                        self.sorted_allele_counts[3] = 1;
                        self.sorted_allele_counts[0] = 0;
                        self.sorted_allele_counts[1] = ploidy - 1;
                        self.distinct_allele_count = 2
                    }
                } else {
                    // Now, all the following ifs are just the way to avoid working with dynamically sizing Vec<Vec<i32>>
                    // as the final size of the resulting new sorted-allele-counts array varies depending on the situation.
                    // this is considerably faster and the logic complexity would not be that different actually so it is worth
                    // the if indentations.
                    //
                    // Notice that at this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
                    //
                    // We only need to look at the two lowest allele indices to decide what to do.
                    let allele0 = self.sorted_allele_counts[0];
                    let freq0 = self.sorted_allele_counts[1];
                    let allele1 = self.sorted_allele_counts[2];
                    let allele0_plus1 = allele0 + 1;

                    let allele0_and_1_are_consecutive = allele0_plus1 == allele1;

                    // The rest of the sorted allele counts array contains junk
                    let sorted_allele_counts_length = self.distinct_allele_count << 1;
                    if freq0 == 1 {
                        // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
                        self.sorted_allele_counts[0..(sorted_allele_counts_length - 2) as usize] = self.sorted_allele_counts[2..(sorted_allele_counts_length - 2) as usize];
                        self.sorted_allele_counts[1] += 1;
                        self.distinct_allele_count -= 1;
                    } else {
                        if allele0_and_1_are_consecutive { // we don't need to add a component for allele0 + 1 since it already exists.
                            self.sorted_allele_counts[0] = 0;
                            self.sorted_allele_counts[1] = freq0 - 1;
                            self.sorted_allele_counts[3] += 1;
                        } else { // we need to insert allele0 + 1 in the sorted-allele-counts array and give it frequency 1.
                            if self.sorted_allele_counts.len() < (sorted_allele_counts_length + 2) as usize {
                                // make room for new component
                                self.sorted_allele_counts.resize((sorted_allele_counts_length + 2) as usize, 0)
                            }
                            self.sorted_allele_counts[4..(sorted_allele_counts_length - 2) as usize] = self.sorted_allele_counts[2..(sorted_allele_counts_length - 2) as usize];
                            self.sorted_allele_counts[0] = 0;
                            self.sorted_allele_counts[1] = freq0 - 1;
                            self.sorted_allele_counts[2] = allele0_plus1;
                            self.sorted_allele_counts[3] = 1;
                            self.distinct_allele_count += 1;
                        }
                    }
                }
                self.index += 1;
                self.log10_combination_count = -1.;
            }
        }
    }


    /**
     * Calculates the next genotype in likelihood indexing order.
     * @return never null.
     */
    pub fn next(&mut self) -> GenotypeAlleleCounts {
        // a few cases worth explicitly optimizing
        if self.distinctAlleleCount == 0 {
            return self.clone();    // only one possible genotype with zero ploidy
        } else if distinctAlleleCount == 1 && ploidy == 1 {
            return GenotypeAlleleCounts::build(1, self.index + 1, &[self.sorted_allele_counts[0] + 1, 1]);    // A -> B , D -> E etc...
        } else if distinctAlleleCount == 1 {
            return GenotypeAlleleCounts::build(ploidy, index + 1, &[0, ploidy - 1, sortedAlleleCounts[0] + 1, 1]);    // AAAAA -> AAAAB, DDD -> AAE etc...
        }

        // The following logic avoids dynamically sizing the new sorted-allele-counts array, which would be very slow
        // At this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
        // We only need to look at the two lowest allele indices to decide what to do.

        let freq0 = self.sorted_allele_counts[1];
        let allele0_plus1 = self.sorted_allele_counts[0] + 1;
        let allele0_and_1_are_consecutive = allele0_plus1 == self.sorted_allele_counts[2];
        let mut new_sorted_allele_counts = Vec::new();

        // The rest of the sorted allele counts array contains junk
        let sorted_allele_count_lengths = self.distinct_allele_count << 1;

        if freq0 == 1 {
            // in this case allele0 won't be present in the result and all its frequency should go to allele0 + 1.
            if allele0_and_1_are_consecutive {  // need just to remove the first allele and 1 to the frequency of the second (freq1 += 1).
                new_sorted_allele_counts = self.sorted_allele_counts[2..sorted_allele_count_lengths as usize].to_vec();
                new_sorted_allele_counts[1] += 1;
            } else {
                new_sorted_allele_counts = self.sorted_allele_counts.to_vec();
                new_sorted_allele_counts.resize(sorted_allele_count_lengths as usize, 0);
                new_sorted_allele_counts[0] = allele0_plus1;
            }
        } else {
            if allele0_and_1_are_consecutive { // we don't need to add a component for allele0 + 1 since it already exists.
                new_sorted_allele_counts = self.sorted_allele_counts.clone();
                new_sorted_allele_counts[0] = 0;
                new_sorted_allele_counts[1] = freq0 - 1;
                new_sorted_allele_counts[3] += 1;
            } else { // we need to insert allele0 + 1 in the sorted-allele-counts array.
                new_sorted_allele_counts = vec![0; (sorted_allele_count_lengths + 2) as usize];
                new_sorted_allele_counts[0] = 0;
                new_sorted_allele_counts[1] = freq0 - 1;
                new_sorted_allele_counts[2] = allele0_plus1;
                new_sorted_allele_counts[3] += 1; // = 1 as the array was freshly created with 0s.
                new_sorted_allele_counts[4..(sorted_allele_counts_length - 2) as usize] = self.sorted_allele_counts[2..(sorted_allele_counts_length - 2) as usize];
            }
        }
        return GenotypeAlleleCounts::build(ploidy, index + 1, &new_sorted_allele_counts[..])
    }

    /**
     * Instantiates the first genotype possible provided a total ploidy.
     * @param ploidy the ploidy of the genotype.
     *
     * @throws Error if ploidy is less than 0.
     *
     * @return never {@code null}.
     */
    pub fn first(ploidy: usize) -> GenotypeAlleleCounts {
        if ploidy < 0 {
            panic!("Ploidy must be >= 0");
        }

        if ploidy == 0 {
            return GenotypeAlleleCounts::build(0,0, &[])
        } else {
            return GenotypeAlleleCounts::build(ploidy, 0, &[0, ploidy])
        };
    }

    /**
     * Returns the allele counts for each allele index to maximum.
     * @param maximumAlleleIndex the maximum allele index required.
     * @throws IllegalArgumentException if {@code maximumAlleleIndex} is less than 0.
     * @return never {@code null}, an array of exactly {@code maximumAlleleIndex + 1} positions with the counts
     * of each allele where the position in the array is equal to its index.
     */
    pub fn allele_counts_by_index(&self, maximum_allele_index: usize) -> Vec<i32> {
        if maximum_allele_index < 0 {
            panic!("The requested allele count cannot be less than 0")
        } else {
            let mut result = vec![0; maximum_allele_index + 1];
            self.copy_allele_counts_by_index(&mut result, 0, 0, maximum_allele_index);

            return result
        }
    }

    /**
     * Returns the index of the allele from its rank in the genotype.
     *
     * @param rank the query rank.
     *
     * @throws IllegalArgumentException if the {@code rank} provided is outside the valid range [0,{@link #distinctAlleleCount()}).
     *
     * @return 0 or greater.
     */
    pub fn allele_index_at(&self, rank: usize) -> usize {
        if rank >= self.distinct_allele_count {
            panic!("The requested rank {} is out of range [0, {})", rank, self.distinct_allele_count);
        }

        self.sorted_allele_counts[rank << 1]
    }

    fn copy_allele_counts_by_index(
        &self,
        dest: &mut Vec<i32>,
        offset: i32,
        minimum_allele_index: usize,
        maximum_allele_index: usize,
    ) {
        // First we determine what section of the sortedAlleleCounts array contains the counts of interest,
        // By the present allele rank range of interest.
        let minimum_allele_rank = self.allele_rank_for(minimum_allele_index);
        let maximum_allele_rank = self.allele_rank_for(maximum_allele_index);

        // If the min or max allele index are absent (returned rank < 0) we note where the would be inserted; that
        // way we avoid going through the rest of positions in the sortedAlleleCounts array.
        // The range of interest is then [startRank,endRank].
        let start_rank = if minimum_allele_index < 0 { -minimum_allele_rank - 1 } else { minimum_allele_rank };
        let end_rank = if maximum_allele_rank < 0 { -maximum_allele_rank - 2 } else { maximum_allele_rank };

        let mut next_index = minimum_allele_index; // next index that we want to output the count for.
        let mut next_rank = start_rank; // next rank to query in sortedAlleleCounts.
        let mut next_sorted_allele_counts_offset = next_rank << 1; // offset in sortedAlleleCounts where the info is present for the next rank.
        let mut next_dest_offset = offset; // next offset in destination array where to set the count for the nextIndex.

        while next_rank <= end_rank {
            next_rank += 1;
            let allele_index = self.sorted_allele_counts[next_sorted_allele_counts_offset];
            next_sorted_allele_counts_offset += 1;
            // fill non-present allele counts with 0s.
            while allele_index > next_index {
                dest[next_dest_offset] = 0;
                next_dest_offset += 1;
                next_index += 1;
            }

            // It is guaranteed that at this point alleleIndex == nextIndex
            // thanks to the condition of the enclosing while: there must be at least one index of interest that
            // is present in the remaining (nextRank,endRank] interval as otherwise endRank would be less than nextRank.
            dest[next_dest_offset] = self.sorted_allele_counts[next_sorted_allele_counts_offset];
            next_dest_offset += 1;
            next_sorted_allele_counts_offset += 1;
            next_index += 1;
        }

        // Finally we take care of trailing requested allele indices.
        while next_index <= maximum_allele_index {
            next_index += 1;
            dest[next_dest_offset] = 0;
            next_dest_offset += 1;
        }
    }

    /**
     * Returns the rank of an allele in the genotype by its index.
     *
     * @param index the target index.
     *
     * @throws IllegalArgumentException if {@code index} is less that 0. Indices can be arbitrarily large.
     *
     * @return -1 or less if the allele index is not present in the genotype, 0 to {@link #distinctAlleleCount()} - 1 otherwise.
     *   If negative, the absolute value can be used to determine where would be that index inserted within {@code [0,{@link #distinctAlleleCount()}]} as
     *   {@code - result - 1}.
     *
     */
    pub fn allele_rank_for(&self, index: usize) -> i64 {
        return self.allele_index_to_rank(index, 0, self.distinct_allele_count as usize)
    }

    /**
     * Returns the count of an allele in the genotype given is rank in the genotype (not the allele index itself).
     *
     * @param rank of the requested allele within the genotype.
     *
     * @throws IllegalArgumentException if {@code rank} is out the the valid range [0,{@link #distinctAlleleCount})
     *
     * @return 1 or greater.
     */
    pub fn allele_count_at(&self, rank: usize) -> usize {
        if rank >= self.distinct_allele_count {
            panic!("The rank is out of range")
        }
        return self.sorted_allele_counts[(rank << 1) + 1]
    }

    /**
     * Returns the count of an allele in the genotype given it index.
     *
     * @return 0 if the allele is not present in the genotype, 1 or more otherwise.
     */
    pub fn allele_count_for(&self, index: usize) -> usize {
        let rank = self.allele_rank_for(index);
        if rank < 0 {
            0
        } else {
            self.allele_count_at(rank)
        }
    }

    /**
     * Implements binary search across allele indexes.
     * @param index the target index.
     * @param from first inclusive possible rank.
     * @param to last exclusive possible rank.
     * @return -1 or less if the allele index is not in the genotype false otherwise. You can obtain
     *  the potential insertion point (within the interval [from,to]) as {@code -result - 1}
     */
    fn allele_index_to_rank(&self, index: usize, from: usize, to: usize) -> i64 {
        if to <= from {
            return -(from as i64) - 1
        }
        if from == to - 1 {
            let only_index = self.sorted_allele_counts[from << 1] as usize;
            if only_index == index { from } else {
                if only_index > index { -(from as i64) - 1 } else { -(to as i64) - 1 }
            }
        }

        let mid = (to + from) >> 1;
        let mid_index = self.sorted_allele_counts[mid << 1] as usize;

        if mid_index == index {
            return mid as i64
        } else if mid_index < index {
            return self.allele_index_to_rank(index, mid + 1, to);
        } else {
            return self.allele_index_to_rank(index, 0, mid);
        }
    }

    pub fn distinct_allele_count(&self) -> usize {
        self.distinct_allele_count
    }
}