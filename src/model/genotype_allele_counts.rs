/**
The following code is adapted from the broadinstitute GATK HaplotypeCaller program
*/

use ordered_float::{NotNan, OrderedFloat};
use ndarray::{Array, Array2, ArrayBase, OwnedRepr};


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
pub struct GenotypeAlleleCounts {
    ploidy: i32,
    distinct_allele_count: i32,
    sorted_allele_counts: Vec<i32>,
    index: usize,
    log10_combination_count: f64,
}

impl GenotypeAlleleCounts {
    pub fn build_empty() -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy: 0,
            index: 0,
            sorted_allele_counts: Vec::new(),
            distinct_allele_count: 0,
            log10_combination_count: -1.,
        }
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
    pub fn build(ploidy: i32, index: usize, sorted_allele_counts: &[i32]) -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy,
            index,
            sorted_allele_counts: Vec::from(sorted_allele_counts),
            distinct_allele_count: sorted_allele_counts.len() as i32 >> 1,
            log10_combination_count: -1.,
        }
    }

    pub fn build_with_count(
        ploidy: i32,
        index: usize,
        sorted_allele_counts: Vec<i32>,
        distinct_allele_count: i32
    ) -> GenotypeAlleleCounts {
        GenotypeAlleleCounts {
            ploidy,
            index,
            sorted_allele_counts,
            distinct_allele_count,
            log10_combination_count: -1.,
        }
    }

    pub fn ploidy(&self) -> i32 {
        self.ploidy
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
    pub fn first(ploidy: i32) -> GenotypeAlleleCounts {
        if ploidy < 0 {
            panic!("Ploidy must be >= 0");
        }

        if ploidy == 0 {
            return GenotypeAlleleCounts::build(0,0, &[])
        } else {
            return GenotypeAlleleCounts::build(ploidy, 0, &[0, ploidy])
        };
    }
}