use ndarray::{Array2, Array};
use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use genotype::genotype_allele_counts::GenotypeAlleleCounts;

pub struct GenotypeLikelihoodCalculators {
    pub ploidy: usize,
    pub allele_count: usize,
    pub allele_first_genotype_offset_by_ploidy: Array2<i32>,
    pub genotype_table_by_ploidy: Vec<Vec<GenotypeAlleleCounts>>,
}

impl GenotypeLikelihoodCalculators {

    pub const MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY: i32 = 1000;

    pub fn build_empty() -> GenotypeLikelihoodCalculators {
        GenotypeLikelihoodCalculators {
            ploidy: 0,
            allele_count: 0,
            allele_first_genotype_offset_by_ploidy: Array2::zeros([0, 0]),
            genotype_table_by_ploidy: Vec::new(),
        }
    }

    pub fn get_instance(ploidy: usize, allele_count: usize) -> GenotypeLikelihoodCalculator {
        let allele_first_offset_by_ploidy =
            GenotypeLikelihoodCalculators::calculate_genotype_counts_using_table_and_validate(ploidy, allele_count);
        let genotype_table_by_ploidy: Vec<Vec<GenotypeAlleleCounts>> = GenotypeLikelihoodCalculators::build_genotype_allele_counts_table(
            ploidy,
            allele_count,
            &allele_first_offset_by_ploidy,
        );
        return GenotypeLikelihoodCalculator::new(
            ploidy,
            allele_count,
            allele_first_offset_by_ploidy,
            genotype_table_by_ploidy
        )
    }

    fn calculate_genotype_counts_using_table_and_validate(ploidy: usize, allele_count: usize) -> Array2<i32> {
        GenotypeLikelihoodCalculators::check_ploidy_and_allele(ploidy, allele_count);

        let result = GenotypeLikelihoodCalculators::build_allele_first_genotype_offset_table(ploidy, allele_count);

        if result[[ploidy, allele_count]] == -1 {
            panic!("The number of genotypes is too large for ploidy {} and allele {}", ploidy, allele_count);
        }

        return result

    }

    fn check_ploidy_and_allele(ploidy: usize, allele_count: usize) {
        if ploidy < 0 {
            panic!("Ploidy cannot be negative");
        }

        if allele_count < 0 {
            panic!("The allele count cannot be negative")
        }
    }

    fn calculate_genotype_count_using_tables(ploidy: usize, allele_count: usize) -> Array2<i32> {
        let result =
            GenotypeLikelihoodCalculators::build_allele_first_genotype_offset_table(ploidy, allele_count);

        return result
    }


    /**
     * Build the table with the genotype offsets based on ploidy and the maximum allele index with representation
     * in the genotype.
     * <p>
     * The result is a matrix containing the offset of the first genotype that contain a particular allele
     * stratified by ploidy.
     * <p>
     *     Row (first dimension) represent the ploidy, whereas
     *     the second dimension represents the allele.
     * </p>
     *
     * <p>
     *     Thus the value a position <i>[p][a]</i> indicates how many genotypes of ploidy <i>p</i> there are before the first
     *     one that contains allele <i>a</i>. <br/>
     *
     *     For example, considering ploidy 3 and alleles A, B, C, D, etc ... (indexed 0, 1, 2, ... respectively):
     *     <br/>
     *     [3][A] == [3][0] == 0 as the first genotype AAA contains A.
     *     <br/>
     *     [3][C] == [3][2] == 4 as the first genotype that contains C, AAC follows: AAA AAB ABB BBB
     *     <br/>
     *     [4][D] == [4][3] == 14  as the first genotype that contains D, AAAD follows: AAAA AAAB AABB ABBB BBBB AAAC
     *     AABC ABBC BBBC AACC ABCC BBCC ACCC BCCC CCCC.
     *
     * </p>
     *
     * <p>
     *     This value are calculated recursively as follows:
     * </p>
     * <pre>
     *
     *     Offset[p][a] := Offset[p-1][a] + Offset[p][a-1] when a > 0, p > 0
     *                     0                               when a == 0
     *                     1                               otherwise
     *
     *
     *         0 1 1  1  1  1   1 ...
     *         0 1 2  3  4  5   6 ...
     *         0 1 3  6 10 15  21 ...
     *         0 1 4 10 20 35  56 ...
     *         0 1 5 15 35 70 126 ...
     *         0 ..................
     * </pre>
     *
     * <p>
     *    Note: if someone can come with a close form computable 0(1) (respect to ploidy and allele count)
     *     please let the author know.
     * </p>
     *
     * <p>
     *     The matrix is guaranteed to have as many rows as indicated by {@code maximumPloidy} + 1; the first
     *     row refers to the special case of ploidy == 0, the second row to ploidy 1 and so forth. Thus the ploidy
     *     matches the index.
     * </p>
     * <p>
     *     The matrix is guaranteed to have as many columns as indicate by {@code maximumAllele} + 1. In this case however
     *     the first allele index 0 is a sense allele (typically the reference allele). The reason to have at least the total
     *     genotype count up to allele count {@link @alleleCapacity} that is equal to the offset of the first genotype
     *     of the following allele; thus we need an extra one.
     * </p>
     *
     * <p>
     *     Although it might seem non-sense to have genotypes of ploidy 0. The values in the first row are used when
     *     filling up values in row 1 and so forth so it is present for programmatic convenience.
     *     Offsets in this row are 0 for the first column and 1 for any others.
     * </p>
     *
     * @param maximumPloidy maximum supported ploidy.
     * @param maximumAllele maximum supported allele index.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative.
     *
     * @return never {@code null}, the matrix described with enough information to address
     *       problems concerning up to the requested maximum allele index and ploidy.
     */
    pub fn build_allele_first_genotype_offset_table(
        ploidy: usize, allele_count: usize
    ) -> Array2<i32> {
        let row_count = ploidy + 1;
        let col_count = allele_count + 1;

        let mut result = Array::zeros([row_count, col_count]);

        result.slice_mut(s![0, 1..col_count]).fill(1);

        for p in 1..row_count {
            for a in 1..col_count {
                result[[p, a]] = result[[p, a - 1]] + result[[p - 1, a]];
                if result[[p, a]] < result[[p, a - 1]] {
                    result[[p, a]] = -1
                }
            }
        }

        return result
    }

    /**
     * Composes a table with the lists of all possible genotype allele counts given the the ploidy and maximum allele index.
     * <p>
     *     The resulting matrix has at least as many rows as {@code maximumPloidy } + 1 as the first row with index 0 correspond
     *     to ploidy == 0. Each row array has as many positions as necessary to contain all possible genotype-allele-counts in increasing order.
     *     This quantity varies with the ploidy.
     * </p>
     *
     * <p>
     *     Therefore <code>result[3][4]</code> would contain the 5th genotype with ploidy 3, and <code>result[4].length</code>
     *     would be equal to the count of possible genotypes for ploidy 4.
     * </p>
     *
     * @param maximumPloidy maximum ploidy to use in queries to the resulting table.
     * @param maximumAllele maximum allele index to use in queries to the resulting table.
     * @param offsetTable an allele first genotype offset table as constructed using {@link #buildAlleleFirstGenotypeOffsetTable(int, int)}
     *                    that supports at least up to {@code maximumAllele} and {@code maximumPloidy}.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative, or {@code offsetTable} is {@code null},
     *   or it does not have the capacity to handle the requested maximum ploidy or allele index.
     *
     * @return never {@code null}.
     */
    fn build_genotype_allele_counts_table(
        ploidy: usize, allele_count: usize, offset_table: &Array2<i32>
    ) -> Vec<Vec<GenotypeAlleleCounts>>{

        GenotypeLikelihoodCalculators::check_ploidy_and_allele(ploidy, allele_count);

        let row_count = ploidy + 1;
        let mut result = vec![Vec::new(); row_count as usize]; // each row has a different number of columns.
        for p in 0..ploidy {
            result[p] = GenotypeLikelihoodCalculators::build_genotype_allele_counts_array(p, allele_count, offset_table);
        }

        return result

    }

    /**
     * Builds a genotype-allele-counts array given the genotype ploidy and how many genotype you need.
     * <p>
     *     The result is guarantee to have exactly {@code length} positions and the elements are sorted
     *     in agreement with the standard way to display genotypes following the VCF standard.
     * </p>
     *
     * <p> Notice that is possible to request ploidy ==0. In that case the resulting array will have repetitions
     * of the empty genotype allele count.
     * </p>
     *
     * <p>
     *     For example,
     *
     *     <pre>
     *         ploidy = 1, length = 5 : [ {A}, {B}, {C}, {D}, {E} ]
     *         ploidy = 2, length = 7 : [ {AA}, {AB}, {BB}, {AC}, {BC}, {CC}, {AD}
     *         ploidy = 3, length = 10 : [ {AAA}, {AAB}, {ABB}, {BBB}, {AAC}, {ABC}, {BBC}, {BCC}, {CCC}, {AAD} ]
     *     </pre>
     * </p>
     *
     * @param ploidy requested ploidy.
     * @param alleleCount number of different alleles that the genotype table must support.
     * @param genotypeOffsetTable table with the offset of the first genotype that contain an allele given
     *                            the ploidy and its index.
     *
     * @return never {@code null}, follows the specification above.
     */
    fn build_genotype_allele_counts_array(
        ploidy: usize, allele_count: usize,
        offset_table: &Array2<i32>
    ) -> Vec<GenotypeAlleleCounts> {
        let length = offset_table[[ploidy, allele_count]];
        let strong_ref_length = if length == -1 {
            GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY
        } else {
            std::cmp::min(length, GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY)
        };

        let mut result = vec![GenotypeAlleleCounts::build_empty()];
        result[0] = GenotypeAlleleCounts::first(ploidy);

        for genotype_index in 0..strong_ref_length as usize {
            result[genotype_index] = result[genotype_index - 1].next();
        }

        return result
    }


}