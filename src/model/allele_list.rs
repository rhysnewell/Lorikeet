use model::variants::Allele;

pub trait AlleleList<T: Sized> {
    fn new(input_list: Vec<T>) -> AlleleList<T>;

    /**
     * Returns the number of alleles in this AlleleList.
     */
    fn number_of_alleles(&self) -> usize;

    /**
     * Returns the index of the given Allele in this AlleleList.
     * Returns a negative number if the given allele is not present in this AlleleList.
     * @throws IllegalArgumentException if allele is null.
     */
    fn index_of_allele(&self, input: T) -> usize;

    /**
     * Returns the allele at the given index in this AlleleList.
     * @throws IllegalArgumentException if index is negative or equal
     * to or higher than the number of elements in this AlleleList {@link AlleleList#numberOfAlleles()}).
     */
    fn get_allele(&self, index: usize) -> T;

    /**
     * Returns <code>true</code> if this AlleleList contains the specified allele
     * and <code>false</code> otherwise.
     */
    fn contains_allele(&self, input: T) -> bool;

    /**
     * Returns an unmodifiable empty allele-list.
     * @param <A> the allele class.
     * @return never {@code null}.
     */
    fn empty_allele_list() -> AlleleList<T>;

    /**
     * Checks whether two allele lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    fn equals(first: AlleleList<T>, second: AlleleList<T>) -> bool;

    /**
     * Resolves the index of the reference allele in an allele-list.
     *
     * <p>
     *     If there is no reference allele, it returns -1. If there is more than one reference allele,
     *     it returns the first occurrence (lowest index).
     * </p>
     *
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return -1 if there is no reference allele, or a values in [0,{@code list.alleleCount()}).
     */
    fn index_of_reference(&self) -> usize;

    /**
     * Returns a {@link List} unmodifiable view of this allele-list
     *
     * @return never {@code null}.
     */
    fn as_list_of_alleles(&self) -> Vec<T>;

    /**
     * Returns a permutation between two allele lists.
     * @param target the target allele list.
     *
     * @throws IllegalArgumentException if {@code target} is {@code null}, or
     * elements in {@code target} is not contained in {@code this}
     *
     * @return never {@code null}
     */
    fn permuation(&self, target: AlleleList<T>) -> AlleleListPermutation<T>;



}

pub trait AlleleListPermutation<T> {
}

pub struct NonPermutation<T> {
    allele_list: Vec<T>
}

pub struct ActualPermuation<T> {
    from: Vec<T>,
    to: Vec<T>,
    from_index: Vec<usize>,
    kept_from_indices: Vec<bool>,
    non_permuted: bool,
    is_partial: bool,
}