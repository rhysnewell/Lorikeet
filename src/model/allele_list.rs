use model::variants::Allele;
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub struct AlleleList<T> {
    list: Vec<T>,
}

impl AlleleList<Allele> {
    pub fn new(input_list: &Vec<Allele>) -> AlleleList<Allele> {
        AlleleList {
            list: input_list.clone(),
        }
    }

    /**
     * Returns the number of alleles in this AlleleList.
     */
    pub fn number_of_alleles(&self) -> usize {
        self.list.len()
    }

    /**
     * Returns the index of the given Allele in this AlleleList.
     * Returns a negative number if the given allele is not present in this AlleleList.
     * @throws IllegalArgumentException if allele is null.
     */
    pub fn index_of_allele(&self, input: &Allele) -> usize {
        self.list.par_iter().position_first(|p| p == input).unwrap()
    }

    /**
     * Returns the allele at the given index in this AlleleList.
     * @throws IllegalArgumentException if index is negative or equal
     * to or higher than the number of elements in this AlleleList {@link AlleleList#numberOfAlleles()}).
     */
    pub fn get_allele(&self, index: usize) -> &Allele {
        &self.list[index]
    }

    /**
     * Returns <code>true</code> if this AlleleList contains the specified allele
     * and <code>false</code> otherwise.
     */
    pub fn contains_allele(&self, input: &Allele) -> bool {
        self.list.contains(input)
    }

    /**
     * Returns an unmodifiable empty allele-list.
     * @param <A> the allele class.
     * @return never {@code null}.
     */
    pub fn empty_allele_list() -> AlleleList<Allele> {
        AlleleList { list: Vec::new() }
    }

    /**
     * Checks whether two allele lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    pub fn equals(first: &AlleleList<Allele>, second: &AlleleList<Allele>) -> bool {
        first == second
    }

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
    pub fn index_of_reference(&self) -> usize {
        self.list
            .par_iter()
            .position_first(|p| p.is_reference())
            .unwrap()
    }

    /**
     * Returns a {@link List} unmodifiable view of this allele-list
     *
     * @return never {@code null}.
     */
    pub fn as_list_of_alleles(&self) -> &Vec<Allele> {
        &self.list
    }

    /**
     * Returns a permutation between two allele lists.
     * @param target the target allele list.
     *
     * @throws IllegalArgumentException if {@code target} is {@code null}, or
     * elements in {@code target} is not contained in {@code this}
     *
     * @return never {@code null}
     */
    pub fn permuation(self, target: AlleleList<Allele>) -> Permutation<Allele> {
        Permutation::new(self, target)
    }
}

pub trait AlleleListPermutation<T> {
    fn is_partial(&self) -> bool;

    fn is_non_permuted(&self) -> bool;

    fn to_index(&self, from_index: usize) -> usize;

    fn from_index(&self, to_index: usize) -> usize;

    fn is_kept(&self, from_index: usize) -> bool;

    fn from_size(&self) -> usize;

    fn to_size(&self) -> usize;

    fn from_list(&self) -> &Vec<T>;

    fn to_list(&self) -> &Vec<T>;

    fn number_of_alleles(&self) -> usize;

    fn index_of_allele(&self, allele: &T) -> usize;

    fn get_allele(&self, index: usize) -> &T;
}

pub enum Permutation<T> {
    NonPermutation {
        allele_list: AlleleList<T>,
    },
    ActualPermutation {
        from: AlleleList<T>,
        to: AlleleList<T>,
        from_index: Vec<usize>,
        kept_from_indices: Vec<bool>,
        non_permuted: bool,
        is_partial: bool,
    },
}

impl Permutation<Allele> {
    pub fn new(original: AlleleList<Allele>, target: AlleleList<Allele>) -> Permutation<Allele> {
        if AlleleList::equals(&original, &target) {
            return Permutation::NonPermutation {
                allele_list: original,
            };
        } else {
            let mut kept_from_indices = vec![false; original.number_of_alleles()];
            let to_size = target.number_of_alleles();
            let from_size = original.number_of_alleles();

            if from_size < to_size {
                panic!("Target allele list is not a permutation of the original allele list")
            }

            let mut from_index = vec![0; to_size];
            let mut non_permuted = from_size == to_size;
            let is_partial = !non_permuted;

            for i in 0..to_size {
                let original_index = original.index_of_allele(target.get_allele(i));
                if original_index < 0 {
                    panic!("Target allele is not a permutation of the original allele list")
                }
                kept_from_indices[original_index] = true;
                from_index[i] = original_index;
                non_permuted = non_permuted & (original_index == i);
            }

            return Permutation::ActualPermutation {
                from: original,
                to: target,
                from_index,
                kept_from_indices,
                is_partial,
                non_permuted,
            };
        }
    }
}

impl AlleleListPermutation<Allele> for Permutation<Allele> {
    fn is_partial(&self) -> bool {
        match self {
            Permutation::NonPermutation { .. } => false,
            Permutation::ActualPermutation { is_partial, .. } => *is_partial,
        }
    }

    fn is_non_permuted(&self) -> bool {
        match self {
            Permutation::NonPermutation { .. } => true,
            Permutation::ActualPermutation { non_permuted, .. } => *non_permuted,
        }
    }

    fn to_index(&self, from_index: usize) -> usize {
        match self {
            Permutation::NonPermutation { .. } => from_index,
            Permutation::ActualPermutation { to, from, .. } => {
                to.index_of_allele(from.get_allele(from_index))
            }
        }
    }

    fn from_index(&self, to_index: usize) -> usize {
        match self {
            Permutation::NonPermutation { .. } => to_index,
            Permutation::ActualPermutation { from_index, .. } => from_index[to_index],
        }
    }

    fn is_kept(&self, from_index: usize) -> bool {
        match self {
            Permutation::NonPermutation { .. } => true,
            Permutation::ActualPermutation {
                kept_from_indices, ..
            } => kept_from_indices[from_index],
        }
    }

    fn from_size(&self) -> usize {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.number_of_alleles(),
            Permutation::ActualPermutation { from, .. } => from.number_of_alleles(),
        }
    }

    fn to_size(&self) -> usize {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.number_of_alleles(),
            Permutation::ActualPermutation { to, .. } => to.number_of_alleles(),
        }
    }

    fn from_list(&self) -> &Vec<Allele> {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.as_list_of_alleles(),
            Permutation::ActualPermutation { from, .. } => from.as_list_of_alleles(),
        }
    }

    fn to_list(&self) -> &Vec<Allele> {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.as_list_of_alleles(),
            Permutation::ActualPermutation { to, .. } => to.as_list_of_alleles(),
        }
    }

    fn number_of_alleles(&self) -> usize {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.number_of_alleles(),
            Permutation::ActualPermutation { to, .. } => to.number_of_alleles(),
        }
    }

    fn index_of_allele(&self, allele: &Allele) -> usize {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.index_of_allele(allele),
            Permutation::ActualPermutation { to, .. } => to.index_of_allele(allele),
        }
    }

    fn get_allele(&self, index: usize) -> &Allele {
        match self {
            Permutation::NonPermutation { allele_list } => allele_list.get_allele(index),
            Permutation::ActualPermutation { to, .. } => to.get_allele(index),
        }
    }
}
