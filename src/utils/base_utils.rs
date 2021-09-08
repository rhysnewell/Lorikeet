use bio::alphabets::dna::alphabet;
use rand::{rngs::ThreadRng, Rng};
use rayon::prelude::*;
use std::cmp::{min, Ordering};
use std::collections::HashSet;

lazy_static! {
    pub static ref BASES: HashSet<u8> = vec!['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8]
        .into_iter()
        .collect::<HashSet<u8>>();
    pub static ref BASE_CHARS: HashSet<char> = vec!['A', 'C', 'G', 'T']
        .into_iter()
        .collect::<HashSet<char>>();
    pub static ref BASES_EXTENDED: HashSet<u8> =
        vec!['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8, 'N' as u8, 'D' as u8]
            .into_iter()
            .collect::<HashSet<u8>>();
    pub static ref BASE_CHARS_EXTENDED: HashSet<char> = vec!['A', 'C', 'G', 'T', 'N', 'D']
        .into_iter()
        .collect::<HashSet<char>>();
}

pub enum Base {
    A(char),
    C(char),
    G(char),
    T(char),
    N(char),
    D(char),
}

impl Base {
    pub fn new(base: char) -> Base {
        match base {
            'A' => Base::A(base),
            'C' => Base::C(base),
            'G' => Base::G(base),
            'T' => Base::T(base),
            'N' => Base::N(base),
            'D' => Base::D(base),
            _ => panic!("Incompatible base char {}", base),
        }
    }
}

pub struct BaseUtils {}

impl BaseUtils {
    /**
     * Lexicographical sorting of base arrays {@link Comparator}.
     */
    pub fn bases_comparator(o1: &[u8], o2: &[u8]) -> Ordering {
        let min_length = min(o1.len(), o2.len());
        for i in 0..min_length {
            let comparison = o1[i].cmp(&o2[i]);
            if comparison != Ordering::Equal {
                return comparison;
            }
        }
        return o1.len().cmp(&o2.len());
    }

    /**
     * Converts a simple base to a base index
     *
     * @param base [AaCcGgTt]
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    // pub fn simple_base_to_base_index(base: u8) -> i32 {
    //
    // }

    /**
     * Returns true iff the base represented by the byte is a 'regular' base (ACGT or *).
     */
    pub fn is_regular_base(base: u8) -> bool {
        return alphabet().is_word(&[base]);
    }

    pub fn is_all_regular_base(bases: &[u8]) -> bool {
        bases.par_iter().all(|base| Self::is_regular_base(*base))
    }

    /**
     * Fill an array section with random bases.
     *
     * @param dest array to fill.
     * @param fromIndex first index to be filled (inclusive).
     * @param toIndex index after last to be filled (exclusive).
     *
     * @throws IllegalArgumentException if {@code dest} is {@code null},
     *              {@code fromIndex} or {@code toIndex} is negative,
     *              {@code fromIndex} or {@code toIndex} are greater than {@code dest} length,
     *              or {@code fromIndex} greater than {@code toIndex}.
     */
    pub fn fill_with_random_bases(dest: &mut Vec<u8>, from_index: usize, to_index: usize) {
        let mut rnd = ThreadRng::default();
        assert!(
            to_index <= dest.len(),
            "both indices must be less or equal to dest length from {} to {} dest length {}",
            from_index,
            to_index,
            dest.len()
        );
        (from_index..to_index)
            .for_each(|i| dest[i] = Self::base_index_to_simple_base(rnd.gen_range(0, 4)))
    }

    /**
     * Converts a base index to a simple base
     *
     * @param baseIndex 0, 1, 2, 3
     * @return A, C, G, T, or '.' if the index can't be understood
     */
    pub fn base_index_to_simple_base(base_index: usize) -> u8 {
        match base_index {
            0 => return b'A',
            1 => return b'C',
            2 => return b'G',
            3 => return b'T',
            _ => return b'.',
        }
    }
}
