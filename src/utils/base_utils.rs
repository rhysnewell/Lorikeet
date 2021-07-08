use std::cmp::{Ordering, min};
use std::collections::HashSet;
use rayon::prelude::*;

lazy_static! {
    pub static ref BASES: HashSet<u8> = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8].into_par_iter().cloned().collect::<HashSet<u8>>();
    pub static ref BASE_CHARS: HashSet<char> = ['A', 'C', 'G', 'T'].into_par_iter().cloned().collect::<HashSet<char>>();

    pub static ref BASES_EXTENDED: HashSet<u8> = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8, 'N' as u8, 'D' as u8].into_par_iter().cloned().collect::<HashSet<u8>>();
    pub static ref BASE_CHARS_EXTENDED: HashSet<char> = ['A', 'C', 'G', 'T', 'N', 'D'].into_par_iter().cloned().collect::<HashSet<char>>();
}

pub enum Base {
    A(char),
    C(char),
    G(char),
    T(char),
    N(char),
    D(char)
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
                return comparison
            }
        }
        return o1.len().cmp(&o2.len())
    }
}