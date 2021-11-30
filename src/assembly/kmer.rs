use std::cmp::min;
use std::hash::{Hash, Hasher};

/**
 * Fast wrapper for byte[] kmers
 *
 * This objects has several important features that make it better than using a raw byte[] for a kmer:
 *
 * -- Can create kmer from a range of a larger byte[], allowing us to avoid Array.copyOfRange
 * -- Fast equals and hashcode methods
 * -- can get actual byte[] of the kmer, even if it's from a larger byte[], and this operation
 *    only does the work of that operation once, updating its internal state
 */
#[derive(Debug, Clone)]
pub struct Kmer {
    // this values may be updated in the course of interacting with this kmer
    pub bases: Vec<u8>,
    start: usize,
    // two constants
    length: usize,
}

// TODO: Change Kmer to take a reference to a sequence and have a lifetime
impl Kmer {
    /**
     * Create a new kmer using all bases in kmer
     * @param kmer a non-null byte[]. The input array must not be modified by the caller.
     */
    pub fn new(kmer: Vec<u8>) -> Kmer {
        Kmer {
            start: 0,
            length: kmer.len(),
            bases: kmer,
        }
    }

    /**
     * Create a new kmer backed by the bases in bases, spanning start -> start + length
     *
     * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
     * of bases for performance reasons
     *
     * @param bases an array of bases
     * @param start the start of the kmer in bases, must be >= 0 and < bases.length
     * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
     */
    pub fn new_with_start_and_length(bases: Vec<u8>, start: usize, length: usize) -> Kmer {
        Kmer {
            bases,
            start,
            length,
        }
    }

    /**
     * Create a derived shallow kmer that starts at newStart and has newLength bases
     * @param newStart the new start of kmer, where 0 means that start of the kmer, 1 means skip the first base
     * @param newLength the new length
     * @return a new kmer based on the data in this kmer.  Does not make a copy, so shares most of the data
     */
    pub fn sub_kmer(&self, new_start: usize, new_length: usize) -> Kmer {
        Kmer {
            bases: self.bases.to_vec(),
            start: self.start + new_start,
            length: new_length,
        }
    }

    /**
     * Get the bases of this kmer.  May create a copy of the bases, depending on how this kmer was constructed.
     *
     * Note that this function is efficient in that if it needs to copy the bases this only occurs once.
     *
     * @return a non-null byte[] containing length() bases of this kmer, regardless of how this kmer was created
     */
    pub fn bases(&mut self) -> &[u8] {
        if self.start != 0 || self.bases.len() != self.length {
            // update operation.  Rip out the exact byte[] and update start so we don't ever do this again
            self.bases =
                self.bases[self.start..min((self.start + self.length), self.bases.len())].to_vec();
            self.start = 0;
        }

        return &self.bases;
    }

    pub fn len(&self) -> usize {
        self.length
    }

    /**
     * Gets a set of differing positions and bases from another k-mer, limiting up to a max distance.
     * For example, if this = "ACATT" and other = "ACGGT":
     * - if maxDistance < 2 then -1 will be returned, since distance between kmers is 2.
     * - If maxDistance >=2, then 2 will be returned, and arrays will be filled as follows:
     * differingIndeces = {2,3}
     * differingBases = {'G','G'}
     * @param other                 Other k-mer to test
     * @param maxDistance           Maximum distance to search. If this and other k-mers are beyond this Hamming distance,
     *                              search is aborted and -1 is returned
     * @param differingIndeces      Array with indices of differing bytes in array
     * @param differingBases        Actual differing bases
     * @return                      Set of mappings of form (int->byte), where each elements represents index
     *                              of k-mer array where bases mismatch, and the byte is the base from other kmer.
     *                              If both k-mers differ by more than maxDistance, returns null
     */
    pub fn get_differing_positions(
        &self,
        other: &Self,
        max_distance: usize,
        differing_indices: &mut Vec<usize>,
        differing_bases: &mut Vec<u8>,
    ) -> i32 {
        let mut dist = 0;
        if self.length == other.length {
            let f2 = &other.bases;
            for i in 0..self.length {
                if self.bases[self.start + i] != f2[i] {
                    differing_indices[dist] = i;
                    differing_bases[dist + 1] = f2[i];
                    dist += 1;
                    if dist > max_distance {
                        return -1;
                    }
                }
            }
        }
        return dist as i32;
    }

    pub fn to_string(&self) -> String {
        return format!(
            "Kmer{{{}}}",
            format!(
                "{}{}{}",
                std::str::from_utf8(&self.bases).unwrap(),
                self.start,
                self.length
            )
        );
    }
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.bases[self.start..min((self.start + self.length), self.bases.len())]
            == other.bases[other.start..min((other.start + other.length), other.bases.len())]
    }
}

impl Eq for Kmer {}

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bases[self.start..min((self.start + self.length), self.bases.len())].hash(state)
    }
}
