use graphs::base_vertex::BaseVertex;
use itertools::Itertools;
use rand::rngs::ThreadRng;
use rand::Rng;
use rayon::prelude::*;
use std::hash::{Hash, Hasher};

lazy_static! {
    static ref suffices_as_byte_array: Vec<Vec<u8>> = (0..=std::u8::MAX)
        .into_par_iter()
        .map(|i| { vec![i] })
        .collect::<Vec<Vec<u8>>>();
}

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 */
#[derive(Debug, Clone)]
pub struct MultiDeBruijnVertex {
    pub(crate) reads: Vec<String>,
    pub(crate) merge_identical_nodes: bool,
    pub(crate) additonal_info: String,
    pub sequence: Vec<u8>,
    identity_code: usize,
}

/**
 * Create a new MultiDeBruijnVertex with kmer sequence
 * @param mergeIdenticalNodes should nodes with the same sequence be treated as equal?
 * @param sequence the kmer sequence
 */
impl MultiDeBruijnVertex {
    pub const KEEP_TRACK_OF_READS: bool = false;

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param mergeIdenticalNodes should nodes with the same sequence be treated as equal?
     * @param sequence the kmer sequence
     */
    pub fn new(sequence: Vec<u8>, merge_identical_nodes: bool) -> MultiDeBruijnVertex {
        let mut rng = ThreadRng::default();

        MultiDeBruijnVertex {
            sequence,
            merge_identical_nodes,
            reads: Vec::new(),
            additonal_info: format!(""),
            identity_code: rng.gen(),
        }
    }

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param sequence the kmer sequence
     */
    pub fn new_with_sequence(seqeunce: Vec<u8>) -> MultiDeBruijnVertex {
        Self::new(seqeunce, false)
    }

    /**
     * Add name information to this vertex for debugging
     *
     * This information will be captured as a list of strings, and displayed in DOT if this
     * graph is written out to disk
     *
     * This functionality is only enabled when KEEP_TRACK_OF_READS is true
     *
     * @param name a non-null string
     */
    pub fn add_read(&mut self, name: String) {
        if Self::KEEP_TRACK_OF_READS {
            self.reads.push(name)
        }
    }

    pub fn get_additional_info(&self) -> String {
        if self.reads.contains(&format!("ref")) {
            return self.additonal_info.clone();
        } else {
            return format!(
                "{}{}",
                self.additonal_info,
                if Self::KEEP_TRACK_OF_READS {
                    format!("__{}", self.reads.iter().join(","))
                } else {
                    format!("")
                }
            );
        }
    }

    /**
     * Get the kmer size for this DeBruijnVertex
     * @return integer >= 1
     */
    pub fn get_kmer_size(&self) -> usize {
        return self.get_sequence().len();
    }

    /**
     * Get the string representation of the suffix of this DeBruijnVertex
     * @return a non-null non-empty string
     */
    pub fn get_suffix_string(&self) -> String {
        return std::str::from_utf8(self.get_suffix_as_array())
            .unwrap()
            .to_string();
    }

    /**
     * Get the suffix byte of this DeBruijnVertex
     *
     * The suffix byte is simply the last byte of the kmer sequence, so if this is holding sequence ACT
     * getSuffix would return T
     *
     * @return a byte
     */
    pub fn get_suffix(&self) -> u8 {
        return self.sequence[self.get_kmer_size() - 1];
    }

    /**
     * Optimized version that returns a byte[] for the single byte suffix of this graph without allocating memory.
     *
     * Should not be modified
     *
     * @return a byte[] that contains 1 byte == getSuffix()
     */
    pub fn get_suffix_as_array(&self) -> &[u8] {
        return &*suffices_as_byte_array[self.get_suffix() as usize];
    }
}

impl Hash for MultiDeBruijnVertex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        if self.merge_identical_nodes {
            self.sequence.hash(state);
        } else {
            self.sequence.hash(state);
            self.identity_code.hash(state);
        }
    }
}

impl PartialEq for MultiDeBruijnVertex {
    fn eq(&self, other: &Self) -> bool {
        if self.merge_identical_nodes {
            self.sequence == other.sequence
        } else {
            false
        }
    }
}

impl Eq for MultiDeBruijnVertex {}
