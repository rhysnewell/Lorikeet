use itertools::Itertools;
use rayon::prelude::*;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use std::collections::hash_map::DefaultHasher;
use std::fmt::Debug;
use std::hash::{Hash, Hasher};

/**
 * A graph vertex that holds some sequence information
 */
pub trait BaseVertex: Debug + Clone + Send + Sync + Eq + PartialEq + Hash {
    fn is_empty(&self) -> bool;

    fn len(&self) -> usize;

    fn to_string(&self) -> String;

    fn get_sequence(&self) -> &[u8];

    fn get_sequence_string(&self) -> &str;

    fn get_additional_sequence(&self, source: bool) -> &[u8];

    fn set_additional_info(&mut self, info: String);

    fn get_additional_info(&self) -> String;

    fn has_ambiguous_sequence(&self) -> bool;

    fn merge_identical_nodes(&self) -> bool {
        false
    }
}

impl BaseVertex for MultiDeBruijnVertex {
    // /**
    //  * Create a new sequence vertex with sequence
    //  *
    //  * This code doesn't copy sequence for efficiency reasons, so sequence must absolutely not be modified
    //  * in any way after passing this sequence to the BaseVertex
    //  *
    //  * @param sequence a non-null sequence of bases contained in this vertex
    //  */
    // pub fn new<'a>(sequence: &'a [u8]) -> Self {
    //     Self {
    //         sequence,
    //         additonal_info: format!(""),
    //         // node_index: 0,
    //     }
    // }

    fn merge_identical_nodes(&self) -> bool {
        self.merge_identical_nodes
    }

    /**
     * Does this vertex have an empty sequence?
     *
     * That is, is it a dummy node that's only present for structural reasons but doesn't actually
     * contribute to the sequence of the graph?
     *
     * @return true if sequence is empty, false otherwise
     */
    fn is_empty(&self) -> bool {
        self.sequence.len() == 0
    }

    /**
     * Get the length of this sequence
     * @return a positive integer >= 1
     */
    fn len(&self) -> usize {
        self.sequence.len()
    }

    fn to_string(&self) -> String {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        return format!(
            "MultiDeBruijnVertex_id_{}_seq_{}",
            hasher.finish(),
            std::str::from_utf8(self.sequence.as_slice()).unwrap()
        );
    }

    /**
     * Get the sequence of bases contained in this vertex
     *
     * Do not modify these bytes in any way!
     *
     * @return a non-null pointer to the bases contained in this vertex
     */
    fn get_sequence(&self) -> &[u8] {
        self.sequence.as_slice()
    }

    fn get_sequence_string(&self) -> &str {
        std::str::from_utf8(self.sequence.as_slice()).unwrap()
    }

    /**
     * Get the sequence unique to this vertex
     *
     * This function may not return the entire sequence stored in the vertex, as kmer graphs
     * really only provide 1 base of additional sequence (the last base of the kmer).
     *
     * The base implementation simply returns the sequence.
     *
     * @param source is this vertex a source vertex (i.e., no in nodes) in the graph
     * @return a byte[] of the sequence added by this vertex to the overall sequence
     */
    fn get_additional_sequence(&self, source: bool) -> &[u8] {
        if source {
            return self.sequence.as_slice();
        } else {
            return self.get_suffix_as_array();
        }
    }

    /**
     * Set additional debugging information for this vertex
     * @param info the new info value.
     */
    fn set_additional_info(&mut self, info: String) {
        self.additonal_info = info
    }

    /**
     * @return the additional information for display about this vertex
     */
    fn get_additional_info(&self) -> String {
        if self.reads.contains(&format!("ref")) {
            return self.additonal_info.to_string();
        } else {
            return format!(
                "{}{}",
                self.additonal_info,
                if Self::KEEP_TRACK_OF_READS {
                    format!("__{}", self.reads.iter().join(","))
                } else {
                    "".to_string()
                }
            );
        }
    }

    /**
     * Checks whether the vertex sequence is ambiguous or not.
     *
     * <p>
     *     Ambiguity may come about as a result of either:
     *     <ul>
     *        <li>by construction as the generating sequence (read or haplotype) had ambiguous bases</li>
     *        <li>or because this vertex is the result of merging two or more vertices with some variation upstream
     *        no more than kmerSize bases away.</li>
     *     </ul>
     * </p>
     *
     * @return {@code true} iff so.
     */
    fn has_ambiguous_sequence(&self) -> bool {
        return self.sequence.par_iter().any(|base| {
            let base = base.to_ascii_uppercase();
            let base = base as char;
            match base {
                'A' | 'T' | 'C' | 'G' => false,
                _ => true,
            }
        });
    }
}
