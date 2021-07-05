use rayon::prelude::*;
use std::hash::Hasher;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;

/**
 * A graph vertex that holds some sequence information
 */
pub trait BaseVertex {
    fn is_empty(&self) -> bool;

    fn len(&self) -> usize;

    fn to_string(&self) -> String;

    fn get_sequence(&self) -> &[u8];

    fn get_sequence_string(&self) -> String;

    fn get_additional_sequence(&self, source: bool) -> &[u8];

    fn set_additional_info(&mut self, info: String);

    fn get_additional_info(&self) -> &String;

    fn has_ambiguous_sequence(&self) -> bool;
}

impl BaseVertex for MultiDeBruijnVertex {
    /**
     * Create a new sequence vertex with sequence
     *
     * This code doesn't copy sequence for efficiency reasons, so sequence must absolutely not be modified
     * in any way after passing this sequence to the BaseVertex
     *
     * @param sequence a non-null sequence of bases contained in this vertex
     */
    // pub fn new<'a>(sequence: &'a [u8]) -> Self {
    //     Self {
    //         sequence,
    //         additonal_info: format!(""),
    //         // node_index: 0,
    //     }
    // }

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
        return format!("MultiDeBruijnVertex_is_{}_seq_{}", self.hash().into(), String::from_utf8_lossy(self.sequence).into_string())
    }

    /**
     * Get the sequence of bases contained in this vertex
     *
     * Do not modify these bytes in any way!
     *
     * @return a non-null pointer to the bases contained in this vertex
     */
    fn get_sequence(&self) -> &[u8] {
        self.sequence
    }

    fn get_sequence_string(&self) -> String {
        std::str::from_utf8(self.sequence).unwrap().to_string()
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
        self.sequence
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
    fn get_additional_info(&self) -> &String {
        if self.reads.contains(&format!("ref")) {
            return &self.additonal_info
        } else {
            return format!("{}{}", self.additonal_info,
                           if Self::KEEP_TRACK_OF_READS { format!("__{}", self.reads.iter().join(",")) } else { "" } )
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
        return self.sequence.par_iter().any(|base|{
            match base.to_ascii_uppercase() {
                'A' | 'T' | 'C' | 'G' => false,
                _ => true
            }
        })
    }
}