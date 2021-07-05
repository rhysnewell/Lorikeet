use graphs::base_vertex::BaseVertex;

/**
 * A graph vertex containing a sequence of bases and a unique ID that
 * allows multiple distinct nodes in the graph to have the same sequence.
 *
 * This is essential when thinking about representing the actual sequence of a haplotype
 * in a graph.  There can be many parts of the sequence that have the same sequence, but
 * are distinct elements in the graph because they have a different position in the graph.  For example:
 *
 * A -> C -> G -> A -> T
 *
 * The two As are not the same, because they occur with different connections.  In a kmer graph equals()
 * is based on the sequence itself, as each distinct kmer can only be represented once.  But the transformation
 * of the kmer graph into a graph of base sequences, without their kmer prefixes, means that nodes that
 * where once unique including their prefix can become equal after shedding the prefix.  So we need to
 * use some mechanism -- here a unique ID per node -- to separate nodes that have the same sequence
 * but are distinct elements of the graph.
 *
 */
pub struct SeqVertex {
    sequence: &[u8],
    additional_info: String,
}

impl SeqVertex {
    pub fn new(sequence: &[u8]) -> Self {
        Self {
            sequence,
            additional_info: format!("")
        }
    }

    /**
     * Get the unique ID for this SeqVertex
     * @return a positive integer >= 0
     */
    pub fn get_id(&self) -> usize {
        self.hash()
    }

    pub fn to_string(&self) -> String {
        format!("SeqVertex_id_{}_seq_{}", self.hash().into(), std::str::from_utf8(self.sequence).unwrap())
    }

    /**
     * Return a new SeqVertex derived from this one but not including the suffix bases
     *
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex with appropriate prefix, or null if suffix removes all bases from this node
     */
    pub fn without_suffix(&self, suffix: &[u8]) -> Option<SeqVertex> {
        let prefix_size = self.sequence.len() - suffix.len();

        return if prefix_size > 0 { Self::new(self.sequence[0..suffix.len()]) } else { None }
    }

    /**
     * Return a new SeqVertex derived from this one but not including prefix or suffix bases
     *
     * @param prefix the previx bases to remove
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex
     */
    pub fn without_prefix_and_suffix(&self, prefix: &[u8], suffix: &[u8]) -> Option<SeqVertex> {
        let start = prefix.len();
        let length = self.sequence.len() - suffix.len() - prefix.len();
        let stop = start + length;

        return if length > 0 { SeqVertex::new(self.sequence[start..stop]) } else { None }
    }

}

impl BaseVertex for SeqVertex {
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

impl Hash for SeqVertex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.hash(state)
    }
}

impl PartialEq for SeqVertex {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

impl Eq for MultiDeBruijnVertex {}