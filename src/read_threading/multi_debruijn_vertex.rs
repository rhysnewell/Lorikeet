use graphs::base_vertex::BaseVertex;

lazy_static! {
    static mut suffices_as_byte_array: Vec<Vec<u8>> = vec![Vec::new(); std::u8::MAX - std::u8::MIN + 1];
    static {
        for i in 0..(*suffices_as_byte_array.len() as u16) {
            suffices_as_byte_array[i] = i.to_be_bytes();
        }
    }
}

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 */
pub struct MultiDeBruijnVertex {
    pub(crate) reads: Vec<String>,
    merge_identical_nodes: bool,
    pub(crate) additonal_info: String,
    pub sequence: &[u8],
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
    pub fn new(sequence: &[u8], merge_identical_nodes: bool) -> Self {
        let base_vertex = BaseVertex::new(sequence);

        Self {
            sequence,
            merge_identical_nodes,
            reads: Vec::new(),
            additonal_info: format!(""),
        }
    }

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param sequence the kmer sequence
     */
    pub fn new_with_sequence(seqeunce: &[u8]) -> Self {
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

    pub fn get_additional_info(&self) {
        if self.reads.contains(&format!("ref")) {
            return self.base_vertex.get_additional_info()
        } else {
            return format!("{}{}", self.base_vertex.get_additional_info(),
                           if Self::KEEP_TRACK_OF_READS { format!("__{}", self.reads.iter().join(",")) } else { "" } )
        }
    }

    /**
     * Get the kmer size for this DeBruijnVertex
     * @return integer >= 1
     */
    pub fn get_kmer_size(&self) -> usize {
        return self.base_vertex.get_sequence().len()
    }

    /**
     * Get the string representation of the suffix of this DeBruijnVertex
     * @return a non-null non-empty string
     */
    pub fn get_suffix_string(&self) -> String {
        return format!("{}", std::str::from_utf8(self.get_suffix_as_array()).unwrap())
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
        return self.base_vertex.sequence[self.get_kmer_size() - 1]
    }

    /**
     * Optimized version that returns a byte[] for the single byte suffix of this graph without allocating memory.
     *
     * Should not be modified
     *
     * @return a byte[] that contains 1 byte == getSuffix()
     */
    fn get_suffix_as_array(&self) -> &[u8] {
        return &suffices_as_byte_array[self.get_suffix()][..]
    }
}

impl Hash for MultiDeBruijnVertex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.hash(state)
    }
}

impl PartialEq for MultiDeBruijnVertex {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

impl Eq for MultiDeBruijnVertex {}