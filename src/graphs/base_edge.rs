use rayon::prelude::*;

/**
 * Simple edge class for connecting nodes in the graph.
 *
 * Works equally well for all graph types (kmer or sequence)
 */
#[derive(Debug, Clone)]
pub struct BaseEdge {
    multiplicity: usize,
    is_ref: bool,
}

impl BaseEdge {
    /**
     * Create a new BaseEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge
     */
    pub fn new(is_ref: bool, multiplicity: usize) -> Self {
        Self {
            multiplicity,
            is_ref
        }
    }

    /**
     * Get the number of observations of paths connecting two vertices
     * @return a positive integer >= 0
     */
    pub fn get_multiplicity(&self) -> usize {
        self.multiplicity
    }

    /**
     * Get the DOT format label for this edge, to be displayed when printing this edge to a DOT file
     * @return a non-null string
     */
    pub fn get_dot_label(&self) -> String {
        self.multiplicity.to_string()
    }

    /**
     * Increase the multiplicity of this edge by incr
     * @param incr the change in this multiplicity, must be >= 0
     */
    pub fn inc_multiplicity(&mut self, incr: usize) {
        self.multiplicity += incr
    }

    /**
     * A special assessor that returns the multiplicity that should be used by pruning algorithm
     *
     * @return the multiplicity value that should be used for pruning
     */
    pub fn get_pruning_multiplicity(&self) -> usize {
        self.multiplicity
    }

    /**
     * Set the multiplicity of this edge to value
     * @param value an integer >= 0
     */
    pub fn set_multiplicity(&mut self, value: usize) {
        self.multiplicity = value
    }

    /**
     * Does this edge indicate a path through the reference graph?
     * @return true if so
     */
    pub fn is_ref(&self) -> bool {
        return self.is_ref
    }

    /**
     * Indicate that this edge follows the reference sequence, or not
     * @param isRef true if this is a reference edge
     */
    pub fn set_is_ref(&mut self, is_ref: bool) {
        self.is_ref = is_ref
    }

    /**
     * Add edge to this edge, updating isRef and multiplicity as appropriate
     *
     * isRef is simply the or of this and edge
     * multiplicity is the sum
     *
     * @param edge the edge to add
     * @return this
     */
    pub fn add(self, edge: Self) -> Self {
        self.multiplicity += edge.multiplicity;
        self.is_ref = self.is_ref || edge.is_ref;

        return self
    }

    /**
     * Create a new BaseEdge with the given multiplicity.
     * The resulting edge is a reference edge if any of the argument edges are reference.
     *
     * @param edges a collection of edges to or their isRef values
     * @param multiplicity our desired multiplicity
     * @return a newly allocated BaseEdge
     */
    pub fn make_o_r_edge(edges: Vec<Self>, multiplicity: bool) -> Self {
        assert!(!edges.is_empty(), "Edges cannot be empty");
        let is_ref = edges.par_iter().any(|e| e.is_ref());

        Self {
            multiplicity,
            is_ref
        }
    }

    pub fn to_string(&self) -> String {
        return format!("BaseEdge{{multiplicity={}, isRef={}}}", self.multiplicity, self.is_ref)
    }
}