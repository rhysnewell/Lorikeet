use rayon::prelude::*;
use std::fmt::Debug;

/**
 * Simple edge class for connecting nodes in the graph.
 *
 * Works equally well for all graph types (kmer or sequence)
 */
pub trait BaseEdge: Debug + Clone + Send + Sync + Eq + PartialEq {
    fn get_multiplicity(&self) -> usize;

    fn get_dot_label(&self) -> String;

    fn inc_multiplicity(&mut self, incr: usize);

    fn get_pruning_multiplicity(&self) -> usize;

    fn set_multiplicity(&mut self, value: usize);

    fn is_ref(&self) -> bool;

    fn set_is_ref(&mut self, is_ref: bool);

    fn add(&mut self, edge: Self);

    fn make_o_r_edge(edges: Vec<Self>, multiplicity: usize, single_sample_capacity: usize) -> Self
    where
        Self: std::marker::Sized;

    fn to_string(&self) -> String;
}

/**
* The most basic implementation of a BaseEdge like object. Only meant as a placeholder for certain
* functions
*/
#[derive(Debug, Clone, Copy, Eq, Ord, PartialOrd, PartialEq)]
pub struct BaseEdgeStruct {
    pub(crate) multiplicity: usize,
    pub(crate) is_ref: bool,
}

impl BaseEdgeStruct {
    pub fn new(is_ref: bool, multiplicity: usize) -> Self {
        Self {
            is_ref,
            multiplicity,
        }
    }

    pub fn set(&mut self, is_ref: bool, multiplicity: usize) {
        self.is_ref = is_ref;
        self.multiplicity = multiplicity;
    }
}

impl BaseEdge for BaseEdgeStruct {
    /**
     * Get the number of observations of paths connecting two vertices
     * @return a positive integer >= 0
     */
    fn get_multiplicity(&self) -> usize {
        self.multiplicity
    }

    /**
     * Get the DOT format label for this edge, to be displayed when printing this edge to a DOT file
     * @return a non-null string
     */
    fn get_dot_label(&self) -> String {
        return self.multiplicity.to_string();
    }

    /**
     * Increase the multiplicity of this edge by incr
     * @param incr the change in this multiplicity, must be >= 0
     */
    fn inc_multiplicity(&mut self, incr: usize) {
        self.multiplicity += incr
    }

    /**
     * A special assessor that returns the multiplicity that should be used by pruning algorithm
     *
     * @return the multiplicity value that should be used for pruning
     */
    fn get_pruning_multiplicity(&self) -> usize {
        self.multiplicity
    }

    /**
     * Set the multiplicity of this edge to value
     * @param value an integer >= 0
     */
    fn set_multiplicity(&mut self, value: usize) {
        self.multiplicity = value
    }

    /**
     * Does this edge indicate a path through the reference graph?
     * @return true if so
     */
    fn is_ref(&self) -> bool {
        return self.is_ref;
    }

    /**
     * Indicate that this edge follows the reference sequence, or not
     * @param isRef true if this is a reference edge
     */
    fn set_is_ref(&mut self, is_ref: bool) {
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
    fn add(&mut self, edge: Self) {
        self.multiplicity += edge.multiplicity;
        self.is_ref = self.is_ref || edge.is_ref;
    }

    /**
     * Create a new BaseEdge with the given multiplicity.
     * The resulting edge is a reference edge if any of the argument edges are reference.
     *
     * @param edges a collection of edges to or their isRef values
     * @param multiplicity our desired multiplicity
     * @return a newly allocated BaseEdge
     */
    fn make_o_r_edge(edges: Vec<Self>, multiplicity: usize, single_sample_capacity: usize) -> Self {
        assert!(!edges.is_empty(), "Edges cannot be empty");
        let is_ref = edges.par_iter().any(|e| e.is_ref());

        Self::new(is_ref, multiplicity)
    }

    fn to_string(&self) -> String {
        return format!(
            "BaseEdge{{multiplicity={}, isRef={}}}",
            self.multiplicity, self.is_ref
        );
    }
}
