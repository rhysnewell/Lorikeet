use std::collections::BinaryHeap;
use graphs::base_edge::BaseEdge;
use rayon::prelude::*;

/**
 * Edge class for connecting nodes in the graph that tracks some per-sample information.
 * <p>
 * This class extends BaseEdge with the additional functionality of tracking the maximum
 * multiplicity seen within any single sample.  The workflow for using this class is:
 * </p>
 * <pre>
 * {@code
 *      MultiSampleEdge e = new MultiSampleEdge(ref, 1)
 *      e.incMultiplicity(1)              // total is 2, per sample is 2, max per sample is 1
 *      e.getPruningMultiplicity()        // = 1
 *      e.flushSingleSampleMultiplicity() // total is 2, per sample is 0, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.incMultiplicity(3)              // total is 5, per sample is 3, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.flushSingleSampleMultiplicity() // total is 5, per sample is 0, max per sample is 3
 *      e.getPruningMultiplicity()        // = 3
 * }
 * </pre>
 */
#[derive(Debug, Clone)]
pub struct MultiSampleEdge {
    current_single_sample_multiplicity: usize,
    single_sample_capacity: usize,
    single_sample_multiplicities: BinaryHeap<usize>,
    reference_path_indexes: Vec<usize>,
    pub(crate) multiplicity: usize,
    pub(crate) is_ref: bool,
}

impl MultiSampleEdge {

    /**
     * Create a new MultiSampleEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge in this sample
     * @param singleSampleCapacity the max number of samples to track edge multiplicities
     */
    pub fn new(is_ref: bool, multiplicity: usize, single_sample_capacity: usize) -> MultiSampleEdge {
        let mut single_sample_multiplicities = BinaryHeap::with_capacity(single_sample_capacity);
        single_sample_multiplicities.push(multiplicity);

        MultiSampleEdge {
            multiplicity,
            is_ref,
            single_sample_multiplicities,
            single_sample_capacity,
            current_single_sample_multiplicity: multiplicity,
            reference_path_indexes: Vec::with_capacity(2),
        }
    }

    pub fn set(&mut self, is_ref: bool, multiplicity: usize, single_sample_capacity: usize) {
        let mut single_sample_multiplicities = BinaryHeap::with_capacity(single_sample_capacity);
        single_sample_multiplicities.push(multiplicity);
        self.multiplicity = multiplicity;
        self.is_ref = is_ref;
        self.single_sample_capacity = single_sample_capacity;
        self.single_sample_multiplicities = single_sample_multiplicities;
        self.current_single_sample_multiplicity = multiplicity;
        self.reference_path_indexes = Vec::with_capacity(2);
    }

    /**
     * update the single sample multiplicities by adding the current single sample multiplicity to the priority queue, and
     * reset the current single sample multiplicity to 0.
     */
    pub fn flush_single_sample_multiplicity(&mut self) {
        self.single_sample_multiplicities.push(self.current_single_sample_multiplicity);
        if self.single_sample_multiplicities.len() == self.single_sample_capacity + 1 {
            self.single_sample_multiplicities.pop();
        } else if self.single_sample_multiplicities.len() > self.single_sample_capacity + 1 {
            panic!("Somehow the per sample multiplicity list has grown too big: {:?}", self.single_sample_multiplicities);
        }

        self.current_single_sample_multiplicity = 0;
    }

    pub fn inc_multiplicity(&mut self, incr: usize) {
        self.multiplicity += incr;
        self.current_single_sample_multiplicity += incr;
    }

    pub fn get_pruning_multiplicity(&self) -> usize {
        *self.single_sample_multiplicities.peek().unwrap()
    }

    pub fn add_reference_index(&mut self, i: usize) {
        self.reference_path_indexes.push(i)
    }

    pub fn get_reference_path_indexes(&self) -> &Vec<usize> {
        &self.reference_path_indexes
    }
}

impl BaseEdge for MultiSampleEdge {

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
        return format!("{}/{}", self.multiplicity.to_string(), self.get_pruning_multiplicity())
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
        return self.is_ref
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

        Self::new(is_ref, multiplicity, single_sample_capacity)
    }

    fn to_string(&self) -> String {
        return format!("BaseEdge{{multiplicity={}, isRef={}}}", self.multiplicity, self.is_ref)
    }
}