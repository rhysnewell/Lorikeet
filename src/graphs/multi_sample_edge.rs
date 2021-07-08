use std::collections::BinaryHeap;
use graphs::base_edge::BaseEdge;

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