use rayon::prelude::*;
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use petgraph::csr::NodeIndex;
use petgraph::graph::{EdgeIndex, Edge, EdgeReference};
use utils::smith_waterman_aligner::SmithWatermanAligner;
use rust_htslib::bam::record::CigarString;
use reads::cigar_utils::CigarUtils;
use std::cmp::Ordering;
use ordered_float::OrderedFloat;
use utils::base_utils::BaseUtils;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 */
#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq, Copy)]
pub struct Path {
    last_vertex: NodeIndex,
    edges_in_order: Vec<EdgeReference<BaseEdge, u32>>,
    // graph: BaseGraph
}

impl Path {
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    pub fn new(initial_vertex: NodeIndex, edges_in_order: Vec<EdgeReference<BaseEdge, u32>>) -> Self {
        Self {
            last_vertex: initial_vertex,
            edges_in_order: edges_in_order,
        }
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
    // pub fn extend_path(p: &mut Self, edge: BaseEdge) {
    //     assert!(p.graph.graph.)
    // }

    fn paths_are_the_same(&self, path: &Self) -> bool {
        self.edges_in_order == path.edges_in_order
    }

    /**
     * Does this path contain the given vertex?
     *
     * @param v a non-null vertex
     * @return true if v occurs within this path, false otherwise
     */
    pub fn contains_vertex(&self, v: NodeIndex) -> bool {
        v == self.get_first_vertex() || self.edges_in_order.par_iter().map(|e| {
            e.target()
        }).any(|edge_target| edge_target == v)
    }

    pub fn to_string(&self) -> String {

    }

    /**
     * Get the edges of this path in order.
     * Returns an unmodifiable view of the underlying list
     * @return a non-null list of edges
     */
    pub fn get_edges(&self) -> &Vec<EdgeIndex> {
        &self.edges_in_order
    }

    pub fn get_last_edge(&self) -> EdgeIndex {
        self.edges_in_order.last().unwrap();
    }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    pub fn get_vertices(&self) -> Vec<NodeIndex> {
        let mut result = Vec::with_capacity(self.edges_in_order.len() + 1);
        result.add(self.get_first_vertex());
        result.par_extend(self.edges_in_order.par_iter().map(|e| e.target()).collect::<Vec<NodeIndex>>());
        return result
    }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    pub fn get_first_vertex(&self) -> NodeIndex {
        if self.edges_in_order.is_empty() {
            return self.last_vertex
        } else {
            return self.edges_in_order[0].source()
        }
    }

    /**
     * Get the final vertex of the path
     * @return a non-null vertex
     */
    pub fn get_last_vertex(&self) -> NodeIndex {
        self.last_vertex
    }

    /**
     * The base sequence for this path. Pull the full sequence for source nodes and then the suffix for all subsequent nodes
     * @return  non-null sequence of bases corresponding to this path
     */
    pub fn get_bases(&self, graph: &BaseGraph) -> &[u8] {
        if self.edges_in_order.is_empty() {
            return graph.graph[self.last_vertex].unwrap().get_additional_sequence(true)
        }

        let mut bases = graph.graph[self.edges_in_order[0].source()].unwrap().get_additional_sequence(true);
        for e in self.edges_in_order {
            bases.par_extend(graph.graph[e].unwrap().get_additional_sequence(true));
        }

        return bases
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @param aligner
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    pub fn calculate_cigar(&self, ref_seq: &[u8], grasph: &BaseGraph, aligner: SmithWatermanAligner) -> CigarString {
        return CigarUtils::calculate_cigar(ref_seq, self.get_bases(graph), aligner)
    }

    /**
     * Length of the path in edges.
     *
     * @return {@code 0} or greater.
     */
    pub fn len(&self) -> usize {
        self.edges_in_order.len()
    }
}

pub struct PathWithGraph

#[derive(Debug, Eq, PartialEq)]
pub struct Chain {
    pub log_odds: OrderedFloat<f64>,
    pub path: &Path,
    graph: &BaseGraph,
}

impl Chain {
    pub fn new(log_odds: OrderedFloat<f64>, path: &Path, graph: &BaseGraph) -> Chain {
        Chain {
            log_odds,
            path,
            graph
        }
    }
}

// The priority queue depends on `Ord`.
// Explicitly implement the trait so the queue becomes a min-heap
// instead of a max-heap.
impl Ord for Chain {
    fn cmp(&self, other: &Self) -> Ordering {
        // In case of a tie we compare positions - this step is necessary
        // to make implementations of `PartialEq` and `Ord` consistent.
        (-other.log_odds).cmp(&(-self.log_odds))
            .then_with(|| BaseUtils::bases_comparator(self.path.get_bases(self.graph), other.path.get_bases(other.graph)))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for Chain {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}