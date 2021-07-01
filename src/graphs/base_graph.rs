use petgraph::graph::{Graph, NodeIndex};
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use rayon::prelude::*;
use petgraph::Direction;
use std::collections::{BTreeSet, BinaryHeap};

/**
 * Common code for graphs used for local assembly.
 */
pub struct BaseGraph {
    kmer_size: usize,
    pub graph: Graph<BaseVertex, BaseEdge>,
}

impl BaseGraph {

    pub fn new(kmer_size: usize) -> BaseGraph {
        BaseGraph {
            kmer_size,
            graph: Graph::<BaseVertex, BaseEdge>::new(),
        }
    }

    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    pub fn is_reference_node(&self, vertex_index: NodeIndex) -> bool {
        self.graph.edges(vertex_index).into_par_iter().any(|e| e.is_ref())
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    pub fn is_source(&self, vertex_index: NodeIndex) -> bool {
        return self.in_degree_of(vertex_index) == 0
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    pub fn is_sink(&self, veretex_index: NodeIndex) -> bool {
        return self.out_degree_of(veretex_index) == 0
    }

    /**
     * @param v the vertex to test
     * @return number of incoming edges
     */
    pub fn in_degree_of(&self, vertex_index: NodeIndex) -> usize {
        self.graph.edges_directed(vertex_index, Direction::Incoming).count()
    }

    /**
     * @param v the vertex to test
     * @return number of outgoing edges
     */
    pub fn out_degree_of(&self, vertex_index: NodeIndex) -> usize {
        self.graph.edges_directed(vertex_index, Direction::Outgoing).count()
    }

    /**
     * Get the set of source vertices of this graph
     * NOTE: We return a BTreeSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being BTreeSet
     * @return a non-null set
     */
    pub fn get_sources(&self) -> BinaryHeap<NodeIndex> {
        return self.graph.node_indices().into_par_iter().filter(|v_index| {
            self.is_source(v_index)
        }).collect::<BinaryHeap<NodeIndex>>()
    }

    /**
     * Get the set of sink vertices of this graph
     * NOTE: We return a BTreeSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being BTreeSet
     * @return a non-null set
     */
    pub fn get_sinks(&self) -> BinaryHeap<NodeIndex> {
        return self.graph.node_indices().into_par_iter().filter(|v_index| {
            self.is_sink(v_index)
        }).collect::<BinaryHeap<NodeIndex>>()
    }
}