use petgraph::graph::{Graph, NodeIndex};
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use rayon::prelude::*;
use petgraph::Direction;
use std::collections::{BTreeSet, BinaryHeap};
use graphs::path::Path;
use std::cmp::Ordering;
use utils::base_utils::BaseUtils;
use linked_hash_set::LinkedHashSet;
use petgraph::csr::EdgeIndex;
use std::fs::File;
use std::io::Write;

/**
 * Common code for graphs used for local assembly.
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseGraph<V: BaseVertex, E: BaseEdge> {
    kmer_size: usize,
    pub graph: Graph<V, E>,
}

impl BaseGraph {

    pub fn new<V: BaseVertex, E: BaseEdge>(kmer_size: usize) -> BaseGraph<V, E> {
        BaseGraph {
            kmer_size,
            graph: Graph::<V, E>::new(),
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

    pub fn compare_paths(&self, first_path: &Path, second_path: &Path) -> Ordering {
        return BaseUtils::bases_comparator(
            first_path.get_bases(self.graph), second_path.get_bases(self.graph)
        )
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    pub fn incoming_vertices_of(&self, v: NodeIndex) -> LinkedHashSet<NodeIndex> {
        return self.graph.neighbors_directed(vertex_index, Direction::Incoming).par_iter().collect::<LinkedHashSet<NodeIndex>>()
    }

    /**
     * Get the set of vertices connected to v by outgoing edges
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    pub fn outgoing_vertices_of(&self, v: NodeIndex) -> LinkedHashSet<NodeIndex> {
        return self.graph.neighbors_directed(vertex_index, Direction::Outgoing).par_iter().collect::<LinkedHashSet<NodeIndex>>()
    }

    pub fn get_edge_target(&self, edge: EdgeIndex) -> NodeIndex {
        self.graph.edge_endpoints(edge).unwrap().1
    }

    pub fn get_edge_source(&self, edge: EdgeIndex) -> NodeIndex {
        self.graph.edge_endpoints(edge).unwrap().0
    }

    /**
    * Removes all provided vertices from the graph
    */
    pub fn remove_all_vertices(&mut self, vertices: &Vec<NodeIndex>) {
        self.graph.retain_nodes(|v| !vertices.contains(&v));
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    pub fn print_graph(&self, destination: &str, write_header: bool, prune_factor: usize) {
        let mut graph_writer = File::create(destination).unwrap();

        if write_header {
            graph_writer.write(b"digraph assemblyGraphs {")
        }

        for edge in self.graph.edge_indices() {

        }
    }
}