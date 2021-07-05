use petgraph::graph::{Graph, NodeIndex};
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use rayon::prelude::*;
use petgraph::Direction;
use std::collections::{BTreeSet, BinaryHeap, HashSet};
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

impl<V: BaseVertex, E: BaseEdge> BaseGraph<V, E> {

    pub fn new(kmer_size: usize) -> BaseGraph<V, E> {
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
            graph_writer.write(b"digraph assemblyGraphs {\n")
        }

        for edge in self.graph.edge_indices() {
            let edge_string = format!("\t{} -> {} ", self.get_edge_source(edge).into(), self.get_edge_target(edge).into());
            let mut edge_label_string;
            let edge_weight = self.graph.edge_weight(edge).unwrap();
            if edge_weight.get_multiplicity() > 0 && edge.get_multiplicity() < prune_factor {
                edge_label_string = format!("[style=dotted,color=grey,label=\"{}\"];", edge.get_dot_label());
            } else {
                edge_label_string = format!("[label=\"{}\"];", edge.get_dot_label());
            }
            graph_writer.write(edge_string.as_bytes());
            graph_writer.write(edge_label_string.as_bytes());
            if edge_weight.is_ref() {
                graph_writer.write(format!("{} [color=red];\n", edge_string).as_bytes());
            }
        }

        for v in self.graph.node_indices() {
            let node_weight = self.graph.node_weight(v).unwrap();
            graph_writer.write(format!("\t{} [label=\"{}\",shape=box]\n", node_weight.to_string(),
                                       format!(
                                           "{}{}",
                                           std::str::from_utf8(self.get_additional_sequence(v, node_weight)).unwrap(),
                                           node_weight.get_additional_info()
                                       )));
        }

        if write_header {
            graph_writer.print("}\n")
        }
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    pub fn get_additional_sequence(&self, index: NodeIndex, v: &V) -> &[u8] {
        return v.get_additional_sequence(self.is_source(index))
    }

    /**
     * Get the set of vertices within distance edges of source, regardless of edge direction
     *
     * @param source the source vertex to consider
     * @param distance the distance
     * @return a set of vertices within distance of source
     */
    fn vertices_within_distance(&self, source: NodeIndex, distance: usize) -> HashSet<NodeIndex> {
        if distance == 0 {
            return vec![source]
        }

        let mut found = HashSet::new();
        found.insert(source);
        for v in self.graph.neighbors(source) {
            found.par_extend(self.vertices_within_distance(v, distance - 1))
        }

        return found
    }

    /**
     * Get a graph containing only the vertices within distance edges of target
     * @param target a vertex in graph
     * @param distance the max distance
     * @return a non-null graph
     */
    pub fn subset_to_neighbors(&self, target: NodeIndex, distance: usize) -> BaseGraph<V, E> {
        let to_keep = self.vertices_within_distance(target, distance);
        let mut to_remove = self.graph.node_indices().collect_vec();
        to_remove.re
    }
}