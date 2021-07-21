use rayon::prelude::*;
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use rust_htslib::bam::record::CigarString;
use reads::cigar_utils::CigarUtils;
use ordered_float::OrderedFloat;
use std::hash::{Hash, Hasher};
use smith_waterman::bindings::SWOverhangStrategy;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 */
#[derive(Debug, Clone)]
pub struct Path<'a, V: BaseVertex, E: BaseEdge> {
    pub(crate) last_vertex: NodeIndex,
    pub(crate) edges_in_order: Vec<EdgeIndex>,
    pub(crate) graph: &'a BaseGraph<V, E>
}

impl<'a, V: BaseVertex + std::marker::Sync, E: BaseEdge + std::marker::Sync> Path<'a, V, E> {
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    pub fn new(
        initial_vertex: NodeIndex, edges_in_order: Vec<EdgeIndex>, graph: &'a BaseGraph<V, E>
    ) -> Self {
        Self {
            last_vertex: initial_vertex,
            edges_in_order,
            graph,
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
    pub fn new_add_edge(&self, edge: EdgeIndex) -> Path<'a, V, E> {
        assert!(self.graph.contains_edge(edge), "Graph must contain edge {:?} but it doesn't", edge);
        assert_eq!(self.graph.get_edge_source(edge), self.last_vertex, "Edges added to path must be contiguous.");
        let mut edges_in_order = self.edges_in_order.clone();
        edges_in_order.push(edge);
        Path {
            graph: self.graph,
            last_vertex: self.graph.get_edge_target(edge),
            edges_in_order,
        }
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edges list of edges to extend. Does not check arguments' validity i.e. doesn't check that edges are in order
     *
     * @throws IllegalArgumentException if {@code p} or {@code edges} are {@code null} or empty, or {@code edges} is
     * not part of {@code p}'s graph, or {@code edges} does not have as a source the last vertex in {@code p}.
     */
    pub fn new_add_edges(&self, edges: Vec<EdgeIndex>) -> Path<'a, V, E> {
        assert!(edges.par_iter().all(|edge| self.graph.contains_edge(*edge)),
                "Graph does not contain an edge that attempted to be added to the path");

        let tmp_vertex = self.last_vertex;
        for edge in edges.iter() {
            if self.graph.get_edge_source(*edge) != tmp_vertex {
                panic!("Edges added to the path must be contiguous")
            };
            tmp_vertex = self.graph.get_edge_target(*edge)
        }
        let mut edges_in_order = self.edges_in_order.clone();
        edges_in_order.par_extend(edges);
        Path {
            graph: self.graph,
            last_vertex: tmp_vertex,
            edges_in_order,
        }
    }

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
            self.graph.graph.edge_endpoints(*e).unwrap().1
        }).any(|edge_target| edge_target == v)
    }

    pub fn to_string(&self) -> String {
        let mut joined_path = self.get_vertices().par_iter().map(|v| {
            self.graph.graph.node_weight(*v).unwrap().get_sequence_string()
        }).collect::<Vec<String>>().join("->");

        return format!("Path{{path={}}}", joined_path)
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
        *self.edges_in_order.last().unwrap()
    }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    pub fn get_vertices(&self) -> Vec<NodeIndex> {
        let mut result = Vec::with_capacity(self.edges_in_order.len() + 1);
        result.push(self.get_first_vertex());
        result.par_extend(self.edges_in_order.par_iter().map(|e| {
            self.graph.graph.edge_endpoints(*e).unwrap().1
        }).collect::<Vec<NodeIndex>>());
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
            return self.graph.graph.edge_endpoints(self.edges_in_order[0]).unwrap().0
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
    pub fn get_bases(&self) -> &[u8] {
        if self.edges_in_order.is_empty() {
            return self.graph.graph.node_weight(self.last_vertex).unwrap().get_additional_sequence(true)
        }

        let mut bases =
            self.graph.graph.node_weight(self.graph.graph.edge_endpoints(self.edges_in_order[0]).unwrap().0).unwrap().get_additional_sequence(true).to_vec();
        for e in self.edges_in_order {
            bases.par_extend(self.graph.graph.node_weight(self.graph.graph.edge_endpoints(e).unwrap().1).unwrap().get_additional_sequence(true));
        }

        return &bases[..]
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @param aligner
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    pub fn calculate_cigar(&self, ref_seq: &[u8],) -> CigarString {
        return CigarUtils::calculate_cigar(ref_seq, self.get_bases(), SWOverhangStrategy::SoftClip).unwrap()
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

impl<'a, V: BaseVertex, E: BaseEdge> PartialEq for Path<'a, V, E> {
    fn eq(&self, other: &Self) -> bool {
        self.edges_in_order == other.edges_in_order
    }
}

// `Eq` needs to be implemented as well.
impl<'a, V: BaseVertex, E: BaseEdge> Eq for Path<'a, V, E> {}

impl<V: BaseVertex, E: BaseEdge> Hash for Path<'_, V, E> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.edges_in_order.hash(state)
    }
}

#[derive(Debug)]
pub struct Chain<'a, V: BaseVertex, E: BaseEdge> {
    pub log_odds: OrderedFloat<f64>,
    pub path: &'a Path<'a, V, E>,
}

impl<'a, V: BaseVertex, E: BaseEdge> Chain<'a, V, E> {
    pub fn new(log_odds: OrderedFloat<f64>, path: &'a Path<'a, V, E>) -> Chain<'a, V, E> {
        Chain {
            log_odds,
            path,
        }
    }
}