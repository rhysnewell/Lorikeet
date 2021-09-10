use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use ordered_float::OrderedFloat;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use rayon::prelude::*;
use reads::cigar_utils::CigarUtils;
use rust_htslib::bam::record::CigarString;
use smith_waterman::bindings::SWOverhangStrategy;
use smith_waterman::smith_waterman_aligner::NEW_SW_PARAMETERS;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::hash::{Hash, Hasher};
use utils::base_utils::BaseUtils;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 */
#[derive(Debug, Clone)]
pub struct Path {
    pub(crate) last_vertex: NodeIndex,
    pub(crate) edges_in_order: Vec<EdgeIndex>,
}

impl Path {
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    pub fn new(initial_vertex: NodeIndex, edges_in_order: Vec<EdgeIndex>) -> Self {
        Self {
            last_vertex: initial_vertex,
            edges_in_order,
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
    pub fn new_add_edge<V: BaseVertex, E: BaseEdge>(
        &self,
        edge: EdgeIndex,
        graph: &BaseGraph<V, E>,
    ) -> Path {
        assert!(
            graph.contains_edge(edge),
            "Graph must contain edge {:?} but it doesn't",
            edge
        );
        assert_eq!(
            graph.get_edge_source(edge),
            self.last_vertex,
            "Edges added to path must be contiguous."
        );
        let mut edges_in_order = self.edges_in_order.clone();
        edges_in_order.push(edge);
        Path {
            last_vertex: graph.get_edge_target(edge),
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
    pub fn new_add_edges<V: BaseVertex, E: BaseEdge>(
        &self,
        edges: Vec<EdgeIndex>,
        graph: &BaseGraph<V, E>,
    ) -> Path {
        assert!(
            edges.par_iter().all(|edge| graph.contains_edge(*edge)),
            "Graph does not contain an edge that attempted to be added to the path"
        );

        let mut tmp_vertex = self.last_vertex;
        for edge in edges.iter() {
            if graph.get_edge_source(*edge) != tmp_vertex {
                panic!("Edges added to the path must be contiguous")
            };
            tmp_vertex = graph.get_edge_target(*edge)
        }
        let mut edges_in_order = self.edges_in_order.clone();
        edges_in_order.par_extend(edges);
        Path {
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
    pub fn contains_vertex<V: BaseVertex, E: BaseEdge>(
        &self,
        v: NodeIndex,
        graph: &BaseGraph<V, E>,
    ) -> bool {
        v == self.get_first_vertex(graph)
            || self
                .edges_in_order
                .par_iter()
                .map(|e| graph.graph.edge_endpoints(*e).unwrap().1)
                .any(|edge_target| edge_target == v)
    }

    pub fn to_string<V: BaseVertex, E: BaseEdge>(&self, graph: &BaseGraph<V, E>) -> String {
        let mut joined_path = self
            .get_vertices(graph)
            .par_iter()
            .map(|v| {
                graph
                    .graph
                    .node_weight(*v)
                    .unwrap()
                    .get_sequence_string()
                    .to_string()
            })
            .collect::<Vec<String>>()
            .join("->");

        return format!("Path{{path={}}}", joined_path);
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
    pub fn get_vertices<V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &BaseGraph<V, E>,
    ) -> Vec<NodeIndex> {
        let mut result = Vec::with_capacity(self.edges_in_order.len() + 1);
        result.push(self.get_first_vertex(graph));
        result.par_extend(
            self.edges_in_order
                .par_iter()
                .map(|e| graph.graph.edge_endpoints(*e).unwrap().1)
                .collect::<Vec<NodeIndex>>(),
        );
        return result;
    }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    pub fn get_first_vertex<V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &BaseGraph<V, E>,
    ) -> NodeIndex {
        if self.edges_in_order.is_empty() {
            return self.last_vertex;
        } else {
            return graph
                .graph
                .edge_endpoints(self.edges_in_order[0])
                .unwrap()
                .0;
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
    pub fn get_bases<V: BaseVertex, E: BaseEdge>(&self, graph: &BaseGraph<V, E>) -> Vec<u8> {
        if self.edges_in_order.is_empty() {
            return graph
                .graph
                .node_weight(self.last_vertex)
                .unwrap()
                .get_additional_sequence(true)
                .to_vec();
        }

        let mut bases = graph
            .graph
            .node_weight(
                graph
                    .graph
                    .edge_endpoints(self.edges_in_order[0])
                    .unwrap()
                    .0,
            )
            .unwrap()
            .get_additional_sequence(true)
            .to_vec();
        for e in self.edges_in_order.iter() {
            bases.par_extend(
                graph
                    .graph
                    .node_weight(graph.graph.edge_endpoints(*e).unwrap().1)
                    .unwrap()
                    .get_additional_sequence(true),
            );
        }

        return bases;
    }

    pub fn get_max_multiplicity<V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &BaseGraph<V, E>,
    ) -> usize {
        self.get_edges()
            .iter()
            .map(|e| graph.graph.edge_weight(*e).unwrap().get_multiplicity())
            .max()
            .unwrap_or(0)
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @param aligner
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    pub fn calculate_cigar<V: BaseVertex, E: BaseEdge>(
        &self,
        ref_seq: &[u8],
        graph: &BaseGraph<V, E>,
    ) -> CigarString {
        return CigarUtils::calculate_cigar(
            ref_seq,
            &self.get_bases(graph),
            SWOverhangStrategy::SoftClip,
            &*NEW_SW_PARAMETERS,
        )
        .unwrap();
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

impl PartialEq for Path {
    fn eq(&self, other: &Self) -> bool {
        self.edges_in_order == other.edges_in_order
    }
}

// `Eq` needs to be implemented as well.
impl Eq for Path {}

impl Hash for Path {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.edges_in_order.hash(state)
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct Chain<'a, V: BaseVertex + Sync + std::fmt::Debug, E: BaseEdge + Sync + std::fmt::Debug> {
    pub log_odds: OrderedFloat<f64>,
    pub path: &'a Path,
    pub graph: &'a BaseGraph<V, E>,
}

impl<'a, V: BaseVertex + Sync + std::fmt::Debug, E: BaseEdge + Sync + std::fmt::Debug>
    Chain<'a, V, E>
{
    pub fn new(
        log_odds: OrderedFloat<f64>,
        path: &'a Path,
        graph: &'a BaseGraph<V, E>,
    ) -> Chain<'a, V, E> {
        Chain {
            log_odds,
            path,
            graph,
        }
    }
}

impl<'a, V: BaseVertex + Sync + std::fmt::Debug, E: BaseEdge + Sync + std::fmt::Debug> Ord
    for Chain<'a, V, E>
{
    fn cmp(&self, other: &Self) -> Ordering {
        (-self.log_odds).cmp(&-other.log_odds).then_with(|| {
            BaseUtils::bases_comparator(
                self.graph
                    .get_sequence_from_index(self.path.get_first_vertex(self.graph)),
                other
                    .graph
                    .get_sequence_from_index(other.path.get_first_vertex(other.graph)),
            )
        })
    }
}

impl<'a, V: BaseVertex, E: BaseEdge> PartialOrd for Chain<'a, V, E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
