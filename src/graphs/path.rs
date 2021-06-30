use rayon::prelude::*;
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 */
pub struct Path {
    last_vertex: BaseVertex,
    edges_in_order: Vec<BaseEdge>,
    graph: BaseGraph
}

impl Path {
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    pub fn new(initial_vertex: BaseVertex, graph: BaseGraph) -> Self {
        Self {
            last_vertex: initial_vertex,
            graph,
            edges_in_order: Vec::new()
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
    pub fn extend_path(p: &mut Self, edge: BaseEdge) {
        assert!(p.graph.graph.)
    }

}