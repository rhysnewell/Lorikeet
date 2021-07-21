use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use petgraph::prelude::NodeIndex;
use std::collections::HashSet;
use petgraph::algo::is_cyclic_directed;

/**
 * A common interface for the different KBestHaplotypeFinder implementations to conform to
 */
pub struct KBestHaplotypeFinder<'a, V: BaseVertex, E: BaseEdge> {
    pub(crate) graph: &'a BaseGraph<V, E>,
    pub(crate) sinks: HashSet<NodeIndex>,
    pub(crate) sources: HashSet<NodeIndex>,
}

impl<'a, V: BaseVertex, E: BaseEdge> KBestHaplotypeFinder<'a, V, E> {
    pub fn new(
        sinks: HashSet<NodeIndex>,
        sources: HashSet<NodeIndex>,
        graph: &'a BaseGraph<V, E>,
    ) -> KBestHaplotypeFinder<'a, V, E> {
        assert!(graph.contains_all_vertices(&sinks), "sink does not belong to the graph");
        assert!(graph.contains_all_vertices(&sources), "source does not belong to the graph");

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        // TODO: Lorikeet shouldn't have the above mentioned bugs, but will test and then implement
        //       above points if needed.
        if is_cyclic_directed(&graph.graph) {
            panic!("Input graph contains cycles")
        };

        KBestHaplotypeFinder {
            sinks,
            sources,
            graph
        }
    }
}