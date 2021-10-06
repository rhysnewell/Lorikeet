use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use hashlink::LinkedHashSet;
use petgraph::algo::is_cyclic_directed;
use petgraph::prelude::EdgeIndex;
use petgraph::prelude::NodeIndex;
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::collections::{BTreeSet, HashSet};

/**
 * A common interface for the different KBestHaplotypeFinder implementations to conform to
 */
pub struct KBestHaplotypeFinder {
    pub sinks: HashSet<NodeIndex>,
    pub sources: HashSet<NodeIndex>,
}

impl KBestHaplotypeFinder {
    pub fn new<V: BaseVertex, E: BaseEdge>(
        sinks: HashSet<NodeIndex>,
        sources: HashSet<NodeIndex>,
        graph: &mut BaseGraph<V, E>,
    ) -> KBestHaplotypeFinder {
        assert!(
            graph.contains_all_vertices(&sinks),
            "sink does not belong to the graph"
        );
        assert!(
            graph.contains_all_vertices(&sources),
            "source does not belong to the graph"
        );

        let mut result = KBestHaplotypeFinder { sinks, sources };

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //     Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //     just the line commented out below if you want to trade early-bug-fail for speed.
        debug!("Graph before cycles removed {:?}", graph);
        result.remove_cycles_if_necessary(graph, false);
        debug!("Graph after cycles removed {:?}", graph);
        result
    }

    pub fn new_with_keep_cycles<V: BaseVertex, E: BaseEdge>(
        sinks: HashSet<NodeIndex>,
        sources: HashSet<NodeIndex>,
        graph: &mut BaseGraph<V, E>,
        keep_cycles: bool,
    ) -> Self {
        assert!(
            graph.contains_all_vertices(&sinks),
            "sink does not belong to the graph"
        );
        assert!(
            graph.contains_all_vertices(&sources),
            "source does not belong to the graph"
        );

        let mut result = KBestHaplotypeFinder { sinks, sources };

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        result.remove_cycles_if_necessary(graph, keep_cycles);
        result
    }

    fn remove_cycles_if_necessary<V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &mut BaseGraph<V, E>,
        keep_cycles: bool,
    ) {
        if keep_cycles {
            // pass
        } else {
            if is_cyclic_directed(&graph.graph) {
                self.remove_cycles_and_vertices_that_dont_lead_to_sinks(graph)
            }
        }
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    fn remove_cycles_and_vertices_that_dont_lead_to_sinks<V: BaseVertex, E: BaseEdge>(
        &self,
        original: &mut BaseGraph<V, E>,
    ) {
        let mut edges_to_remove = HashSet::with_capacity(original.graph.edge_count());
        let mut vertices_to_remove = HashSet::with_capacity(original.graph.node_count());

        let mut found_some_path = false;
        for source in self.sources.iter() {
            let mut parent_vertices = HashSet::with_capacity(original.graph.node_count());
            found_some_path = self.find_guilty_vertices_and_edges_to_remove_cycles(
                original,
                *source,
                &mut edges_to_remove,
                &mut vertices_to_remove,
                &mut parent_vertices,
            ) || found_some_path;
        }

        assert!(
            found_some_path,
            "Could not find any path from the source to the sink after removing cycles"
        );
        assert!(
            !(edges_to_remove.is_empty() && vertices_to_remove.is_empty()),
            "Cannot find a way to remove the cycles"
        );

        debug!("Nodes to remove {:?}", &vertices_to_remove);
        debug!("Edges to remove {:?}", &edges_to_remove);
        original.remove_all_edges(&edges_to_remove);
        original.remove_all_vertices(&vertices_to_remove);
    }

    /**
     * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
     *
     * @param graph the original graph.
     * @param currentVertex current search vertex.
     * @param sinks considered sink vertices.
     * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
     * @param verticesToRemove collection of vertices that can be removed.
     * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
     *                       reached from those vertices using edges existing in {@code graph}.
     *
     * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
     *  {@code false} otherwise.
     */
    fn find_guilty_vertices_and_edges_to_remove_cycles<V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &mut BaseGraph<V, E>,
        current_vertex: NodeIndex,
        edges_to_remove: &mut HashSet<EdgeIndex>,
        vertices_to_remove: &mut HashSet<NodeIndex>,
        parent_vertices: &mut HashSet<NodeIndex>,
    ) -> bool {
        if self.sinks.contains(&current_vertex) {
            return true;
        }

        let outgoing_edges = graph
            .graph
            .edges_directed(current_vertex, Direction::Outgoing)
            .map(|e| e.id())
            .collect::<LinkedHashSet<EdgeIndex>>();

        debug!(
            "Current vertex {:?} outgoing edges {:?}",
            &current_vertex, &outgoing_edges
        );
        parent_vertices.insert(current_vertex);

        let mut reaches_sink = false;
        for e in outgoing_edges {
            let child = graph.get_edge_target(e);
            debug!("Current edge {:?} child {:?}", &e, &child);
            if parent_vertices.contains(&child) {
                edges_to_remove.insert(e);
            } else {
                let child_reach_sink = self.find_guilty_vertices_and_edges_to_remove_cycles(
                    graph,
                    child,
                    edges_to_remove,
                    vertices_to_remove,
                    parent_vertices,
                );
                reaches_sink = reaches_sink || child_reach_sink;
            };
        }

        if !reaches_sink {
            vertices_to_remove.insert(current_vertex);
        };

        return reaches_sink;
    }
}
