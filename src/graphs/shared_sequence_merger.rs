use graphs::base_edge::BaseEdge;
use graphs::seq_graph::SeqGraph;
use graphs::seq_vertex::SeqVertex;
use hashlink::LinkedHashSet;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::collections::HashSet;

/**
 * Merges the incoming vertices of a vertex V of a graph
 *
 * Looks at the vertices that are incoming to V (i.e., have an outgoing edge connecting to V).  If
 * they all have the same sequence, merges them into the sequence of V, and updates the graph
 * as appropriate
 */
pub struct SharedSequenceMerger {}

impl SharedSequenceMerger {
    /**
     * Attempt to merge the incoming vertices of v
     *
     * @param graph the graph containing the vertex v
     * @param v the vertex whose incoming vertices we want to merge
     * @return true if some useful merging was done, false otherwise
     */
    pub fn merge<E: BaseEdge>(graph: &mut SeqGraph<E>, v: NodeIndex) -> bool {
        assert!(
            graph.base_graph.graph.contains_node(v),
            "Graph doesn't contain vertex {:?}",
            v
        );

        let prevs = graph.base_graph.incoming_vertices_of(v);
        let can_merge = Self::can_merge(&graph, v, &prevs);
        debug!("Can Merge {}: {}", can_merge, prevs.is_empty());
        if !can_merge {
            return false;
        } else {
            let mut edges_to_remove = Vec::new();
            let mut prev_seq = graph
                .base_graph
                .graph
                .node_weight(*prevs.iter().next().unwrap())
                .unwrap()
                .sequence
                .clone();
            prev_seq.extend_from_slice(&graph.base_graph.graph.node_weight(v).unwrap().sequence);
            let new_v = SeqVertex::new(prev_seq);
            let new_v_index = graph.base_graph.add_node(&new_v);

            for prev in prevs.iter() {
                let incoming_edges = graph
                    .base_graph
                    .graph
                    .edges_directed(*prev, Direction::Incoming)
                    .map(|e| (e.id(), e.weight().clone()))
                    .collect::<Vec<(EdgeIndex, E)>>();
                for (prev_in_index, prev_in_weight) in incoming_edges {
                    graph.base_graph.graph.add_edge(
                        graph.base_graph.get_edge_source(prev_in_index),
                        new_v_index,
                        prev_in_weight,
                    );
                    edges_to_remove.push(prev_in_index);
                }
            }

            let outgoing_edges = graph
                .base_graph
                .graph
                .edges_directed(v, Direction::Outgoing)
                .map(|e| (e.id(), e.weight().clone()))
                .collect::<Vec<(EdgeIndex, E)>>();
            for (e_index, e_weight) in outgoing_edges {
                graph.base_graph.graph.add_edge(
                    new_v_index,
                    graph.base_graph.get_edge_target(e_index),
                    e_weight,
                );
            }

            debug!("Vertices to remove {} + {:?} edges to remove {}", prevs.len(), &v, edges_to_remove.len());
            graph.base_graph.remove_all_vertices(&prevs);
            graph.base_graph.graph.remove_node(v);
            graph.base_graph.remove_all_edges(&edges_to_remove);

            return true;
        }
    }

    /**
     * Can we safely merge the incoming vertices of v
     *
     * @param graph the graph containing v and incomingVertices
     * @param v the vertex we want to merge into
     * @param incomingVertices the incoming vertices of v
     * @return true if we can safely merge incomingVertices
     */
    fn can_merge<E: BaseEdge>(
        graph: &SeqGraph<E>,
        v: NodeIndex,
        incoming_vertices: &LinkedHashSet<NodeIndex>,
    ) -> bool {
        if incoming_vertices.is_empty() {
            return false;
        };

        let first = incoming_vertices.iter().next().unwrap();
        let mut count = 0;
        for prev in incoming_vertices {
            count += 1;
            debug!("{count} {} -> {}",
                   std::str::from_utf8(graph.base_graph.graph.node_weight(*first).unwrap().sequence.as_slice()).unwrap(),
                   std::str::from_utf8(graph.base_graph.graph.node_weight(*prev).unwrap().sequence.as_slice()).unwrap());
            if graph.base_graph.graph.node_weight(*prev).unwrap().sequence
                != graph.base_graph.graph.node_weight(*first).unwrap().sequence
            {
                debug!(
                    "1 {}: {} -> {}",
                    count,
                    std::str::from_utf8(graph.base_graph.graph.node_weight(*first).unwrap().sequence.as_slice()).unwrap(),
                    std::str::from_utf8(graph.base_graph.graph.node_weight(*prev).unwrap().sequence.as_slice()).unwrap()
                );
                return false;
            }

            let prev_outs = graph.base_graph.outgoing_vertices_of(*prev);
            if prev_outs.len() != 1 {
                // prev -> v must be the only edge from prev
                debug!("2 {} outgoing {}", count, prev_outs.len());
                return false;
            }

            if prev_outs.iter().next().unwrap() != &v {
                // don't allow cyles
                debug!("3 {}", count);
                return false;
            }

            if graph.base_graph.in_degree_of(*prev) == 0 {
                // cannot merge when any of the incoming nodes are sources
                debug!("4 {}", count);
                return false;
            }
        }

        return true;
    }
}
