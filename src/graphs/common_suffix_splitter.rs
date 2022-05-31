use graphs::base_edge::BaseEdge;
use graphs::base_vertex::BaseVertex;
use graphs::graph_utils::GraphUtils;
use graphs::seq_graph::SeqGraph;
use graphs::seq_vertex::SeqVertex;
use hashlink::LinkedHashSet;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use std::collections::HashSet;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V.  This replaces each
 * Vi with three nodes prefix, seq_i, and suffix connected in a simple chain.
 *
 * This operation can be performed in a very general case, without too much worry about the incoming
 * and outgoing edge structure of each Vi.  The partner algorithm SharedSequenceMerger can
 * put these pieces back together in a smart way that maximizes the sharing of nodes
 * while respecting complex connectivity.
 *
 */
pub struct CommonSuffixSplitter {}

impl CommonSuffixSplitter {
    /**
     * Simple single-function interface to split and then update a graph
     *
     * @param graph the graph containing the vertices in toMerge
     * @param v The bottom node whose incoming vertices we'd like to split
     * @return true if some useful splitting was done, false otherwise
     */
    pub fn split<E: BaseEdge>(graph: &mut SeqGraph<E>, v: NodeIndex) -> bool {
        assert!(
            graph.base_graph.graph.contains_node(v),
            "Graph doesn't contain the vertex {:?}",
            v
        );

        let to_split = graph.base_graph.incoming_vertices_of(v);

        let suffix_v_template = Self::common_suffix(&graph, v, &to_split);
        match suffix_v_template {
            None => return false,
            Some(suffix_v_template) => {
                let mut edges_to_remove = Vec::new();

                for mid in to_split.iter() {
                    let suffix_v = graph.base_graph.add_node(suffix_v_template.clone());
                    let prefix_v = graph
                        .base_graph
                        .graph
                        .node_weight(*mid)
                        .unwrap()
                        .without_suffix(suffix_v_template.get_sequence());
                    let out = graph.base_graph.outgoing_edge_of(*mid).unwrap();
                    // debug!("Mid {:?} out {:?}", &mid, &out);
                    let mut incoming_target;
                    match prefix_v {
                        None => {
                            // this node is entirely explained by suffix
                            incoming_target = suffix_v;
                        }
                        Some(prefix_v) => {
                            let prefix_v_index = graph.base_graph.add_node(prefix_v);
                            incoming_target = prefix_v_index;
                            let mut out_weight =
                                graph.base_graph.graph.edge_weight(out).unwrap().clone();
                            out_weight.set_multiplicity(1);
                            graph
                                .base_graph
                                .graph
                                .add_edge(prefix_v_index, suffix_v, out_weight);
                            edges_to_remove.push(out);
                        }
                    }
                    // debug!("Incoming target {:?}", &incoming_target);
                    graph.base_graph.graph.add_edge(
                        suffix_v,
                        graph.base_graph.get_edge_target(out),
                        graph.base_graph.graph.edge_weight(out).unwrap().clone(),
                    );
                    let incoming_edges = graph
                        .base_graph
                        .graph
                        .edges_directed(*mid, Direction::Incoming)
                        .map(|e| (e.id(), e.weight().clone()))
                        .collect::<Vec<(EdgeIndex, E)>>();
                    for (incoming_index, incoming_weight) in incoming_edges {
                        graph.base_graph.graph.add_edge(
                            graph.base_graph.get_edge_source(incoming_index),
                            incoming_target,
                            incoming_weight,
                        );
                        edges_to_remove.push(incoming_index);
                    }
                }
                // debug!("Node to remove {:?}", &to_split);
                // debug!("Edges to remove {:?}", &edges_to_remove);
                graph.base_graph.remove_all_vertices(&to_split);
                graph.base_graph.remove_all_edges(&edges_to_remove);
                return true;
            }
        }
    }

    fn common_suffix<E: BaseEdge>(
        graph: &SeqGraph<E>,
        v: NodeIndex,
        to_split: &LinkedHashSet<NodeIndex>,
    ) -> Option<SeqVertex> {
        if to_split.len() < 2 {
            // Can only split at least 2 vertices
            return None;
        } else if !Self::safe_to_split(graph, v, to_split) {
            return None;
        }

        let suffix_v_template = Self::common_suffix_of_set(to_split, graph);

        if suffix_v_template.is_empty() {
            return None;
        } else if Self::would_eliminate_ref_source(graph, &suffix_v_template, to_split) {
            return None;
        } else if Self::all_vertices_are_the_common_suffix(&suffix_v_template, to_split, graph) {
            return None;
        } else {
            return Some(suffix_v_template);
        }
    }

    /**
     * Would all vertices that we'd split just result in the common suffix?
     *
     * That is, suppose we have prefix nodes ABC and ABC.  After splitting all of the vertices would
     * just be ABC again, and we'd enter into an infinite loop.
     *
     * @param commonSuffix the common suffix of all vertices in toSplits
     * @param toSplits the collection of vertices we want to split
     * @return true if all of the vertices are equal to the common suffix
     */
    fn all_vertices_are_the_common_suffix<E: BaseEdge>(
        common_suffix: &SeqVertex,
        to_splits: &LinkedHashSet<NodeIndex>,
        graph: &SeqGraph<E>,
    ) -> bool {
        to_splits
            .iter()
            .all(|v| graph.base_graph.graph.node_weight(*v).unwrap().len() == common_suffix.len())
    }

    /**
     * Would factoring out this suffix result in elimating the reference source vertex?
     * @param graph the graph
     * @param commonSuffix the common suffix of all toSplits
     * @param toSplits the list of vertices we're are trying to split
     * @return true if toSplit contains the reference source and this ref source has all and only the bases of commonSuffix
     */
    fn would_eliminate_ref_source<E: BaseEdge>(
        graph: &SeqGraph<E>,
        common_suffix: &SeqVertex,
        to_splits: &LinkedHashSet<NodeIndex>,
    ) -> bool {
        for to_split in to_splits {
            if graph.base_graph.is_ref_source(*to_split) {
                return graph.base_graph.graph.node_weight(*to_split).unwrap().len()
                    == common_suffix.len();
            }
        }
        return false;
    }

    /**
     * Can we safely split up the vertices in toMerge?
     *
     * @param graph a graph
     * @param bot a vertex whose incoming vertices we want to split
     * @param toMerge the set of vertices we'd be splitting up
     * @return true if we can safely split up toMerge
     */
    fn safe_to_split<E: BaseEdge>(
        graph: &SeqGraph<E>,
        bot: NodeIndex,
        to_merge: &LinkedHashSet<NodeIndex>,
    ) -> bool {
        let outgoing_of_bot = graph.base_graph.outgoing_vertices_of(bot);
        for m in to_merge {
            let outs = graph
                .base_graph
                .graph
                .edges_directed(*m, Direction::Outgoing)
                .map(|e| e.id())
                .collect::<HashSet<EdgeIndex>>();
            if m == &bot
                || outs.len() != 1
                || !graph.base_graph.outgoing_vertices_of(*m).contains(&bot)
            {
                // m == bot => don't allow self cycles in the graph
                return false;
            }

            if outgoing_of_bot.contains(m) {
                // forbid cycles from bottom -> mid
                return false;
            }
        }

        return true;
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return a single vertex that contains the common suffix of all middle vertices
     */
    fn common_suffix_of_set<E: BaseEdge>(
        middle_vertices: &LinkedHashSet<NodeIndex>,
        graph: &SeqGraph<E>,
    ) -> SeqVertex {
        let kmers = GraphUtils::get_kmers(middle_vertices, &graph.base_graph);
        let min = GraphUtils::min_kmer_length(&kmers);
        let suffix_len = GraphUtils::common_maximum_suffix_length(&kmers, min);
        let kmer = &kmers[0];
        let suffix = kmer[kmer.len() - suffix_len..kmer.len()].to_vec();
        // debug!(
        //     "Choosing suffix {:?}",
        //     std::str::from_utf8(&suffix).unwrap()
        // );
        return SeqVertex::new(suffix);
    }
}
