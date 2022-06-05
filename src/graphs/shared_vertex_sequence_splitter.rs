use graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use graphs::base_vertex::BaseVertex;
use graphs::graph_utils::GraphUtils;
use graphs::seq_graph::SeqGraph;
use graphs::seq_vertex::SeqVertex;
use hashlink::LinkedHashSet;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
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
 * where prefix and suffix are shared sequences across all vertices V
 *
 * This algorithm creates a new SeqGraph with the following configuration
 *
 * prefix -> has outgoing edges to all seq_i
 * suffix -> has incoming edges for all seq_i
 *
 * There are a few special cases that must be handled.  First, Vi could be simply
 * == to the prefix or the suffix.  These generate direct connections between
 * the prefix and suffix nodes, and they are handled internally by the algorithm.
 *
 * Note that for convenience, we will always create newTop and newBottom nodes, but
 * these may be empty node (i.e., they contain no sequence).  That allows them to be
 * trivially merged, if desired, when the graph is incorporated into an overall
 * graph.
 *
 * The product of this operation is a SeqGraph that contains the split.  There's a
 * function to merge reconnect this graph into the graph that contains the middle nodes
 *
 * The process guarentees a few things about the output:
 *
 * -- Preserves the paths and weights among all vertices
 *
 * It produces a graph that has some unusual properties
 *
 * -- May add nodes with no sequence (isEmpty() == true) to preserve connectivity among the graph
 * -- May introduce edges with no multiplicity to preserve paths through the graph
 *
 * The overall workflow of using this class is simple:
 *
 * find vertices V in graph that you want to split out
 * s = new SharedVertexSequenceSplitter(graph, V)
 * s.updateGraph(graph)
 *
 * to update the graph with the modifications created by this splitter
 */
pub struct SharedVertexSequenceSplitter<'a, E: BaseEdge> {
    outer: &'a mut SeqGraph<E>,
    prefix_v: SeqVertex,
    suffix_v: SeqVertex,
    prefix_v_index: NodeIndex,
    suffix_v_index: NodeIndex,
    pub to_splits: LinkedHashSet<NodeIndex>,
    split_graph: SeqGraph<E>,
    new_middles: Vec<NodeIndex>,
    edges_to_remove: Vec<EdgeIndex>,
}

impl<'a, E: BaseEdge> SharedVertexSequenceSplitter<'a, E> {
    /**
     * Create a utility that will change the given graph so that the vertices in toSplitsArg (with their shared suffix and prefix
     * sequences) are extracted out.
     *
     * @param graph the graph containing the vertices in toSplitsArg
     * @param toSplitsArg a collection of vertices to split.  Must be contained within graph, and have only connections
     *                    from a single shared top and/or bottom node
     */
    pub fn new(
        graph: &'a mut SeqGraph<E>,
        to_split_args: LinkedHashSet<NodeIndex>,
    ) -> SharedVertexSequenceSplitter<'a, E> {
        assert!(
            to_split_args.len() > 1,
            "Can only split at least 2 vertices but only got {}",
            to_split_args.len()
        );
        assert!(
            graph.base_graph.contains_all_vertices(&to_split_args),
            "Graph doesn't contain all the provided vertices"
        );

        let (prefix, suffix) = Self::common_prefix_and_suffix_of_vertices(&to_split_args, &graph);

        Self {
            outer: graph,
            prefix_v: prefix,
            suffix_v: suffix,
            prefix_v_index: NodeIndex::new(0),
            suffix_v_index: NodeIndex::new(1),
            to_splits: to_split_args,
            split_graph: SeqGraph::new(11),
            new_middles: Vec::new(),
            edges_to_remove: Vec::new(),
        }
    }

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @see #updateGraph(SeqVertex, SeqVertex) for a full description of top and bottom
     *
     * @param top the top vertex, may be null
     * @param bottom the bottom vertex, may be null
     * @return true if some useful splitting was done, false otherwise
     */
    pub fn split_and_update(mut self, top: Option<NodeIndex>, bottom: Option<NodeIndex>) -> bool {
        self.split();
        self.update_graph(top, bottom);
        return true;
    }

    /**
     * Does either the common suffix or prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if either suffix or prefix length >= minCommonSequence
     */
    pub fn meets_min_mergable_sequence_for_either_prefix_or_suffix(
        &self,
        min_common_sequence: usize,
    ) -> bool {
        return self.meets_min_mergable_sequence_for_prefix(min_common_sequence)
            || self.meets_min_mergable_sequence_for_suffix(min_common_sequence);
    }

    /**
     * Does the common prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if prefix length >= minCommonSequence
     */
    pub fn meets_min_mergable_sequence_for_prefix(&self, min_common_sequence: usize) -> bool {
        return self.prefix_v.len() >= min_common_sequence;
    }

    /**
     * Does the common suffix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if suffix length >= minCommonSequence
     */
    pub fn meets_min_mergable_sequence_for_suffix(&self, min_common_sequence: usize) -> bool {
        return self.suffix_v.len() >= min_common_sequence;
    }

    /**
     * Actually do the splitting up of the vertices
     *
     * Must be called before calling updateGraph
     */
    pub fn split(&mut self) {
        self.split_graph = SeqGraph::new(self.outer.base_graph.get_kmer_size());
        self.edges_to_remove = Vec::new();
        self.new_middles = Vec::new();

        let fix_indices = self
            .split_graph
            .base_graph
            .add_vertices(vec![&self.prefix_v, &self.suffix_v]);
        self.prefix_v_index = fix_indices[0];
        self.suffix_v_index = fix_indices[1];
        let to_splits = self.to_splits.clone();
        for mid in to_splits {
            let mut to_mid =
                self.process_edge_to_remove(mid, self.outer.base_graph.incoming_edge_of(mid));
            let from_mid =
                self.process_edge_to_remove(mid, self.outer.base_graph.outgoing_edge_of(mid));

            let remaining = self
                .outer
                .base_graph
                .graph
                .node_weight(mid)
                .unwrap()
                .without_prefix_and_suffix(
                    self.prefix_v.get_sequence(),
                    self.suffix_v.get_sequence(),
                );
            match remaining {
                Some(remaining) => {
                    // there's some sequence prefix + seq + suffix, so add the node and make edges
                    let remaining_index = self.split_graph.base_graph.add_node(&remaining);
                    self.new_middles.push(remaining_index);

                    // update edge from top -> middle to be top -> without suffix
                    self.split_graph.base_graph.graph.add_edge(
                        self.prefix_v_index,
                        remaining_index,
                        to_mid,
                    );
                    self.split_graph.base_graph.graph.add_edge(
                        remaining_index,
                        self.suffix_v_index,
                        from_mid,
                    );
                }
                None => {
                    // prefix + suffix completely explain this node
                    to_mid.add(from_mid);
                    self.split_graph.base_graph.add_or_update_edge(
                        self.prefix_v_index,
                        self.suffix_v_index,
                        to_mid,
                    );
                }
            };
        }
    }

    /**
     * Update graph outer, replacing the previous middle vertices that were split out with the new
     * graph structure of the split, linking this subgraph into the graph at top and bot (the
     * vertex connecting the middle nodes and the vertex outgoing of all middle node)
     *
     * @param top an optional top node that must have outgoing edges to all split vertices.  If null, this subgraph
     *            will be added without any incoming edges
     * @param bot an optional bottom node that must have incoming edges to all split vertices.  If null, this subgraph
     *            will be added without any outgoing edges to the rest of the graph
     */
    pub fn update_graph(&mut self, top: Option<NodeIndex>, bot: Option<NodeIndex>) {
        assert!(
            self.outer.base_graph.contains_all_vertices(&self.to_splits),
            "Graph does not contain all of the original vertices to split"
        );
        assert!(
            top.is_some() || bot.is_some(),
            "Cannot update graph without at least one top or bot vertex, but both were None"
        );

        self.outer.base_graph.remove_all_vertices(&self.to_splits);
        self.outer
            .base_graph
            .remove_all_edges(&self.edges_to_remove);

        let new_nodes = self.outer.base_graph.add_vertices(
            self.split_graph
                .base_graph
                .get_node_weights(&self.new_middles)
                .into_iter()
                .map(|v| v.unwrap())
                .collect::<Vec<&SeqVertex>>(),
        );

        let has_prefix_suffix_edge = self
            .split_graph
            .base_graph
            .graph
            .contains_edge(self.prefix_v_index, self.suffix_v_index);
        let has_only_prefix_suffix_edges = has_prefix_suffix_edge
            && self
                .split_graph
                .base_graph
                .out_degree_of(self.prefix_v_index)
                == 1;
        let need_prefix_node =
            !self.prefix_v.is_empty() || (top.is_none() && !has_only_prefix_suffix_edges);
        let need_suffix_node =
            !self.suffix_v.is_empty() || (bot.is_none() && !has_only_prefix_suffix_edges);


        let mut outer_prefix_index = None;
        if need_prefix_node {
            outer_prefix_index = Some(self.add_prefix_node_and_edges(top));
        }

        let mut outer_suffix_index = None;
        if need_suffix_node {
            outer_suffix_index = Some(self.add_suffix_node_and_edges(bot));
        }

        // if prefix / suffix are needed, keep them
        // This needs to be run after we add the prefix and suffix nodes
        let top_for_connect = if need_prefix_node {
            outer_prefix_index
        } else {
            top
        };

        let bot_for_connect = if need_suffix_node {
            outer_suffix_index
        } else {
            bot
        };

        if top_for_connect.is_some() {
            self.add_edges_from_top_node(top_for_connect.unwrap(), bot_for_connect);
        };

        if bot_for_connect.is_some() {
            self.add_edges_to_bottom_node(bot_for_connect.unwrap());
        }
    }

    fn add_edges_to_bottom_node(&mut self, bot_for_connect: NodeIndex) {
        let incoming_edges = self
            .split_graph
            .base_graph
            .graph
            .edges_directed(self.suffix_v_index, Direction::Incoming)
            .map(|e| (e.id(), e.weight().clone()))
            .collect::<Vec<(EdgeIndex, E)>>();
        for (e_index, e_weight) in incoming_edges {
            let outer_source_index = self.outer.base_graph.index_from_vertex(
                &self.split_graph.base_graph.graph
                    [self.split_graph.base_graph.get_edge_source(e_index)],
            );
            match outer_source_index {
                Some(outer_source_index) => {
                    self.outer.base_graph.graph.add_edge(
                        outer_source_index,
                        bot_for_connect,
                        e_weight,
                    );
                }
                None => {
                    // pass
                }
            }
        }
    }

    fn add_edges_from_top_node(
        &mut self,
        top_for_connect: NodeIndex,
        bot_for_connect: Option<NodeIndex>,
    ) {
        // clone and collect because we can't mutate self while holding any references
        let outgoing = self
            .split_graph
            .base_graph
            .graph
            .edges_directed(self.prefix_v_index, Direction::Outgoing)
            .map(|e| (e.id(), e.weight().clone()))
            .collect::<Vec<(EdgeIndex, E)>>();

        for (e_index, e_weight) in outgoing {
            let target = self.split_graph.base_graph.get_edge_target(e_index);
            if target == self.suffix_v_index {
                // going straight from prefix -> suffix
                match bot_for_connect {
                    Some(bot_for_connect) => {
                        self.outer.base_graph.graph.add_edge(
                            top_for_connect,
                            bot_for_connect,
                            e_weight,
                        );
                    }
                    None => {
                        // pass
                    }
                }
            } else {
                let outer_target = self
                    .outer
                    .base_graph
                    .index_from_vertex(&self.split_graph.base_graph.graph[target]);

                match outer_target {
                    Some(outer_target) => {
                        self.outer.base_graph.graph.add_edge(
                            top_for_connect,
                            outer_target,
                            e_weight,
                        );
                    }
                    None => {
                        // pass
                    }
                }
            }
        }
    }

    fn add_suffix_node_and_edges(&mut self, bot: Option<NodeIndex>) -> NodeIndex {
        let suffix_v_index = self.outer.base_graph.add_node(&self.suffix_v);
        match bot {
            Some(bot) => {
                self.outer.base_graph.graph.add_edge(
                    suffix_v_index,
                    bot,
                    BaseEdge::make_o_r_edge(
                        self.split_graph
                            .base_graph
                            .graph
                            .edges_directed(self.suffix_v_index, Direction::Incoming)
                            .map(|e| e.weight().clone())
                            .collect::<Vec<E>>(),
                        1,
                        0,
                    ),
                );
            }
            None => {
                // pass
            }
        }

        suffix_v_index
    }

    fn add_prefix_node_and_edges(&mut self, top: Option<NodeIndex>) -> NodeIndex {
        let prefix_v_index = self.outer.base_graph.add_node(&self.prefix_v);
        match top {
            Some(top) => {
                self.outer.base_graph.graph.add_edge(
                    top,
                    prefix_v_index,
                    BaseEdge::make_o_r_edge(
                        self.split_graph
                            .base_graph
                            .graph
                            .edges_directed(self.prefix_v_index, Direction::Outgoing)
                            .map(|e| e.weight().clone())
                            .collect::<Vec<E>>(),
                        1,
                        0,
                    ),
                );
            }
            None => {
                // pass
            }
        }

        prefix_v_index
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return
     */
    pub fn common_prefix_and_suffix_of_vertices(
        middle_vertices: &LinkedHashSet<NodeIndex>,
        graph: &SeqGraph<E>,
    ) -> (SeqVertex, SeqVertex) {
        let mut kmers = Vec::new();
        let mut min = std::usize::MAX;
        for v in middle_vertices {
            let sequence = graph.base_graph.get_sequence_from_index(*v);
            min = std::cmp::min(min, sequence.len());
            kmers.push(sequence);
        }

        let prefix_len = GraphUtils::common_maximum_prefix_length(&kmers);
        let suffix_len = GraphUtils::common_maximum_suffix_length(&kmers, min - prefix_len);

        let kmer = &kmers[0];
        let prefix = kmer[0..prefix_len].to_vec();
        let suffix = kmer[kmer.len() - suffix_len..kmer.len()].to_vec();
        return (SeqVertex::new(prefix), SeqVertex::new(suffix));
    }

    /**
     * Helper function that returns an edge that we should use for splitting
     *
     * If e is null, creates a new 0 multiplicity edge, set to ref is any edges to V are ref
     * If e is not null, returns a new copy of e, and schedules e for removal
     *
     * @param e a non-null edge
     * @return a non-null edge
     */
    fn process_edge_to_remove(&mut self, v: NodeIndex, e: Option<EdgeIndex>) -> E {
        match e {
            None => {
                // there's no edge, so we return a newly allocated one and don't schedule e for removal
                // the weight must be 0 to preserve sum through the diamond

                let mut edge_template = E::new(self.outer.base_graph.is_reference_node(v), 0, 0);
                // debug!("Edge template {:?}", &e);
                return edge_template;
            }
            Some(e) => {
                // schedule edge for removal, and return a freshly allocated one for our graph to use
                self.edges_to_remove.push(e);
                return self.outer.base_graph.graph.edge_weight(e).cloned().unwrap();
            }
        }
    }

    pub fn get_prefix(&self) -> &SeqVertex {
        &self.prefix_v
    }

    pub fn get_prefix_index(&self) -> NodeIndex {
        self.prefix_v_index
    }

    pub fn get_suffix(&self) -> &SeqVertex {
        &self.suffix_v
    }

    pub fn get_suffix_index(&self) -> NodeIndex {
        self.suffix_v_index
    }

    pub fn get_split_graph(&self) -> &SeqGraph<E> {
        &self.split_graph
    }

    pub fn get_new_middles(&self) -> &Vec<NodeIndex> {
        &self.new_middles
    }
}
