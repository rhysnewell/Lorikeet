use graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use graphs::common_suffix_splitter::CommonSuffixSplitter;
use graphs::seq_graph::SeqGraph;
use graphs::shared_sequence_merger::SharedSequenceMerger;
use graphs::shared_vertex_sequence_splitter::SharedVertexSequenceSplitter;
use graphs::vertex_based_transformer::VertexBasedTransformer::{
    MergeCommonSuffices, MergeDiamonds, MergeTails, SplitCommonSuffices,
};
use hashlink::LinkedHashSet;
use petgraph::stable_graph::NodeIndex;
use petgraph::Direction;
use std::collections::HashSet;

/**
 * Base enum for transformation operations that need to iterate over proposed vertices, where
 * each proposed vertex is a seed vertex for a potential transformation.
 *
 * transformUntilComplete will iteratively apply the tryToTransform function on each vertex in the graph
 * until no vertex can be found that can be transformed.
 *
 * Note that in order to eventually terminate tryToTransform must transform the graph such that eventually
 * no vertices are candidates for further transformations.
 */
pub enum VertexBasedTransformer<'a, E: BaseEdge> {
    /**
     * Merge diamond configurations:
     *
     * Performs the transformation:
     *
     * { A -> x + S_i + y -> Z }
     *
     * goes to:
     *
     * { A -> x -> S_i -> y -> Z }
     *
     * for all nodes that match this configuration.
     */
    MergeDiamonds {
        graph: &'a mut SeqGraph<E>,
        dont_modify_graph_even_if_possible: bool,
    },
    /**
     * Merge tail configurations:
     *
     * Performs the transformation:
     *
     * { A -> x + S_i + y }
     *
     * goes to:
     *
     * { A -> x -> S_i -> y }
     *
     * for all nodes that match this configuration.
     *
     * Differs from the diamond transform in that no bottom node is required
     */
    MergeTails {
        graph: &'a mut SeqGraph<E>,
        dont_modify_graph_even_if_possible: bool,
    },
    /**
     * Performs the transformation:
     *
     * { x + S_i + y -> Z }
     *
     * goes to:
     *
     * { x -> S_i -> y -> Z }
     *
     * for all nodes that match this configuration.
     *
     * Differs from the diamond transform in that no top node is required
     */
    SplitCommonSuffices {
        graph: &'a mut SeqGraph<E>,
        dont_modify_graph_even_if_possible: bool,
        already_split: &'a mut LinkedHashSet<NodeIndex>,
    },
    /**
     * Merge headless configurations:
     *
     * Performs the transformation:
     *
     * { x + S_i -> y -> Z }
     *
     * goes to:
     *
     * { x -> S_i -> y + Z }
     *
     * for all nodes that match this configuration.
     */
    MergeCommonSuffices {
        graph: &'a mut SeqGraph<E>,
        dont_modify_graph_even_if_possible: bool,
    },
}

pub enum VertexBasedTransformerOptions {
    MergeDiamonds,
    MergeCommonSuffices,
    MergeTails,
    SplitCommonSuffices,
}

impl VertexBasedTransformerOptions {
    /**
     * Merge until the graph has no vertices that are candidates for merging
     */
    pub fn transform_until_complete<E: BaseEdge>(&mut self, graph: &mut SeqGraph<E>) -> bool {
        match self {
            Self::MergeDiamonds => {
                let mut did_at_least_one_transform = false;
                let mut found_nodes_to_merge = true;
                let mut transformed_count = 0;
                while found_nodes_to_merge {
                    found_nodes_to_merge = false;
                    let node_indices = graph
                        .base_graph
                        .graph
                        .node_indices()
                        .collect::<LinkedHashSet<NodeIndex>>();
                    for v in node_indices {
                        found_nodes_to_merge = MergeDiamonds {
                            graph,
                            dont_modify_graph_even_if_possible: false,
                        }
                        .try_to_transform(v);
                        if found_nodes_to_merge {
                            did_at_least_one_transform = true;
                            transformed_count += 1;
                            debug!(
                                "Transformed {} graph: V {} E {}",
                                transformed_count,
                                graph.base_graph.graph.node_count(),
                                graph.base_graph.graph.edge_count()
                            );
                            break;
                        }
                    }
                }

                return did_at_least_one_transform;
            }
            Self::SplitCommonSuffices => {
                let mut did_at_least_one_transform = false;
                let mut found_nodes_to_merge = true;
                let mut already_split = LinkedHashSet::new();
                let mut transformed_count = 0;
                while found_nodes_to_merge {
                    found_nodes_to_merge = false;
                    let node_indices = graph
                        .base_graph
                        .graph
                        .node_indices()
                        .collect::<LinkedHashSet<NodeIndex>>();
                    for v in node_indices {
                        found_nodes_to_merge = SplitCommonSuffices {
                            graph,
                            dont_modify_graph_even_if_possible: false,
                            already_split: &mut already_split,
                        }
                        .try_to_transform(v);

                        if found_nodes_to_merge {
                            did_at_least_one_transform = true;
                            transformed_count += 1;
                            debug!(
                                "Transformed {} graph: V {} E {}",
                                transformed_count,
                                graph.base_graph.graph.node_count(),
                                graph.base_graph.graph.edge_count()
                            );
                            break;
                        }
                    }
                }

                return did_at_least_one_transform;
            }
            Self::MergeCommonSuffices => {
                let mut did_at_least_one_transform = false;
                let mut found_nodes_to_merge = true;
                let mut transformed_count = 0;
                while found_nodes_to_merge {
                    found_nodes_to_merge = false;
                    let node_indices = graph
                        .base_graph
                        .graph
                        .node_indices()
                        .collect::<LinkedHashSet<NodeIndex>>();
                    for v in node_indices {
                        found_nodes_to_merge = MergeCommonSuffices {
                            graph,
                            dont_modify_graph_even_if_possible: false,
                        }
                        .try_to_transform(v);
                        if found_nodes_to_merge {
                            did_at_least_one_transform = true;
                            transformed_count += 1;
                            debug!(
                                "Transformed {} graph: V {} E {}",
                                transformed_count,
                                graph.base_graph.graph.node_count(),
                                graph.base_graph.graph.edge_count()
                            );
                            break;
                        }
                    }
                }

                return did_at_least_one_transform;
            }
            Self::MergeTails => {
                let mut did_at_least_one_transform = false;
                let mut found_nodes_to_merge = true;
                let mut transformed_count = 0;
                while found_nodes_to_merge {
                    found_nodes_to_merge = false;
                    let node_indices = graph
                        .base_graph
                        .graph
                        .node_indices()
                        .collect::<LinkedHashSet<NodeIndex>>();
                    for v in node_indices {
                        found_nodes_to_merge = MergeTails {
                            graph,
                            dont_modify_graph_even_if_possible: false,
                        }
                        .try_to_transform(v);
                        if found_nodes_to_merge {
                            did_at_least_one_transform = true;
                            transformed_count += 1;
                            debug!(
                                "Transformed {} graph: V {} E {}",
                                transformed_count,
                                graph.base_graph.graph.node_count(),
                                graph.base_graph.graph.edge_count()
                            );
                            break;
                        }
                    }
                }

                return did_at_least_one_transform;
            }
        }
    }
}

impl<'a, E: BaseEdge> VertexBasedTransformer<'a, E> {
    // The minimum number of common bp from the prefix (head merging) or suffix (tail merging)
    // required before we'll merge in such configurations.  A large value here is critical to avoid
    // merging inappropriate head or tail nodes, which introduces large insertion / deletion events
    // as the merge operation creates a link among the non-linked sink / source vertices
    pub const MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES: usize = 10;

    /**
     * Merge, if possible, seeded on the vertex v
     * @param v the proposed seed vertex to merge
     * @return true if some useful merging happened, false otherwise
     */
    pub fn try_to_transform(&mut self, v: NodeIndex) -> bool {
        match self {
            Self::MergeDiamonds {
                graph,
                dont_modify_graph_even_if_possible,
            } => {
                let middles = graph.base_graph.outgoing_vertices_of(v);
                if middles.len() <= 1 {
                    // we can only merge if there's at least two middle nodes
                    return false;
                }

                let mut bottom = None;
                for mi in middles.iter() {
                    // all nodes must have at least 1 connection
                    if graph.base_graph.out_degree_of(*mi) < 1 {
                        return false;
                    }

                    // can only have 1 incoming node, the root vertex
                    if graph.base_graph.in_degree_of(*mi) != 1 {
                        return false;
                    }

                    // make sure that all outgoing vertices of mi go only to the bottom node
                    for mt in graph.base_graph.outgoing_vertices_of(*mi) {
                        bottom = match bottom {
                            None => Some(mt),
                            Some(bottom) => {
                                if bottom != mt {
                                    return false;
                                } else {
                                    Some(bottom)
                                }
                            }
                        }
                    }
                }

                let bottom = bottom.unwrap();
                // bottom has some connections coming in from other nodes, don't allow
                if graph.base_graph.in_degree_of(bottom) != middles.len() {
                    return false;
                }

                if *dont_modify_graph_even_if_possible {
                    return true;
                }

                // actually do the merging, returning true if at least 1 base was successfully split
                let mut splitter = SharedVertexSequenceSplitter::new(graph, middles);
                return splitter.meets_min_mergable_sequence_for_either_prefix_or_suffix(1)
                    && splitter.split_and_update(Some(v), Some(bottom));
            }
            Self::MergeTails {
                graph,
                dont_modify_graph_even_if_possible,
            } => {
                let tails = graph.base_graph.outgoing_vertices_of(v);
                if tails.len() <= 1 {
                    return false;
                };

                for t in tails.iter() {
                    if !graph.base_graph.is_sink(*t) || graph.base_graph.in_degree_of(*t) > 1 {
                        return false;
                    }
                }

                if *dont_modify_graph_even_if_possible {
                    return true;
                };

                let mut splitter = SharedVertexSequenceSplitter::new(graph, tails);
                return splitter.meets_min_mergable_sequence_for_either_prefix_or_suffix(
                    Self::MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES,
                ) && splitter.split_and_update(Some(v), None);
            }
            Self::SplitCommonSuffices {
                graph,
                dont_modify_graph_even_if_possible,
                already_split,
            } => {
                if already_split.contains(&v) {
                    return false;
                } else {
                    already_split.insert(v);
                    return CommonSuffixSplitter::split(graph, v);
                }
                // return CommonSuffixSplitter::split(graph, v);
            }
            Self::MergeCommonSuffices {
                graph,
                dont_modify_graph_even_if_possible,
            } => return SharedSequenceMerger::merge(graph, v),
        }
    }
}
