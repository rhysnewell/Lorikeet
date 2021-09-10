use graphs::seq_graph::SeqGraph;
use graphs::base_edge::BaseEdgeStruct;
use petgraph::stable_graph::NodeIndex;

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
pub enum VertexBasedTransformer<'a> {
    MergeDiamonds {
        graph: &'a mut SeqGraph<BaseEdgeStruct>,
        dont_modify_graph_even_if_possible: bool,
    },
    MergeTails {
        graph: &'a mut SeqGraph<BaseEdgeStruct>,
        dont_modify_graph_even_if_possible: bool,
    },
    SplitCommonSuffices {
        graph: &'a mut SeqGraph<BaseEdgeStruct>,
        dont_modify_graph_even_if_possible: bool,
    },
    MergeCommonSuffices {
        graph: &'a mut SeqGraph<BaseEdgeStruct>,
        dont_modify_graph_even_if_possible: bool,
    },
}

impl<'a> VertexBasedTransformer<'a> {

    /**
     * Merge until the graph has no vertices that are candidates for merging
     */
    pub fn transform_until_complete(&mut self) -> bool {
        match self {
            Self::MergeDiamonds {
                graph,
                ..
            } | Self::MergeCommonSuffices {
                graph,
                ..
            } | Self::MergeTails {
                graph,
                ..
            } | Self::SplitCommonSuffices {
                graph,
                ..
            } => {
                let mut did_at_least_one_transform = false;
                let mut found_nodes_to_merge = true;
                while found_nodes_to_merge {
                    found_nodes_to_merge = false;
                    for v in graph.base_graph.graph.node_indices() {
                        found_nodes_to_merge = self.try_to_transform(v);
                        if found_nodes_to_merge {
                            did_at_least_one_transform = true;
                            break
                        }
                    }
                }
            }
        }
    }

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
                if middles.len() <= 1 { // we can only merge if there's at least two middle nodes
                    return false
                }

                let mut bottom = None;
                for mi in middles.iter() {
                    // all nodes must have at least 1 connection
                    if graph.base_graph.out_degree_of(*mi) < 1 {
                        return false
                    }

                    // can only have 1 incoming node, the root vertex
                    if graph.base_graph.in_degree_of(*mi) != 1 {
                        return false
                    }

                    // make sure that all outgoing vertices of mi go only to the bottom node
                    for mt in graph.base_graph.outgoing_vertices_of(*mi) {
                        bottom = match bottom {
                            None => Some(mt),
                            Some(bottom) => {
                                if bottom != mt {
                                    return false
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
                    return false
                }

                if dont_modify_graph_even_if_possible {
                    return true
                }

                // actually do the merging, returning true if at least 1 base was successfully split
                let
            }
        }
    }
}