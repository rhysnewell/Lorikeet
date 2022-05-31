use graphs::path::Path;
use graphs::base_graph::BaseGraph;
use graphs::base_edge::BaseEdge;
use graphs::base_vertex::BaseVertex;

/**
 * Prune all chains from this graph where all edges in the path have multiplicity < pruneFactor
 *
 *
 * For A -[1]> B -[1]> C -[1]> D would be removed with pruneFactor 2
 * but A -[1]> B -[2]> C -[1]> D would not be because the linear chain includes an edge with weight >= 2
 *
 */
#[derive(Debug, Clone)]
pub struct LowWeightChainPruner {
    pub(crate) prune_factor: usize
}

impl LowWeightChainPruner {
    pub fn new(prune_factor: usize) -> Self {
        Self {
            prune_factor
        }
    }

    pub fn needs_pruning<
        'a,
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync,
    >(&self, graph: &BaseGraph<V, E>, chain: &Path) -> bool {
        chain.get_edges().iter().all(|e| match graph.graph.edge_weight(*e) {
            None => panic!("Edge index not in graph"),
            Some(edge) => {
                edge.get_pruning_multiplicity() < self.prune_factor && !edge.is_ref()
            }
        })
    }
}