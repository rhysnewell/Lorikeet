use graphs::base_graph::BaseGraph;
use graphs::base_edge::BaseEdge;
use graphs::base_vertex::BaseVertex;
use graphs::path::Path;
use petgraph::stable_graph::EdgeIndex;

pub trait ChainPruner<V: BaseVertex, E: BaseEdge> {
    fn prune_low_weight_chains(&self, graph: &mut BaseGraph<V, E>);

    fn find_all_chains(graph: &BaseGraph<V, E>) -> Vec<Path<'_, V, E>>;

    fn find_chain<'a>(start_edge: EdgeIndex, graph: &'a BaseGraph<V, E>) -> Path<'a, V, E>;

    fn chains_to_remove<'a>(&self, chains: &'a Vec<Path<'a, V, E>>, graph: &'a BaseGraph<V, E>) -> Vec<&'a Path<'a, V, E>>;
}