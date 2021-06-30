use graphs::base_graph::BaseGraph;
use rayon::prelude::*;
use std::collections::HashSet;

pub struct ChainPruner {

}

impl ChainPruner {
    pub fn prune_low_weight_chain(graph: &BaseGraph) {
        let chains = Self::find_all_chains(graph);
    }

    pub fn find_all_chains(graph: &BaseGraph) {
        let mut chain_starts = graph.get_sources();
        let mut chains = Vec::new();
        let mut already_seen = HashSet::new();

        while !chain_starts.is_empty() {
            let chain_start = chain_starts.pop_first()
        }
    }
}