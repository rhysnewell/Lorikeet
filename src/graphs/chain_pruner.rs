use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::path::Path;
use petgraph::stable_graph::EdgeIndex;
use std::collections::VecDeque;
use graphs::adaptive_chain_pruner::AdaptiveChainPruner;
use graphs::low_weight_chain_pruner::LowWeightChainPruner;
use rayon::prelude::*;
use petgraph::visit::EdgeRef;
use std::collections::HashSet;
use petgraph::stable_graph::NodeIndex;
use petgraph::Direction;


#[derive(Debug, Clone)]
pub enum ChainPruner {
    AdaptiveChainPruner(AdaptiveChainPruner),
    LowWeightChainPruner(LowWeightChainPruner)
}

impl ChainPruner {
    pub fn prune_low_weight_chains<
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync
    >(&self, graph: &mut BaseGraph<V, E>) {
        let chains = Self::find_all_chains(&graph);
        debug!("Chains {}", chains.len());

        let chains_to_remove = self.chains_to_remove(&chains, &graph);
        debug!("Chains to remove {}", chains_to_remove.len());
        chains_to_remove
            .into_iter()
            .for_each(|chain| graph.remove_all_edges(chain.get_edges()));

        graph.remove_singleton_orphan_vertices();
    }

    pub fn find_all_chains<
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync
    >(graph: &BaseGraph<V, E>) -> VecDeque<Path> {
        let mut chain_starts = graph.get_sources();
        let mut chains = VecDeque::new();
        let mut already_seen = chain_starts
            .clone()
            .into_iter()
            .collect::<HashSet<NodeIndex>>();
        while !chain_starts.is_empty() {
            let chain_start = chain_starts.pop_front().unwrap();
            for out_edge in graph.graph.edges_directed(chain_start, Direction::Outgoing) {
                let chain = Self::find_chain(out_edge.id(), graph);
                let chain_end = chain.get_last_vertex();
                chains.push_back(chain);
                if !already_seen.contains(&chain_end) {
                    chain_starts.push_back(chain_end);
                    already_seen.insert(chain_end);
                }
            }
        }
        chains
    }

    /**
     *
     * @return a fully extended linear path
     */
    pub fn find_chain<
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync
    >(start_edge: EdgeIndex, graph: &BaseGraph<V, E>) -> Path {
        let mut edges = vec![start_edge];
        let start_edge_endpoints = graph.graph.edge_endpoints(start_edge).unwrap();
        let first_vertex_id = start_edge_endpoints.0;
        let mut last_vertex_id = start_edge_endpoints.1;

        loop {
            // chain ends if:
            //  1) no out edges;
            //  2) multiple out edges;
            //  3) multiple in edges;
            //  4) cycle back to start of chain
            let out_edges = graph
                .graph
                .edges_directed(last_vertex_id, Direction::Outgoing)
                .map(|e| e.id())
                .collect::<HashSet<EdgeIndex>>();
            if out_edges.len() != 1
                || graph.in_degree_of(last_vertex_id) > 1
                || last_vertex_id == first_vertex_id
            {
                break;
            }
            let next_edge = out_edges.into_iter().next().unwrap();
            edges.push(next_edge);
            last_vertex_id = graph.graph.edge_endpoints(next_edge).unwrap().1;
        }

        Path::new(last_vertex_id, edges)
    }

    pub fn chains_to_remove<'a,
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync
    >(
        &self,
        chains: &'a VecDeque<Path>,
        graph: &BaseGraph<V, E>,
    ) -> Vec<&'a Path> {
        match self {
            ChainPruner::AdaptiveChainPruner(adaptive) => {
                if chains.is_empty() {
                    return Vec::new();
                }
                let probable_error_chains =
                    adaptive.likely_error_chains(&chains, graph, adaptive.initial_error_probability);
                debug!("Probable error chains {}", probable_error_chains.len());

                let error_count = probable_error_chains
                    .into_iter()
                    .map(|chain| {
                        // chain
                        //     .get_edges()
                        //     .par_iter()
                        //     .map(|e| graph.graph.edge_weight(*e).unwrap().get_multiplicity())
                        //     .sum::<usize>()
                        graph
                            .graph
                            .edge_weight(chain.get_last_edge())
                            .unwrap()
                            .get_multiplicity()
                    })
                    .sum::<usize>();
                debug!("Error count {}", error_count);

                let total_bases = chains
                    .iter()
                    .map(|chain| {
                        chain
                            .get_edges()
                            .iter()
                            .map(|e| graph.graph.edge_weight(*e).unwrap().get_multiplicity())
                            .sum::<usize>()
                    })
                    .sum::<usize>();

                debug!("Total bases {}", total_bases);
                let error_rate = error_count as f64 / total_bases as f64;
                debug!("Error rate {}", error_rate);
                adaptive.likely_error_chains(&chains, graph, error_rate)
                    .into_par_iter()
                    .filter(|c| {
                        !c.get_edges()
                            .iter()
                            .any(|e| graph.graph.edge_weight(*e).unwrap().is_ref())
                    })
                    .collect::<Vec<&Path>>()
            },
            ChainPruner::LowWeightChainPruner(low_weight) => {
                chains.into_par_iter().filter(|chain| low_weight.needs_pruning(graph, chain)).collect()
            }
        }

    }
}