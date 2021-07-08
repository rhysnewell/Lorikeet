use graphs::base_graph::BaseGraph;
use rayon::prelude::*;
use std::collections::HashSet;
use petgraph::Direction;
use petgraph::graph::EdgeReference;
use graphs::base_edge::BaseEdge;
use graphs::base_vertex::BaseVertex;
use graphs::path::Path;
use itertools::Itertools;
use graphs::adaptive_chain_pruner::AdaptiveChainPruner;
use petgraph::graph::EdgeIndex;
use petgraph::visit::EdgeRef;

pub trait ChainPruner<V: BaseVertex, E: BaseEdge> {
    fn prune_low_weight_chain(&self, graph: &mut BaseGraph<V, E>);

    fn find_all_chains(graph: &BaseGraph<V, E>) -> Vec<Path<'_, V, E>>;

    fn find_chain<'a>(start_edge: EdgeIndex, graph: &'a BaseGraph<V, E>) -> Path<'a, V, E>;

    fn chains_to_remove<'a>(&self, chains: &'a Vec<Path<'a, V, E>>, graph: &'a BaseGraph<V, E>) -> Vec<&'a Path<'a, V, E>>;
}

impl<V: BaseVertex + std::marker::Sync, E: BaseEdge + std::marker::Sync> ChainPruner<V, E> for AdaptiveChainPruner {
    fn prune_low_weight_chain(&self, graph: &mut BaseGraph<V, E>) {
        let chains = Self::find_all_chains(graph);
        let chains_to_remove = self.chains_to_remove(&chains, graph);
        chains_to_remove.iter().for_each(|chain| graph.remove_all_edges(chain.get_edges()));
        graph.remove_singleton_orphan_vertices();
    }

    fn find_all_chains(graph: &BaseGraph<V, E>) -> Vec<Path<'_, V, E>> {
        let mut chain_starts = graph.get_sources();
        let mut chains = Vec::new();
        let mut already_seen = HashSet::new();

        while !chain_starts.is_empty() {
            let chain_start = chain_starts.pop().unwrap();
            for out_edge in graph.graph.edges_directed(chain_start, Direction::Outgoing) {
                let chain = Self::find_chain(out_edge.id(), graph);
                let chain_end = chain.get_last_vertex();
                chains.push(chain);
                if !already_seen.contains(&chain_end) {
                    chain_starts.push(chain_end);
                    already_seen.insert(chain_end);
                }
            }
        }
        return chains
    }

    fn find_chain<'a>(start_edge: EdgeIndex, graph: &'a BaseGraph<V, E>) -> Path<'a, V, E> {
        let mut edges = Vec::new();
        edges.push(start_edge);
        let start_edge_endpoints = graph.graph.edge_endpoints(start_edge).unwrap();
        let first_vertex_id = start_edge_endpoints.0;
        let mut last_vertex_id = start_edge_endpoints.1;

        loop {
            let out_edges = graph.graph
                .edges_directed(last_vertex_id, Direction::Outgoing).collect_vec();
            if out_edges.len() != 1 || graph.in_degree_of(last_vertex_id) > 1 || last_vertex_id == first_vertex_id {
                break
            }
            let next_edge = out_edges[0];
            edges.push(next_edge.id());
            last_vertex_id = graph.graph.edge_endpoints(next_edge.id()).unwrap().1;
        }

        return Path::new(last_vertex_id, edges, graph)
    }

    fn chains_to_remove<'a>(&self, chains: &'a Vec<Path<'a, V, E>>, graph: &'a BaseGraph<V, E>) -> Vec<&'a Path<'a, V, E>> {
        if chains.is_empty() {
            return Vec::new()
        }
        let probable_error_chains = self.likely_error_chains(&chains, graph, self.initial_error_probability);
        let error_count = probable_error_chains.into_par_iter().map(|chain| {
            chain.get_edges().par_iter().map(|e| chain.graph.graph.edge_weight(*e).unwrap().get_multiplicity()).sum::<usize>()
        }).sum::<usize>();
        let total_bases = chains.par_iter().map(|chain| {
            chain.get_edges().par_iter().map(|e| chain.graph.graph.edge_weight(*e).unwrap().get_multiplicity()).sum::<usize>()
        }).sum::<usize>();

        let error_rate = error_count as f64 / total_bases as f64;

        return self.likely_error_chains(&chains, graph, error_rate).into_iter().par_bridge().filter(|c| {
            !c.get_edges().par_iter().any(|e| c.graph.graph.edge_weight(*e).unwrap().is_ref())
        }).collect::<Vec<&Path<V, E>>>()
    }
}