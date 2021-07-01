use graphs::base_graph::BaseGraph;
use rayon::prelude::*;
use std::collections::HashSet;
use petgraph::Direction;
use petgraph::graph::EdgeReference;
use graphs::base_edge::BaseEdge;
use graphs::path::Path;
use itertools::Itertools;

pub trait ChainPruner {

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
            let chain_start = chain_starts.pop().unwrap();
            for out_edge in graph.graph.edges_directed(chain_start, Direction::Outgoing) {
                let chain = Self::find_chain(out_edge, graph);
                let chain_end = chain
            }
        }
    }

    pub fn find_chain(start_edge: EdgeReference<BaseEdge, u32>, graph: &BaseGraph) -> Path {
        let mut edges = Vec::new();
        edges.push(start_edge.id());
        let first_vertex_id = start_edge.source();
        let mut last_vertex_id = start_edge.target();

        loop {
            let out_edges = graph.graph
                .edges_directed(last_vertex_id, Direction::Outgoing).collect_vec();
            if out_edges.len() != 1 || graph.in_degree_of(last_vertex_id) > 1 || last_vertex_id == first_vertex_id {
                break
            }
            let next_edge = out_edges[0];
            edges.push(next_edge);
            last_vertex_id = next_edge.target();
        }

        return Path::new(last_vertex_id, edges)
    }
}