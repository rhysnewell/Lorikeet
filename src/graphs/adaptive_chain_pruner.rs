use rayon::prelude::*;
use graphs::base_graph::BaseGraph;
use graphs::path::Path;
use utils::math_utils::MathUtils;
use petgraph::Direction;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;

pub struct AdaptiveChainPruner {
    initial_error_probability: f64,
    log_odds_threshold: f64,
    seeding_log_odds_threshold: f64, // threshold for seeding subgraph of good vertices
    max_unpruned_variants: usize
}

impl AdaptiveChainPruner {
    pub fn new(
        initial_error_probability: f64,
        log_odds_threshold: f64,
        seeding_log_odds_threshold: f64,
        max_unpruned_variants: usize
    ) -> AdaptiveChainPruner {

        AdaptiveChainPruner {
            initial_error_probability,
            log_odds_threshold,
            seeding_log_odds_threshold,
            max_unpruned_variants
        }
    }

    fn likely_error_chains(&self, chains: Vec<Path>, graph: &BaseGraph, error_rate: f64) -> Vec<Path> {
        let chain_log_odds = chains.par_iter().map(|chain| {
            self.chain_log_odds(chain, graph, error_rate)
        }).collect::<Vec<(f64, f64)>>();

        //TODO: Finish this function
    }

    // left and right chain log odds
    fn chain_log_odds(&self, chain: &Path, graph: &BaseGraph, error_rate: f64) -> (f64, f64) {
        let left_total_multiplicity = graph.graph.edges_directed(chain.get_first_vertex(), Direction::Outgoing).map(|e| e.weight().get_multiplicity()).sum();
        let right_total_multiplicity = graph.graph.edges_directed(chain.get_last_vertex(), Direction::Incoming).map(|e| e.weight().get_multiplicity()).sum();

        let left_multiplicity = graph.graph.edge_weight(chain.get_edges()[0]).unwrap().get_multiplicity();
        let right_multiplicity = graph.graph.edge_weight(chain.get_last_edge()).unwrap().get_multiplicity();

        let left_log_odds = if graph.is_source(chain.get_first_vertex()) { 0.0 } else {
            HaplotypeCallerEngine::log_likelihood_ratio_constant_error(left_total_multiplicity - left_multiplicity, left_multiplicity, error_rate)
        };
        let right_log_odds = if graph.is_source(chain.get_last_vertex()) { 0.0 } else {
            HaplotypeCallerEngine::log_likelihood_ratio_constant_error(right_total_multiplicity - right_multiplicity, right_multiplicity, error_rate)
        };

        return (left_log_odds, right_log_odds)
    }
}