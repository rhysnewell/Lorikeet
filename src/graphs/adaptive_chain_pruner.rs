use rayon::prelude::*;
use graphs::base_graph::BaseGraph;
use graphs::path::{Path, Chain};
use utils::math_utils::MathUtils;
use petgraph::Direction;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use multimap::MultiMap;
use std::collections::{BinaryHeap, HashMap, HashSet};
use comparator::{as_fn, comparing, Comparator};
use utils::base_utils::BaseUtils;
use linked_hash_set::LinkedHashSet;
use ordered_float::OrderedFloat;
use graphs::chain_pruner::ChainPruner;

pub struct AdaptiveChainPruner {
    pub initial_error_probability: f64,
    log_odds_threshold: f64,
    seeding_log_odds_threshold: f64, // threshold for seeding subgraph of good vertices
    max_unpruned_variants: usize
}

// TODO: This whole thing needs some unit tests, I don't trust parts of this code as there might be
//       some issues with refrence/ownership and also chained comparators.
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

    pub fn likely_error_chains(&self, chains: &Vec<Path>, graph: &BaseGraph, error_rate: f64) -> HashSet<Path> {
        // pre-compute the left and right log odds of each chain
        let chain_log_odds = chains.par_iter().map(|chain| {
            let result = self.chain_log_odds(chain, graph, error_rate);
            (chain, result)
        }).collect::<HashMap<Path, (f64, f64)>>();

        // compute correspondence of vertices to incident chains with log odds above the seeding and extending thresholds
        let mut vertex_to_seedable_chains = MultiMap::new();
        let mut vertex_to_good_incoming_chains = MultiMap::new();
        let mut vertex_to_good_outgoing_chains = MultiMap::new();

        for chain in chains.iter() {
            if chain_log_odds.get(chain).unwrap().1 >= self.log_odds_threshold {
               vertex_to_good_incoming_chains.insert(chain.get_last_vertex(), chain);
            }

            if chain_log_odds.get(chain).unwrap().0 >= self.log_odds_threshold {
                vertex_to_good_outgoing_chains.insert(chain.get_first_vertex(), chain);
            }

            // seed-worthy chains must pass the more stringent seeding log odds threshold on both sides
            // in addition to that, we only seed from vertices with multiple such chains incoming or outgoing (see below)
            if chain_log_odds.get(chain).unwrap().1 >= self.seeding_log_odds_threshold
                && chain_log_odds.get(chain).unwrap().0 >= self.seeding_log_odds_threshold {
                vertex_to_seedable_chains.insert(chain.get_first_vertex(), chain);
                vertex_to_seedable_chains.insert(chain.get_last_vertex(), chain);
            }
        }

        // We have a priority queue of chains to add to the graph, with priority given by the log odds (higher first)
        // for determinism we have a tie breaker based on chains' first vertex
        // Note that chains can we added twice to the queue, once for each side
        // The comparator for this heap is implemented in (Chain)[Chain]
        let mut chains_to_add = BinaryHeap::new();

        // seed the subgraph of good chains by starting with some definitely-good chains.  These include the max-weight chain
        // and chains emanating from vertices with two incoming or two outgoing chains (plus one outgoing or incoming for a total of 3 or more) with good log odds
        // The idea is that a high-multiplicity error chain A that branches into a second error chain B and a continuation-of-the-original-error chain A'
        // may have a high log odds for A'.  However, only in the case of true variation will multiple branches leaving the same vertex have good log odds.
        let max_weight_chain = Self::get_max_weight_chains(&chains, graph);
        chains_to_add.push(Chain::new(OrderedFloat::new(std::f64::INFINITY), &max_weight_chain, graph));
        let mut processed_vertices = LinkedHashSet::new();
        for (vertex, paths) in vertex_to_seedable_chains.iter() {
            if paths.len() > 2 {
                vertex_to_good_outgoing_chains.get(vertex).iter().for_each(|chain| {
                    chains_to_add.push(Chain::new(OrderedFloat::new(chain_log_odds.get(chain).unwrap().0), chain, graph))
                });
                vertex_to_good_incoming_chains.get(vertex).iter().for_each(|chain| {
                    chains_to_add.push(Chain::new(OrderedFloat::new(chain_log_odds.get(chain).unwrap().1), chain, graph))
                });
                processed_vertices.insert(*vertex);
            }
        }

        let mut good_chains = LinkedHashSet::new();
        let mut vertices_that_already_have_outgoing_good_chains = HashSet::new();
        let mut variant_count = 0;

        // starting from the high-confidence seed vertices, grow the "good" subgraph along chains with above-threshold log odds,
        // discovering good chains as we go.
        while !chains_to_add.is_empty() && variant_count <= self.max_unpruned_variants {
            let chain = *chains_to_add.pop().unwrap().path.clone();

            if !good_chains.insert(chain) {
                continue
            }

            // When we add an outgoing chain starting from a vertex with other good outgoing chains, we add a variant
            let new_variant = !vertices_that_already_have_outgoing_good_chains.insert(chain.get_first_vertex());
            if new_variant {
                variant_count += 1;
            }
            // check whether we've already added this chain from the other side or we've exceeded the variant count limit
            if new_variant && variant_count > self.max_unpruned_variants {
                continue
            }

            for vertex in vec![chain.get_first_vertex(), chain.get_last_vertex()] {
                if !processed_vertices.contains(vertex) {
                    vertex_to_good_outgoing_chains.get(vertex).iter().for_each(|chain| {
                        chains_to_add.push(Chain::new(OrderedFloat::new(chain_log_odds.get(chain).unwrap().0), chain, graph))
                    });
                    vertex_to_good_incoming_chains.get(vertex).iter().for_each(|chain| {
                        chains_to_add.push(Chain::new(OrderedFloat::new(chain_log_odds.get(chain).unwrap().1), chain, graph))
                    });
                    processed_vertices.insert(*vertex);
                }
            }
        }

        return chains.into_par_iter().filter(|c| {
            !good_chains.contains(c)
        }).collect::<HashSet<Path>>()
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

    // find the chain containing the edge of greatest weight, taking care to break ties deterministically
    fn get_max_weight_chains(chains: &Vec<Path>, graph: &BaseGraph) -> Path {

        chains.par_iter().max_by(
            as_fn(comparing(|c1, c2| {
                c1.get_edges().par_iter().map(|edge| {
                    edge.weight().get_multiplicity()
                }).max().cmp(c2.get_edges().par_iter().map(|edge| {
                    edge.weight().get_multiplicity()
                }).max())
            }).then_compare_by_key(|c1, c2| c1.len().cmp(c2.len()))
                .then_compare_by_key(|c1, c2| graph.compare_paths(c1, c2)))
        );
    }
}