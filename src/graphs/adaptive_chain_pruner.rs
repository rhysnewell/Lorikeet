use compare::{Compare, Extract};
use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::chain_pruner::ChainPruner;
use graphs::path::{Chain, Path};
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use hashlink::LinkedHashSet;
use multimap::MultiMap;
use ordered_float::OrderedFloat;
use petgraph::prelude::{EdgeIndex, EdgeRef};
use petgraph::stable_graph::NodeIndex;
use petgraph::Direction;
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};
use utils::base_utils::BaseUtils;

#[derive(Debug, Clone)]
pub struct AdaptiveChainPruner {
    pub initial_error_probability: f64,
    log_odds_threshold: f64,
    seeding_log_odds_threshold: f64, // threshold for seeding subgraph of good vertices
    max_unpruned_variants: usize,
}

// TODO: This whole thing needs some unit tests, I don't trust parts of this code as there might be
//       some issues with refrence/ownership and also chained comparators.
impl AdaptiveChainPruner {
    pub fn new(
        initial_error_probability: f64,
        log_odds_threshold: f64,
        seeding_log_odds_threshold: f64,
        max_unpruned_variants: usize,
    ) -> AdaptiveChainPruner {
        AdaptiveChainPruner {
            initial_error_probability,
            log_odds_threshold,
            seeding_log_odds_threshold,
            max_unpruned_variants,
        }
    }

    pub fn likely_error_chains<
        'a,
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync,
    >(
        &self,
        chains: &'a VecDeque<Path>,
        graph: &BaseGraph<V, E>,
        error_rate: f64,
    ) -> HashSet<&'a Path> {
        // pre-compute the left and right log odds of each chain
        let chain_log_odds = chains
            .iter()
            .map(|chain| {
                let result = self.chain_log_odds(chain, graph, error_rate);
                (chain, result)
            })
            .collect::<HashMap<&Path, (f64, f64)>>();

        // compute correspondence of vertices to incident chains with log odds above the seeding and extending thresholds
        let mut vertex_to_seedable_chains = MultiMap::new();
        let mut vertex_to_good_incoming_chains = MultiMap::new();
        let mut vertex_to_good_outgoing_chains = MultiMap::new();
        let mut right = 0;
        let mut left = 0;
        let mut seed = 0;

        for (index, chain) in chains.iter().enumerate() {
            if chain_log_odds.get(chain).unwrap().1 >= self.log_odds_threshold
                || graph.edge_is_ref(chain.edges_in_order[0])
            {
                vertex_to_good_incoming_chains.insert(chain.get_last_vertex(), chain);
                right += 1;
            }

            if chain_log_odds.get(chain).unwrap().0 >= self.log_odds_threshold
                || graph.edge_is_ref(chain.edges_in_order[0])
            {
                vertex_to_good_outgoing_chains.insert(chain.get_first_vertex(graph), chain);
                left += 1;
            }

            // seed-worthy chains must pass the more stringent seeding log odds threshold on both sides
            // in addition to that, we only seed from vertices with multiple such chains incoming or outgoing (see below)
            if chain_log_odds.get(chain).unwrap().1 >= self.seeding_log_odds_threshold
                && chain_log_odds.get(chain).unwrap().0 >= self.seeding_log_odds_threshold
            {
                vertex_to_seedable_chains.insert(chain.get_first_vertex(graph), chain);
                vertex_to_seedable_chains.insert(chain.get_last_vertex(), chain);
                seed += 2;
            }
        }

        // We have a priority queue of chains to add to the graph, with priority given by the log odds (higher first)
        // for determinism we have a tie breaker based on chains' first vertex
        // Note that chains can be added twice to the queue, once for each side
        // The comparator for this heap is implemented in (Chain)[Chain]
        let mut chains_to_add = BinaryHeap::new();

        // seed the subgraph of good chains by starting with some definitely-good chains.  These include the max-weight chain
        // and chains emanating from vertices with two incoming or two outgoing chains (plus one outgoing or incoming for a total of 3 or more) with good log odds
        // The idea is that a high-multiplicity error chain A that branches into a second error chain B and a continuation-of-the-original-error chain A'
        // may have a high log odds for A'.  However, only in the case of true variation will multiple branches leaving the same vertex have good log odds.
        let max_weight_chain = Self::get_max_weight_chain(chains, graph);
        chains_to_add.push(Reverse(Chain::new(
            OrderedFloat::from(std::f64::INFINITY),
            &max_weight_chain,
            graph,
        )));
        let mut processed_vertices = HashSet::new();
        for (vertex, paths) in vertex_to_seedable_chains.iter_all() {
            if paths.len() > 2 {
                vertex_to_good_outgoing_chains
                    .get_vec(vertex)
                    .unwrap_or(&vec![])
                    .into_iter()
                    .for_each(|chain| {
                        chains_to_add.push(Reverse(Chain::new(
                            OrderedFloat::from(chain_log_odds.get(chain).unwrap().0),
                            chain,
                            graph,
                        )))
                    });
                vertex_to_good_incoming_chains
                    .get_vec(vertex)
                    .unwrap_or(&vec![])
                    .into_iter()
                    .for_each(|chain| {
                        chains_to_add.push(Reverse(Chain::new(
                            OrderedFloat::from(chain_log_odds.get(chain).unwrap().1),
                            chain,
                            graph,
                        )))
                    });
                processed_vertices.insert(*vertex);
            }
        }

        let mut good_chains = HashSet::new();
        let mut vertices_that_already_have_outgoing_good_chains = HashSet::new();
        let mut variant_count = 0;
        // starting from the high-confidence seed vertices, grow the "good" subgraph along chains with above-threshold log odds,
        // discovering good chains as we go.
        let mut checked = 0;
        while !chains_to_add.is_empty() && variant_count <= self.max_unpruned_variants {
            let chain = chains_to_add.pop().unwrap().0;

            if !good_chains.insert(chain.path) {
                continue;
            }

            // When we add an outgoing chain starting from a vertex with other good outgoing chains, we add a variant
            let new_variant = !vertices_that_already_have_outgoing_good_chains
                .insert(chain.path.get_first_vertex(graph));
            if new_variant {
                variant_count += 1;
            }
            // check whether we've already added this chain from the other side or we've exceeded the variant count limit
            if new_variant && variant_count > self.max_unpruned_variants {
                continue;
            }

            for vertex in &[
                chain.path.get_first_vertex(graph),
                chain.path.get_last_vertex(),
            ] {
                if !processed_vertices.contains(vertex) {
                    vertex_to_good_outgoing_chains
                        .get_vec(vertex)
                        .unwrap_or(&vec![])
                        .into_iter()
                        .for_each(|chain| {
                            chains_to_add.push(Reverse(Chain::new(
                                OrderedFloat::from(chain_log_odds.get(chain).unwrap().0),
                                chain,
                                graph,
                            )))
                        });
                    vertex_to_good_incoming_chains
                        .get_vec(vertex)
                        .unwrap_or(&vec![])
                        .into_iter()
                        .for_each(|chain| {
                            chains_to_add.push(Reverse(Chain::new(
                                OrderedFloat::from(chain_log_odds.get(chain).unwrap().1),
                                chain,
                                graph,
                            )))
                        });
                    processed_vertices.insert(*vertex);
                }
            }
        }

        return chains
            .iter()
            .filter(|c| !good_chains.contains(c))
            .collect::<HashSet<&Path>>();
    }

    // left and right chain log odds
    fn chain_log_odds<V: BaseVertex + std::marker::Sync, E: BaseEdge + std::marker::Sync>(
        &self,
        chain: &Path,
        graph: &BaseGraph<V, E>,
        error_rate: f64,
    ) -> (f64, f64) {
        let left_total_multiplicity = graph
            .graph
            .edges_directed(chain.get_first_vertex(graph), Direction::Outgoing)
            .map(|e| e.weight().get_multiplicity())
            .sum::<usize>();
        let right_total_multiplicity = graph
            .graph
            .edges_directed(chain.get_last_vertex(), Direction::Incoming)
            .map(|e| e.weight().get_multiplicity())
            .sum::<usize>();

        let left_multiplicity = graph
            .graph
            .edge_weight(chain.get_edges()[0])
            .unwrap()
            .get_multiplicity();
        let right_multiplicity = graph
            .graph
            .edge_weight(chain.get_last_edge())
            .unwrap()
            .get_multiplicity();

        let left_log_odds = if graph.is_source(chain.get_first_vertex(graph)) {
            0.0
        } else {
            HaplotypeCallerEngine::log_likelihood_ratio_constant_error(
                left_total_multiplicity - left_multiplicity,
                left_multiplicity,
                error_rate,
            )
        };
        let right_log_odds = if graph.is_sink(chain.get_last_vertex()) {
            0.0
        } else {
            HaplotypeCallerEngine::log_likelihood_ratio_constant_error(
                right_total_multiplicity - right_multiplicity,
                right_multiplicity,
                error_rate,
            )
        };

        (left_log_odds, right_log_odds)
    }

    // find the chain containing the edge of greatest weight, taking care to break ties deterministically
    fn get_max_weight_chain<
        'a,
        V: BaseVertex + std::marker::Sync,
        E: BaseEdge + std::marker::Sync,
    >(
        chains: &'a VecDeque<Path>,
        graph: &BaseGraph<V, E>,
    ) -> &'a Path {
        // let comparator = Extract::new(|chain: &Path| {
        //     chain
        //         .get_edges()
        //         .par_iter()
        //         .map(|edge| graph.graph.edge_weight(*edge).unwrap().get_multiplicity())
        //         .max().unwrap_or(0)
        // })
        // .then(Extract::new(|l: &Path| l.len()))
        // .then(Extract::new(|l: &Path| l.get_bases(graph)));

        chains
            .iter()
            // .max_by(|l, r| comparator.compare(l, r))
            .max_by(|l, r| {
                l.get_max_multiplicity(graph)
                    .cmp(&r.get_max_multiplicity(graph))
                    .then_with(|| l.len().cmp(&r.len()))
                    .then_with(|| {
                        BaseUtils::bases_comparator(&l.get_bases(graph), &r.get_bases(graph))
                    })
            })
            .unwrap()
    }
}

impl<V: BaseVertex + std::marker::Sync, E: BaseEdge + std::marker::Sync> ChainPruner<V, E>
    for AdaptiveChainPruner
{
    fn prune_low_weight_chains(&self, graph: &mut BaseGraph<V, E>) {
        let chains = Self::find_all_chains(&graph);

        let chains_to_remove = self.chains_to_remove(&chains, &graph);
        chains_to_remove
            .into_iter()
            .for_each(|chain| graph.remove_all_edges(chain.get_edges()));

        graph.remove_singleton_orphan_vertices();
    }

    fn find_all_chains(graph: &BaseGraph<V, E>) -> VecDeque<Path> {
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
    fn find_chain(start_edge: EdgeIndex, graph: &BaseGraph<V, E>) -> Path {
        let mut edges = vec![start_edge];
        let start_edge_endpoints = graph.graph.edge_endpoints(start_edge).unwrap();
        let first_vertex_id = start_edge_endpoints.0;
        let mut last_vertex_id = start_edge_endpoints.1;

        loop {
            // chain ends if: 1) no out edges; 2) multiple out edges;
            // 3) multiple in edges; 4) cycle back to start of chain
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

    fn chains_to_remove<'a>(
        &self,
        chains: &'a VecDeque<Path>,
        graph: &BaseGraph<V, E>,
    ) -> Vec<&'a Path> {
        if chains.is_empty() {
            return Vec::new();
        }
        let probable_error_chains =
            self.likely_error_chains(&chains, graph, self.initial_error_probability);

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
        let error_rate = error_count as f64 / total_bases as f64;

        self.likely_error_chains(&chains, graph, error_rate)
            .into_iter()
            .filter(|c| {
                !c.get_edges()
                    .iter()
                    .any(|e| graph.graph.edge_weight(*e).unwrap().is_ref())
            })
            .collect::<Vec<&Path>>()
    }
}
