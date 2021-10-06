use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::k_best_haplotype::KBestHaplotype;
use graphs::k_best_haplotype_finder::KBestHaplotypeFinder;
use petgraph::prelude::NodeIndex;
use petgraph::Direction;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link BaseGraph instance}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * For Rust translation:
 * @author Rhys Newell &lt;rhys.newell@hdr.qut.edu.au&gt;
 */
pub struct GraphBasedKBestHaplotypeFinder {
    pub k_best_haplotype_finder: KBestHaplotypeFinder,
}

impl GraphBasedKBestHaplotypeFinder {
    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the graph to search.
     * @param sources source vertices for all haplotypes.
     * @param sinks sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code sources} or {@code sinks} is {@code null} or</li>
     *     <li>any of {@code sources}' or any {@code sinks}' member is not a vertex in {@code graph}.</li>
     * </ul>
     */
    pub fn new<V: BaseVertex, E: BaseEdge>(
        graph: &mut BaseGraph<V, E>,
        sources: HashSet<NodeIndex>,
        sinks: HashSet<NodeIndex>,
    ) -> GraphBasedKBestHaplotypeFinder {
        GraphBasedKBestHaplotypeFinder {
            k_best_haplotype_finder: KBestHaplotypeFinder::new(sinks, sources, graph),
        }
    }

    pub fn new_from_singletons<V: BaseVertex, E: BaseEdge>(
        graph: &mut BaseGraph<V, E>,
        source: NodeIndex,
        sink: NodeIndex,
    ) -> GraphBasedKBestHaplotypeFinder {
        let mut sources = HashSet::new();
        sources.insert(source);

        let mut sinks = HashSet::new();
        sinks.insert(sink);

        Self::new(graph, sources, sinks)
    }

    /**
     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
     */
    pub fn find_best_haplotypes<V: BaseVertex, E: BaseEdge>(
        &self,
        max_number_of_haplotypes: usize,
        graph: &BaseGraph<V, E>,
    ) -> Vec<KBestHaplotype> {
        let mut result = Vec::new();
        debug!(
            "Sources {:?} Sinks {:?}",
            &self.k_best_haplotype_finder.sources, &self.k_best_haplotype_finder.sinks
        );

        let mut queue: BinaryHeap<KBestHaplotype> = self
            .k_best_haplotype_finder
            .sources
            .par_iter()
            .map(|source| KBestHaplotype::new(*source, graph))
            .collect::<BinaryHeap<KBestHaplotype>>();

        debug!("Graph {:?}", &graph);
        let mut vertex_counts = graph
            .graph
            .node_indices()
            .par_bridge()
            .map(|v| (v, 0))
            .collect::<HashMap<NodeIndex, usize>>();

        debug!(
            "Sinks {:?} queue {:?} vertex counts {:?}",
            &self.k_best_haplotype_finder.sinks, &queue, &vertex_counts
        );

        while !queue.is_empty() && result.len() < max_number_of_haplotypes {
            let mut path_to_extend = queue.pop().unwrap();
            debug!("Path to extend {:?}", &path_to_extend);
            let vertex_to_extend = path_to_extend.path.last_vertex;
            if self
                .k_best_haplotype_finder
                .sinks
                .contains(&vertex_to_extend)
            {
                result.push(path_to_extend);
            } else {
                if vertex_counts.contains_key(&vertex_to_extend) {
                    let vertex_count = vertex_counts.entry(vertex_to_extend).or_insert(0);
                    *vertex_count += 1;

                    if *vertex_count < max_number_of_haplotypes {
                        let outgoing_edges =
                            graph.edges_directed(vertex_to_extend, Direction::Outgoing);
                        debug!("Outgoing edges {:?}", &outgoing_edges);
                        let mut total_outgoing_multiplicity = 0;
                        for edge in outgoing_edges.iter() {
                            total_outgoing_multiplicity +=
                                graph.graph.edge_weight(*edge).unwrap().get_multiplicity();
                        }

                        for edge in outgoing_edges.iter() {
                            queue.push(path_to_extend.new_from_edge(
                                *edge,
                                total_outgoing_multiplicity,
                                graph,
                            ));
                        }
                    }
                }
            }
        }

        debug!("Results {:?}", &result);

        return result;
    }
}
