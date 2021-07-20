use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use graphs::k_best_haplotype_finder::KBestHaplotypeFinder;
use graphs::base_graph::BaseGraph;
use petgraph::prelude::NodeIndex;
use graphs::k_best_haplotype::KBestHaplotype;
use std::collections::{BinaryHeap, HashMap, HashSet};
use petgraph::Direction;

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link BaseGraph instance}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * For Rust translation:
 * @author Rhys Newell &lt;rhys.newell@hdr.qut.edu.au&gt;
 */
pub struct GraphBasedKBestHaplotypeFinder<'a, V: BaseVertex, E: BaseEdge> {
    k_best_haplotype_finder: KBestHaplotypeFinder<'a, V, E>
}

impl<'a, V: BaseVertex, E: BaseEdge> GraphBasedKBestHaplotypeFinder<'a, V, E> {

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
    pub fn new(
        graph: &'a BaseGraph<V, E>,
        sources: HashSet<NodeIndex>,
        sinks: HashSet<NodeIndex>
    ) -> GraphBasedKBestHaplotypeFinder<'a, V, E> {
        GraphBasedKBestHaplotypeFinder {
            k_best_haplotype_finder: KBestHaplotypeFinder::new(sinks, sources, graph)
        }
    }

    pub fn new_from_singletons(
        graph: &'a BaseGraph<V, E>,
        source: NodeIndex,
        sink: NodeIndex
    ) -> GraphBasedKBestHaplotypeFinder<'a, V, E> {
        let mut sources = HashSet::new();
        sources.insert(source);

        let mut sinks = HashSet::new();
        sinks.insert(sink);

        Self::new(graph, sources, sinks)
    }

    /**
     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
     */
    pub fn find_best_haplotypes(&self, max_number_of_haplotypes: usize) -> Vec<KBestHaplotype<'a, V, E>> {
        let mut result = Vec::new();

        let mut queue: BinaryHeap<KBestHaplotype<V, E>> = self.k_best_haplotype_finder.sources.par_iter().map(|source| {
            KBestHaplotype::new(source, self.k_best_haplotype_finder.graph)
        }).collect::<BinaryHeap<KBestHaplotype<V, E>>>();

        let mut vertex_counts = self.k_best_haplotype_finder.graph.graph.node_indices()
            .par_bridge().map(|v| (v, 0)).collect::<HashMap<NodeIndex, usize>>();

        while !queue.is_empty() && result.len() < max_number_of_haplotypes {
            let mut path_to_extend = queue.pop().unwrap();
            let vertex_to_extend = path_to_extend.path.last_vertex;
            if self.k_best_haplotype_finder.sinks.contains(&vertex_to_extend) {
                result.push(path_to_extend);
            } else {
                if vertex_counts.contains_key(&vertex_to_extend) {
                    let vertex_count = vertex_counts.entry(&vertex_to_extend).or_insert(0);
                    vertex_count += 1;
                    if vertex_count < max_number_of_haplotypes {
                        let outgoing_edges = self.k_best_haplotype_finder.graph.graph
                            .edges_directed(vertex_to_extend, Direction::Outgoing);
                        let total_outgoing_multiplicity = 0;
                        for edge in outgoing_edges {
                            total_outgoing_multiplicity +=
                                self.k_best_haplotype_finder.graph.graph
                                    .edge_weight(edge.id()).unwrap().get_multiplicity();
                        };

                        for edge in outgoing_edges {
                            queue.push(path_to_extend.new_from_edge(edge.id(), total_outgoing_multiplicity));
                        };
                    }
                }
            }
        }

        return result
    }
}