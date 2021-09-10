use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::path::Path;
use haplotype::haplotype::Haplotype;
use ordered_float::OrderedFloat;
use petgraph::prelude::{EdgeIndex, NodeIndex};
use std::cmp::Ordering;
use utils::simple_interval::Locatable;

/**
 * Represents a result from a K-best haplotype search.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
#[derive(Debug, Clone)]
pub struct KBestHaplotype {
    pub score: f64,
    pub is_reference: bool,
    pub path: Path,
    pub kmer_size: usize,
}

impl KBestHaplotype {
    pub fn new<V: BaseVertex, E: BaseEdge>(
        initial_vertex: NodeIndex,
        graph: &BaseGraph<V, E>,
    ) -> KBestHaplotype {
        KBestHaplotype {
            score: 0.0,
            is_reference: true,
            path: Path::new(initial_vertex, Vec::new()),
            kmer_size: graph.get_kmer_size(),
        }
    }

    pub fn new_from_edge<V: BaseVertex, E: BaseEdge>(
        &self,
        edge: EdgeIndex,
        total_outgoing_multiplicity: usize,
        graph: &BaseGraph<V, E>,
    ) -> KBestHaplotype {
        let edge_weight = graph.graph.edge_weight(edge).unwrap();

        let path = self.path.new_add_edge(edge, graph);
        let score = self.score
            + Self::compute_log_penalty_score(
                edge_weight.get_multiplicity(),
                total_outgoing_multiplicity,
            );
        let is_reference = self.is_reference & edge_weight.is_ref();

        KBestHaplotype {
            path,
            score,
            is_reference,
            kmer_size: graph.get_kmer_size(),
        }
    }

    pub fn new_from_edges<V: BaseVertex, E: BaseEdge>(
        &self,
        edges_to_extend: Vec<EdgeIndex>,
        edge_penalty: f64,
        graph: &BaseGraph<V, E>,
    ) -> KBestHaplotype {
        let edge_weight = graph
            .graph
            .edge_weight(*edges_to_extend.last().unwrap())
            .unwrap();
        let score = self.score + edge_penalty;
        let is_reference = self.is_reference & edge_weight.is_ref();
        let path = self.path.new_add_edges(edges_to_extend, graph);

        KBestHaplotype {
            score,
            is_reference,
            path,
            kmer_size: graph.get_kmer_size(),
        }
    }

    pub fn compute_log_penalty_score(
        edge_multiplicity: usize,
        total_outgoing_multiplicity: usize,
    ) -> f64 {
        return (edge_multiplicity as f64).log10() - (total_outgoing_multiplicity as f64).log10();
    }

    pub fn haplotype<L: Locatable, V: BaseVertex, E: BaseEdge>(
        &self,
        graph: &BaseGraph<V, E>,
    ) -> Haplotype<'static, L> {
        let mut haplotype = Haplotype::new(&self.path.get_bases(graph), self.is_reference);
        haplotype.score = OrderedFloat(self.score);
        return haplotype;
    }
}

impl Ord for KBestHaplotype {
    fn cmp(&self, other: &Self) -> Ordering {
        let result = OrderedFloat::from(self.score)
            .cmp(&OrderedFloat::from(other.score))
            .reverse();
        if result == Ordering::Equal {
            return other
                .path
                .edges_in_order
                .cmp(&self.path.edges_in_order)
                .reverse();
        } else {
            return result;
        }
    }
}

impl PartialOrd for KBestHaplotype {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for KBestHaplotype {
    fn eq(&self, other: &Self) -> bool {
        self == other
    }
}

impl Eq for KBestHaplotype {}
