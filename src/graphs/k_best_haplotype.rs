use graphs::path::Path;
use graphs::base_vertex::BaseVertex;
use graphs::base_edge::BaseEdge;
use petgraph::prelude::{NodeIndex, EdgeIndex};
use graphs::base_graph::BaseGraph;
use haplotype::haplotype::Haplotype;
use utils::simple_interval::Locatable;
use std::cmp::Ordering;
use ordered_float::OrderedFloat;
use utils::base_utils::BaseUtils;

/**
 * Represents a result from a K-best haplotype search.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
#[derive(Debug, Clone)]
pub struct KBestHaplotype<'a, V: BaseVertex, E: BaseEdge> {
    pub(crate) score: f64,
    pub(crate) is_reference: bool,
    pub(crate) path: Path<'a, V, E>
}

impl<'a, V: BaseVertex, E: BaseEdge> KBestHaplotype<'a, V, E> {
    pub fn new(initial_vertex: NodeIndex, graph: &'a BaseGraph<V, E>) -> KBestHaplotype<'a, V, E> {
        KBestHaplotype {
            score: 0.0,
            is_reference: true,
            path: Path::new(initial_vertex, Vec::new(), graph)
        }
    }

    pub fn new_from_edge(&self, edge: EdgeIndex, total_outgoing_multiplicity: usize) -> KBestHaplotype<'a, V, E> {
        let edge_weight = self.path.graph.graph.edge_weight(edge).unwrap();

        let path = self.path.new_add_edge(edge);
        let score = self.score + Self::compute_log_penalty_score(
            edge_weight.get_multiplicity(),
            total_outgoing_multiplicity
        );
        let is_reference = self.is_reference & edge_weight.is_ref();

        KBestHaplotype {
            path,
            score,
            is_reference
        }
    }

    pub fn new_from_edges(
        &self,
        edges_to_extend: Vec<EdgeIndex>,
        edge_penalty: f64
    ) -> KBestHaplotype<'a, V, E> {
        let edge_weight = self.path.graph.graph.edge_weight(edges_to_extend.last().unwrap()).unwrap();
        let score = self.score + edge_penalty;
        let is_reference = self.is_reference & edge_weight.is_ref();
        let path = self.path.new_add_edges(edges_to_extend);

        KBestHaplotype {
            score,
            is_reference,
            path,
        }
    }

    pub fn compute_log_penalty_score(
        edge_multiplicity: usize,
        total_outgoing_multiplicity: usize
    ) -> f64 {
        return (edge_multiplicity as f64).log10() - (total_outgoing_multiplicity as f64).log10()
    }

    pub fn haplotype<L: Locatable>(
        &self
    ) -> Haplotype<'a, L> {
        let mut haplotype = Haplotype::new(self.path.get_bases(), self.is_reference);
        haplotype.score = self.score;
        return haplotype
    }
}

impl<'a, V: BaseVertex, E: BaseEdge> Ord for KBestHaplotype<'a, V, E> {
    fn cmp(&self, other: &Self) -> Ordering {
        let result = OrderedFloat::from(self.score).cmp(OrderedFloat::from(other.score)).reverse();
        if result == Ordering::Equal {
            return BaseUtils::bases_comparator(self.path.get_bases(), other.path.get_bases()).reverse()
        } else {
            return result
        }
    }
}

impl<'a, V: BaseVertex, E: BaseEdge> PartialOrd for KBestHaplotype<'a, V, E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, V: BaseVertex, E: BaseEdge> PartialEq for KBestHaplotype<'a, V, E> {
    fn eq(&self, other: &Self) -> bool {
        self == other
    }
}

impl<'a, V: BaseVertex, E: BaseEdge> Eq for KBestHaplotype<'a, V, E> {}
