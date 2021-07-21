use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use graphs::seq_graph::SeqGraph;
use graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use std::collections::HashSet;
use haplotype::haplotype::Haplotype;
use utils::simple_interval::Locatable;
use linked_hash_set::LinkedHashSet;

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum Status {
    Failed,
    JustAssembledReference,
    AssembledSomeVariation,
}

/**
 * Result of assembling, with the resulting graph and status
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AssemblyResult<'a, L: Locatable, A: AbstractReadThreadingGraph<'a>> {
    pub(crate) status: Status,
    pub(crate) threading_graph: Option<A>,
    pub(crate) graph: Option<SeqGraph<'a, BaseEdgeStruct>>,
    pub(crate) discovered_haplotypes: LinkedHashSet<Haplotype<'a, L>>,
    pub(crate) contains_suspect_haploptypes: bool,
}

impl<'a, L: Locatable, A: AbstractReadThreadingGraph<'a>> AssemblyResult<'a, L, A> {
    /**
     * Create a new assembly result
     * @param status the status, cannot be null
     * @param graph the resulting graph of the assembly, can only be null if result is failed
     */
    pub fn new(
        status: Status,
        graph: Option<SeqGraph<'a, BaseEdgeStruct>>,
        threading_graph: Option<A>
    ) -> AssemblyResult<'a, L, A> {
        AssemblyResult {
            status,
            graph,
            threading_graph,
            discovered_haplotypes: LinkedHashSet::new(),
            contains_suspect_haploptypes: false,
        }
    }

    pub fn set_discovered_haplotypes(&mut self, discovered_haplotypes: LinkedHashSet<Haplotype<'a, L>>) {
        self.discovered_haplotypes = discovered_haplotypes
    }

    pub fn set_contains_suspect_haplotypes(&mut self, contains_suspect_haplotypes: bool) {
        self.contains_suspect_haploptypes = contains_suspect_haplotypes
    }

    pub fn get_kmer_size(&self) -> usize {
        match self.graph {
            None => {
                match self.threading_graph {
                    None => panic!("No established kmer size for assembly result"),
                    Some(threading_graph) => threading_graph.get_kmer_size()
                }
            },
            Some(graph) => {
                graph.base_graph.get_kmer_size()
            }
        }
    }
}