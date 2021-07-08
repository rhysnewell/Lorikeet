use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;

pub enum Status {
    Failed,
    JustAssembledReference,
    AssembledSomeVariation,
}

/**
 * Result of assembling, with the resulting graph and status
 */
pub struct AssemblyResult<'a> {
    status: Status,
    threading_graph: AbstractReadThreadingGraph<'a>,
    // graph:
}