use graphs::seq_graph::SeqGraph;
use graphs::base_edge::BaseEdgeStruct;
use petgraph::stable_graph::NodeIndex;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V
 *
 * This algorithm creates a new SeqGraph with the following configuration
 *
 * prefix -> has outgoing edges to all seq_i
 * suffix -> has incoming edges for all seq_i
 *
 * There are a few special cases that must be handled.  First, Vi could be simply
 * == to the prefix or the suffix.  These generate direct connections between
 * the prefix and suffix nodes, and they are handled internally by the algorithm.
 *
 * Note that for convenience, we will always create newTop and newBottom nodes, but
 * these may be empty node (i.e., they contain no sequence).  That allows them to be
 * trivially merged, if desired, when the graph is incorporated into an overall
 * graph.
 *
 * The product of this operation is a SeqGraph that contains the split.  There's a
 * function to merge reconnect this graph into the graph that contains the middle nodes
 *
 * The process guarentees a few things about the output:
 *
 * -- Preserves the paths and weights among all vertices
 *
 * It produces a graph that has some unusual properties
 *
 * -- May add nodes with no sequence (isEmpty() == true) to preserve connectivity among the graph
 * -- May introduce edges with no multiplicity to preserve paths through the graph
 *
 * The overall workflow of using this class is simple:
 *
 * find vertices V in graph that you want to split out
 * s = new SharedVertexSequenceSplitter(graph, V)
 * s.updateGraph(graph)
 *
 * to update the graph with the modifications created by this splitter
 */
pub struct SharedVertexSequenceSplitter<'a> {
    outer: &'a mut SeqGraph<BaseEdgeStruct>,
    prefix_v: NodeIndex,
    suffix_v: NodeIndex,
    to_splits: HashSet<NodeIndex>,
    split_graph: Option<SeqGraph<BaseEdgeStruct>>,
    new_middles: Option<HashSet<NodeIndex>>,

}

impl<'a> SharedVertexSequenceSplitter<'a> {

}