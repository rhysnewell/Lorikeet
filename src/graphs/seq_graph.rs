use graphs::base_graph::BaseGraph;
use graphs::base_edge::{BaseEdgeStruct, BaseEdge};
use graphs::seq_vertex::SeqVertex;
use petgraph::graph::NodeIndex;
use itertools::zip;
use petgraph::Direction;

/**
 * A graph that contains base sequence at each node
 */
#[derive(Debug, Clone)]
pub struct SeqGraph {
    base_graph: BaseGraph<SeqVertex, BaseEdgeStruct>,
}

impl SeqGraph {
    const PRINT_SIMPLIFY_GRAPHS: bool = false;
    /**
     * How many cycles of the graph simplifications algorithms will we run before
     * thinking something has gone wrong and throw an exception?
     */
    const MAX_REASONABLE_SIMPLIFICATION_CYCLES: usize = 100;

    /**
     * Construct an empty SeqGraph where we'll add nodes based on a kmer size of kmer
     *
     * The kmer size is purely information.  It is useful when converting a Debruijn graph -> SeqGraph
     * for us to track the kmer used to make the transformation.
     *
     * @param kmer kmer
     */
    pub fn new(kmer: usize) -> SeqGraph {
        SeqGraph {
            base_graph: BaseGraph::new(kmer)
        }
    }

    /**
     * Simplify this graph, merging vertices together and restructuring the graph in an
     * effort to minimize the number of overall vertices in the graph without changing
     * in any way the sequences implied by a complex enumeration of all paths through the graph.
     */
    pub fn simplify_graph(&mut self) {
        self.simplify_graph_with_cycles(std::usize::MAX)
    }

    fn simplify_graph_with_cycles(&mut self, max_cycles: usize) {
        // start off with one round of zipping of chains for performance reasons
        self.zip_linear_chains();

        let prev_graph = None;
        for i in 0..max_cycles {
            if i > Self::MAX_REASONABLE_SIMPLIFICATION_CYCLES {
                warn!("Infinite loop detected in simpliciation routines.  Writing current graph to debugMeRhys.dot");
                self.
            }
        }
    }

    /**
     * Zip up all of the simple linear chains present in this graph.
     *
     * Merges together all pairs of vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
     *
     * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
     *
     * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
     *
     * @return true if any such pair of vertices could be found, false otherwise
     */
    pub fn zip_linear_chains(&self) -> bool {
        let mut zip_starts = Vec::new();
        for source in self.base_graph.graph.node_indices() {
            if self.is_linear_chain_start(source) {
                zip_starts.push(source)
            }
        }

        if zip_starts.is_empty() {
            return false
        }

        // At this point, zipStarts contains all of the vertices in this graph that might start some linear
        // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
        // zipping them up with mergeLinearChain, if possible
        let mut merged_one = false;
        for zip_start in zip_starts.into_iter() {
            let linear_chain = self.trace_linear_chain(zip_start);

            // merge the linearized chain, recording if we actually did some useful work
            merged_one = merged_one | self.merge_linear_chain(&linear_chain);
        }

        return merged_one
    }

    /**
     * Is source vertex potentially a start of a linear chain of vertices?
     *
     * We are a start of a zip chain if our out degree is 1 and either the
     * the vertex has no incoming connections or 2 or more (we must start a chain) or
     * we have exactly one incoming vertex and that one has out-degree > 1 (i.e., source's incoming
     * vertex couldn't be a start itself
     *
     * @param source a non-null vertex
     * @return true if source might start a linear chain
     */
    fn is_linear_chain_start(&self, source: NodeIndex) -> bool {
        return self.base_graph.out_degree_of(source) == 1 &&
            (self.base_graph.in_degree_of(source) != 1 ||
                self.base_graph.out_degree_of(self.base_graph.incoming_vertices_of(source).iter().next()) > 1)
    }

    /**
     * Get all of the vertices in a linear chain of vertices starting at zipStart
     *
     * Build a list of vertices (in order) starting from zipStart such that each sequential pair of vertices
     * in the chain A and B can be zipped together.
     *
     * @param zipStart a vertex that starts a linear chain
     * @return a list of vertices that comprise a linear chain starting with zipStart.  The resulting
     *         list will always contain at least zipStart as the first element.
     */
    fn trace_linear_chain(&self, zip_start: NodeIndex) -> Vec<NodeIndex> {
        let mut linear_chain = Vec::new();
        linear_chain.push(zip_start);

        let last_is_ref = self.base_graph.is_reference_node(zip_start); // remember because this calculation is expensive
        let mut last = zip_start;

        loop {
            if self.base_graph.out_degree_of(last) != 1 {
                // cannot extend a chain from last if last has multiple outgoing branches
                break
            }

            // there can only be one (outgoing edge of last) by contract
            let target = self.base_graph.get_edge_target(
                self.base_graph.graph.first_edge(last, Direction::Outgoing).unwrap_or_else(break)
            );

            if self.base_graph.in_degree_of(target) != 1 || last == target {
                // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
                break
            }

            let target_is_ref = self.base_graph.is_reference_node(target);
            if last_is_ref != target_is_ref { // both our isRef states must be equal
                break
            }

            linear_chain.push(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            last_is_ref = target_is_ref;
        }

        return linear_chain
    }

    /**
     * Merge a linear chain of vertices into a single combined vertex, and update this graph to such that
     * the incoming edges into the first element of the linearChain and the outgoing edges from linearChain.getLast()
     * all point to this new combined vertex.
     *
     * @param linearChain a non-empty chain of vertices that can be zipped up into a single vertex
     * @return true if we actually merged at least two vertices together
     */
    fn merge_linear_chain(&self, linear_chain: &Vec<NodeIndex>) -> bool {
        match self.merge_linear_chain_vertex(linear_chain) {
            Some(index) => true,
            None => false
        }
    }

    fn merge_linear_chain_vertex(&self, linear_chain: &Vec<NodeIndex>) -> Option<NodeIndex> {
        assert!(!linear_chain.is_empty(), "Cannot have linear chain with 0 elements but got {:?}", linear_chain);

        let first = linear_chain.first().unwrap();
        let last = linear_chain.last().unwrap();

        if first == last {
            return None // only one element in the chain, cannot be extended
        }

        // create the combined vertex, and add it to the graph
        // TODO -- performance problem -- can be optimized if we want
        let added_vertex = self.merge_linear_chain_vertices(linear_chain);

        let added_node_index = self.base_graph.graph.add_node(added_vertex);

        // update the incoming and outgoing edges to point to the new vertex
        for edge in self.base_graph.graph.edges_directed(last, Direction::Outgoing) {
            self.base_graph.graph.add_edge(added_node_index, self.base_graph.get_edge_target(edge.id()), edge.weight().clone())
        }
        for edge in self.base_graph.graph.edges_directed(last, Direction::Incoming) {
            self.base_graph.graph.add_edge(self.base_graph.get_edge_source(edge.id()), added_node_index, edge.weight().clone())
        }

        self.base_graph.graph.retain_nodes(|v| !linear_chain.contains(&v));

        return Some(added_node_index)
    }

    fn merge_linear_chain_vertices(&self, vertices: &Vec<NodeIndex>) -> SeqVertex {
        let seqs = vertices.par_iter().map(|v| {
            self.base_graph.graph[v].unwrap().get_sequence()
        }).collect::<Vec<&[u8]>>();

        let seqs_joined = seqs.concat();

        return SeqVertex::new(seqs_joined);
    }
}