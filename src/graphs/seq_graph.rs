use graphs::base_edge::BaseEdge;
use graphs::base_graph::BaseGraph;
use graphs::base_vertex::BaseVertex;
use graphs::seq_vertex::SeqVertex;
use graphs::vertex_based_transformer::VertexBasedTransformer::{
    MergeCommonSuffices, MergeDiamonds, MergeTails, SplitCommonSuffices,
};
use graphs::vertex_based_transformer::VertexBasedTransformerOptions;
use petgraph::algo::is_cyclic_directed;
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::collections::HashSet;

/**
 * A graph that contains base sequence at each node
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SeqGraph<E: BaseEdge> {
    pub base_graph: BaseGraph<SeqVertex, E>,
}

impl<E: BaseEdge + std::marker::Sync> SeqGraph<E> {
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
    pub fn new(kmer: usize) -> SeqGraph<E> {
        SeqGraph {
            base_graph: BaseGraph::<SeqVertex, E>::new(kmer),
        }
    }

    /**
     * Simplify this graph, merging vertices together and restructuring the graph in an
     * effort to minimize the number of overall vertices in the graph without changing
     * in any way the sequences implied by a complex enumeration of all paths through the graph.
     */
    pub fn simplify_graph(&mut self, name: &str) {
        self.simplify_graph_with_cycles(std::usize::MAX, name)
    }

    pub fn simplify_graph_with_cycles(&mut self, max_cycles: usize, name: &str) {
        // start off with one round of zipping of chains for performance reasons
        self.zip_linear_chains();
        let mut prev_graph = None;
        for i in 0..max_cycles {
            if i > Self::MAX_REASONABLE_SIMPLIFICATION_CYCLES {
                warn!("Infinite loop detected in simpliciation routines.  Writing current graph to debugMeRhys.dot");
                self.base_graph.print_graph("debugMeRhys.dot", true, 0);
                panic!(
                    "Infinite loop detected in simplification routines for kmer graph {}",
                    self.base_graph.get_kmer_size()
                );
            }

            let did_some_work = self.simplify_graph_once(i, name);
            if !did_some_work {
                // no simplification algorithm could be run so break
                break;
            }

            // we get five cycles before we start looking for changes in the graph
            // by cloning ourselves and then checking for any changes
            if i > 5 {
                // the previous graph and this graph have the same structure, so the simplification
                // algorithms are looping endless between states.  Just break and consider ourselves done
                prev_graph = match prev_graph {
                    None => Some(self.base_graph.clone()),
                    Some(prev_graph) => {
                        if prev_graph == self.base_graph {
                            break;
                        } else {
                            Some(self.base_graph.clone())
                        }
                    }
                }
            }
        }
    }

    /**
     * Run one full cycle of the graph simplification algorithms
     * @return true if any algorithms said they did some simplification
     */
    fn simplify_graph_once(&mut self, iteration: usize, name: &str) -> bool {
        // iterate until we haven't don't anything useful
        debug!(
            "Before anything Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );

        let mut did_some_work = false;
        did_some_work |=
            VertexBasedTransformerOptions::MergeDiamonds.transform_until_complete(self);
        debug!(
            "After diamonds Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );
        self.print_graph_simplification(&format!(
            "{}_simplify_graph.{}.1.diamonds.dot",
            name, iteration
        ));
        did_some_work |= VertexBasedTransformerOptions::MergeTails.transform_until_complete(self);
        debug!(
            "After tails Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );
        self.print_graph_simplification(&format!(
            "{}_simplify_graph.{}.2.tails.dot",
            name, iteration
        ));

        did_some_work |=
            VertexBasedTransformerOptions::SplitCommonSuffices.transform_until_complete(self);
        debug!(
            "After split common Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );

        self.print_graph_simplification(&format!(
            "{}_simplify_graph.{}.3.split_suffix.dot",
            name, iteration
        ));

        did_some_work |=
            VertexBasedTransformerOptions::MergeCommonSuffices.transform_until_complete(self);
        debug!(
            "After merge common Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );

        self.print_graph_simplification(&format!(
            "{}_simplify_graph.{}.4.merge_suffix.dot",
            name, iteration
        ));

        did_some_work |= self.zip_linear_chains();
        debug!(
            "After zip Edges {}, Nodes {}, Cyclic? {}",
            self.base_graph.edge_set().len(),
            self.base_graph.vertex_set().len(),
            is_cyclic_directed(&self.base_graph.graph)
        );

        return did_some_work;
    }

    /**
     * Print simplication step of this graph, if PRINT_SIMPLIFY_GRAPHS is enabled
     * @param file the destination for the graph DOT file
     */
    fn print_graph_simplification(&self, path: &str) {
        if Self::PRINT_SIMPLIFY_GRAPHS {
            self.base_graph
                .subset_to_neighbours(self.base_graph.get_reference_source_vertex().unwrap(), 5)
                .print_graph(path, true, 0)
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
    pub fn zip_linear_chains(&mut self) -> bool {
        // create the list of start sites [doesn't modify graph yet]
        let mut zip_starts = Vec::new();
        for source in self.base_graph.graph.node_indices() {
            if self.is_linear_chain_start(source) {
                zip_starts.push(source)
            }
        }

        if zip_starts.is_empty() {
            return false;
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

        return merged_one;
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
        return self.base_graph.out_degree_of(source) == 1
            && (self.base_graph.in_degree_of(source) != 1
                || self.base_graph.out_degree_of(
                    *self
                        .base_graph
                        .incoming_vertices_of(source)
                        .iter()
                        .next()
                        .unwrap(),
                ) > 1);
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
        let mut linear_chain = vec![zip_start];

        let mut last_is_ref = self.base_graph.is_reference_node(zip_start); // remember because this calculation is expensive
        let mut last = zip_start;

        loop {
            if self.base_graph.out_degree_of(last) != 1 {
                // cannot extend a chain from last if last has multiple outgoing branches
                break;
            }

            // there can only be one (outgoing edge of last) by contract
            let target = self.base_graph.get_edge_target(
                match self
                    .base_graph
                    .graph
                    .edges_directed(last, Direction::Outgoing)
                    .next()
                {
                    Some(edge_index) => edge_index.id(),
                    None => break,
                },
            );

            if self.base_graph.in_degree_of(target) != 1 || last == target {
                // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
                break;
            }

            let target_is_ref = self.base_graph.is_reference_node(target);
            if last_is_ref != target_is_ref {
                // both our isRef states must be equal
                break;
            }

            linear_chain.push(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            last_is_ref = target_is_ref;
        }

        return linear_chain;
    }

    /**
     * Merge a linear chain of vertices into a single combined vertex, and update this graph to such that
     * the incoming edges into the first element of the linearChain and the outgoing edges from linearChain.getLast()
     * all point to this new combined vertex.
     *
     * @param linearChain a non-empty chain of vertices that can be zipped up into a single vertex
     * @return true if we actually merged at least two vertices together
     */
    fn merge_linear_chain(&mut self, linear_chain: &Vec<NodeIndex>) -> bool {
        match self.merge_linear_chain_vertex(linear_chain) {
            Some(index) => true,
            None => false,
        }
    }

    fn merge_linear_chain_vertex(&mut self, linear_chain: &Vec<NodeIndex>) -> Option<NodeIndex> {
        assert!(
            !linear_chain.is_empty(),
            "Cannot have linear chain with 0 elements but got {:?}",
            linear_chain
        );

        let first = linear_chain.first().unwrap();
        let last = linear_chain.last().unwrap();

        if first == last {
            return None; // only one element in the chain, cannot be extended
        }

        // create the combined vertex, and add it to the graph
        // TODO -- performance problem -- can be optimized if we want
        let added_vertex = self.merge_linear_chain_vertices(linear_chain);

        let added_node_index = self.base_graph.add_node(added_vertex);

        // update the incoming and outgoing edges to point to the new vertex
        for edge in self.base_graph.edges_directed(*last, Direction::Outgoing) {
            self.base_graph.graph.add_edge(
                added_node_index,
                self.base_graph.get_edge_target(edge),
                self.base_graph.graph.edge_weight(edge).unwrap().clone(),
            );
        }
        for edge in self.base_graph.edges_directed(*first, Direction::Incoming) {
            self.base_graph.graph.add_edge(
                self.base_graph.get_edge_source(edge),
                added_node_index,
                self.base_graph.graph.edge_weight(edge).unwrap().clone(),
            );
        }

        self.base_graph
            .graph
            .retain_nodes(|gr, v| !linear_chain.contains(&v));

        return Some(added_node_index);
    }

    fn merge_linear_chain_vertices(&self, vertices: &Vec<NodeIndex>) -> SeqVertex {
        let seqs = vertices
            .iter()
            .map(|v| {
                self.base_graph
                    .graph
                    .node_weight(*v)
                    .unwrap()
                    .get_sequence()
            })
            .collect::<Vec<&[u8]>>();

        let seqs_joined = seqs.concat();

        return SeqVertex::new(seqs_joined);
    }
}
