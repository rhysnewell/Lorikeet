use graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use graphs::base_vertex::BaseVertex;
use graphs::path::Path;
use graphs::seq_graph::SeqGraph;
use graphs::seq_vertex::SeqVertex;
use hashlink::LinkedHashSet;
use petgraph::prelude::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::{Externals, IndexType, StableDiGraph};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use petgraph::{algo, Directed, EdgeType};
use rayon::prelude::*;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::hash::Hash;
use std::io::Write;
use utils::base_utils::BaseUtils;

/**
 * Common code for graphs used for local assembly.
 */
#[derive(Debug, Clone)]
pub struct BaseGraph<V: BaseVertex + Hash, E: BaseEdge> {
    kmer_size: usize,
    pub graph: StableDiGraph<V, E>,
}

impl<V: BaseVertex + Hash, E: BaseEdge> Eq for BaseGraph<V, E> {}

impl<V: BaseVertex + Hash, E: BaseEdge> PartialEq for BaseGraph<V, E> {
    fn eq(&self, other: &Self) -> bool {
        self.graph_equals(other)
    }
}

impl<V: BaseVertex + Hash, E: BaseEdge> BaseGraph<V, E> {
    pub fn new(kmer_size: usize) -> BaseGraph<V, E> {
        BaseGraph {
            kmer_size,
            graph: StableDiGraph::<V, E>::new(),
        }
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    pub fn to_sequence_graph(&self) -> SeqGraph<BaseEdgeStruct> {
        let mut seq_graph = SeqGraph::new(self.kmer_size);
        let mut vertex_map = HashMap::new();

        // create all of the equivalent seq graph vertices
        for dv in self.graph.node_indices() {
            let v_weight = self.graph.node_weight(dv).unwrap();
            let mut sv = SeqVertex::new(
                v_weight
                    .get_additional_sequence(self.is_source(dv))
                    .to_vec(),
            );
            sv.set_additional_info(v_weight.get_additional_info());
            let sv_ind = seq_graph.base_graph.add_node(&sv);
            vertex_map.insert(dv, sv_ind);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for e in self.graph.edge_indices() {
            let seq_in_v = *vertex_map.get(&self.get_edge_source(e)).unwrap();
            let seq_out_v = *vertex_map.get(&self.get_edge_target(e)).unwrap();
            let e_weight = self.graph.edge_weight(e).unwrap();
            seq_graph.base_graph.graph.add_edge(
                seq_in_v,
                seq_out_v,
                BaseEdgeStruct::new(e_weight.is_ref(), e_weight.get_multiplicity(), 0),
            );
        }

        return seq_graph;
    }

    pub fn vertex_set(&self) -> HashSet<&V> {
        self.graph.node_weights().collect::<HashSet<&V>>()
    }

    pub fn edge_set(&self) -> HashSet<&E> {
        self.graph.edge_weights().collect::<HashSet<&E>>()
    }

    pub fn edges_directed(&self, node: NodeIndex, direction: Direction) -> Vec<EdgeIndex> {
        self.graph
            .edges_directed(node, direction)
            .map(|e| e.id())
            .collect::<Vec<EdgeIndex>>()
    }

    /**
     * Checks whether the graph contains all the vertices in a collection.
     *
     * @param vertices the vertices to check. Must not be null and must not contain a null.
     *
     * @throws IllegalArgumentException if {@code vertices} is {@code null}.
     *
     * @return {@code true} if all the vertices in the input collection are present in this graph.
     * Also if the input collection is empty. Otherwise it returns {@code false}.
     */
    pub fn contains_all_vertices<'a, I: IntoIterator<Item = &'a NodeIndex>>(
        &self,
        nodes: I,
    ) -> bool {
        return nodes.into_iter().all(|v| self.graph.contains_node(*v));
    }

    pub fn contains_edge(&self, edge: EdgeIndex) -> bool {
        self.graph.edge_endpoints(edge).is_some()
    }

    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    pub fn is_reference_node(&self, vertex_index: NodeIndex) -> bool {
        if self
            .graph
            .edges_directed(vertex_index, Direction::Incoming)
            .into_iter()
            .any(|e| e.weight().is_ref())
            || self
                .graph
                .edges_directed(vertex_index, Direction::Outgoing)
                .into_iter()
                .any(|e| e.weight().is_ref())
        {
            return true;
        } else {
            self.vertex_set().len() == 1
        }
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    pub fn is_source(&self, vertex_index: NodeIndex) -> bool {
        return self.in_degree_of(vertex_index) == 0;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    pub fn is_sink(&self, veretex_index: NodeIndex) -> bool {
        return self.out_degree_of(veretex_index) == 0;
    }

    /**
     * @param v the vertex to test
     * @return number of incoming edges
     */
    pub fn in_degree_of(&self, vertex_index: NodeIndex) -> usize {
        self.graph
            .edges_directed(vertex_index, Direction::Incoming)
            .count()
    }

    /**
     * @param v the vertex to test
     * @return number of outgoing edges
     */
    pub fn out_degree_of(&self, vertex_index: NodeIndex) -> usize {
        self.graph
            .edges_directed(vertex_index, Direction::Outgoing)
            .count()
    }

    /**
     * Get the set of source vertices of this graph
     * NOTE: We return a BTreeSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being BTreeSet
     * @return a non-null set
     */
    pub fn get_sources(&self) -> VecDeque<NodeIndex> {
        return self
            .graph
            .externals(Direction::Incoming)
            .collect::<VecDeque<NodeIndex>>();
    }

    pub fn get_sources_generic(&self) -> Externals<'_, V, Directed, u32> {
        return self.graph.externals(Direction::Incoming);
    }

    /**
     * Get the set of sink vertices of this graph
     * NOTE: We return a BTreeSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being BTreeSet
     * @return a non-null set
     */
    pub fn get_sinks(&self) -> BinaryHeap<NodeIndex> {
        return self
            .graph
            .externals(Direction::Outgoing)
            .collect::<BinaryHeap<NodeIndex>>();
    }

    pub fn get_sinks_generic(&self) -> Externals<'_, V, Directed, u32> {
        return self.graph.externals(Direction::Outgoing);
    }

    pub fn compare_paths(&self, first_path: &Path, second_path: &Path) -> Ordering {
        return BaseUtils::bases_comparator(
            &first_path.get_bases(&self),
            &second_path.get_bases(&self),
        );
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * NOTE: We return a HashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being HashSets.
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    pub fn incoming_vertices_of(&self, v: NodeIndex) -> LinkedHashSet<NodeIndex> {
        return self
            .graph
            .neighbors_directed(v, Direction::Incoming)
            .collect::<LinkedHashSet<NodeIndex>>();
    }

    /**
     * Get the set of vertices connected to v by outgoing edges
     * NOTE: We return a HashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being HashSets.
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    pub fn outgoing_vertices_of(&self, v: NodeIndex) -> LinkedHashSet<NodeIndex> {
        return self
            .graph
            .neighbors_directed(v, Direction::Outgoing)
            .collect::<LinkedHashSet<NodeIndex>>();
    }

    /**
     * Get the incoming edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single incoming edge to v, or null if none exists
     */
    pub fn incoming_edge_of(&self, v: NodeIndex) -> Option<EdgeIndex> {
        self.get_directed_singleton_edge(v, Direction::Incoming)
    }

    /**
     * Get the incoming edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single incoming edge to v, or null if none exists
     */
    pub fn outgoing_edge_of(&self, v: NodeIndex) -> Option<EdgeIndex> {
        self.get_directed_singleton_edge(v, Direction::Outgoing)
    }

    fn get_directed_singleton_edge(&self, v: NodeIndex, direction: Direction) -> Option<EdgeIndex> {
        self.get_singleton_edge(
            self.graph
                .edges_directed(v, direction)
                .map(|e| e.id())
                .collect::<Vec<EdgeIndex>>(),
        )
    }

    /**
     * Helper function that gets the a single edge from edges, null if edges is empty, or
     * throws an error is edges has more than 1 element
     * @param edges a set of edges
     * @return a edge
     */
    fn get_singleton_edge(&self, edges: Vec<EdgeIndex>) -> Option<EdgeIndex> {
        assert!(
            edges.len() <= 1,
            "Cannot get a single incoming edge for a vertex with multiple incoming edges {:?}",
            edges
        );
        return if edges.is_empty() {
            None
        } else {
            Some(edges[0])
        };
    }

    pub fn get_edge_target(&self, edge: EdgeIndex) -> NodeIndex {
        self.graph.edge_endpoints(edge).unwrap().1
    }

    pub fn get_edge_source(&self, edge: EdgeIndex) -> NodeIndex {
        self.graph.edge_endpoints(edge).unwrap().0
    }

    pub fn edge_is_ref(&self, edge: EdgeIndex) -> bool {
        self.graph.edge_weight(edge).unwrap().is_ref()
    }

    /**
     * Removes all provided vertices from the graph
     */
    pub fn remove_all_vertices<'a, I>(&'a mut self, vertices: I)
    where
        I: IntoIterator<Item = &'a NodeIndex>,
    {
        for vertex in vertices.into_iter() {
            self.graph.remove_node(*vertex);
        }
    }

    /**
     * Removes all provided edges from the graph
     */
    pub fn remove_all_edges<'a, I>(&'a mut self, edges: I)
    where
        I: IntoIterator<Item = &'a EdgeIndex>,
    {
        for edge in edges.into_iter() {
            self.graph.remove_edge(*edge);
        }
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    pub fn is_ref_source(&self, v: NodeIndex) -> bool {
        // confirm that no incoming edges are reference edges
        if self
            .graph
            .edges_directed(v, Direction::Incoming)
            .any(|e| e.weight().is_ref())
        {
            return false;
        }

        // confirm that there is an outgoing reference edge
        if self
            .graph
            .edges_directed(v, Direction::Outgoing)
            .any(|e| e.weight().is_ref())
        {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref source, otherwise it's not
        return self.graph.node_indices().count() == 1;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    pub fn is_ref_sink(&self, v: NodeIndex) -> bool {
        // confirm that no incoming edges are reference edges
        if self
            .graph
            .edges_directed(v, Direction::Outgoing)
            .any(|e| e.weight().is_ref())
        {
            return false;
        }

        // confirm that there is an outgoing reference edge
        if self
            .graph
            .edges_directed(v, Direction::Incoming)
            .any(|e| e.weight().is_ref())
        {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref source, otherwise it's not
        return self.graph.node_indices().count() == 1;
    }

    /**
     * Remove all vertices in the graph that have in and out degree of 0
     */
    pub fn remove_singleton_orphan_vertices(&mut self) {
        // Run through the graph and clean up singular orphaned nodes
        //Note: need to collect nodes to remove first because we can't directly modify the list we're iterating over
        let to_remove = self
            .graph
            .node_indices()
            .filter(|v| self.is_singleton_orphan(*v))
            .collect::<HashSet<NodeIndex>>();
        debug!("Orphans {}", to_remove.len());
        self.remove_all_vertices(&to_remove)
    }

    fn is_singleton_orphan(&self, v: NodeIndex) -> bool {
        return self.in_degree_of(v) == 0 && self.out_degree_of(v) == 0 && !self.is_ref_source(v);
    }

    /**
     * @return the reference source vertex pulled from the graph, can be None if it doesn't exist in the graph
     */
    pub fn get_reference_source_vertex(&self) -> Option<NodeIndex> {
        return self
            .graph
            .node_indices()
            .filter(|v| self.is_ref_source(*v))
            .next();
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be None if it doesn't exist in the graph
     */
    pub fn get_reference_sink_vertex(&self) -> Option<NodeIndex> {
        return self
            .graph
            .node_indices()
            .filter(|v| self.is_ref_sink(*v))
            .next();
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @param allowNonRefPaths if true, allow sub-paths that are non-reference if there is only a single outgoing edge
     * @param blacklistedEdge optional edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
     */
    pub fn get_next_reference_vertex(
        &self,
        v: Option<NodeIndex>,
        allow_non_ref_paths: bool,
        blacklisted_edge: Option<EdgeIndex>,
    ) -> Option<NodeIndex> {
        if v.is_none() {
            return None;
        }

        let v = v.unwrap();
        let outgoing_edges = self.graph.edges_directed(v, Direction::Outgoing);

        for edge_to_test in outgoing_edges {
            if edge_to_test.weight().is_ref() {
                return Some(self.get_edge_target(edge_to_test.id()));
            }
        }

        if !allow_non_ref_paths {
            return None;
        }

        //singleton or empty set
        let outgoing_edges = self.graph.edges_directed(v, Direction::Outgoing);
        let edges = match blacklisted_edge {
            Some(blacklisted_edge) => {
                let edges = outgoing_edges
                    .filter(|e| e.id() != blacklisted_edge)
                    .map(|e| e.id())
                    .collect::<Vec<EdgeIndex>>();
                edges
            }
            None => {
                let edges = outgoing_edges.map(|e| e.id()).collect::<Vec<EdgeIndex>>();
                edges
            }
        };

        return if edges.len() == 1 {
            Some(self.get_edge_target(edges[0]))
        } else {
            None
        };
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists or null otherwise.
     */
    pub fn get_prev_reference_vertex(&self, v: Option<NodeIndex>) -> Option<NodeIndex> {
        match v {
            None => return None,
            Some(v) => self
                .graph
                .edges_directed(v, Direction::Incoming)
                .map(|e| self.get_edge_source(e.id()))
                .filter(|v| self.is_reference_node(*v))
                .take(1)
                .next(),
        }
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    pub fn print_graph(&self, destination: &str, write_header: bool, prune_factor: usize) {
        let mut graph_writer = File::create(destination).unwrap();

        if write_header {
            graph_writer.write(b"digraph assemblyGraphs {\n");
        };

        for edge in self.graph.edge_indices() {
            let edge_string = format!(
                "\t{} -> {} ",
                self.get_edge_source(edge).index(),
                self.get_edge_target(edge).index()
            );
            let mut edge_label_string;
            let edge_weight = self.graph.edge_weight(edge).unwrap();
            if edge_weight.get_multiplicity() > 0 && edge_weight.get_multiplicity() < prune_factor {
                edge_label_string = format!(
                    "[style=dotted,color=grey,label=\"{}\"];",
                    edge_weight.get_dot_label()
                );
            } else {
                edge_label_string = format!("[label=\"{}\"];", edge_weight.get_dot_label());
            }
            graph_writer.write(edge_string.as_bytes());
            graph_writer.write(edge_label_string.as_bytes());
            if edge_weight.is_ref() {
                graph_writer.write(format!("{} [color=red];\n", edge_string).as_bytes());
            }
        }

        for v in self.graph.node_indices() {
            let node_weight = self.graph.node_weight(v).unwrap();
            graph_writer.write(
                format!(
                    "\t{} [label=\"{}\",shape=box]\n",
                    node_weight.to_string(),
                    format!(
                        "{}{}",
                        std::str::from_utf8(self.get_additional_sequence(v, node_weight)).unwrap(),
                        node_weight.get_additional_info()
                    )
                )
                .as_bytes(),
            );
        }

        if write_header {
            graph_writer.write(b"}\n");
        };
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    pub fn get_additional_sequence<'a>(&self, index: NodeIndex, v: &'a V) -> &'a [u8] {
        return v.get_additional_sequence(self.is_source(index));
    }

    pub fn get_sequence_from_index<'a>(&'a self, index: NodeIndex) -> &'a [u8] {
        return self
            .graph
            .node_weight(index)
            .unwrap()
            .get_additional_sequence(self.is_source(index));
    }

    /**
     * Get the set of vertices within distance edges of source, regardless of edge direction
     *
     * @param source the source vertex to consider
     * @param distance the distance
     * @return a set of vertices within distance of source
     */
    fn vertices_within_distance(&self, source: NodeIndex, distance: usize) -> HashSet<NodeIndex> {
        if distance == 0 {
            return vec![source].into_iter().collect::<HashSet<NodeIndex>>();
        }

        let mut found = HashSet::new();
        found.insert(source);
        for v in self.graph.neighbors(source) {
            found.extend(self.vertices_within_distance(v, distance - 1))
        }

        return found;
    }

    /**
     * Get a graph containing only the vertices within distance edges of target
     * @param target a vertex in graph
     * @param distance the max distance
     * @return a non-null graph
     */
    pub fn subset_to_neighbours(&self, target: NodeIndex, distance: usize) -> BaseGraph<V, E> {
        let to_keep = self.vertices_within_distance(target, distance);
        let mut result = self.graph.clone();
        result.retain_nodes(|gr, v| to_keep.contains(&v));
        BaseGraph {
            kmer_size: self.kmer_size,
            graph: result,
        }
    }

    /**
     * Checks for the presence of directed cycles in the graph.
     *
     * @return {@code true} if the graph has cycles, {@code false} otherwise.
     */
    pub fn has_cycles(&self) -> bool {
        algo::is_cyclic_directed(&self.graph)
    }

    /**
     * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
     *
     * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
     * as it requires vertices to not only be connected by a series of directed edges but also prunes away
     * paths that do not also meet eventually with the reference sink vertex
     */
    pub fn remove_paths_not_connected_to_ref(&mut self) {
        if self.get_reference_source_vertex().is_none()
            || self.get_reference_sink_vertex().is_none()
        {
            panic!("Graph must have ref source and sink vertices");
        }

        // get the set of vertices we can reach by going forward from the ref source
        let mut on_path_from_ref_source = self.get_all_connected_nodes_from_node(
            self.get_reference_source_vertex().unwrap(),
            Direction::Outgoing,
        );

        // get the set of vertices we can reach by going backward from the ref sink
        let mut on_path_from_ref_sink = self.get_all_connected_nodes_from_node(
            self.get_reference_sink_vertex().unwrap(),
            Direction::Incoming,
        );

        // we want to remove anything that's not in both the sink and source sets
        let mut vertices_to_remove = self.graph.node_indices().collect::<HashSet<NodeIndex>>();
        on_path_from_ref_source.retain(|v| on_path_from_ref_sink.contains(v));
        vertices_to_remove.retain(|v| !on_path_from_ref_source.contains(v));
        // let not_in_both = on_path_from_ref_source
        //     .symmetric_difference(&on_path_from_ref_sink)
        //     .map(|v| *v)
        //     .collect::<HashSet<NodeIndex>>();

        // vertices_to_remove.retain(|v| not_in_both.contains(v));
        self.remove_all_vertices(&vertices_to_remove);

        // simple sanity checks that this algorithm is working.
        if self.get_sinks().len() > 1 {
            panic!(
                "Should have eliminated all but the reference sink, but found {:?}",
                self.get_sinks()
            );
        };

        if self.get_sources().len() > 1 {
            panic!(
                "Should have eliminated all but the reference sink, but found {:?}",
                self.get_sources()
            );
        };
    }

    /**
     * Iterate through the BaseGraph and return all vertices reachable by specified vertex in given
     * direction
     *
     * Note that if both followIncomingEdges and followOutgoingEdges are false, we simply return the
     * start vertex
     *
     * @param graph the graph to iterator over.  Cannot be null
     * @param start the vertex to start at.  Cannot be null
     * @param followIncomingEdges should we follow incoming edges during our
     *                            traversal? (goes backward through the graph)
     * @param followOutgoingEdges should we follow outgoing edges during out traversal?
     */
    fn get_all_connected_nodes_from_node(
        &self,
        starting_node: NodeIndex,
        direction: Direction,
    ) -> HashSet<NodeIndex> {
        let mut to_visit = BinaryHeap::new();
        let mut visited = HashSet::new();
        to_visit.push(starting_node);

        let mut v = None;
        while !to_visit.is_empty() {
            v = to_visit.pop();
            if !visited.contains(&v.unwrap()) {
                visited.insert(v.unwrap());
                to_visit.extend(
                    self.graph
                        .neighbors_directed(v.unwrap(), direction)
                        .collect::<Vec<NodeIndex>>(),
                );
            }
        }

        return visited;
    }

    /**
     * Remove edges that are connected before the reference source and after the reference sink
     *
     * Also removes all vertices that are orphaned by this process
     */
    pub fn clean_non_ref_paths(&mut self) {
        if self.get_reference_source_vertex().is_none()
            || self.get_reference_sink_vertex().is_none()
        {
            // pass
        } else {
            // Remove non-ref edges connected before and after the reference path
            let mut edges_to_check: BinaryHeap<EdgeIndex> = BinaryHeap::new();
            edges_to_check.extend(
                self.graph
                    .edges_directed(
                        self.get_reference_source_vertex().unwrap(),
                        Direction::Incoming,
                    )
                    .map(|edge_ref| edge_ref.id()),
            );

            while !edges_to_check.is_empty() {
                let e = edges_to_check.pop().unwrap();
                if !self.graph.edge_weight(e).unwrap().is_ref() {
                    edges_to_check.extend(
                        self.graph
                            .edges_directed(self.get_edge_source(e), Direction::Incoming)
                            .map(|edge_ref| edge_ref.id()),
                    );
                    self.graph.remove_edge(e);
                }
            }

            edges_to_check.extend(
                self.graph
                    .edges_directed(
                        self.get_reference_sink_vertex().unwrap(),
                        Direction::Outgoing,
                    )
                    .map(|edge_ref| edge_ref.id()),
            );

            while !edges_to_check.is_empty() {
                let e = edges_to_check.pop().unwrap();
                if !self.graph.edge_weight(e).unwrap().is_ref() {
                    edges_to_check.extend(
                        self.graph
                            .edges_directed(self.get_edge_target(e), Direction::Outgoing)
                            .map(|edge_ref| edge_ref.id()),
                    );
                    self.graph.remove_edge(e);
                }
            }

            self.remove_singleton_orphan_vertices();
        }
    }

    /**
     * Remove all vertices on the graph that cannot be accessed by following any edge,
     * regardless of its direction, from the reference source vertex
     */
    pub fn remove_vertices_not_connected_to_ref_regardless_of_edge_direction(&mut self) {
        let mut to_remove = self.graph.node_indices().collect::<HashSet<NodeIndex>>();

        let ref_v = self.get_reference_source_vertex();
        match ref_v {
            None => self.remove_all_vertices(&to_remove),
            Some(ref_v) => {
                for node in self.get_all_connected_nodes_from_node(ref_v, Direction::Incoming) {
                    to_remove.remove(&node);
                }

                for node in self.get_all_connected_nodes_from_node(ref_v, Direction::Outgoing) {
                    to_remove.remove(&node);
                }
                to_remove.remove(&ref_v);

                self.remove_all_vertices(&to_remove);
            }
        }
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    pub fn get_reference_bytes(
        &self,
        from_vertex: NodeIndex,
        to_vertex: Option<NodeIndex>,
        include_start: bool,
        include_stop: bool,
    ) -> Vec<u8> {
        let mut bytes = Vec::new();
        let mut v = Some(from_vertex);
        if include_start {
            bytes.extend(
                self.graph
                    .node_weight(v.unwrap())
                    .unwrap()
                    .get_additional_sequence(true),
            );
        };

        v = self.get_next_reference_vertex(v, false, None);
        while v.is_some() && v != to_vertex {
            bytes.extend(
                self.graph
                    .node_weight(v.unwrap())
                    .unwrap()
                    .get_additional_sequence(true),
            );
            v = self.get_next_reference_vertex(v, false, None);
        }

        if include_stop && v.is_some() && v == to_vertex {
            bytes.extend(
                self.graph
                    .node_weight(v.unwrap())
                    .unwrap()
                    .get_additional_sequence(true),
            );
        };

        return bytes;
    }

    pub fn add_vertices<'a, I>(&mut self, vertices: I) -> Vec<NodeIndex>
    where V: BaseVertex + 'a, I: IntoIterator<Item=&'a V>
    {
        let node_indices = vertices
            .into_iter()
            .map(|v| self.add_node(v))
            .collect::<Vec<NodeIndex>>();

        return node_indices;
    }

    pub fn add_edges(&mut self, start: NodeIndex, remaining: Vec<NodeIndex>, template: E) {
        let mut prev = start;
        for next in remaining {
            self.graph.add_edge(prev, next, template.clone());
            prev = next;
        }
    }

    /**
     * Gets all node weights in vector of node indices. Returns a vector of Options pointing to references
     */
    pub fn get_node_weights(&self, node_indices: &Vec<NodeIndex>) -> Vec<Option<&V>> {
        node_indices
            .into_iter()
            .map(|n| self.graph.node_weight(*n))
            .collect::<Vec<Option<&V>>>()
    }

    /**
     * Semi-lenient comparison of two graphs, truing true if g1 and g2 have similar structure
     *
     * By similar this means that both graphs have the same number of vertices, where each vertex can find
     * a vertex in the other graph that's seqEqual to it.  A similar constraint applies to the edges,
     * where all edges in g1 must have a corresponding edge in g2 where both source and target vertices are
     * seqEqual
     *
     * @param g1 the first graph to compare
     * @param g2 the second graph to compare
     * @param <T> the type of the nodes in those graphs
     * @return true if g1 and g2 are equals
     */
    pub fn graph_equals(&self, other: &Self) -> bool {
        let vertices_1 = self.graph.node_indices().collect::<HashSet<NodeIndex>>();
        let vertices_2 = other.graph.node_indices().collect::<HashSet<NodeIndex>>();
        let edges_1 = self.graph.edge_indices().collect::<HashSet<EdgeIndex>>();
        let edges_2 = other.graph.edge_indices().collect::<HashSet<EdgeIndex>>();

        if vertices_1.len() != vertices_2.len() || edges_1.len() != edges_2.len() {
            return false;
        }

        if vertices_1.len() == 0
            && vertices_2.len() == 0
            && edges_1.len() == 0
            && edges_2.len() == 0
        {
            return true;
        }

        //for every vertex in g1 there is a vertex in g2 with an equal getSequenceString
        let ok = vertices_1
            .iter()
            .map(|v1| self.graph.node_weight(*v1).unwrap().get_sequence())
            .all(|v1_seq_string| {
                vertices_2
                    .iter()
                    .any(|v2| v1_seq_string == other.graph.node_weight(*v2).unwrap().get_sequence())
            });
        if !ok {
            return false;
        }

        //for every edge in g1 there is an equal edge in g2
        let ok_g1 = edges_1
            .iter()
            .all(|e1| edges_2.iter().any(|e2| self.seq_equals(*e1, *e2, &other)));
        if !ok_g1 {
            return false;
        }

        return edges_2
            .iter()
            .all(|e2| edges_1.iter().any(|e1| other.seq_equals(*e2, *e1, self)));
    }

    fn seq_equals(&self, edge_1: EdgeIndex, edge_2: EdgeIndex, other: &Self) -> bool {
        return self
            .graph
            .node_weight(self.get_edge_source(edge_1))
            .unwrap()
            .get_sequence()
            == other
                .graph
                .node_weight(other.get_edge_source(edge_2))
                .unwrap()
                .get_sequence()
            && self
                .graph
                .node_weight(self.get_edge_target(edge_1))
                .unwrap()
                .get_sequence()
                == other
                    .graph
                    .node_weight(other.get_edge_target(edge_2))
                    .unwrap()
                    .get_sequence();
    }

    /**
     * Add nodes to the graph, if the node is already present in the graph then retrieve and return
     * its current index. Use this instead of petgraph's usual add_node as we do not want to be able
     * have identical nodes especially when dealing with k-mers
     * **Returns** the NodeIndex of the provided vertex/node
     */
    pub fn add_node(&mut self, v: &V) -> NodeIndex {
        if v.merge_identical_nodes() {
            let index = self.graph.node_indices().find(|i| &self.graph[*i] == v);
            match index {
                None => return self.graph.add_node(v.clone()),
                Some(index) => return index,
            }
        } else {
            return self.graph.add_node(v.clone());
        }
    }

    /**
     * Add edge between source -> target if none exists, or add e to an already existing one if present
     *
     * @param source source vertex
     * @param target vertex
     * @param e edge to add
     */
    pub fn add_or_update_edge(&mut self, source: NodeIndex, target: NodeIndex, e: E) {
        let prev = self.graph.find_edge(source, target);
        match prev {
            None => {
                self.graph.add_edge(source, target, e);
            }
            Some(prev) => {
                let mut prev_weight = self.graph.edge_weight_mut(prev).unwrap();
                prev_weight.add(e);
            }
        }
    }

    /**
     * Returns the index of a given node if it is found in the graph. None if not found.
     */
    pub fn index_from_vertex(&self, v: &V) -> Option<NodeIndex> {
        let index = self.graph.node_indices().find(|i| &self.graph[*i] == v);
        return index;
    }
}

pub struct TestGraph {
    pub graph: BaseGraph<MultiDeBruijnVertex, BaseEdgeStruct>,
}

impl TestGraph {
    pub fn new(kmer_size: usize) -> Self {
        Self {
            graph: BaseGraph::new(kmer_size),
        }
    }

    /**
     * Only used for testing purposes
     * Add edge to assembly graph connecting the two kmers
     * @param kmer1 the source kmer for the edge
     * @param kmer2 the target kmer for the edge
     * @param isRef true if the added edge is a reference edge
     */
    pub fn add_kmers_to_graph(
        &mut self,
        kmer_1: &[u8],
        kmer_2: &[u8],
        is_ref: bool,
        multiplicity: usize,
    ) {
        assert_eq!(
            kmer_1.len(),
            kmer_2.len(),
            "Attempting to add kmers of differing lengths"
        );
        let v1 = MultiDeBruijnVertex::new(kmer_1.to_vec(), true);
        let v2 = MultiDeBruijnVertex::new(kmer_2.to_vec(), true);
        let to_add = BaseEdgeStruct::new(is_ref, multiplicity, 0);

        let node_indices = self.graph.add_vertices(vec![&v1, &v2]);
        self.graph
            .add_or_update_edge(node_indices[0], node_indices[1], to_add);
    }

    // /**
    //  * Convert this kmer graph to a simple sequence graph.
    //  *
    //  * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
    //  * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
    //  *
    //  * @return a newly allocated SequenceGraph
    //  */
    // pub fn to_sequence_graph(&self) -> SeqGraph<BaseEdgeStruct>
}
