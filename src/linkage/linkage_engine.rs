use coverm::bam_generator::generate_indexed_named_bam_readers_from_bam_files;
use coverm::bam_generator::IndexedNamedBamReader;
use coverm::bam_generator::NamedBamReaderGenerator;
use hashlink::{LinkedHashMap, LinkedHashSet};
use itertools::Itertools;
use model::byte_array_allele::Allele;
use model::variant_context::VariantContext;
use ndarray::{Array1, Array2};
use ordered_float::OrderedFloat;
use petgraph::algo::{all_simple_paths, min_spanning_tree, tarjan_scc};
use petgraph::data::{Element, FromElements};
use petgraph::prelude::{NodeIndex, UnGraph};
use petgraph::{Direction, Undirected};
use rayon::prelude::*;
use rust_htslib::bam::Record;
use std::cmp::min;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// LinkageEngine aims to take a set of variant clusters and link them back together into likely
/// strain genomes. It does this by taking all of the reads that mapped to all of the variants in a
/// variant grouping and counting how many reads are shared between groups. To do this effectively,
/// It builds an undirected graph where each node is a variant group and each edge represents a
/// shared read, where the multiplicity of the edge represents how many reads were shared.
pub struct LinkageEngine<'a> {
    grouped_contexts: LinkedHashMap<i32, Vec<&'a VariantContext>>,
    grouped_mean_read_depth: LinkedHashMap<i32, f64>,
    samples: &'a Vec<String>,
    cluster_separations: &'a Array2<f64>,
}

impl<'a> LinkageEngine<'a> {
    const MIN_DETECTABLE_DEPTH_EPSILON: f64 = 0.25;

    pub fn new(
        grouped_contexts: LinkedHashMap<i32, Vec<&'a VariantContext>>,
        samples: &'a Vec<String>,
        cluster_separations: &'a Array2<f64>,
    ) -> LinkageEngine<'a> {
        Self {
            grouped_contexts,
            grouped_mean_read_depth: LinkedHashMap::new(),
            samples,
            cluster_separations,
        }
    }

    pub fn retrieve_grouped_contexts(mut self) -> LinkedHashMap<i32, Vec<&'a VariantContext>> {
        self.grouped_contexts
    }

    /// Runs the linkage algorithm on the grouped variant contexts and returns
    /// an ordered vector of length n where n is number of variant contexts. Each index of the vector
    /// corresponds to the index of the variant contexts present in the HaplotypeClusteringEngine.
    /// Each element of the vector is another vector containing the strain indices
    /// that are associated with the corresponding variant context
    pub fn run_linkage(
        mut self,
        indexed_bam_readers: &Vec<String>,
        n_threads: usize,
    ) -> Vec<Vec<i32>> {
        let read_ids_in_groups = self.get_reads_for_groups(indexed_bam_readers, n_threads);
        debug!("group mean read depths {:?}", &self.grouped_mean_read_depth);
        let graph = self.build_graph(read_ids_in_groups);
        debug!("Graph {:?}", &graph);
        if graph.edge_count() == 0 {
            // no connection formed, so each variant group is its own strain
            return graph.node_weights().map(|n| vec![*n]).collect();
        }
        let connected_components = self.extract_connected_components(graph);

        self.compute_strain_denominations(connected_components)
    }

    /// Compute the different denominations of strain groupings from a given component.
    /// The minimum spanning tree for a component is calculated. The tree is then rooted using
    /// the tip node with the highest read count/depth. The tip nodes of the tree are then visited
    /// in ascending order of read depth. Low coverage nodes are visited first, tracing path from
    /// this node back to the rooting node. If the starting node has not been visited before or if
    /// the current cumulative depth of this node (i.e. the read depth accumulated from previous visits)
    /// divided by the nodes total depth is greater than some threshold then trace the path from this
    /// node back to the rooting node, incrementing each node by coverage difference of the current node,
    /// and treat that path as a strain.
    /// Cumulative depth can be thought of as a rising water table, and if the start of path sits above
    /// the water table then save that path and then raise the water table to the current nodes read depth.
    /// TODO: Test how this performs on more complex MSTs
    fn compute_strain_denominations(
        &self,
        connected_components: Vec<UnGraph<i32, f64>>,
    ) -> Vec<Vec<i32>> {
        let mut all_strains = Vec::new();
        for component_graph in connected_components {
            // compute the minimum spanning tree
            let mut mst = UnGraph::from_elements(min_spanning_tree(&component_graph));
            debug!("MST {:?}", &mst);

            // sorted list of nodes in ascending order of read depth, only checking external nodes
            // i.e. nodes without incoming edges or outgoing edges. i.e. edge count <= 1
            let mut starting_nodes_vec = mst
                .node_indices()
                .filter(|node| mst.edges(*node).count() <= 1)
                .map(|node| {
                    (
                        Reverse(OrderedFloat(
                            *self
                                .grouped_mean_read_depth
                                .get(mst.node_weight(node).unwrap())
                                .unwrap(),
                        )),
                        node,
                    )
                })
                .collect::<Vec<(Reverse<OrderedFloat<f64>>, NodeIndex)>>();

            starting_nodes_vec.par_sort_unstable();

            // summit - all paths must lead to here
            // This is the highest depth terminating node. i.e. This node is a tip of the tree
            // with the highest depth compared to all other tips
            let highest_depth_node = starting_nodes_vec.first().unwrap().clone().1;

            // Turn the vec into a BinaryHeap
            let mut starting_nodes = BinaryHeap::from(starting_nodes_vec);
            debug!(
                "Starting nodes {:?} highest depth node {:?} VG {}",
                &starting_nodes,
                highest_depth_node,
                mst.node_weight(highest_depth_node).unwrap()
            );

            // Can't have more strains than nodes in the tree
            let mut strains = Vec::with_capacity(mst.node_count());
            let mut seen_nodes = HashSet::with_capacity(mst.node_count());
            // Container holding the cumulative depth of each node in this MST
            // The value for a node will increment by the last starting node if the path
            // from the lowest depth node passed through this node. The value will increment by the
            // depth of the starting node
            let mut nodes_cumulative_depth = HashMap::with_capacity(mst.node_count());
            while !starting_nodes.is_empty() {
                let current_node_info = starting_nodes.pop().unwrap();
                let current_depth = current_node_info.0 .0 .0; // I hate this??? wot
                let current_node = current_node_info.1;

                let current_node_cumulative_depth =
                    *nodes_cumulative_depth.entry(current_node).or_insert(0.0);

                if (1.0 - (current_node_cumulative_depth / current_depth))
                    >= Self::MIN_DETECTABLE_DEPTH_EPSILON
                    || !seen_nodes.contains(mst.node_weight(current_node).unwrap())
                {
                    let paths = all_simple_paths::<LinkedHashSet<NodeIndex>, _>(
                        &mst,
                        current_node,
                        highest_depth_node,
                        0,
                        None,
                    )
                    .collect::<Vec<LinkedHashSet<NodeIndex>>>();
                    debug!(
                        "Paths {:?} current node {:?} depth {} cumulative depth {} value {}",
                        &paths,
                        current_node,
                        current_depth,
                        current_node_cumulative_depth,
                        (1.0 - (current_node_cumulative_depth / current_depth))
                    );

                    if paths.len() == 1 {
                        strains.push(
                            paths
                                .into_iter()
                                .next()
                                .unwrap()
                                .into_iter()
                                .enumerate()
                                .map(|(idx, node)| {
                                    let variant_group = *mst.node_weight(node).unwrap();
                                    seen_nodes.insert(variant_group);
                                    let node_cumulative_depth =
                                        nodes_cumulative_depth.entry(node).or_insert(0.0);
                                    // add the depth of the current node to this node
                                    // subtract the cumulative depth as it has already
                                    // been seen that many times
                                    *node_cumulative_depth +=
                                        (current_depth - current_node_cumulative_depth);

                                    // check if this the next node in the path, if so append to
                                    // binary heap as it will act as the next branch tip once current
                                    // tips run out
                                    if idx == 1 {
                                        starting_nodes
                                            .push((Reverse(OrderedFloat(current_depth)), node));
                                    };
                                    variant_group
                                })
                                .collect::<Vec<i32>>(),
                        );
                    }
                }
            }

            // check if highest depth node is still above the water table after all other paths
            // have been visited. This also should catch component graphs which contained only
            // a single node
            let highest_depth = *self
                .grouped_mean_read_depth
                .get(mst.node_weight(highest_depth_node).unwrap())
                .unwrap();

            let highest_depth_node_cumulative_depth = *nodes_cumulative_depth
                .entry(highest_depth_node)
                .or_insert(0.0);
            if (1.0 - (highest_depth_node_cumulative_depth / highest_depth))
                >= Self::MIN_DETECTABLE_DEPTH_EPSILON
                || !seen_nodes.contains(mst.node_weight(highest_depth_node).unwrap())
            {
                let variant_group = *mst.node_weight(highest_depth_node).unwrap();
                seen_nodes.insert(variant_group);
                strains.push(vec![variant_group]);
            }

            all_strains.extend(strains)
        }

        all_strains
    }

    /// Extracts the connect components of the built graph and turns them into new independent graphs
    /// a connected component is defined as an induced subgraph in which any two vertices are
    /// connected to each other by paths, and which is connected to no additional vertices in the rest of the graph
    /// By contract, the input graph must have only one edge connected any two nodes.
    /// For example, for a directed graph:
    /// G:
    /// a ----> b       e ----> f
    /// ^       |       ^       |
    /// |       v       |       v
    /// d <---- c       h <---- g
    /// represents a graph with two components
    ///
    /// This function would return:
    /// G1:
    /// a ----> b
    /// ^       |
    /// |       v
    /// d <---- c
    ///
    /// G2:
    /// e ----> f
    /// ^       |
    /// |       v
    /// h <---- g
    fn extract_connected_components(&self, graph: UnGraph<i32, f64>) -> Vec<UnGraph<i32, f64>> {
        let tarjan_components = tarjan_scc(&graph);

        // shortcut for when only one component
        if tarjan_components.len() == 1 {
            return vec![graph];
        }

        let mut result = Vec::with_capacity(tarjan_components.len());
        for component in tarjan_components {
            // container for previously seen nodes and their updated value
            let mut already_seen = LinkedHashMap::with_capacity(component.len());
            if component.len() == 1 {
                // Single node graph
                let mut sub_graph = UnGraph::new_undirected();
                let new_node1 =
                    Self::extract_node(&graph, &mut sub_graph, &mut already_seen, component[0]);
                result.push(sub_graph);
                continue;
            }

            let mut sub_graph = UnGraph::new_undirected();
            for nodes in component.into_iter().combinations(2) {
                let node1 = nodes[0];
                let node2 = nodes[1];
                if graph.contains_edge(node1, node2) {
                    let new_node1 =
                        Self::extract_node(&graph, &mut sub_graph, &mut already_seen, node1);

                    let new_node2 =
                        Self::extract_node(&graph, &mut sub_graph, &mut already_seen, node2);

                    // Should be only one edge
                    let weight = graph.edges_connecting(node1, node2).next().unwrap();
                    sub_graph.add_edge(new_node1, new_node2, *weight.weight());
                }
            }

            result.push(sub_graph)
        }

        return result;
    }

    /// Get the ids of reads that map to the padded area around the variant locations in a variant
    /// group
    fn get_reads_for_groups(
        &mut self,
        indexed_bam_readers: &Vec<String>,
        n_threads: usize,
    ) -> LinkedHashMap<i32, HashSet<String>> {
        let mut all_grouped_reads = LinkedHashMap::with_capacity(self.grouped_contexts.len());
        indexed_bam_readers
            .par_iter()
            .enumerate()
            .map(|(sample_idx, bam_generator)| {
                // Get the appropriate sample index based on how many references we are using
                let bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                    vec![&bam_generator],
                    n_threads as u32,
                )
                .into_iter()
                .next()
                .unwrap();

                let mut bam_generated = bam_generator.start();

                bam_generated.set_threads(n_threads);

                let mut grouped_reads = LinkedHashMap::with_capacity(self.grouped_contexts.len());
                let mut record = Record::new();
                for (group, variants) in self.grouped_contexts.iter() {
                    for variant in variants {
                        bam_generated
                            .fetch((
                                variant.loc.tid as i32,
                                variant.loc.start as i64,
                                variant.loc.end as i64 + 1,
                            ))
                            .expect(&format!(
                                "Failed to fetch interval {}:{}-{}",
                                variant.loc.tid, variant.loc.start, variant.loc.end
                            ));

                        let records = grouped_reads.entry(group.clone()).or_insert(HashSet::new()); // container for the records to be collected

                        while bam_generated.read(&mut record) == true {
                            // be very lenient with filtering
                            if record.is_unmapped() {
                                continue;
                            }

                            // TODO: Filter for only reads that contain the variant in question
                            //       not sure how to this besides aligning or other expensive method?
                            //       currently just check against reference
                            let mut read_index = record.pos() - variant.loc.start as i64;
                            if read_index < 0 {
                                read_index = variant.loc.start as i64 - record.pos();
                            }
                            if variant.get_alternate_alleles()[0].get_bases()
                                == &record.seq().as_bytes()[read_index as usize
                                    ..(read_index as usize
                                        + variant.get_alternate_alleles()[0].get_bases().len())]
                            {
                                // Read containing potential alternate allele
                                let read_id = format!(
                                    "{}_{}",
                                    sample_idx,
                                    std::str::from_utf8(record.qname()).unwrap()
                                );
                                records.insert(read_id);
                            }
                        }
                    }
                }

                grouped_reads
            })
            .collect::<Vec<LinkedHashMap<i32, HashSet<String>>>>()
            .into_iter()
            .for_each(|sample_grouping| {
                for (vg, reads) in sample_grouping {
                    let all_result = all_grouped_reads.entry(vg).or_insert(HashSet::new());
                    all_result.par_extend(reads);
                }
            });

        let grouped_mean_read_depth = all_grouped_reads
            .iter()
            .map(|(group, reads)| {
                (
                    group.clone(),
                    (reads.len() as f64) / self.grouped_contexts.get(group).unwrap().len() as f64,
                )
            })
            .collect::<LinkedHashMap<i32, f64>>();

        self.grouped_mean_read_depth = grouped_mean_read_depth;
        all_grouped_reads
    }

    fn build_graph(
        &mut self,
        grouped_reads: LinkedHashMap<i32, HashSet<String>>,
    ) -> UnGraph<i32, f64> {
        let mut graph = UnGraph::new_undirected();
        let mut node_indices = LinkedHashMap::with_capacity(grouped_reads.len());
        for (group1, reads1) in grouped_reads.iter() {
            let node1 = if node_indices.contains_key(group1) {
                *node_indices.get(group1).unwrap()
            } else {
                let node = graph.add_node(*group1);
                node_indices.insert(*group1, node);
                node
            };
            for (group2, reads2) in grouped_reads.iter() {
                if group1 == group2 {
                    continue;
                }
                let node2 = if node_indices.contains_key(group2) {
                    *node_indices.get(group2).unwrap()
                } else {
                    let node = graph.add_node(*group2);
                    node_indices.insert(*group2, node);
                    node
                };

                // Don't count twice
                if !graph.contains_edge(node1, node2) {
                    // How many read ids are shared
                    let intersection = reads1.intersection(reads2).count() as f64;
                    let mut under_sep_thresh = false;
                    if *group1 != 0 && *group2 != 0 {
                        under_sep_thresh = self.cluster_separations
                            [[*group1 as usize - 1, *group2 as usize - 1]]
                            < 3.0;
                    }
                    if intersection > 0.0 || under_sep_thresh {
                        let union = min(reads1.len(), reads2.len()) as f64;

                        // The weight needs to be low for highly connected nodes and
                        // high for poorly connected nodes. This is because minimum spanning trees
                        // generate the tree with lowest edge weights
                        let weight = 1.0 - (intersection / union);

                        if weight < 1.0 {
                            // variant groups connected by reads are favoured
                            debug!(
                                "{}:{} weight {} intersection {} union {}",
                                group1, group2, weight, intersection, union
                            );
                            graph.add_edge(node1, node2, weight);
                        } else if under_sep_thresh {
                            // We will form an edge here but it will essentially
                            // have to be directly on top of the other variant group
                            // for the edge to be favoured in the MST
                            let weight = self.cluster_separations
                                [[*group1 as usize - 1, *group2 as usize - 1]]
                                / 1.5;
                            debug!(
                                "{}:{} weight {} intersection {} union {}",
                                group1, group2, weight, intersection, union
                            );
                            graph.add_edge(node1, node2, weight);
                        }
                    }
                }
            }
        }

        graph
    }

    fn extract_node(
        graph: &UnGraph<i32, f64>,
        sub_graph: &mut UnGraph<i32, f64>,
        already_seen: &mut LinkedHashMap<NodeIndex, NodeIndex>,
        node: NodeIndex,
    ) -> NodeIndex {
        if already_seen.contains_key(&node) {
            *already_seen.get(&node).unwrap()
        } else {
            // extract variant group from og graph
            let vg = graph.node_weight(node).unwrap();
            let new_node = sub_graph.add_node(*vg);
            already_seen.insert(node, new_node);
            new_node
        }
    }
}
