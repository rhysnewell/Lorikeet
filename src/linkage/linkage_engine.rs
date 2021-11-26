use bstr::ByteSlice;
use coverm::bam_generator::generate_indexed_named_bam_readers_from_bam_files;
use coverm::bam_generator::IndexedNamedBamReader;
use coverm::bam_generator::NamedBamReaderGenerator;
use coverm::FlagFilter;
use hashlink::{LinkedHashMap, LinkedHashSet};
use itertools::Itertools;
use log::{log_enabled, Level};
use model::byte_array_allele::Allele;
use model::variant_context::VariantContext;
use ndarray::{Array1, Array2};
use ordered_float::OrderedFloat;
use petgraph::algo::{all_simple_paths, min_spanning_tree, tarjan_scc};
use petgraph::data::{Element, FromElements};
use petgraph::dot::Dot;
use petgraph::prelude::{EdgeRef, NodeIndex, UnGraph};
use petgraph::{Direction, Undirected};
use rayon::prelude::*;
use rust_htslib::bam::Record;
use std::cmp::{min, max};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::{read, File};
use std::io::Write;
use std::path::Path;
use petgraph::graph::EdgeIndex;

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
    previous_groups: &'a HashMap<i32, i32>,
    exclusive_groups: &'a HashMap<i32, HashSet<i32>>,
}

impl<'a> LinkageEngine<'a> {
    const MIN_DETECTABLE_DEPTH_EPSILON: f64 = 0.35;

    pub fn new(
        grouped_contexts: LinkedHashMap<i32, Vec<&'a VariantContext>>,
        samples: &'a Vec<String>,
        cluster_separations: &'a Array2<f64>,
        previous_groups: &'a HashMap<i32, i32>,
        exclusive_groups: &'a HashMap<i32, HashSet<i32>>,
    ) -> LinkageEngine<'a> {
        Self {
            grouped_contexts,
            grouped_mean_read_depth: LinkedHashMap::new(),
            samples,
            cluster_separations,
            previous_groups,
            exclusive_groups,
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
        output_path: &str,
        flag_filters: &FlagFilter,
    ) -> Vec<LinkedHashSet<i32>> {
        let read_ids_in_groups =
            self.get_reads_for_groups(indexed_bam_readers, flag_filters, n_threads);
        debug!("group mean read depths {:?}", &self.grouped_mean_read_depth);
        let graph = self.build_graph(read_ids_in_groups);
        debug!("Graph {:?}", &graph);
        if log_enabled!(Level::Debug) {
            let output_dot = format!("{}_vg_graph.dot", output_path);
            let file_path = Path::new(&output_dot);

            let mut file_open = File::create(file_path)
                .unwrap_or_else(|_| panic!("Unable to create dot file: {}", output_dot));

            writeln!(file_open, "{:?}", Dot::new(&graph));
        }
        if graph.edge_count() == 0 {
            // no connection formed, so each variant group is its own strain
            return graph
                .node_weights()
                .map(|n| {
                    let mut new_strain = LinkedHashSet::with_capacity(1);
                    new_strain.insert(*n);
                    new_strain
                })
                .collect();
        }
        let connected_components = self.extract_connected_components(graph);

        self.compute_strain_denominations(connected_components, output_path)
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
        output_path: &str,
    ) -> Vec<LinkedHashSet<i32>> {
        let mut all_strains = Vec::new();
        for (idx, mut component_graph) in connected_components.into_iter().enumerate() {
            // let max_depth_of_component = self.max_depth_of_graph(&component_graph);
            // let min_depth_of_component = self.min_depth_of_graph(&component_graph);
            // (0..component_graph.edge_count()).for_each(|ei| {
            //
            //     // Alter the edge weights by the depth of the connected nodes in order to push
            //     // low weight nodes to externals of MST
            //     let node1 = component_graph.raw_edges()[ei].source();
            //     let node2 = component_graph.raw_edges()[ei].target();
            //     let mut mean_depth =
            //         *self.grouped_mean_read_depth.get(component_graph.node_weight(node1).unwrap()).unwrap() +
            //         *self.grouped_mean_read_depth.get(component_graph.node_weight(node2).unwrap()).unwrap();
            //     mean_depth = mean_depth / 2.0;
            //     let depth_factor = min_depth_of_component / mean_depth; // closer to the minimum means the edge weight gets close to be doubles
            //     let ew = component_graph.edge_weight_mut(EdgeIndex::new(ei)).unwrap();
            //     *ew = *ew + *ew * depth_factor;
            // });

            // compute the minimum spanning tree
            let mut mst = UnGraph::from_elements(min_spanning_tree(&component_graph));
            debug!("MST {:?}", &mst);
            if log_enabled!(Level::Debug) {
                let output_dot = format!("{}_mst_{}.dot", output_path, idx);
                let file_path = Path::new(&output_dot);

                let mut file_open = File::create(file_path)
                    .unwrap_or_else(|_| panic!("Unable to create dot file: {}", output_dot));

                writeln!(file_open, "{:?}", Dot::new(&mst));
            }

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
            // This is the highest depth node. Can be internal or external.
            // let highest_depth_node = mst
            //     .node_indices()
            //     .max_by(|node_1, node_2| {
            //         self.grouped_mean_read_depth
            //             .get(mst.node_weight(*node_1).unwrap())
            //             .unwrap()
            //             .partial_cmp(
            //                 self.grouped_mean_read_depth
            //                     .get(mst.node_weight(*node_2).unwrap())
            //                     .unwrap(),
            //             )
            //             .unwrap()
            //     })
            //     .unwrap();

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

                let mut depth_being_added_to_other_nodes =
                    current_depth - current_node_cumulative_depth;

                let paths = all_simple_paths::<LinkedHashSet<NodeIndex>, _>(
                    &mst,
                    current_node,
                    highest_depth_node,
                    0,
                    None,
                )
                    .collect::<Vec<LinkedHashSet<NodeIndex>>>();
                debug!(
                    "Paths {:?} current node {:?} depth {} cumulative depth {} detection {}",
                    &paths, current_node, current_depth, current_node_cumulative_depth, 1.0 - (current_node_cumulative_depth / current_depth)
                );
                if ((1.0 - (current_node_cumulative_depth / current_depth))
                    >= Self::MIN_DETECTABLE_DEPTH_EPSILON
                    && depth_being_added_to_other_nodes > 0.0)
                    || !seen_nodes.contains(mst.node_weight(current_node).unwrap())
                {
                    debug!("Potential new strain");
                    if paths.len() == 1 {
                        let mut path = paths.into_iter().next().unwrap();
                        // check if any of the nodes in this path have been consumed entirely
                        let consumed_nodes = self.check_nodes_in_path(
                            &mst,
                            &mut path,
                            &mut nodes_cumulative_depth,
                            depth_being_added_to_other_nodes,
                            current_node,
                        );

                        match consumed_nodes {
                            None => {
                                // Nodes in path above the water table so proceed
                                self.make_strain_and_update(
                                    path,
                                    &mut seen_nodes,
                                    &mut nodes_cumulative_depth,
                                    &mut starting_nodes,
                                    &mut strains,
                                    depth_being_added_to_other_nodes,
                                    current_node_cumulative_depth,
                                    &mst,
                                )
                            }
                            Some(consumed_nodes) => {
                                if consumed_nodes.len() > 1 {
                                    // This path has multiple nodes that are already at capacity. So we need to
                                    // merge it with it's next closest relative and stop it from
                                    // flooding
                                    let groups_at_capacity = consumed_nodes
                                        .iter()
                                        .map(|consumed_node| {
                                            *mst.node_weight(*consumed_node).unwrap()
                                        })
                                        .collect::<Vec<i32>>();
                                    self.merge_paths(
                                        &mut strains,
                                        path,
                                        &mst,
                                        &component_graph,
                                        &mut seen_nodes,
                                        &mut nodes_cumulative_depth,
                                        groups_at_capacity,
                                        depth_being_added_to_other_nodes,
                                        &mut starting_nodes,
                                        current_node_cumulative_depth,
                                    );
                                } else {
                                    // This path had one node at capacity. Attempt to bridge the gap
                                    // over this node
                                    path.remove(consumed_nodes.get(0).unwrap());
                                    self.make_strain_and_update(
                                        path,
                                        &mut seen_nodes,
                                        &mut nodes_cumulative_depth,
                                        &mut starting_nodes,
                                        &mut strains,
                                        depth_being_added_to_other_nodes,
                                        current_node_cumulative_depth,
                                        &mst,
                                    )
                                }
                            }
                        }
                    }
                } else {
                    debug!("Node below the water. Update cumulative depths...");
                    if paths.len() == 1 && current_node != highest_depth_node {
                        let mut path = paths.into_iter().next().unwrap();
                        path.into_iter()
                            .enumerate()
                            .for_each(|(idx, node)| {
                                let variant_group = *mst.node_weight(node).unwrap();
                                seen_nodes.insert(variant_group);
                                let node_cumulative_depth = nodes_cumulative_depth.entry(node).or_insert(0.0);
                                // add the depth of the current node to this node
                                *node_cumulative_depth += depth_being_added_to_other_nodes;

                                // check if this the next node in the path, if so append to
                                // binary heap as it will act as the next branch tip once current
                                // tips run out
                                if idx == 1 {
                                    starting_nodes.push((
                                        Reverse(OrderedFloat(
                                            *self
                                                .grouped_mean_read_depth
                                                .get(mst.node_weight(node).unwrap())
                                                .unwrap(),
                                        )),
                                        node,
                                    ));
                                };
                            });
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
                let mut new_strain = LinkedHashSet::with_capacity(1);
                new_strain.insert(variant_group);
                strains.push(new_strain);
            }

            all_strains.extend(strains)
        }

        all_strains
    }

    fn make_strain_and_update(
        &self,
        path: LinkedHashSet<NodeIndex>,
        seen_nodes: &mut HashSet<i32>,
        nodes_cumulative_depth: &mut HashMap<NodeIndex, f64>,
        starting_nodes: &mut BinaryHeap<(Reverse<OrderedFloat<f64>>, NodeIndex)>,
        strains: &mut Vec<LinkedHashSet<i32>>,
        current_depth: f64,
        current_node_cumulative_depth: f64,
        mst: &UnGraph<i32, f64>,
    ) {
        strains.push(
            path.into_iter()
                .enumerate()
                .map(|(idx, node)| {
                    let variant_group = *mst.node_weight(node).unwrap();
                    seen_nodes.insert(variant_group);
                    let node_cumulative_depth = nodes_cumulative_depth.entry(node).or_insert(0.0);
                    // add the depth of the current node to this node
                    *node_cumulative_depth += current_depth;

                    // check if this the next node in the path, if so append to
                    // binary heap as it will act as the next branch tip once current
                    // tips run out
                    if idx == 1 {
                        starting_nodes.push((
                            Reverse(OrderedFloat(
                                *self
                                    .grouped_mean_read_depth
                                    .get(mst.node_weight(node).unwrap())
                                    .unwrap(),
                            )),
                            node,
                        ));
                    };
                    variant_group
                })
                .collect::<LinkedHashSet<i32>>(),
        );
    }

    /// Takes potential strain path that contains a node whose variant group is at capacity, i.e.
    /// cannot be allocated to another unique strain, and finds the best already existing strain to
    /// merge it with. This is dictated by a few things, first, the strain to merged into must contain
    /// the node that is currently at capacity, and second, the strain to ber merged into must shared the
    /// most nodes with the current path. If there are more than one strain that meets this criteria,
    /// then we look back at the original component graph to see if there are any connections to between
    /// nodes in the new strain path and previous strain paths. We then merge based on connectivity.
    /// If they are again equal, we choose the longest path.
    fn merge_paths(
        &self,
        strains: &mut Vec<LinkedHashSet<i32>>,
        path: LinkedHashSet<NodeIndex>,
        mst: &UnGraph<i32, f64>,
        component_graph: &UnGraph<i32, f64>,
        seen_nodes: &mut HashSet<i32>,
        nodes_cumulative_depth: &mut HashMap<NodeIndex, f64>,
        groups_at_capacity: Vec<i32>,
        current_depth: f64,
        starting_nodes: &mut BinaryHeap<(Reverse<OrderedFloat<f64>>, NodeIndex)>,
        current_node_cumulative_depth: f64,
    ) {
        // take the path and turn it into a linked hash set of variant group ids.
        let groups_in_path = path
            .iter()
            .map(|node| *mst.node_weight(*node).unwrap())
            .collect::<LinkedHashSet<i32>>();

        // container with all the strains that are currently the best fit for this path
        let mut current_closest_strain_indices = Vec::with_capacity(strains.len());

        let mut current_max_shared_nodes = 0;
        for (index, strain) in strains.iter().enumerate() {
            // Check if this strain should be excluded
            let mut exclude_this_option = false;
            for group in groups_in_path.iter() {
                match self.exclusive_groups.get(group) {
                    None => continue,
                    Some(exclusive_groups) => {
                        for group_in_strain in strain {
                            if exclusive_groups.contains(group_in_strain) {
                                exclude_this_option = true;
                                break;
                            }
                        }
                    }
                }
            }

            if strain.iter().any(|vg| groups_at_capacity.contains(vg)) && !exclude_this_option {
                let shared_nodes = groups_in_path.intersection(strain).count();
                if shared_nodes == current_max_shared_nodes {
                    // same amount of nodes as a previous strain so extend the possible options
                    current_closest_strain_indices.push(index);
                } else if shared_nodes > current_max_shared_nodes {
                    // Update and reset the max shared nodes and current closest strain container
                    current_max_shared_nodes = shared_nodes;
                    current_closest_strain_indices.clear();
                    current_closest_strain_indices.push(index);
                }
            }
        }

        if current_closest_strain_indices.len() == 1 {
            // Only one strain shared the most nodes with this path, so easy merge
            let strain_to_merge_into = current_closest_strain_indices.into_iter().next().unwrap();

            debug!(
                "Merging vgs {:?} into strain {}",
                &groups_in_path, strain_to_merge_into
            );

            // merge and update cumulative depths for unseen nodes
            self.merge_path_and_update_unseen(
                path,
                seen_nodes,
                nodes_cumulative_depth,
                &mut strains[strain_to_merge_into],
                current_depth,
                current_node_cumulative_depth,
                mst,
            );
        } else if current_closest_strain_indices.len() > 1 {
            debug!(
                "Multiple close strains: {}",
                current_closest_strain_indices.len()
            );
            let original_node_indices_for_new_path = groups_in_path
                .iter()
                .map(|group| self.node_weight_to_node_index(component_graph, group))
                .collect::<Vec<NodeIndex>>();

            // Need to compare edges in original component graph to see which strain this new path
            // is most closely related to
            let mut index_of_max = 0;
            let mut index_of_previous = 0;
            let mut max_edge_count = 0;
            let mut previous_max_edge_count = 0;
            let mut cumulative_edge_weights_max = 0.0; // smaller cumulative edge weights are better
            let mut cumulative_edge_weights_previous = 0.0; // smaller cumualtive edge weights are better

            for strain_index in current_closest_strain_indices.into_iter() {
                let strain_variant_groups = &strains[strain_index];
                let mut edge_count = 0;
                let mut cumulative_edge_weights = 0.0;
                for old_strain_group in strain_variant_groups {
                    let old_group_node_index =
                        self.node_weight_to_node_index(component_graph, old_strain_group);
                    for new_group_node_index in &original_node_indices_for_new_path {
                        // how many edges connect the old node and this node (with weights)
                        let edges_connecting = component_graph
                            .edges_connecting(old_group_node_index, *new_group_node_index);
                        edges_connecting.into_iter().for_each(|edge| {
                            edge_count += 1;
                            cumulative_edge_weights += *edge.weight();
                        });
                    }
                }

                // update the best values
                if edge_count >= max_edge_count {
                    if edge_count == max_edge_count {
                        if cumulative_edge_weights < cumulative_edge_weights_max {
                            // update previous best
                            index_of_previous = index_of_max;
                            previous_max_edge_count = max_edge_count;
                            cumulative_edge_weights_previous = cumulative_edge_weights_max;

                            // update current best
                            index_of_max = strain_index;
                            max_edge_count = edge_count;
                            cumulative_edge_weights_max = cumulative_edge_weights;
                        } else {
                            // update previous with current
                            index_of_previous = strain_index;
                            previous_max_edge_count = edge_count;
                            cumulative_edge_weights_previous = cumulative_edge_weights;
                        }
                    } else {
                        // update previous best
                        index_of_previous = index_of_max;
                        previous_max_edge_count = max_edge_count;
                        cumulative_edge_weights_previous = cumulative_edge_weights_max;

                        // update current best
                        index_of_max = strain_index;
                        max_edge_count = edge_count;
                        cumulative_edge_weights_max = cumulative_edge_weights;
                    }
                } else if edge_count >= previous_max_edge_count {
                    if edge_count == previous_max_edge_count {
                        if cumulative_edge_weights < cumulative_edge_weights_previous {
                            // update previous with current
                            index_of_previous = strain_index;
                            previous_max_edge_count = edge_count;
                            cumulative_edge_weights_previous = cumulative_edge_weights;
                        }
                    } else {
                        // update previous with current
                        index_of_previous = strain_index;
                        previous_max_edge_count = edge_count;
                        cumulative_edge_weights_previous = cumulative_edge_weights;
                    }
                }
            }

            if max_edge_count >= previous_max_edge_count {
                if max_edge_count == previous_max_edge_count {
                    if cumulative_edge_weights_max <= cumulative_edge_weights_previous {
                        debug!(
                            "using max edge count {} weights {}: Merging vgs {:?} into strain {}",
                            max_edge_count,
                            cumulative_edge_weights_max,
                            &groups_in_path,
                            index_of_max
                        );

                        // merge into this max
                        self.merge_path_and_update_unseen(
                            path,
                            seen_nodes,
                            nodes_cumulative_depth,
                            &mut strains[index_of_max],
                            current_depth,
                            current_node_cumulative_depth,
                            mst,
                        );
                    } else {
                        debug!(
                            "using max edge count {} weights {}: Merging vgs {:?} into strain {}",
                            previous_max_edge_count,
                            cumulative_edge_weights_previous,
                            &groups_in_path,
                            index_of_previous
                        );
                        // merge into previous max
                        self.merge_path_and_update_unseen(
                            path,
                            seen_nodes,
                            nodes_cumulative_depth,
                            &mut strains[index_of_previous],
                            current_depth,
                            current_node_cumulative_depth,
                            mst,
                        );
                    }
                } else {
                    debug!(
                        "using max edge count {} weights {}: Merging vgs {:?} into strain {}",
                        max_edge_count, cumulative_edge_weights_max, &groups_in_path, index_of_max
                    );
                    self.merge_path_and_update_unseen(
                        path,
                        seen_nodes,
                        nodes_cumulative_depth,
                        &mut strains[index_of_max],
                        current_depth,
                        current_node_cumulative_depth,
                        mst,
                    );
                }
            } else {
                debug!(
                    "using max edge count {} weights {}: Merging vgs {:?} into strain {}",
                    previous_max_edge_count,
                    cumulative_edge_weights_previous,
                    &groups_in_path,
                    index_of_previous
                );
                // should never reach here?
                self.merge_path_and_update_unseen(
                    path,
                    seen_nodes,
                    nodes_cumulative_depth,
                    &mut strains[index_of_previous],
                    current_depth,
                    current_node_cumulative_depth,
                    mst,
                );
            }
        } else {
            // next closest strain had conflicting variant groups so just make this a new strain
            self.make_strain_and_update(
                path,
                seen_nodes,
                nodes_cumulative_depth,
                starting_nodes,
                strains,
                current_depth,
                current_node_cumulative_depth,
                mst,
            )
        }
    }

    fn merge_path_and_update_unseen(
        &self,
        path: LinkedHashSet<NodeIndex>,
        seen_nodes: &mut HashSet<i32>,
        nodes_cumulative_depth: &mut HashMap<NodeIndex, f64>,
        strain_to_merge_into: &mut LinkedHashSet<i32>,
        depth_being_added_to_nodes: f64,
        current_node_cumulative_depth: f64,
        mst: &UnGraph<i32, f64>,
    ) {
        // merge and update cumulative depths for unseen nodes
        path.into_iter().enumerate().for_each(|(idx, node)| {
            let variant_group = *mst.node_weight(node).unwrap();
            if !seen_nodes.contains(&variant_group) {
                seen_nodes.insert(variant_group);
            }
            let node_cumulative_depth = nodes_cumulative_depth.entry(node).or_insert(0.0);
            // add the depth of the current node to this node
            // subtract the cumulative depth as it has already
            // been seen that many times
            *node_cumulative_depth += depth_being_added_to_nodes;

            strain_to_merge_into.insert(variant_group);
        });
    }

    fn node_weight_to_node_index(&self, graph: &UnGraph<i32, f64>, node_weight: &i32) -> NodeIndex {
        NodeIndex::new(
            graph
                .node_weights()
                .position(|node| node == node_weight)
                .unwrap(),
        )
    }

    /// Checks all nodes in a path. If a node in a path has gone underneath the water table then it
    /// is likely that this path was meant to be merged with another path. Thus, this function checks
    /// that all nodes in the path are above the water table, if not then it returns the first node that
    /// is below the water table. Otherwise None
    fn check_nodes_in_path(
        &self,
        mst: &UnGraph<i32, f64>,
        path: &mut LinkedHashSet<NodeIndex>,
        nodes_cumulative_depth: &mut HashMap<NodeIndex, f64>,
        depth_being_added_to_nodes: f64,
        current_node: NodeIndex,
    ) -> Option<Vec<NodeIndex>> {
        let mut nodes_at_capacity = None;
        let current_group = *mst.node_weight(current_node).unwrap();

        let excluded_groups = self.exclusive_groups.get(&current_group);
        let mut to_remove = HashSet::new();
        for (idx, node) in path.iter().enumerate() {
            let variant_group = *mst.node_weight(*node).unwrap();
            match &excluded_groups {
                None => {
                    // good to go
                }
                Some(excluded_groups) => {
                    if excluded_groups.contains(&variant_group) {
                        to_remove.insert(*node);
                        continue;
                    }
                }
            }
            let node_cumulative_depth = nodes_cumulative_depth.entry(*node).or_insert(0.0); // This nodes current capacity
            let threshold = self.grouped_mean_read_depth.get(&variant_group).unwrap(); // This nodes maximum capcity
            let updated_depth = *node_cumulative_depth + depth_being_added_to_nodes;

            if (*node_cumulative_depth - *threshold).abs() <= f64::EPSILON
                || (&updated_depth > threshold && *node_cumulative_depth > 0.0)
            {
                // more efficient than division
                debug!(
                    "Node at capacity {:?} vg {} cumulative depth {} capacity {}",
                    node, variant_group, *node_cumulative_depth, *threshold
                );
                match &mut nodes_at_capacity {
                    None => nodes_at_capacity = Some(vec![*node]),
                    Some(nodes_at_capacity) => nodes_at_capacity.push(*node),
                }
            }
        }

        to_remove.into_iter().for_each(|remove_node| {
            path.remove(&remove_node);
        });

        nodes_at_capacity
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
        flag_filters: &FlagFilter,
        n_threads: usize,
    ) -> LinkedHashMap<i32, HashSet<String>> {
        let mut all_grouped_reads = LinkedHashMap::with_capacity(self.grouped_contexts.len());
        let mut all_grouped_read_counts = LinkedHashMap::with_capacity(self.grouped_contexts.len());

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

                // bam_generated.set_threads(n_threads);

                let mut grouped_reads = LinkedHashMap::with_capacity(self.grouped_contexts.len());
                let mut grouped_read_counts =
                    LinkedHashMap::with_capacity(self.grouped_contexts.len());
                let mut record = Record::new();
                for (group, variants) in self.grouped_contexts.iter() {
                    for variant in variants {
                        bam_generated
                            .fetch((
                                variant.loc.tid as i32,
                                variant.loc.start as i64,
                                variant.loc.end as i64 + 1,
                            ))
                            .unwrap_or_else(|_| {
                                panic!(
                                    "Failed to fetch interval {}:{}-{}",
                                    variant.loc.tid, variant.loc.start, variant.loc.end
                                )
                            });

                        let records = grouped_reads.entry(*group).or_insert(HashSet::new()); // container for the records to be collected
                        let counts = grouped_read_counts.entry(*group).or_insert(0.0);
                        let allele_depth = variant.genotypes.genotypes()[sample_idx].ad[1] as f64;
                        let mut read_count = 0.0;
                        while bam_generated.read(&mut record) == true {
                            // be very lenient with filtering
                            if record.is_unmapped() || record.seq_len() == 0 {
                                continue;
                            }

                            // TODO: Filter for only reads that contain the variant in question
                            //       not sure how to this besides aligning or other expensive method?
                            //       currently just check against reference
                            // let mut read_index = record.pos() - variant.loc.start as i64;
                            let variant_start = variant.loc.start as i64;
                            let mut partial_match = false;
                            let record_seq = record.seq().as_bytes();

                            let mut read_index = variant_start - record.pos();
                            let mut stop = 0;
                            if read_index < 0 {
                                partial_match = true;
                                read_index = 0;
                            } else if read_index as usize >= record_seq.len() {
                                read_index = record_seq.len() as i64 - 1;
                                partial_match = true;
                            }

                            let alternate_allele = variant.get_alternate_alleles()[0];
                            if read_index as usize + alternate_allele.get_bases().len()
                                <= record_seq.len()
                                && !partial_match
                            {
                                if alternate_allele.get_bases()
                                    == &record_seq[read_index as usize
                                        ..(read_index as usize
                                            + alternate_allele.get_bases().len())]
                                {
                                    // Read containing potential alternate allele
                                    let read_id = format!(
                                        "{}_{}",
                                        sample_idx,
                                        std::str::from_utf8(record.qname()).unwrap()
                                    );
                                    records.insert(read_id);
                                    read_count += 1.0;
                                }
                            } else if partial_match {
                                // substring match
                                let record_bases = &record_seq[read_index as usize
                                    ..min(
                                        record_seq.len(),
                                        read_index as usize + alternate_allele.bases.len(),
                                    )];
                                if alternate_allele.get_bases().contains_str(record_bases) {
                                    // Read containing potential alternate allele
                                    let read_id = format!(
                                        "{}_{}",
                                        sample_idx,
                                        std::str::from_utf8(record.qname()).unwrap()
                                    );
                                    records.insert(read_id);
                                    read_count += 1.0;
                                }
                            }
                        }

                        if read_count > allele_depth {
                            *counts += read_count;
                        } else {
                            *counts += allele_depth;
                        }
                    }
                }

                (grouped_reads, grouped_read_counts)
            })
            .collect::<Vec<(LinkedHashMap<i32, HashSet<String>>, LinkedHashMap<i32, f64>)>>()
            .into_iter()
            .for_each(|(sample_grouping, sample_counts)| {
                for (vg, reads) in sample_grouping {
                    let all_result = all_grouped_reads.entry(vg).or_insert(HashSet::new());
                    all_result.par_extend(reads);

                    let all_counts = all_grouped_read_counts.entry(vg).or_insert(0.0);
                    *all_counts += *sample_counts.get(&vg).unwrap();
                }
            });

        let grouped_mean_read_depth = all_grouped_read_counts
            .iter()
            .map(|(group, counts)| {
                (
                    *group,
                    counts / self.grouped_contexts.get(group).unwrap().len() as f64,
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
            if *group1 < 0 {
                continue;
            }
            let node1 = if node_indices.contains_key(group1) {
                *node_indices.get(group1).unwrap()
            } else {
                let node = graph.add_node(*group1);
                node_indices.insert(*group1, node);
                node
            };

            // Use the previous grouping values for the purpose of linking by distance
            let mut ind1 = *group1 as usize;
            if self.previous_groups.contains_key(group1) {
                ind1 = *self.previous_groups.get(group1).unwrap() as usize
            }

            for (group2, reads2) in grouped_reads.iter() {
                if group1 == group2 || *group2 < 0 || self.check_exclusion(group1, group2) {
                    continue;
                }

                let node2 = if node_indices.contains_key(group2) {
                    *node_indices.get(group2).unwrap()
                } else {
                    let node = graph.add_node(*group2);
                    node_indices.insert(*group2, node);
                    node
                };

                let mut ind2 = *group2 as usize;
                if self.previous_groups.contains_key(group2) {
                    ind2 = *self.previous_groups.get(group2).unwrap() as usize
                }

                if ind1 == ind2 {
                    // don't form edges between identical groups/indices
                    // the previous group of one of the current groups is the same as other group
                    continue;
                }

                // Don't count twice
                if !graph.contains_edge(node1, node2) {
                    // How many read ids are shared
                    let intersection = reads1.intersection(reads2).count() as f64;

                    let mut under_sep_thresh = false;
                    under_sep_thresh = self.cluster_separations[[ind1, ind2]] < 2.5;
                    if intersection > 0.0 || under_sep_thresh {
                        let union = reads1.union(reads2).count() as f64;

                        // The weight needs to be low for highly connected nodes and
                        // high for poorly connected nodes. This is because minimum spanning trees
                        // generate the tree with lowest edge weights
                        let mut weight = 1.0 - (intersection / union);

                        // We will form an edge here but it will essentially
                        // have to be directly on top of the other variant group
                        // for the edge to be favoured in the MST
                        let depth_1 = *self.grouped_mean_read_depth.get(group1).unwrap();
                        let depth_2 = *self.grouped_mean_read_depth.get(group2).unwrap();

                        let depth_factor = 1.0 - min(OrderedFloat(depth_1), OrderedFloat(depth_2)).0.ln() / max(OrderedFloat(depth_1), OrderedFloat(depth_2)).0.ln();

                        if weight < 1.00 || intersection >= 50.0 {
                            // variant groups connected by reads are favoured

                            debug!(
                                "{}:{} weight {} intersection {} union {}",
                                group1, group2, weight, intersection, union
                            );
                            weight = weight + weight * depth_factor;
                            debug!("updated weight {}", weight);
                            graph.add_edge(node1, node2, weight);
                        } else if under_sep_thresh {

                            let mut weight = self.cluster_separations[[ind1, ind2]];
                            debug!(
                                "{}:{} weight {} intersection {} union {}",
                                group1, group2, weight, intersection, union
                            );
                            weight = weight + weight * depth_factor;
                            debug!("updated weight {}", weight);
                            graph.add_edge(node1, node2, weight);
                        }
                    }
                }
            }
        }

        graph
    }

    fn check_exclusion(&self, group1: &i32, group2: &i32) -> bool {
        // check that neither group1 or group2 is exclusive of each other
        if self.exclusive_groups.contains_key(group1) {
            if self.exclusive_groups.get(group1).unwrap().contains(group2) {
                return true;
            }
        }

        if self.exclusive_groups.contains_key(group2) {
            if self.exclusive_groups.get(group2).unwrap().contains(group1) {
                return true;
            }
        }

        false
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

    fn max_depth_of_graph(&self, component_graph: &UnGraph<i32, f64>) -> f64 {
        component_graph.node_weights().into_iter().map(|vg| {
            OrderedFloat(*self.grouped_mean_read_depth.get(vg).unwrap())
        }).max().unwrap().0
    }

    fn min_depth_of_graph(&self, component_graph: &UnGraph<i32, f64>) -> f64 {
        component_graph.node_weights().into_iter().map(|vg| {
            OrderedFloat(*self.grouped_mean_read_depth.get(vg).unwrap())
        }).min().unwrap().0
    }
}
