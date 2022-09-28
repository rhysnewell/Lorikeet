#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate coverm;
extern crate lorikeet_genome;
extern crate rayon;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;
extern crate bio;
extern crate itertools;
extern crate petgraph;
extern crate rand;
extern crate term;

use lorikeet_genome::graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::graphs::common_suffix_splitter::CommonSuffixSplitter;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::graphs::shared_sequence_merger::SharedSequenceMerger;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::utils::simple_interval::SimpleInterval;
use petgraph::prelude::NodeIndex;
use std::cmp::Ordering::Equal;

pub struct SplitMergeData {
    graph: SeqGraph<BaseEdgeStruct>,
    v: NodeIndex,
    common_suffix: String,
}

impl SplitMergeData {
    pub fn new(graph: SeqGraph<BaseEdgeStruct>, v: NodeIndex, common_suffix: String) -> Self {
        Self {
            graph,
            v,
            common_suffix,
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "SplitMergeData [ graph = {:?}, v = {:?}, common_suffix = {} ]",
            self.graph, self.v, self.common_suffix
        )
    }
}

pub fn make_split_merge_data(test_function: Box<dyn Fn(SplitMergeData)>) {
    let bases = vec!["A", "C", "G", "T"];
    for common_suffix in vec!["", "A", "AT"] {
        for n_bots in vec![0, 1, 2] {
            for n_mids in vec![1, 2, 3] {
                for n_tops in 0..n_mids {
                    for n_top_connections in 1..=n_mids {
                        let mut multi = 1;
                        let mut graph = SeqGraph::new(11);
                        let v = SeqVertex::new(b"GGGG".to_vec());
                        let v_i = graph.base_graph.add_node(&v);

                        let mut tops = Vec::new();
                        let mut mids = Vec::new();

                        for i in 0..n_mids {
                            let mid = SeqVertex::new(
                                format!("{}{}", bases[i], common_suffix).as_bytes().to_vec(),
                            );
                            let mid_i = graph.base_graph.add_node(&mid);
                            graph.base_graph.add_edges(
                                mid_i,
                                vec![v_i],
                                BaseEdgeStruct::new(i == 0, multi, 0),
                            );
                            multi += 1;

                            mids.push(mid_i);

                            tops.push(SeqVertex::new(bases[i].as_bytes().to_vec()));
                        }

                        let top_indices = graph.base_graph.add_vertices(tops.iter());
                        for t_i in top_indices {
                            for i in 0..n_top_connections {
                                graph.base_graph.add_or_update_edge(
                                    t_i,
                                    mids[i],
                                    BaseEdgeStruct::new(i == 0, multi, 0),
                                );
                                multi += 1;
                            }
                        }

                        for i in 0..n_bots {
                            let bot = SeqVertex::new(bases[i].as_bytes().to_vec());
                            let bot_i = graph.base_graph.add_node(&bot);
                            graph.base_graph.add_or_update_edge(
                                v_i,
                                bot_i,
                                BaseEdgeStruct::new(i == 0, multi, 0),
                            );
                            multi += 1;
                        }

                        test_function(SplitMergeData::new(graph, v_i, common_suffix.to_string()));
                    }
                }
            }
        }
    }
}

#[test]
fn test_split() {
    make_split_merge_data(test_split_function())
}

fn test_split_function() -> Box<dyn Fn(SplitMergeData)> {
    Box::new(|mut data: SplitMergeData| {
        let expected_merge =
            !data.common_suffix.is_empty() && data.graph.base_graph.in_degree_of(data.v) > 1;

        let mut original = data.graph.clone();
        let succeed = CommonSuffixSplitter::split(&mut data.graph, data.v);
        assert_eq!(succeed, expected_merge, "Not expected result");
        if succeed {
            assert_eq!(
                data.graph
                    .base_graph
                    .graph
                    .node_weight(
                        data.graph
                            .base_graph
                            .incoming_vertices_of(data.v)
                            .into_iter()
                            .next()
                            .unwrap()
                    )
                    .unwrap()
                    .get_sequence_string(),
                &data.common_suffix,
                // (data.graph.base_graph.incoming_vertices_of(data.v).into_iter().next().unwrap())
            )
        }

        assert_same_haplotypes(
            format!(
                "suffixSplit.{}.{}",
                data.common_suffix,
                data.graph.base_graph.vertex_set().len()
            ),
            data.graph,
            original,
        )
    })
}

fn assert_same_haplotypes(
    name: String,
    mut actual: SeqGraph<BaseEdgeStruct>,
    mut original: SeqGraph<BaseEdgeStruct>,
) {
    let o_sources = original
        .base_graph
        .get_sources_generic()
        .map(|v| v)
        .collect();
    let o_sinks = original.base_graph.get_sinks_generic().map(|v| v).collect();
    let mut original_k_best_haplotypes =
        GraphBasedKBestHaplotypeFinder::new(&mut original.base_graph, o_sources, o_sinks)
            .find_best_haplotypes(std::usize::MAX, &original.base_graph);
    original_k_best_haplotypes.sort_by(|o1, o2| {
        let base_cmp = o1
            .path
            .get_bases(&original.base_graph)
            .cmp(&o2.path.get_bases(&original.base_graph));

        if base_cmp != Equal {
            base_cmp
        } else {
            o1.score.partial_cmp(&o2.score).unwrap()
        }
    });
    let sorted_original_k_best_haplotypes = original_k_best_haplotypes
        .into_iter()
        .map(|k| k.haplotype(&original.base_graph))
        .collect::<Vec<Haplotype<SimpleInterval>>>();

    let a_sources = actual.base_graph.get_sources_generic().map(|v| v).collect();
    let a_sinks = actual.base_graph.get_sinks_generic().map(|v| v).collect();
    let mut actual_k_best_haplotypes =
        GraphBasedKBestHaplotypeFinder::new(&mut actual.base_graph, a_sources, a_sinks)
            .find_best_haplotypes(std::usize::MAX, &actual.base_graph);

    actual_k_best_haplotypes.sort_by(|o1, o2| {
        let base_cmp = o1
            .path
            .get_bases(&actual.base_graph)
            .cmp(&o2.path.get_bases(&actual.base_graph));

        if base_cmp != Equal {
            base_cmp
        } else {
            o1.score.partial_cmp(&o2.score).unwrap()
        }
    });
    let sorted_actual_k_best_haplotypes = actual_k_best_haplotypes
        .into_iter()
        .map(|k| k.haplotype(&actual.base_graph))
        .collect::<Vec<Haplotype<SimpleInterval>>>();

    assert_eq!(
        sorted_actual_k_best_haplotypes,
        sorted_original_k_best_haplotypes,
    );
}

#[test]
fn test_split_prev_have_multiple_edges() {
    let mut original = SeqGraph::new(11);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"A".to_vec());
    let v3 = SeqVertex::new(b"A".to_vec());
    let v4 = SeqVertex::new(b"A".to_vec());

    let node_indices = original.base_graph.add_vertices(vec![&v1, &v2, &v3, &v4]);
    original.base_graph.add_or_update_edge(
        node_indices[0],
        node_indices[2],
        BaseEdgeStruct::new(false, 1, 0),
    );

    assert!(
        !CommonSuffixSplitter::split(&mut original, node_indices[2]),
        "Cannot split graph with only one vertex"
    );

    original.base_graph.add_or_update_edge(
        node_indices[1],
        node_indices[2],
        BaseEdgeStruct::new(false, 1, 0),
    );
    original.base_graph.add_or_update_edge(
        node_indices[1],
        node_indices[3],
        BaseEdgeStruct::new(false, 1, 0),
    );

    assert!(
        !CommonSuffixSplitter::split(&mut original, node_indices[3]),
        "Cannot split graph with multiple outgoing edges from middle nodes"
    );
}

#[test]
fn test_split_no_cycles() {
    let mut original = SeqGraph::new(11);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"AC".to_vec());
    let v3 = SeqVertex::new(b"TC".to_vec());
    let v4 = SeqVertex::new(b"G".to_vec());

    let node_indices = original.base_graph.add_vertices(vec![&v1, &v2, &v3, &v4]);
    original.base_graph.add_edges(
        node_indices[0],
        vec![node_indices[2], node_indices[3]],
        BaseEdgeStruct::new(false, 1, 0),
    );
    original.base_graph.add_edges(
        node_indices[0],
        vec![node_indices[1], node_indices[3]],
        BaseEdgeStruct::new(false, 1, 0),
    );

    let mut clone = original.clone();
    assert!(
        CommonSuffixSplitter::split(&mut clone, node_indices[3]),
        "Should be able to split pre-cycle graph"
    );

    original.base_graph.add_or_update_edge(
        node_indices[3],
        node_indices[3],
        BaseEdgeStruct::new(false, 1, 0),
    );
    assert!(
        !CommonSuffixSplitter::split(&mut original, node_indices[3]),
        "Cannot split graph with a cycle of the bottom list"
    );
}

#[test]
fn test_split_complex_cycle() {
    let mut original = SeqGraph::new(11);
    let r1 = SeqVertex::new(b"ACTG".to_vec());
    let r2 = SeqVertex::new(b"ATGC".to_vec());
    let cat1 = SeqVertex::new(b"CAT".to_vec());
    let cat2 = SeqVertex::new(b"CAT".to_vec());
    let c1 = SeqVertex::new(b"C".to_vec());
    let c2 = SeqVertex::new(b"C".to_vec());

    let node_indices = original
        .base_graph
        .add_vertices(vec![&r1, &r2, &cat1, &cat2, &c1, &c2]);
    original.base_graph.add_edges(
        node_indices[0],
        vec![
            node_indices[2],
            node_indices[4],
            node_indices[3],
            node_indices[4],
        ],
        BaseEdgeStruct::new(false, 1, 0),
    );

    original.base_graph.add_edges(
        node_indices[1],
        vec![node_indices[5], node_indices[3]],
        BaseEdgeStruct::new(false, 1, 0),
    );

    for v in node_indices {
        let mut graph = original.clone();
        let success = CommonSuffixSplitter::split(&mut graph, v);
        assert!(!success);
    }
}


#[test]
fn test_split_infinite_cycle_failure() {
    let mut original = SeqGraph::new(11);
    let v1 = SeqVertex::new(b"GC".to_vec());
    let v2 = SeqVertex::new(b"X".to_vec());
    let v3 = SeqVertex::new(b"N".to_vec());
    let v4 = SeqVertex::new(b"C".to_vec());

    let node_indices = original
        .base_graph
        .add_vertices(vec![&v1, &v2, &v3, &v4]);
    original.base_graph.add_edges(
        node_indices[0],
        vec![
            node_indices[1]
        ],
        BaseEdgeStruct::new(false, 12, 0),
    );

    original.base_graph.add_edges(
        node_indices[1],
        vec![
            node_indices[2]
        ],
        BaseEdgeStruct::new(false, 23, 0),
    );

    original.base_graph.add_edges(
        node_indices[2],
        vec![
            node_indices[3]
        ],
        BaseEdgeStruct::new(false, 34, 0),
    );

    original.base_graph.add_edges(
        node_indices[3],
        vec![
            node_indices[1]
        ],
        BaseEdgeStruct::new(false, 42, 0),
    );

    let mut graph = original.clone();
    let success = CommonSuffixSplitter::split(&mut graph, node_indices[1]);
    assert!(success);

    let vertex_set = graph.base_graph.graph.node_indices().map(|node| node.clone()).collect::<Vec<NodeIndex>>();
    for v in vertex_set {
        let success = CommonSuffixSplitter::split(&mut graph, v);
        assert!(!success);
    }
}