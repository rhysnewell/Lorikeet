#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

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
#[macro_use]
extern crate ntest;

use lorikeet_genome::graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use lorikeet_genome::graphs::base_graph::{BaseGraph, TestGraph};
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::graphs::vertex_based_transformer::VertexBasedTransformer;
use lorikeet_genome::graphs::vertex_based_transformer::VertexBasedTransformer::MergeTails;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use rayon::prelude::*;
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};

struct MergeNodesWithNoVariationTestProvider<'a> {
    sequence: &'a [u8],
    kmer_length: usize,
}

impl<'a> MergeNodesWithNoVariationTestProvider<'a> {
    pub fn new(seq: &'a [u8], kmer: usize) -> MergeNodesWithNoVariationTestProvider<'a> {
        Self {
            sequence: seq,
            kmer_length: kmer,
        }
    }

    pub fn calc_graph(&self) -> SeqGraph<BaseEdgeStruct> {
        let mut debruijn_graph = TestGraph::new(11);
        let kmers_in_sequence = self.sequence.len() - self.kmer_length + 1;
        for i in 0..kmers_in_sequence - 1 {
            // get the kmers
            let kmer_1 = &self.sequence[i..i + self.kmer_length];
            let kmer_2 = &self.sequence[i + 1..i + 1 + self.kmer_length];
            println!(
                "{} -> {} adding 1",
                std::str::from_utf8(kmer_1).unwrap(),
                std::str::from_utf8(kmer_2).unwrap()
            );
            debruijn_graph.add_kmers_to_graph(kmer_1, kmer_2, false, 1);
        }
        let mut seq_graph = debruijn_graph.graph.to_sequence_graph();
        seq_graph.simplify_graph("anon");
        return seq_graph;
    }
}

fn test_merge_nodes_with_no_variation(cfg: MergeNodesWithNoVariationTestProvider) {
    let actual = cfg.calc_graph();
    assert_eq!(
        actual.base_graph.graph.node_count(),
        1,
        "graph {:?}",
        &actual
    );
    let actual_v = actual.base_graph.graph.node_indices().next().unwrap();
    assert_eq!(
        actual
            .base_graph
            .graph
            .node_weight(actual_v)
            .unwrap()
            .sequence
            .as_slice(),
        cfg.sequence
    );
}

#[test]
fn make_merge_nodes_with_no_variation_test() {
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"GGTTAACC", 3));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"GGTTAACC", 4));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"GGTTAACC", 5));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"GGTTAACC", 6));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"GGTTAACC", 7));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(
        b"GGTTAACCATGCAGACGGGAGGCTGAGCGAGAGTTTT",
        6,
    ));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 66));
    test_merge_nodes_with_no_variation(MergeNodesWithNoVariationTestProvider::new(b"AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 76));
}

fn test_linear_zip(graph: &SeqGraph<BaseEdgeStruct>, expected: &SeqGraph<BaseEdgeStruct>) {
    let mut merged = graph.clone();
    merged.zip_linear_chains();
    if !merged.base_graph.graph_equals(&expected.base_graph) {
        graph.base_graph.print_graph("graph.dot", true, 0);
        merged.base_graph.print_graph("merged.dot", true, 0);
        expected.base_graph.print_graph("expected.dot", true, 0);
        assert!(false, "Graphs were not equal, check dot files");
    } else {
        assert!(merged.base_graph.graph_equals(&expected.base_graph))
    }
}

#[test]
fn make_linear_zip_data() {
    let mut graph = SeqGraph::new(11);
    let mut expected = SeqGraph::new(11);

    // empty graph => empty graph
    test_linear_zip(&graph, &expected);

    let a1 = SeqVertex::new(b"A".to_vec());
    let c1 = SeqVertex::new(b"C".to_vec());
    let ac1 = SeqVertex::new(b"AC".to_vec());

    // just a single vertex
    let a1_c1 = graph.base_graph.add_vertices(vec![a1.clone(), c1.clone()]);
    expected
        .base_graph
        .add_vertices(vec![a1.clone(), c1.clone()]);
    test_linear_zip(&graph, &expected);

    graph
        .base_graph
        .add_edges(a1_c1[0], vec![a1_c1[1]], BaseEdgeStruct::new(false, 1, 0));
    expected = SeqGraph::new(11);
    expected.base_graph.add_node(ac1.clone());
    test_linear_zip(&graph, &expected);

    // three long chain merged corrected
    let g1 = SeqVertex::new(b"G".to_vec());
    let g1_id = graph.base_graph.add_node(g1.clone());
    graph
        .base_graph
        .graph
        .add_edge(a1_c1[1], g1_id, BaseEdgeStruct::new(false, 1, 0));
    expected = SeqGraph::new(11);
    expected
        .base_graph
        .add_node(SeqVertex::new(b"ACG".to_vec()));
    test_linear_zip(&graph, &expected);

    // adding something that isn't connected isn't a problem
    let t1 = SeqVertex::new(b"T".to_vec());
    graph.base_graph.add_node(t1.clone());
    expected = SeqGraph::new(11);
    expected
        .base_graph
        .add_vertices(vec![SeqVertex::new(b"ACG".to_vec()), t1.clone()]);
    test_linear_zip(&graph, &expected);

    // splitting chain with branch produces the correct zipped subgraphs
    let a2 = SeqVertex::new(b"A".to_vec());
    let c2 = SeqVertex::new(b"C".to_vec());
    graph = SeqGraph::new(11);
    let node_indices = graph.base_graph.add_vertices(vec![
        a1.clone(),
        c1.clone(),
        g1.clone(),
        t1.clone(),
        a2.clone(),
        c2.clone(),
    ]);
    println!(
        "Adding indices {:?}",
        node_indices[1..(node_indices.len() - 1)].to_vec()
    );
    graph.base_graph.add_edges(
        node_indices[0],
        node_indices[1..(node_indices.len() - 1)].to_vec(),
        BaseEdgeStruct::new(false, 1, 0),
    );
    graph.base_graph.graph.add_edge(
        node_indices[2],
        node_indices[5],
        BaseEdgeStruct::new(false, 1, 0),
    );
    expected = SeqGraph::new(11);
    let acg = SeqVertex::new(b"ACG".to_vec());
    let ta = SeqVertex::new(b"TA".to_vec());
    let expected_indices =
        expected
            .base_graph
            .add_vertices(vec![acg.clone(), ta.clone(), c2.clone()]);
    expected.base_graph.graph.add_edge(
        expected_indices[0],
        expected_indices[1],
        BaseEdgeStruct::new(false, 1, 0),
    );
    expected.base_graph.graph.add_edge(
        expected_indices[0],
        expected_indices[2],
        BaseEdgeStruct::new(false, 1, 0),
    );
    test_linear_zip(&graph, &expected);

    // can merge chains with loops in them
    {
        graph = SeqGraph::new(11);
        let node_indices = graph
            .base_graph
            .add_vertices(vec![a1.clone(), c1.clone(), g1.clone()]);
        graph.base_graph.add_edges(
            node_indices[0],
            node_indices[1..].to_vec(),
            BaseEdgeStruct::new(false, 1, 0),
        );
        let a1_a1_index = graph.base_graph.graph.add_edge(
            node_indices[0],
            node_indices[0],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected = SeqGraph::new(11);

        let ac = SeqVertex::new(b"AC".to_vec());
        let cg = SeqVertex::new(b"CG".to_vec());
        let expected_indices = expected
            .base_graph
            .add_vertices(vec![a1.clone(), cg.clone()]);
        expected.base_graph.graph.add_edge(
            expected_indices[0],
            expected_indices[1],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected.base_graph.graph.add_edge(
            expected_indices[0],
            expected_indices[0],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_linear_zip(&graph, &expected);

        graph.base_graph.graph.remove_edge(a1_a1_index);
        let c1_c1_index = graph.base_graph.graph.add_edge(
            node_indices[1],
            node_indices[1],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_linear_zip(&graph, &graph);

        graph.base_graph.graph.remove_edge(c1_c1_index);
        let g1_g1_index = graph.base_graph.graph.add_edge(
            node_indices[2],
            node_indices[2],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected = SeqGraph::new(11);
        let expected_indices = expected
            .base_graph
            .add_vertices(vec![ac.clone(), g1.clone()]);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![expected_indices[1], expected_indices[1]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_linear_zip(&graph, &expected)
    }

    {
        // check build n element long chains
        let bases = vec!["A", "C", "G", "T", "TT", "GG", "CC", "AA"];
        for len in vec![1, 2, 10, 100, 1000] {
            println!("Length {}", len);
            graph = SeqGraph::new(11);
            expected = SeqGraph::new(11);
            let mut last_index = None;
            let mut expected_bases = Vec::new();
            for i in 0..len {
                let seq = bases[i % bases.len()];
                println!("Adding seq {} to {:?}", &seq, &expected_bases);
                expected_bases.extend_from_slice(seq.as_bytes());
                println!("Result {:?}", &expected_bases);
                let a = SeqVertex::new(seq.as_bytes().to_vec());
                let a_index = graph.base_graph.add_node(a);
                match last_index {
                    None => {
                        // pass
                    }
                    Some(last_index) => {
                        graph.base_graph.graph.add_edge(
                            last_index,
                            a_index,
                            BaseEdgeStruct::new(false, 1, 0),
                        );
                    }
                };
                last_index = Some(a_index);
            }
            expected
                .base_graph
                .add_node(SeqVertex::new(expected_bases.clone()));
            test_linear_zip(&graph, &expected);
        }
    }

    // check that edge connections are properly maintained
    {
        let mut edge_weight = 1;
        for n_incoming in vec![0, 2, 5, 10] {
            for n_outgoing in vec![0, 2, 5, 10] {
                graph = SeqGraph::new(11);
                expected = SeqGraph::new(11);
                let node_indices =
                    graph
                        .base_graph
                        .add_vertices(vec![a1.clone(), c1.clone(), g1.clone()]);
                graph.base_graph.add_edges(
                    node_indices[0],
                    node_indices[1..].to_vec(),
                    BaseEdgeStruct::new(false, 1, 0),
                );
                let expected_index = expected.base_graph.add_node(acg.clone());

                for v in make_vertices(n_incoming) {
                    let e = BaseEdgeStruct::new(false, edge_weight, 0);
                    edge_weight += 1;
                    let v_index = graph.base_graph.add_node(v.clone());
                    graph
                        .base_graph
                        .graph
                        .add_edge(v_index, node_indices[0], e.clone());
                    let expected_v_index = expected.base_graph.add_node(v);
                    expected
                        .base_graph
                        .graph
                        .add_edge(expected_v_index, expected_index, e);
                }

                for v in make_vertices(n_outgoing) {
                    let e = BaseEdgeStruct::new(false, edge_weight, 0);
                    edge_weight += 1;
                    let v_index = graph.base_graph.add_node(v.clone());
                    graph
                        .base_graph
                        .graph
                        .add_edge(node_indices[2], v_index, e.clone());
                    let expected_v_index = expected.base_graph.add_node(v);
                    expected
                        .base_graph
                        .graph
                        .add_edge(expected_index, expected_v_index, e);
                }

                test_linear_zip(&graph, &expected);
            }
        }
    }
}

fn make_vertices(n: usize) -> Vec<SeqVertex> {
    let mut vs = Vec::new();
    let mut bases = vec!["A", "C", "G", "T", "TT", "GG", "CC", "AA"];

    for i in 0..n {
        vs.push(SeqVertex::new(bases[i % bases.len()].as_bytes().to_vec()));
    }

    return vs;
}

fn test_merging<E: BaseEdge>(graph: &SeqGraph<E>, expected: &SeqGraph<E>) {
    let mut merged = graph.clone();
    merged.simplify_graph_with_cycles(1, "anon");
    if !merged.base_graph.graph_equals(&expected.base_graph) {
        graph.base_graph.print_graph("graph.dot", true, 0);
        merged.base_graph.print_graph("merged.dot", true, 0);
        expected.base_graph.print_graph("expected.dot", true, 0);
        assert!(false, "Graphs were not equal, check dot files");
    } else {
        assert!(merged.base_graph.graph_equals(&expected.base_graph))
    }
}

#[test]
fn make_merging_data() {
    let mut graph: SeqGraph<BaseEdgeStruct> = SeqGraph::new(11);
    let poly_a = vec![b'A'; VertexBasedTransformer::<BaseEdgeStruct>::MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES];

    let pre1 = SeqVertex::new(vec![poly_a.clone(), b"CT".to_vec()].concat());
    let pre2 = SeqVertex::new(vec![poly_a, b"GT".to_vec()].concat());
    let top = SeqVertex::new(b"A".to_vec());
    let middle1 = SeqVertex::new(b"GC".to_vec());
    let middle2 = SeqVertex::new(b"TC".to_vec());
    let middle3 = SeqVertex::new(b"AC".to_vec());
    let middle4 = SeqVertex::new(b"GCAC".to_vec());
    let bottom = SeqVertex::new(b"AA".to_vec());
    let tail1 = SeqVertex::new(b"GC".to_vec());
    let tail2 = SeqVertex::new(b"GC".to_vec());

    // just a single vertex
    let pre1_index = graph.base_graph.add_node(pre1.clone());
    test_merging(&graph, &graph);

    let mut top_index;
    // pre1 -> top = pre1 + top
    {
        top_index = graph.base_graph.add_node(top.clone());
        graph.base_graph.add_or_update_edge(
            pre1_index,
            top_index,
            BaseEdgeStruct::new(false, 1, 0),
        );
        let pre1_top = SeqVertex::new(vec![pre1.sequence.clone(), top.sequence.clone()].concat());
        let mut expected = SeqGraph::new(11);
        expected.base_graph.add_node(pre1_top);
        test_merging(&graph, &expected);
    }

    let mut middle1_index;
    // pre1 -> top -> middle1 = pre1 + top + middle1
    {
        middle1_index = graph.base_graph.add_node(middle1.clone());
        graph.base_graph.add_or_update_edge(
            top_index,
            middle1_index,
            BaseEdgeStruct::new(false, 1, 0),
        );
        let mut expected = SeqGraph::new(11);
        let pre1_top_middle1 = SeqVertex::new(
            vec![
                pre1.sequence.clone(),
                top.sequence.clone(),
                middle1.sequence.clone(),
            ]
            .concat(),
        );
        expected.base_graph.add_node(pre1_top_middle1);
        test_merging(&graph, &expected);
    }

    let mut middle2_index;
    // pre1 -> top -> middle1 & top -> middle2 = pre1 + top -> middle1 & -> middle2
    {
        middle2_index = graph.base_graph.add_node(middle2.clone());
        graph.base_graph.add_or_update_edge(
            top_index,
            middle2_index,
            BaseEdgeStruct::new(false, 1, 0),
        );
        let pre1_top = SeqVertex::new(vec![pre1.sequence.clone(), top.sequence.clone()].concat());
        let mut expected = SeqGraph::new(11);
        let expected_indices =
            expected
                .base_graph
                .add_vertices(vec![pre1_top, middle1.clone(), middle2.clone()]);
        expected.base_graph.add_or_update_edge(
            expected_indices[0],
            expected_indices[1],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected.base_graph.add_or_update_edge(
            expected_indices[0],
            expected_indices[2],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&graph, &expected);
    }

    let mut bottom_index;
    let mut middle3_index;
    let mut middle4_index;
    // An actual diamond event to merge!
    {
        bottom_index = graph.base_graph.add_node(bottom.clone());
        graph.base_graph.add_or_update_edge(
            middle1_index,
            bottom_index,
            BaseEdgeStruct::new(false, 1, 0),
        );
        graph.base_graph.add_or_update_edge(
            middle2_index,
            bottom_index,
            BaseEdgeStruct::new(false, 1, 0),
        );
        let mut expected = SeqGraph::new(11);
        let pre1_top = SeqVertex::new(vec![pre1.sequence.clone(), top.sequence.clone()].concat());
        let new_middle1 = SeqVertex::new(b"G".to_vec());
        let new_middle2 = SeqVertex::new(b"T".to_vec());
        let new_bottom = SeqVertex::new(vec![b"C".to_vec(), bottom.sequence.clone()].concat());
        let expected_indices =
            expected
                .base_graph
                .add_vertices(vec![pre1_top, new_middle1, new_middle2, new_bottom]);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![expected_indices[1], expected_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![expected_indices[2], expected_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&graph, &expected);

        middle3_index = graph.base_graph.add_node(middle3.clone());
        graph.base_graph.add_edges(
            top_index,
            vec![middle3_index, bottom_index],
            BaseEdgeStruct::new(false, 1, 0),
        );
        let new_middle3 = SeqVertex::new(b"A".to_vec());
        let new_middle3_index = expected.base_graph.add_node(new_middle3);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![new_middle3_index, expected_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&graph, &expected);

        middle4_index = graph.base_graph.add_node(middle4.clone());
        graph.base_graph.add_edges(
            top_index,
            vec![middle4_index, bottom_index],
            BaseEdgeStruct::new(false, 1, 0),
        );
        let new_middle4 = SeqVertex::new(b"GCA".to_vec());
        let new_middle4_index = expected.base_graph.graph.add_node(new_middle4);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![new_middle4_index, expected_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&graph, &expected);
    }

    {
        // all the nodes -> lots of merging and motion of nodes
        let mut all = SeqGraph::new(11);
        let all_indices = all.base_graph.add_vertices(vec![
            pre1.clone(),
            pre2.clone(),
            top.clone(),
            middle1.clone(),
            middle2.clone(),
            bottom.clone(),
            tail1.clone(),
            tail2.clone(),
        ]);
        all.base_graph.add_edges(
            all_indices[0],
            vec![
                all_indices[2],
                all_indices[3],
                all_indices[5],
                all_indices[6],
            ],
            BaseEdgeStruct::new(false, 1, 0),
        );
        all.base_graph.add_edges(
            all_indices[1],
            vec![
                all_indices[2],
                all_indices[4],
                all_indices[5],
                all_indices[7],
            ],
            BaseEdgeStruct::new(false, 1, 0),
        );

        let mut expected = SeqGraph::new(11);
        let new_pre1 = SeqVertex::new(vec![vec![b'A'; VertexBasedTransformer::<BaseEdgeStruct>::MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES], b"C".to_vec()].concat());
        let new_pre2 = SeqVertex::new(vec![vec![b'A'; VertexBasedTransformer::<BaseEdgeStruct>::MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES], b"G".to_vec()].concat());
        let new_top = SeqVertex::new(b"TA".to_vec());
        let new_middle1 = SeqVertex::new(b"G".to_vec());
        let new_middle2 = SeqVertex::new(b"T".to_vec());
        let new_bottom = SeqVertex::new(vec![b"C".to_vec(), bottom.sequence.clone()].concat());
        let expected_indices = expected.base_graph.add_vertices(vec![
            new_pre1.clone(),
            new_pre2.clone(),
            new_top.clone(),
            new_middle1.clone(),
            new_middle2.clone(),
            new_bottom.clone(),
            tail1.clone(),
            tail2.clone(),
        ]);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![
                expected_indices[2],
                expected_indices[3],
                expected_indices[5],
                expected_indices[6],
            ],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected.base_graph.add_edges(
            expected_indices[1],
            vec![
                expected_indices[2],
                expected_indices[4],
                expected_indices[5],
                expected_indices[7],
            ],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&all, &expected);
    }

    // test the case where we delete a middle node away because the common sequence is all of its sequence
    {
        let mut graph2 = SeqGraph::new(11);
        let my_top = SeqVertex::new(b"A".to_vec());
        let mid1 = SeqVertex::new(b"AC".to_vec());
        let mid2 = SeqVertex::new(b"C".to_vec());
        let bot = SeqVertex::new(b"G".to_vec());
        let graph2_indices = graph2.base_graph.add_vertices(vec![
            my_top.clone(),
            mid1.clone(),
            mid2.clone(),
            bot.clone(),
        ]);
        graph2.base_graph.add_edges(
            graph2_indices[0],
            vec![graph2_indices[1], graph2_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        graph2.base_graph.add_edges(
            graph2_indices[0],
            vec![graph2_indices[2], graph2_indices[3]],
            BaseEdgeStruct::new(false, 1, 0),
        );

        let mut expected = SeqGraph::new(11);
        let new_mid1 = SeqVertex::new(b"A".to_vec());
        let new_bottom = SeqVertex::new(b"CG".to_vec());
        let expected_indices = expected.base_graph.add_vertices(vec![
            my_top.clone(),
            new_mid1.clone(),
            new_bottom.clone(),
        ]);
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![expected_indices[1], expected_indices[2]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        expected.base_graph.add_edges(
            expected_indices[0],
            vec![expected_indices[2]],
            BaseEdgeStruct::new(false, 1, 0),
        );
        test_merging(&graph2, &expected);
    }
}

// A -> ACT -> C [non-ref]
// A -> ACT -> C [non-ref]
// A -> ACT -> C [ref]
//
// Should become A -> ACT -> C [ref and non-ref edges]
//
// #[test]
fn test_bubble_same_bases_with_ref() {
    let mut graph = SeqGraph::new(11);
    let top = SeqVertex::new(b"A".to_vec());
    let mid1 = SeqVertex::new(b"ACT".to_vec());
    let mid2 = SeqVertex::new(b"ACT".to_vec());
    let bot = SeqVertex::new(b"C".to_vec());
    let graph_indices =
        graph
            .base_graph
            .add_vertices(vec![top.clone(), mid1.clone(), mid2.clone(), bot.clone()]);
    graph.base_graph.add_edges(
        graph_indices[0],
        vec![graph_indices[2], graph_indices[3]],
        BaseEdgeStruct::new(false, 1, 0),
    );
    graph.base_graph.graph.add_edge(
        graph_indices[0],
        graph_indices[1],
        BaseEdgeStruct::new(true, 1, 0),
    );
    graph.base_graph.graph.add_edge(
        graph_indices[1],
        graph_indices[3],
        BaseEdgeStruct::new(true, 1, 0),
    );

    let mut expected = SeqGraph::new(11);
    expected
        .base_graph
        .add_vertices(vec![SeqVertex::new(b"AACTC".to_vec())]);
    // let expected_indices = expected.base_graph.add_vertices(vec![top.clone(), SeqVertex::new(b"ACTC".to_vec())]);
    // expected.base_graph.graph.add_edge(graph_indices[0], graph_indices[1], BaseEdgeStruct::new(true, 1, 0));
    // expected.base_graph.graph.add_edge(graph_indices[1], graph_indices[3], BaseEdgeStruct::new(true, 1, 0));
    test_bubble_merging(&graph, &expected);
}

fn test_bubble_merging<E: BaseEdge>(graph: &SeqGraph<E>, expected: &SeqGraph<E>) {
    let mut merged = graph.clone();
    merged.simplify_graph("anon");
    if !merged.base_graph.graph_equals(&expected.base_graph) {
        graph.base_graph.print_graph("graph.dot", true, 0);
        merged.base_graph.print_graph("merged.dot", true, 0);
        expected.base_graph.print_graph("expected.dot", true, 0);
        assert!(false, "Graphs were not equal, check dot files");
    } else {
        graph.base_graph.print_graph("graph.dot", true, 0);
        merged.base_graph.print_graph("merged.dot", true, 0);
        expected.base_graph.print_graph("expected.dot", true, 0);
        assert!(merged.base_graph.graph_equals(&expected.base_graph))
    }
}
