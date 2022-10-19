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
extern crate hashlink;
extern crate itertools;
extern crate petgraph;
extern crate rand;
extern crate term;

use hashlink::LinkedHashSet;
use itertools::Itertools;
use lorikeet_genome::graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::graphs::graph_utils::GraphUtils;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::graphs::shared_sequence_merger::SharedSequenceMerger;
use lorikeet_genome::graphs::shared_vertex_sequence_splitter::SharedVertexSequenceSplitter;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::byte_array_allele::Allele;
use lorikeet_genome::utils::simple_interval::SimpleInterval;
use petgraph::prelude::NodeIndex;
use std::cmp::Ordering::Equal;
use std::collections::HashSet;

fn test_prefix_suffix(strings: Vec<&str>, expected_prefix_len: usize, expected_suffix_len: usize) {
    let mut bytes = Vec::new();
    let mut min = std::usize::MAX;
    for s in &strings {
        bytes.push(s.as_bytes());
        min = std::cmp::min(min, s.len());
    }

    let actual_prefix_len = GraphUtils::common_maximum_prefix_length(&bytes);
    assert_eq!(actual_prefix_len, expected_prefix_len, "Failed prefix test");

    let actual_suffix_len =
        GraphUtils::common_maximum_suffix_length(&bytes, min.saturating_sub(actual_prefix_len));
    assert_eq!(
        actual_suffix_len, expected_suffix_len,
        "Failed suffix test, min {} and prefix len {}",
        min, actual_prefix_len
    );
}

#[test]
fn make_prefix_suffix_data() {
    test_prefix_suffix(vec!["A", "C"], 0, 0);
    test_prefix_suffix(vec!["C", "C"], 1, 0);
    test_prefix_suffix(vec!["ACT", "AGT"], 1, 1);
    test_prefix_suffix(vec!["ACCT", "AGT"], 1, 1);
    test_prefix_suffix(vec!["ACT", "ACT"], 3, 0);
    test_prefix_suffix(vec!["ACTA", "ACT"], 3, 0);
    test_prefix_suffix(vec!["ACTA", "ACTG"], 3, 0);
    test_prefix_suffix(vec!["ACTA", "ACTGA"], 3, 1);
    test_prefix_suffix(vec!["GCTGA", "ACTGA"], 0, 4);
    test_prefix_suffix(vec!["A", "C", "A"], 0, 0);
    test_prefix_suffix(vec!["A", "A", "A"], 1, 0);
    test_prefix_suffix(vec!["A", "AA", "A"], 1, 0);
    test_prefix_suffix(vec!["A", "ACA", "A"], 1, 0);
    test_prefix_suffix(vec!["ACT", "ACAT", "ACT"], 2, 1);
    test_prefix_suffix(vec!["ACT", "ACAT", "ACGT"], 2, 1);
    test_prefix_suffix(vec!["AAAT", "AAA", "CAAA"], 0, 0);
    test_prefix_suffix(vec!["AACTTT", "AAGTTT", "AAGCTTT"], 2, 3);
    test_prefix_suffix(vec!["AAA", "AAA", "CAAA"], 0, 3);
    test_prefix_suffix(vec!["AAA", "AAA", "AAA"], 3, 0);
    test_prefix_suffix(vec!["AC", "ACA", "AC"], 2, 0);
}

fn test_splitter(strings: Vec<&str>, expected_prefix_len: usize, expected_suffix_len: usize) {
    let mut graph = SeqGraph::<BaseEdgeStruct>::new(11);

    let mut v = Vec::new();
    for s in &strings {
        v.push(SeqVertex::new(s.as_bytes().to_vec()));
    }

    let nodes = graph.base_graph.add_vertices(v.iter());

    let expected_prefix = &strings[0].as_bytes()[0..expected_prefix_len];
    let expected_suffix = &strings[0].as_bytes()[strings[0].len() - expected_suffix_len..];

    let to_split = nodes.into_iter().collect::<LinkedHashSet<NodeIndex>>();
    let mut splitter = SharedVertexSequenceSplitter::new(&mut graph, to_split);
    splitter.split();

    assert_eq!(splitter.get_prefix().get_sequence(), expected_prefix);
    assert_eq!(splitter.get_suffix().get_sequence(), expected_suffix);

    assert!(
        splitter
            .get_split_graph()
            .base_graph
            .out_degree_of(splitter.get_prefix_index())
            <= strings.len()
    );
    assert_eq!(
        splitter
            .get_split_graph()
            .base_graph
            .in_degree_of(splitter.get_prefix_index()),
        0
    );

    assert!(
        splitter
            .get_split_graph()
            .base_graph
            .in_degree_of(splitter.get_suffix_index())
            <= strings.len()
    );
    assert_eq!(
        splitter
            .get_split_graph()
            .base_graph
            .out_degree_of(splitter.get_suffix_index()),
        0
    );

    for mid in splitter.get_new_middles() {
        assert!(splitter
            .get_split_graph()
            .base_graph
            .graph
            .find_edge(splitter.get_prefix_index(), *mid)
            .is_some());
        assert!(splitter
            .get_split_graph()
            .base_graph
            .graph
            .find_edge(*mid, splitter.get_suffix_index())
            .is_some());
    }
}

#[test]
fn make_splitter_data() {
    test_splitter(vec!["A", "C"], 0, 0);
    test_splitter(vec!["C", "C"], 1, 0);
    test_splitter(vec!["ACT", "AGT"], 1, 1);
    test_splitter(vec!["ACCT", "AGT"], 1, 1);
    test_splitter(vec!["ACT", "ACT"], 3, 0);
    test_splitter(vec!["ACTA", "ACT"], 3, 0);
    test_splitter(vec!["ACTA", "ACTG"], 3, 0);
    test_splitter(vec!["ACTA", "ACTGA"], 3, 1);
    test_splitter(vec!["GCTGA", "ACTGA"], 0, 4);
    test_splitter(vec!["A", "C", "A"], 0, 0);
    test_splitter(vec!["A", "A", "A"], 1, 0);
    test_splitter(vec!["A", "AA", "A"], 1, 0);
    test_splitter(vec!["A", "ACA", "A"], 1, 0);
    test_splitter(vec!["ACT", "ACAT", "ACT"], 2, 1);
    test_splitter(vec!["ACT", "ACAT", "ACGT"], 2, 1);
    test_splitter(vec!["AAAT", "AAA", "CAAA"], 0, 0);
    test_splitter(vec!["AACTTT", "AAGTTT", "AAGCTTT"], 2, 3);
    test_splitter(vec!["AAA", "AAA", "CAAA"], 0, 3);
    test_splitter(vec!["AAA", "AAA", "AAA"], 3, 0);
    test_splitter(vec!["AC", "ACA", "AC"], 2, 0);
}

fn test_splitter_complete_cycle(strings: Vec<&str>, has_top: bool, has_bot: bool) {
    println!(
        "Strings {:?} has top {} has bot {}",
        &strings, has_top, has_bot
    );
    let mut graph = SeqGraph::new(11);

    let mut edge_weight = 1;
    let top = if has_top {
        Some(SeqVertex::new(b"AAAAAAAA".to_vec()))
    } else {
        None
    };

    let bot = if has_bot {
        Some(SeqVertex::new(b"GGGGGGGG".to_vec()))
    } else {
        None
    };

    let mut v = Vec::new();
    for s in &strings {
        v.push(SeqVertex::new(s.as_bytes().to_vec()));
    }
    println!("vertices {:?}", &v);
    let mut nodes = graph.base_graph.add_vertices(v.iter());
    let first = nodes[0];

    let mut top_index = None;
    if has_top {
        top_index = Some(graph.base_graph.add_node(top.as_ref().unwrap()));
        for vi in nodes.iter() {
            graph.base_graph.add_or_update_edge(
                *top_index.as_ref().unwrap(),
                *vi,
                BaseEdgeStruct::new(*vi == first, edge_weight, 0),
            );
            edge_weight += 1;
        }
    }

    let mut bot_index = None;
    if has_bot {
        bot_index = Some(graph.base_graph.add_node(bot.as_ref().unwrap()));
        for vi in nodes.iter() {
            graph.base_graph.add_or_update_edge(
                *vi,
                *bot_index.as_ref().unwrap(),
                BaseEdgeStruct::new(*vi == first, edge_weight, 0),
            );
            edge_weight += 1;
        }
    }

    println!("{has_top} {has_bot}");
    let sources = graph
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = graph
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    println!("{} {}", sources.len(), sinks.len());
    println!(
        "Initial: Nodes {} Edges {}",
        graph.base_graph.graph.node_count(),
        graph.base_graph.graph.edge_count()
    );

    let mut graph_clone = graph.clone();
    let original_paths =
        GraphBasedKBestHaplotypeFinder::new(&mut graph_clone.base_graph, sources, sinks)
            .find_best_haplotypes(std::usize::MAX, &graph_clone.base_graph);

    let haplotypes = original_paths
        .iter()
        .map(|p| p.haplotype(&graph_clone.base_graph))
        .collect::<HashSet<Haplotype<SimpleInterval>>>();
    println!("Haplotypes {}", haplotypes.len());
    let to_split = nodes.into_iter().collect();
    let mut splitter = SharedVertexSequenceSplitter::new(&mut graph, to_split);
    println!(
        "Prefix {} Suffix {}",
        std::str::from_utf8(splitter.get_prefix().sequence.as_slice()).unwrap(),
        std::str::from_utf8(splitter.get_suffix().sequence.as_slice()).unwrap()
    );
    splitter.split();
    println!("Middle {:?}", splitter.get_new_middles());
    println!("To split {:?}", &splitter.to_splits);
    println!(
        "Split: Nodes {} Edges {}",
        splitter.get_split_graph().base_graph.graph.node_count(),
        splitter.get_split_graph().base_graph.graph.edge_count()
    );
    let sources = splitter
        .get_split_graph()
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = splitter
        .get_split_graph()
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    println!("Split sources {} sinks {}", sources.len(), sinks.len());

    splitter.update_graph(top_index, bot_index);
    println!(
        "Update: Nodes {} Edges {}",
        graph.base_graph.graph.node_count(),
        graph.base_graph.graph.edge_count()
    );

    let sources = graph
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = graph
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();

    println!("{} {}", sources.len(), sinks.len());
    let split_paths = GraphBasedKBestHaplotypeFinder::new(&mut graph.base_graph, sources, sinks)
        .find_best_haplotypes(std::usize::MAX, &graph.base_graph);

    split_paths.iter().for_each(|p| {
        assert!(
            haplotypes.contains(&p.haplotype(&graph.base_graph)),
            "Haplotype {:?} not in {:?}",
            p.haplotype::<SimpleInterval, SeqVertex, BaseEdgeStruct>(&graph.base_graph),
            &haplotypes
        )
    });

    let mut sorted_original_paths: Vec<Haplotype<SimpleInterval>> = original_paths
        .iter()
        .map(|p| p.haplotype(&graph_clone.base_graph))
        .unique()
        .collect_vec();
    sorted_original_paths.sort_by(|h1, h2| h1.get_bases().cmp(&h2.get_bases()));

    let mut sorted_split_paths: Vec<Haplotype<SimpleInterval>> = split_paths
        .iter()
        .map(|p| p.haplotype(&graph.base_graph))
        .unique()
        .collect_vec();
    sorted_split_paths.sort_by(|h1, h2| h1.get_bases().cmp(&h2.get_bases()));

    assert_eq!(sorted_original_paths, sorted_split_paths)
}

#[test]
fn make_complete_cycle_data() {
    for has_top in vec![true, false] {
        for has_bot in vec![true, false] {
            if !has_top && !has_bot {
                continue;
            }

            test_splitter_complete_cycle(vec!["A", "A"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "C"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "AC"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "CA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["AC", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["AT", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["ATA", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["ATAA", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["ATAACA", "ACA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["CCCAAA", "AAA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["CCCAAAAAA", "AAA"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["CCCAAAAAA", "CCCAAA"], has_top, has_bot);

            test_splitter_complete_cycle(vec!["A", "A", "A"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "A", "C"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["A", "C", "C"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["AC", "C", "C"], has_top, has_bot);
            test_splitter_complete_cycle(vec!["CA", "C", "C"], has_top, has_bot);
            // all merged
            test_splitter_complete_cycle(vec!["AGA", "AGA", "AGA"], has_top, has_bot);
            // prefix and suffix
            test_splitter_complete_cycle(vec!["AGA", "AGA", "ACA"], has_top, has_bot);
            // 2 -> prefix, leave C
            test_splitter_complete_cycle(vec!["AGA", "AGA", "AGAC"], has_top, has_bot);
            // 2 -> prefix, leave CCC
            test_splitter_complete_cycle(vec!["AGA", "AGA", "AGACCC"], has_top, has_bot);
            // 2 -> suffix, leave A/T
            test_splitter_complete_cycle(vec!["TAGA", "TAGA", "AAGA"], has_top, has_bot);
            // 2 -> suffix, leave T, delete 1
            test_splitter_complete_cycle(vec!["TAGA", "TAGA", "AGA"], has_top, has_bot);
        }
    }
}
