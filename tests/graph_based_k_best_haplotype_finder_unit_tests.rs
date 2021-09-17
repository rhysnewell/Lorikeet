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

use lorikeet_genome::graphs::base_edge::BaseEdge;
use lorikeet_genome::graphs::base_edge::BaseEdgeStruct;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::model::byte_array_allele::Allele;
use lorikeet_genome::smith_waterman::bindings::SWParameters;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::utils::simple_interval::SimpleInterval;
use petgraph::prelude::NodeIndex;
use std::collections::HashSet;

lazy_static! {
    static ref PATH_TO_REFERENCE_SW_PARAMETERS: SWParameters = *NEW_SW_PARAMETERS;
}

#[test]
fn test_score() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"C".to_vec());
    let v3 = SeqVertex::new(b"G".to_vec());
    let v4 = SeqVertex::new(b"T".to_vec());
    let v5 = SeqVertex::new(b"A".to_vec());

    let nodes = g.base_graph.add_vertices(vec![v1, v2, v3, v4, v5]);
    g.base_graph
        .graph
        .add_edge(nodes[0], nodes[1], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[2], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[3], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[4], BaseEdgeStruct::new(false, 1, 0));

    let sources = g
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = g
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    let mut finder = GraphBasedKBestHaplotypeFinder::new(&mut g.base_graph, sources, sinks);
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 3);

    let haplotypes = finder.find_best_haplotypes(std::usize::MAX, &g.base_graph);
    let ACG = haplotypes
        .into_iter()
        .filter(|h| {
            h.haplotype::<SimpleInterval, SeqVertex, BaseEdgeStruct>(&g.base_graph)
                .bases_match(b"ACG")
        })
        .next()
        .unwrap();
    assert!(relative_eq!(
        ACG.score,
        -0.47712125471966244,
        epsilon = std::f64::EPSILON
    ));
}

#[test]
fn test_cycle_remove() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"C".to_vec());
    let v3 = SeqVertex::new(b"G".to_vec());
    let v4 = SeqVertex::new(b"T".to_vec());

    let nodes = g.base_graph.add_vertices(vec![v1, v2, v3, v4]);
    g.base_graph
        .graph
        .add_edge(nodes[0], nodes[1], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[2], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[2], nodes[1], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[2], nodes[3], BaseEdgeStruct::new(false, 1, 0));

    let sources = g
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = g
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    let mut finder = GraphBasedKBestHaplotypeFinder::new(&mut g.base_graph, sources, sinks);
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 1);
}

#[test]
fn test_no_source_or_sink() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"C".to_vec());

    let nodes = g.base_graph.add_vertices(vec![v1, v2]);
    g.base_graph
        .graph
        .add_edge(nodes[0], nodes[1], BaseEdgeStruct::new(false, 1, 0));

    let sources = g
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = g
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    let mut finder =
        GraphBasedKBestHaplotypeFinder::new(&mut g.base_graph, sources.clone(), sinks.clone());
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 1);

    let mut finder = GraphBasedKBestHaplotypeFinder::new(&mut g.base_graph, HashSet::new(), sinks);
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 0);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 1);

    let mut finder =
        GraphBasedKBestHaplotypeFinder::new(&mut g.base_graph, sources, HashSet::new());
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 0);
}

#[test]
fn test_dead_node() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"C".to_vec());
    let v3 = SeqVertex::new(b"G".to_vec());
    let v4 = SeqVertex::new(b"T".to_vec());
    let v5 = SeqVertex::new(b"A".to_vec());

    let nodes = g.base_graph.add_vertices(vec![v1, v2, v3, v4, v5]);
    g.base_graph
        .graph
        .add_edge(nodes[0], nodes[1], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[2], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[2], nodes[1], BaseEdgeStruct::new(false, 1, 0)); // cycle
    g.base_graph
        .graph
        .add_edge(nodes[1], nodes[4], BaseEdgeStruct::new(false, 1, 0));
    g.base_graph
        .graph
        .add_edge(nodes[0], nodes[3], BaseEdgeStruct::new(false, 1, 0));

    let sources = g
        .base_graph
        .get_sources_generic()
        .collect::<HashSet<NodeIndex>>();
    let sinks = g
        .base_graph
        .get_sinks_generic()
        .collect::<HashSet<NodeIndex>>();
    let mut g_copy = g.clone();
    let mut finder = GraphBasedKBestHaplotypeFinder::new(&mut g_copy.base_graph, sources, sinks);
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 2);

    let mut finder =
        GraphBasedKBestHaplotypeFinder::new_from_singletons(&mut g.base_graph, nodes[0], nodes[3]);
    assert_eq!(finder.k_best_haplotype_finder.sources.len(), 1);
    assert_eq!(finder.k_best_haplotype_finder.sinks.len(), 1);
}

fn test_basic_path_finding(n_start_nodes: usize, n_branches_per_bubble: usize, n_end_nodes: usize) {
    let mut graph = SeqGraph::new(11);

    let middle_top = SeqVertex::new(b"GTAC".to_vec());
    let middle_bottom = SeqVertex::new(b"ACTG".to_vec());
    let mut weight = 1;
    let source_sink = graph
        .base_graph
        .add_vertices(vec![middle_top, middle_bottom]);
    let starts = create_vertices(
        &mut graph,
        n_start_nodes,
        None,
        Some(source_sink[0]),
        &mut weight,
    );
    let bubbles = create_vertices(
        &mut graph,
        n_branches_per_bubble,
        Some(source_sink[0]),
        Some(source_sink[1]),
        &mut weight,
    );
    let ends = create_vertices(
        &mut graph,
        n_end_nodes,
        Some(source_sink[1]),
        None,
        &mut weight,
    );

    let expected_num_paths = n_start_nodes * n_branches_per_bubble * n_end_nodes;
    let haplotypes = GraphBasedKBestHaplotypeFinder::new(&mut graph.base_graph, starts, ends)
        .find_best_haplotypes(std::usize::MAX, &graph.base_graph);
    assert_eq!(haplotypes.len(), expected_num_paths);
    (1..haplotypes.len())
        .into_iter()
        .for_each(|n| assert!(haplotypes[n - 1].score >= haplotypes[n].score));
}

fn create_vertices(
    graph: &mut SeqGraph<BaseEdgeStruct>,
    n: usize,
    source: Option<NodeIndex>,
    target: Option<NodeIndex>,
    weight: &mut usize,
) -> HashSet<NodeIndex> {
    let mut seqs = vec![b"A".to_vec(), b"C".to_vec(), b"G".to_vec(), b"T".to_vec()];
    let mut vertices = HashSet::new();
    for i in 0..n {
        let mut v = SeqVertex::new(seqs[i].clone());
        let node = graph.base_graph.add_node(v);
        vertices.insert(node);
        match source {
            None => {
                //
            }
            Some(source) => {
                graph.base_graph.graph.add_edge(
                    source,
                    node,
                    BaseEdgeStruct::new(false, *weight, 0),
                );
                *weight += 1;
            }
        }

        match target {
            None => {
                //
            }
            Some(target) => {
                graph.base_graph.graph.add_edge(
                    node,
                    target,
                    BaseEdgeStruct::new(false, *weight, 0),
                );
                *weight += 1;
            }
        }
    }

    return vertices;
}

#[test]
fn make_basic_path_finding_data() {
    for n_start_nodes in vec![1, 2, 3] {
        for n_branches_per_bubble in vec![2, 3] {
            for n_end_nodes in vec![1, 2, 3] {
                test_basic_path_finding(n_start_nodes, n_branches_per_bubble, n_end_nodes)
            }
        }
    }
}
