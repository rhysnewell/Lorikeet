#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use gkl::smithwaterman::Parameters;
use lorikeet_genome::graphs::base_edge::BaseEdge;
use lorikeet_genome::graphs::base_edge::BaseEdgeStruct;
use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::graphs::path::Path;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::model::byte_array_allele::Allele;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::AVXMode;
use lorikeet_genome::reads::cigar_builder::CigarBuilder;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::utils::simple_interval::SimpleInterval;
use petgraph::prelude::NodeIndex;
use rust_htslib::bam::record::Cigar;
use std::collections::HashSet;


lazy_static! {
    static ref PATH_TO_REFERENCE_SW_PARAMETERS: Parameters = *NEW_SW_PARAMETERS;
}

#[test]
fn test_score() {
    let mut g = SeqGraph::new(3);
    let v1 = SeqVertex::new(b"A".to_vec());
    let v2 = SeqVertex::new(b"C".to_vec());
    let v3 = SeqVertex::new(b"G".to_vec());
    let v4 = SeqVertex::new(b"T".to_vec());
    let v5 = SeqVertex::new(b"A".to_vec());

    let nodes = g.base_graph.add_vertices(vec![&v1, &v2, &v3, &v4, &v5]);
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

    let nodes = g.base_graph.add_vertices(vec![&v1, &v2, &v3, &v4]);
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

    let nodes = g.base_graph.add_vertices(vec![&v1, &v2]);
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

    let nodes = g.base_graph.add_vertices(vec![&v1, &v2, &v3, &v4, &v5]);
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
        .add_vertices(vec![&middle_top, &middle_bottom]);
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
        let node = graph.base_graph.add_node(&v);
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

fn test_basic_bubble_data(ref_bubble_length: usize, alt_bubble_length: usize) {
    // Construct the assembly graph
    let mut graph = SeqGraph::new(3);
    let pre_ref = "ATGC";
    let post_ref = "GGGC";

    let v = SeqVertex::new(pre_ref.as_bytes().to_vec());
    let v2_ref = SeqVertex::new(vec![b'A'; ref_bubble_length]);
    let v2_alt = SeqVertex::new(
        ("A".repeat(alt_bubble_length - 1) + "T")
            .as_bytes()
            .to_vec(),
    );
    let v3 = SeqVertex::new(post_ref.as_bytes().to_vec());

    let v_indices = graph
        .base_graph
        .add_vertices(vec![&v, &v2_ref, &v2_alt, &v3]);
    graph.base_graph.add_or_update_edge(
        v_indices[0],
        v_indices[1],
        BaseEdgeStruct::new(true, 10, 0),
    );
    graph.base_graph.add_or_update_edge(
        v_indices[1],
        v_indices[3],
        BaseEdgeStruct::new(true, 10, 0),
    );
    let e_alt_1 = graph.base_graph.graph.add_edge(
        v_indices[0],
        v_indices[2],
        BaseEdgeStruct::new(false, 5, 0),
    );
    let e_alt_2 = graph.base_graph.graph.add_edge(
        v_indices[2],
        v_indices[3],
        BaseEdgeStruct::new(false, 5, 0),
    );

    // Construct the path
    let mut path = Path::new(v_indices[0], Vec::new());
    path = path.new_add_edge(e_alt_1, &graph.base_graph);
    path = path.new_add_edge(e_alt_2, &graph.base_graph);

    // Construct the actual cigar
    let mut expected_cigar = CigarBuilder::new(true);
    expected_cigar.add(Cigar::Match(pre_ref.len() as u32));

    if ref_bubble_length > alt_bubble_length {
        expected_cigar.add(Cigar::Del((ref_bubble_length - alt_bubble_length) as u32));
        expected_cigar.add(Cigar::Match(alt_bubble_length as u32));
    } else if ref_bubble_length < alt_bubble_length {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
        expected_cigar.add(Cigar::Ins((alt_bubble_length - ref_bubble_length) as u32));
    } else {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
    }

    expected_cigar.add(Cigar::Match(post_ref.len() as u32));

    let reference = (pre_ref.to_string() + v2_ref.get_sequence_string()) + post_ref;
    assert_eq!(
        path.calculate_cigar(
            reference.as_bytes(),
            &graph.base_graph,
            AVXMode::detect_mode()
        )
        .to_string(),
        expected_cigar.make(false).unwrap().to_string(),
        "Cigar string mismatch"
    );
}

#[test]
fn make_basic_bubble_data_provider() {
    for ref_bubble_length in vec![1, 5, 10] {
        for alt_bubble_length in vec![1, 5, 10] {
            test_basic_bubble_data(ref_bubble_length, alt_bubble_length);
        }
    }
}

fn test_triple_bubble_data(
    ref_bubble_length: usize,
    alt_bubble_length: usize,
    off_ref_beginning: bool,
    off_ref_ending: bool,
) {
    let mut graph = SeqGraph::new(11);
    let pre_alt_option = "ATCGATCGATCGATCGATCG";
    let post_alt_option = "CCCC";
    let pre_ref = "ATGG";
    let post_ref = "GGCCG";
    let mid_ref1 = "TTCCT";
    let mid_ref2 = "CCCAAAAAAAAAAAA";

    let pre_v = SeqVertex::new(pre_alt_option.as_bytes().to_vec());
    let v = SeqVertex::new(pre_ref.as_bytes().to_vec());
    let v2_ref = SeqVertex::new("A".repeat(ref_bubble_length).as_bytes().to_vec());
    let v2_alt = SeqVertex::new(
        ("A".repeat(alt_bubble_length - 1) + "T")
            .as_bytes()
            .to_vec(),
    );
    let v4_ref = SeqVertex::new("C".repeat(ref_bubble_length).as_bytes().to_vec());
    let v4_alt = SeqVertex::new(
        ("C".repeat(alt_bubble_length - 1) + "T")
            .as_bytes()
            .to_vec(),
    );
    let v6_ref = SeqVertex::new("G".repeat(ref_bubble_length).as_bytes().to_vec());
    let v6_alt = SeqVertex::new(
        ("G".repeat(alt_bubble_length - 1) + "T")
            .as_bytes()
            .to_vec(),
    );
    let v3 = SeqVertex::new(mid_ref1.as_bytes().to_vec());
    let v5 = SeqVertex::new(mid_ref2.as_bytes().to_vec());
    let v7 = SeqVertex::new(post_ref.as_bytes().to_vec());
    let post_v = SeqVertex::new(post_alt_option.as_bytes().to_vec());

    let reference = format!(
        "{}{}{}{}{}{}{}",
        pre_ref,
        v2_ref.get_sequence_string(),
        mid_ref1,
        v4_ref.get_sequence_string(),
        mid_ref2,
        v6_ref.get_sequence_string(),
        post_ref
    );

    let pre_v_i = graph.base_graph.add_node(&pre_v);
    let v_i = graph.base_graph.add_node(&v);
    let v2_ref_i = graph.base_graph.add_node(&v2_ref);
    let v2_alt_i = graph.base_graph.add_node(&v2_alt);
    let v3_i = graph.base_graph.add_node(&v3);
    let v4_ref_i = graph.base_graph.add_node(&v4_ref);
    let v4_alt_i = graph.base_graph.add_node(&v4_alt);
    let v5_i = graph.base_graph.add_node(&v5);
    let v6_ref_i = graph.base_graph.add_node(&v6_ref);
    let v6_alt_i = graph.base_graph.add_node(&v6_alt);
    let v7_i = graph.base_graph.add_node(&v7);
    let post_v_i = graph.base_graph.add_node(&post_v);

    let p1_i = graph
        .base_graph
        .graph
        .add_edge(pre_v_i, v_i, BaseEdgeStruct::new(false, 1, 0));
    let r1_i = graph
        .base_graph
        .graph
        .add_edge(v_i, v2_ref_i, BaseEdgeStruct::new(true, 10, 0));
    let r2_i = graph
        .base_graph
        .graph
        .add_edge(v2_ref_i, v3_i, BaseEdgeStruct::new(true, 10, 0));
    let a1_i = graph
        .base_graph
        .graph
        .add_edge(v_i, v2_alt_i, BaseEdgeStruct::new(false, 5, 0));
    let a2_i = graph
        .base_graph
        .graph
        .add_edge(v2_alt_i, v3_i, BaseEdgeStruct::new(false, 5, 0));
    let r3_i = graph
        .base_graph
        .graph
        .add_edge(v3_i, v4_ref_i, BaseEdgeStruct::new(true, 10, 0));
    let r4_i = graph
        .base_graph
        .graph
        .add_edge(v4_ref_i, v5_i, BaseEdgeStruct::new(true, 10, 0));
    let a3_i = graph
        .base_graph
        .graph
        .add_edge(v3_i, v4_alt_i, BaseEdgeStruct::new(false, 5, 0));
    let a4_i = graph
        .base_graph
        .graph
        .add_edge(v4_alt_i, v5_i, BaseEdgeStruct::new(false, 5, 0));
    let r5_i = graph
        .base_graph
        .graph
        .add_edge(v5_i, v6_ref_i, BaseEdgeStruct::new(true, 11, 0));
    let r6_i = graph
        .base_graph
        .graph
        .add_edge(v6_ref_i, v7_i, BaseEdgeStruct::new(true, 11, 0));
    let a5_i = graph
        .base_graph
        .graph
        .add_edge(v5_i, v6_alt_i, BaseEdgeStruct::new(false, 55, 0));
    let a6_i = graph
        .base_graph
        .graph
        .add_edge(v6_alt_i, v7_i, BaseEdgeStruct::new(false, 55, 0));
    let p2_i = graph
        .base_graph
        .graph
        .add_edge(v7_i, post_v_i, BaseEdgeStruct::new(false, 1, 0));

    let mut path = Path::new(if off_ref_beginning { pre_v_i } else { v_i }, Vec::new());
    if off_ref_beginning {
        path = path.new_add_edge(p1_i, &graph.base_graph);
    }
    path = path.new_add_edge(a1_i, &graph.base_graph);
    path = path.new_add_edge(a2_i, &graph.base_graph);
    path = path.new_add_edge(r3_i, &graph.base_graph);
    path = path.new_add_edge(r4_i, &graph.base_graph);
    path = path.new_add_edge(a5_i, &graph.base_graph);
    path = path.new_add_edge(a6_i, &graph.base_graph);
    if off_ref_ending {
        path = path.new_add_edge(p2_i, &graph.base_graph);
    }

    let mut expected_cigar = CigarBuilder::new(true);
    if off_ref_beginning {
        expected_cigar.add(Cigar::Ins(pre_alt_option.len() as u32));
    }
    expected_cigar.add(Cigar::Match(pre_ref.len() as u32));
    //first bubble
    if ref_bubble_length > alt_bubble_length {
        expected_cigar.add(Cigar::Del((ref_bubble_length - alt_bubble_length) as u32));
        expected_cigar.add(Cigar::Match(alt_bubble_length as u32));
    } else if ref_bubble_length < alt_bubble_length {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
        expected_cigar.add(Cigar::Ins((alt_bubble_length - ref_bubble_length) as u32));
    } else {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
    }

    expected_cigar.add(Cigar::Match(mid_ref1.len() as u32));

    // second bubble is ref path
    expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
    expected_cigar.add(Cigar::Match(mid_ref2.len() as u32));

    // third bubble
    if ref_bubble_length > alt_bubble_length {
        expected_cigar.add(Cigar::Del((ref_bubble_length - alt_bubble_length) as u32));
        expected_cigar.add(Cigar::Match(alt_bubble_length as u32));
    } else if ref_bubble_length < alt_bubble_length {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
        expected_cigar.add(Cigar::Ins((alt_bubble_length - ref_bubble_length) as u32));
    } else {
        expected_cigar.add(Cigar::Match(ref_bubble_length as u32));
    }
    expected_cigar.add(Cigar::Match(post_ref.len() as u32));

    if off_ref_ending {
        expected_cigar.add(Cigar::Ins(post_alt_option.len() as u32));
    };

    assert_eq!(
        path.calculate_cigar(reference.as_bytes(), &graph.base_graph, AVXMode::detect_mode()).to_string(),
        expected_cigar.make(false).unwrap().to_string(),
        "Cigar string mismatch ref bubble len {} alt bubble len {} off ref beginning {} off ref ending {}",
        ref_bubble_length,
        alt_bubble_length,
        off_ref_beginning,
        off_ref_ending
    )
}

#[test]
fn make_triple_bubble_data_provider() {
    for ref_bubble_length in vec![1, 5, 10] {
        for alt_bubble_length in vec![1, 5, 10] {
            for off_ref_ending in vec![true, false] {
                for off_ref_beginning in vec![false] {
                    test_triple_bubble_data(
                        ref_bubble_length,
                        alt_bubble_length,
                        off_ref_beginning,
                        off_ref_ending,
                    );
                }
            }
        }
    }
}

#[test]
fn test_intra_node_insertion_deletion() {
    let mut graph = SeqGraph::new(11);
    let top = SeqVertex::new(b"T".to_vec());
    let bot = SeqVertex::new(b"T".to_vec());
    let alternate = SeqVertex::new(b"AAACCCCC".to_vec());
    let reference = SeqVertex::new(b"CCCCCGGG".to_vec());

    let nodes = graph
        .base_graph
        .add_vertices(vec![&top, &bot, &alternate, &reference]);
    graph.base_graph.add_edges(
        nodes[0],
        vec![nodes[3], nodes[1]],
        BaseEdgeStruct::new(true, 1, 0),
    );
    graph.base_graph.add_edges(
        nodes[0],
        vec![nodes[2], nodes[1]],
        BaseEdgeStruct::new(false, 1, 0),
    );

    let mut sources = HashSet::new();
    sources.insert(nodes[0]);

    let mut sinks = HashSet::new();
    sinks.insert(nodes[1]);
    let best_paths = GraphBasedKBestHaplotypeFinder::new(&mut graph.base_graph, sources, sinks)
        .find_best_haplotypes(std::usize::MAX, &graph.base_graph);
    assert_eq!(best_paths.len(), 2);
    let ref_path = best_paths[0].clone();
    let alt_path = best_paths[1].clone();

    let ref_string = format!(
        "{}{}{}",
        top.get_sequence_string(),
        reference.get_sequence_string(),
        bot.get_sequence_string()
    );
    assert_eq!(
        ref_path
            .path
            .calculate_cigar(
                ref_string.as_bytes(),
                &graph.base_graph,
                AVXMode::detect_mode()
            )
            .to_string()
            .as_str(),
        "10M"
    );
    assert_eq!(
        alt_path
            .path
            .calculate_cigar(
                ref_string.as_bytes(),
                &graph.base_graph,
                AVXMode::detect_mode()
            )
            .to_string()
            .as_str(),
        "1M3I5M3D1M"
    );
}

/*
This is a test of what can go awry if path pruning is based on the number of incoming edges to a given vertex.
An illustration of the graph in this test:
          / ----- top \
         /(33)         \
refStart -(32) -- mid -- midExt ---------------------------- refEnd
         \(34)        / (1)                               /
          \ ---- bot /- (33) botExt - (15) - botExtTop - /
                                   \ (18)               /
                                    \ ------ botExtBot /

 The expected best paths are (refStart->top->midExt->refEnd), and (refStart->mid->midExt->refEnd) because the bottom
 path is penalized for two extra forks despite it greedily looking the best at the start.

 Because the old behavior used to base pruning on the number of incoming edges < K, the edge (bot->midExt) would be created
 first, then (top -> midExt) but (mid -> midExt) would never be created because midExt already has 2 incoming edges.
 */
#[test]
fn test_degenerate_path_pruning_optimization() {
    // Construct an assembly graph demonstrating this issue
    let mut graph = SeqGraph::new(11);
    let top = SeqVertex::new(b"T".to_vec());
    let mid = SeqVertex::new(b"C".to_vec());
    let mid_and_top_ext = SeqVertex::new(b"GGG".to_vec());
    let bot = SeqVertex::new(b"G".to_vec());
    let bot_ext = SeqVertex::new(b"AAA".to_vec());
    let bot_ext_top = SeqVertex::new(b"A".to_vec());
    let bot_ext_bot = SeqVertex::new(b"T".to_vec());

    // Ref sources and sink vertices
    let ref_start = SeqVertex::new(b"CCCCCGGG".to_vec());
    let ref_end = SeqVertex::new(b"TTTT".to_vec());

    let nodes = graph.base_graph.add_vertices(vec![
        &top,
        &bot,
        &mid, // 0 1 2
        &mid_and_top_ext,
        &bot, // 3 4
        &bot_ext,
        &bot_ext_top,
        &bot_ext_bot, // 5 6 7
        &ref_start,
        &ref_end, // 8 9
    ]);

    // First diamon with 3 mostly equivalent cost paths
    graph.base_graph.add_edges(
        nodes[8],
        vec![nodes[1], nodes[5]],
        BaseEdgeStruct::new(false, 34, 0),
    );
    graph.base_graph.add_edges(
        nodes[8],
        vec![nodes[0], nodes[3]],
        BaseEdgeStruct::new(true, 33, 0),
    );
    graph.base_graph.add_edges(
        nodes[8],
        vec![nodes[2], nodes[3]],
        BaseEdgeStruct::new(false, 32, 0),
    );

    // The best looking path reconnects with a very poor edge multiplicity to mid_and_top_ext (this is what causes potential bug)
    graph
        .base_graph
        .add_edges(nodes[1], vec![nodes[3]], BaseEdgeStruct::new(false, 1, 0));
    // There is another diamond at bot and ext that will end up discounting that path from being in the k best
    graph.base_graph.add_edges(
        nodes[5],
        vec![nodes[6], nodes[9]],
        BaseEdgeStruct::new(false, 15, 0),
    );
    graph.base_graph.add_edges(
        nodes[5],
        vec![nodes[7], nodes[9]],
        BaseEdgeStruct::new(false, 18, 0),
    );

    // Whereas this path is smooth sailing from ref end
    graph
        .base_graph
        .add_edges(nodes[3], vec![nodes[9]], BaseEdgeStruct::new(true, 65, 0));

    let best_paths = GraphBasedKBestHaplotypeFinder::new_from_singletons(
        &mut graph.base_graph,
        nodes[8],
        nodes[9],
    )
    .find_best_haplotypes(5, &graph.base_graph);

    assert_eq!(best_paths.len(), 5);
    let ref_path = best_paths[0].clone();
    let alt_path = best_paths[1].clone();

    assert_eq!(
        ref_path.path.get_vertices(&graph.base_graph),
        vec![nodes[8], nodes[0], nodes[3], nodes[9]]
    );
    assert_eq!(
        alt_path.path.get_vertices(&graph.base_graph),
        vec![nodes[8], nodes[2], nodes[3], nodes[9]]
    );
}

#[test]
fn test_hard_sw_path() {
    let mut graph = SeqGraph::new(11);
    let top = SeqVertex::new(b"NNN".to_vec());
    let bot = SeqVertex::new(b"NNN".to_vec());
    let alternate = SeqVertex::new(b"ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA".to_vec());
    let reference = SeqVertex::new(b"TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA".to_vec());

    let nodes = graph.base_graph.add_vertices(vec![
        &top,
        &bot,
        &alternate,
        &reference,
    ]);
    graph.base_graph.add_edges(
        nodes[0],
        vec![nodes[3], nodes[1]],
        BaseEdgeStruct::new(true, 1, 0),
    );
    graph.base_graph.add_edges(
        nodes[0],
        vec![nodes[2], nodes[1]],
        BaseEdgeStruct::new(false, 1, 0),
    );

    let mut sources = HashSet::new();
    sources.insert(nodes[0]);

    let mut sinks = HashSet::new();
    sinks.insert(nodes[1]);
    let best_paths = GraphBasedKBestHaplotypeFinder::new(&mut graph.base_graph, sources, sinks)
        .find_best_haplotypes(std::usize::MAX, &graph.base_graph);
    assert_eq!(best_paths.len(), 2);
    let ref_path = best_paths[0].clone();
    let alt_path = best_paths[1].clone();

    let ref_string = format!(
        "{}{}{}",
        top.get_sequence_string(),
        reference.get_sequence_string(),
        bot.get_sequence_string()
    );
    assert_eq!(
        ref_path
            .path
            .calculate_cigar(
                ref_string.as_bytes(),
                &graph.base_graph,
                AVXMode::detect_mode()
            )
            .to_string()
            .as_str(),
        "51M"
    );
    assert_eq!(
        alt_path
            .path
            .calculate_cigar(
                ref_string.as_bytes(),
                &graph.base_graph,
                AVXMode::detect_mode()
            )
            .to_string()
            .as_str(),
        "3M6I48M"
    );
}
