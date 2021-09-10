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

use bio::io::fasta::IndexedReader;
use itertools::Itertools;
use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::assembly::assembly_result_set::AssemblyResultSet;
use lorikeet_genome::genotype::genotype_builder::Genotype;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::genotype::genotype_likelihoods::GenotypeLikelihoods;
use lorikeet_genome::graphs::adaptive_chain_pruner::AdaptiveChainPruner;
use lorikeet_genome::graphs::base_edge::{BaseEdge, BaseEdgeStruct};
use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::graphs::chain_pruner::ChainPruner;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::graphs::seq_vertex::SeqVertex;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::allele_frequency_calculator::AlleleFrequencyCalculator;
use lorikeet_genome::model::allele_likelihoods::AlleleLikelihoods;
use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::{allele_list::AlleleList, variants::SPAN_DEL_ALLELE};
use lorikeet_genome::pair_hmm::pair_hmm::PairHMM;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::PairHMMInputScoreImputator;
use lorikeet_genome::read_error_corrector::nearby_kmer_error_corrector::{
    CorrectionSet, NearbyKmerErrorCorrector,
};
use lorikeet_genome::read_error_corrector::read_error_corrector::ReadErrorCorrector;
use lorikeet_genome::read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use lorikeet_genome::read_threading::read_threading_assembler::ReadThreadingAssembler;
use lorikeet_genome::read_threading::read_threading_graph::ReadThreadingGraph;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reference::reference_reader_utils::ReferenceReaderUtils;
use lorikeet_genome::smith_waterman::bindings::SWParameters;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::base_utils::BaseUtils;
use lorikeet_genome::utils::math_utils::{MathUtils, LOG10_ONE_HALF};
use lorikeet_genome::utils::quality_utils::QualityUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use rand::rngs::ThreadRng;
use rand::seq::index::sample;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView, Seq};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::fs::File;
use std::ops::Deref;
use std::sync::Mutex;

// fn test_prune_low_weight_chains<E: BaseEdge>(
//     name: &str,
//     mut graph: SeqGraph<E>,
//     prune_factor: usize,
//     remaining_vertices: HashSet<SeqVertex>,
//     chain_count_before_pruning: Option<usize>,
// ) {
//     println!("remaining vertices {:?}", &remaining_vertices);
//     let copy = remaining_vertices.clone();
//     let mut pruner = AdaptiveChainPruner::new(
//         0.001, MathUtils::log10_to_log(1.0), MathUtils::log10_to_log(4.0), 50
//     );
//
//     match chain_count_before_pruning {
//         Some(chain_count_before_pruning) => {
//             assert_eq!(AdaptiveChainPruner::find_all_chains(&graph.base_graph).len(), chain_count_before_pruning);
//         },
//         _ => {
//             // pass
//         }
//     };
//
//     pruner.prune_low_weight_chains(&mut graph.base_graph);
//     assert_eq!(graph.base_graph.graph.node_weights().cloned().collect::<HashSet<SeqVertex>>(), copy);
// }
//
// #[test]
// fn make_prune_chains_data() {
//     let v1 = SeqVertex::new("A".to_string());
//     let v2 = SeqVertex::new("C".to_string());
//     let v3 = SeqVertex::new("G".to_string());
//     let v4 = SeqVertex::new("T".to_string());
//     let v5 = SeqVertex::new("AA".to_string());
//     let v6 = SeqVertex::new("CC".to_string());
//
//     for edge_weight in vec![1, 2, 3] {
//         for prune_factor in vec![1, 2, 3, 4] {
//             for is_ref in vec![true, false] {
//                 // just an isolated chain
//                 let n_expected = if edge_weight < prune_factor && !is_ref {
//                     3
//                 } else {
//                     0
//                 };
//                 println!("edge weight {} prune factor {} is ref {}", edge_weight, prune_factor, is_ref);
//                 let mut graph = SeqGraph::new(11);
//                 let new_nodes = graph.base_graph.add_vertices(vec![v1.clone(), v2.clone(), v3.clone()]);
//                 graph.base_graph.add_edges(new_nodes[0], new_nodes[1..].to_vec(), BaseEdgeStruct::new(is_ref, edge_weight));
//                 test_prune_low_weight_chains(
//                     "combinatorial",
//                     graph.clone(),
//                     prune_factor,
//                     if n_expected > 0 {
//                         HashSet::new()
//                     } else {
//                         graph.base_graph.graph.node_weights().cloned().collect::<HashSet<SeqVertex>>()
//                     },
//                     Some(1)
//                 );
//             }
//         }
//     }
// }

// test that in graph with good path A -> B -> C and bad edges A -> D -> C and D -> B that the adjacency of bad edges --
// such that when bad edges meet the multiplicities do not indicate an error - does not harm pruning.
// we test with and without a true variant path A -> E -> C
#[test]
fn test_adaptive_pruning_with_adjacent_bad_edges() {
    let good_multiplicity = 1000;
    let variant_multiplicity = 50;
    let bad_multiplicity = 5;

    let source = SeqVertex::new("source".to_string());
    let sink = SeqVertex::new("sink".to_string());
    let A = SeqVertex::new("A".to_string());
    let B = SeqVertex::new("B".to_string());
    let C = SeqVertex::new("C".to_string());
    let D = SeqVertex::new("D".to_string());
    let E = SeqVertex::new("E".to_string());

    for variant_present in vec![false, true] {
        let mut graph = SeqGraph::new(20);

        let node_indices = graph.base_graph.add_vertices(vec![
            source.clone(),
            A.clone(),
            B.clone(),
            C.clone(),
            D.clone(),
            sink.clone(),
        ]);
        graph.base_graph.add_edges(
            NodeIndex::new(0),
            vec![
                NodeIndex::new(1),
                NodeIndex::new(2),
                NodeIndex::new(3),
                NodeIndex::new(5),
            ],
            BaseEdgeStruct::new(true, good_multiplicity),
        );
        graph.base_graph.add_edges(
            NodeIndex::new(1),
            vec![NodeIndex::new(4), NodeIndex::new(3)],
            BaseEdgeStruct::new(false, bad_multiplicity),
        );
        graph.base_graph.add_edges(
            NodeIndex::new(4),
            vec![NodeIndex::new(2)],
            BaseEdgeStruct::new(false, bad_multiplicity),
        );

        if variant_present {
            graph.base_graph.add_vertices(vec![E.clone()]);
            graph.base_graph.add_edges(
                NodeIndex::new(1),
                vec![NodeIndex::new(6), NodeIndex::new(2)],
                BaseEdgeStruct::new(false, variant_multiplicity),
            );
        };

        let mut pruner = AdaptiveChainPruner::new(0.01, 2.0, MathUtils::log10_to_log(4.0), 50);
        pruner.prune_low_weight_chains(&mut graph.base_graph);

        assert!(!graph.base_graph.graph.contains_node(NodeIndex::new(4)));
        if variant_present {
            assert!(
                graph.base_graph.graph.contains_node(NodeIndex::new(6)),
                "Node not found {:?}",
                graph.base_graph.graph.node_indices()
            );
        }
    }
}

/**
 * Comprehensive test of adaptive pruning -- given an alt haplotype and a ref haplotype
 * @param kmerSize
 * @param ref           reference haplotype
 * @param alt           alt haplotype
 * @param altFraction   alt allele fraction to simulate somatic, mitochondrial etc variants
 * @param errorRate     substitution error rate of simulated reads
 * @param depthPerAlignmentStart    number of reads starting at each position.  Note that holding this constant yields
 *                                  low coverage at the beginning of the graph and high in the middle and end, simulating
 *                                  the leading edge of an exome target, for example
 * @param logOddsThreshold  log-10 odds threshold for pruning chains
 */
fn test_adaptive_pruning(
    kmer_size: usize,
    reference: &[u8],
    alternate: &[u8],
    alt_fraction: f64,
    error_rate: f64,
    depth_alignment_start: usize,
    log_odds_threshold: f64,
) {
    let mut rng = ThreadRng::default();
    let mut graph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    graph.add_sequence(
        "anonymous".to_string(),
        "XXX_UNNAMED_XXX",
        reference,
        0,
        reference.len(),
        1,
        true,
    );

    let reads = (0..reference.len())
        .map(|start| {
            (0..depth_alignment_start)
                .map(|n| {
                    generate_read_with_errors(
                        if rng.gen_range(0.0, 1.0) < alt_fraction {
                            alternate
                        } else {
                            reference
                        },
                        start,
                        None,
                        error_rate,
                        &mut rng,
                    )
                })
                .collect::<Vec<Vec<u8>>>()
        })
        .flat_map(|s| s)
        .collect::<Vec<Vec<u8>>>();

    reads.into_iter().for_each(|r| {
        graph.add_sequence(
            "anonymous".to_string(),
            "XXX_UNNAMED_XXX",
            r.as_slice(),
            0,
            r.len(),
            1,
            false,
        )
    });
    println!("Reads added");

    // note: these are the steps in ReadThreadingAssembler::createGraph
    graph.build_graph_if_necessary();
    println!("Graph built");

    let mut pruner =
        AdaptiveChainPruner::new(0.001, log_odds_threshold, MathUtils::log10_to_log(4.0), 50);
    pruner.prune_low_weight_chains(graph.get_base_graph_mut());
    println!("Low weight chains pruned");

    graph.recover_dangling_tails(1, 3, false, &*STANDARD_NGS);
    graph.recover_dangling_heads(1, 3, false, &*STANDARD_NGS);
    graph.remove_paths_not_connected_to_ref();

    let mut seq_graph = graph.to_sequence_graph();
    seq_graph.zip_linear_chains();
    seq_graph.base_graph.remove_singleton_orphan_vertices();
    seq_graph
        .base_graph
        .remove_vertices_not_connected_to_ref_regardless_of_edge_direction();
    seq_graph.simplify_graph();
    seq_graph.base_graph.remove_paths_not_connected_to_ref();
    seq_graph.simplify_graph();

    let best_paths = GraphBasedKBestHaplotypeFinder::new(
        &seq_graph.base_graph,
        seq_graph
            .base_graph
            .get_sources_generic()
            .collect::<HashSet<NodeIndex>>(),
        seq_graph
            .base_graph
            .get_sinks_generic()
            .collect::<HashSet<NodeIndex>>(),
    )
    .find_best_haplotypes(10);

    let alt_index = (0..10)
        .filter(|n| {
            best_paths[*n]
                .haplotype::<SimpleInterval, SeqVertex, BaseEdgeStruct>(&seq_graph.base_graph)
                .bases_match(alternate)
        })
        .next();
    assert!(alt_index.is_some(), "Alt index not present in haplotypes");

    // ref path should not be pruned even if all reads are alt
    let ref_index = (0..10)
        .filter(|n| {
            best_paths[*n]
                .haplotype::<SimpleInterval, SeqVertex, BaseEdgeStruct>(&seq_graph.base_graph)
                .bases_match(reference)
        })
        .next();
    assert!(ref_index.is_some(), "Ref index not present in haplotypes");

    // the haplotype score is the sum of the log-10 of all branching fractions, so the alt haplotype score should come out to
    // around the log-10 of the allele fraction up to some fudge factor, assuming we didn't do any dumb pruning
    assert!(relative_eq!(
        best_paths[alt_index.unwrap()].score,
        (alt_fraction).log10(),
        epsilon = 0.5
    ));
    assert!(best_paths.len() < 15);
}

// #[test]
fn get_chain_pruner_data() {
    let mut rng = ThreadRng::default();
    let ref_length = 100;
    let left_snv_position = 15;
    let middle_snv_position = ref_length / 2;
    let right_snv_position = ref_length - left_snv_position;

    let mut reference = vec![0; ref_length];
    reference
        .iter_mut()
        .for_each(|b| *b = BaseUtils::base_index_to_simple_base(rng.gen_range(0, 4)));
    reference[left_snv_position] = b'A';
    reference[middle_snv_position] = b'G';
    reference[right_snv_position] = b'T';

    let mut left_snv = reference.clone();
    left_snv[left_snv_position] = b'G';

    let mut middle_snv = reference.clone();
    middle_snv[middle_snv_position] = b'T';

    let mut right_snv = reference.clone();
    right_snv[right_snv_position] = b'A';

    // kmer size, ref bases, alt bases, alt fraction, base error rate, depth per start, log odds threshold, max unpruned variants
    test_adaptive_pruning(
        10,
        reference.as_slice(),
        left_snv.as_slice(),
        0.5,
        0.001,
        20,
        MathUtils::log10_to_log(1.0),
    );
    test_adaptive_pruning(
        10,
        reference.as_slice(),
        left_snv.as_slice(),
        1.0,
        0.001,
        20,
        MathUtils::log10_to_log(1.0),
    );
    test_adaptive_pruning(
        10,
        reference.as_slice(),
        middle_snv.as_slice(),
        0.1,
        0.001,
        5,
        MathUtils::log10_to_log(1.0),
    );
    test_adaptive_pruning(
        25,
        reference.as_slice(),
        middle_snv.as_slice(),
        0.1,
        0.001,
        5,
        MathUtils::log10_to_log(1.0),
    );
    test_adaptive_pruning(
        25,
        reference.as_slice(),
        middle_snv.as_slice(),
        0.01,
        0.001,
        1000,
        MathUtils::log10_to_log(1.0),
    ); // note the extreme depth -- this would confuse non-adaptive pruning
    test_adaptive_pruning(
        10,
        reference.as_slice(),
        right_snv.as_slice(),
        0.1,
        0.001,
        2,
        MathUtils::log10_to_log(1.0),
    );
}

fn generate_read_with_errors(
    sequence: &[u8],
    start: usize,
    end: Option<usize>,
    error_rate: f64,
    rng: &mut ThreadRng,
) -> Vec<u8> {
    let end = match end {
        None => sequence.len(),
        Some(end) => end,
    };

    let adjusted_error_rate = (error_rate * 4.0) / 3.0; // one time in four a random base won't be an error
    let mut result = vec![0; end - start];
    (start..end).for_each(|n| {
        result[n - start] = if rng.gen_range(0.0, 1.0) > adjusted_error_rate {
            sequence[n]
        } else {
            BaseUtils::base_index_to_simple_base(rng.gen_range(0, 4))
        }
    });

    return result;
}
