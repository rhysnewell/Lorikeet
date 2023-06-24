#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;

use gkl::smithwaterman::Parameters;
use hashlink::LinkedHashMap;
use itertools::Itertools;
use lorikeet_genome::assembly::kmer::Kmer;
use lorikeet_genome::graphs::base_edge::BaseEdge;
use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::graphs::graph_based_k_best_haplotype_finder::GraphBasedKBestHaplotypeFinder;
use lorikeet_genome::read_threading::abstract_read_threading_graph::{
    AbstractReadThreadingGraph, SequenceForKmers,
};
use lorikeet_genome::read_threading::read_threading_graph::ReadThreadingGraph;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reads::read_utils::ReadUtils;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    STANDARD_NGS,
};
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::base_utils::BaseUtils;

use petgraph::Direction;
use rand::prelude::ThreadRng;
use rand::Rng;
use rust_htslib::bam::record::CigarString;
use std::cmp::min;
use std::collections::HashSet;
use std::convert::TryFrom;

lazy_static! {
    static ref DANGLING_END_SW_PARAMETERS: Parameters = *STANDARD_NGS;
}

fn get_bytes(alignment: &str) -> String {
    alignment.replace('-', "")
}

fn assert_non_uniques(assembler: &mut ReadThreadingGraph, non_uniques: HashSet<Kmer>, pending: &mut LinkedHashMap<usize, Vec<SequenceForKmers>>) {
    let mut actual = HashSet::new();
    assembler.build_graph_if_necessary(pending);
    println!("Pending: {:?}", pending);

    let non_uniques_from_graph = assembler
        .get_non_uniques()
        .iter()
        .cloned()
        .collect::<Vec<Kmer>>();
    println!("Non-unique {:?}", &non_uniques_from_graph);
    for kmer in non_uniques_from_graph {
        actual.insert(kmer);
    }

    assert_eq!(non_uniques, actual);
}

#[test]
fn test_simple_haplotype_rethreading() {
    let mut assembler = ReadThreadingGraph::default_with_kmer_size(11);
    let reference = "CATGCACTTTAAAACTTGCCTTTTTAACAAGACTTCCAGATG";
    let alternate = "CATGCACTTTAAAACTTGCCGTTTTAACAAGACTTCCAGATG";
    let mut pending = LinkedHashMap::new();

    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        alternate.as_bytes(),
        0,
        alternate.len(),
        1,
        false,
    );

    assembler.build_graph_if_necessary(&mut pending);
    assert_ne!(
        reference.len() - 11 + 1,
        assembler.get_base_graph().vertex_set().len(),
        "the number of vertex in the graph is the same as if there was no alternative sequence"
    );
    assert_eq!(
        reference.len() - 11 + 1 + 11,
        assembler.get_base_graph().vertex_set().len(),
        "the number of vertex in the graph is not the same as if there an alternative sequence"
    );

    let start_alt = assembler.find_kmer(&Kmer::new_with_start_and_length(
        alternate.as_bytes(),
        20,
        11,
    ));
    assert!(start_alt.is_some());
}

#[test]
fn test_non_unique_middle() {
    let mut assembler = ReadThreadingGraph::default_with_kmer_size(3);
    let reference = "GACACACAGTCA";
    let read1 = "GACACGTCA";
    let read2 = "CACGTCA";
    let mut pending = LinkedHashMap::new();
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read1.as_bytes(),
        0,
        read1.len(),
        1,
        false,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read2.as_bytes(),
        0,
        read2.len(),
        1,
        false,
    );

    assert_non_uniques(&mut assembler, vec![Kmer::new(b"ACA"), Kmer::new(b"CAC")].into_iter().collect(), &mut pending);
}

#[test]
fn test_read_creation_non_unique() {
    let mut assembler = ReadThreadingGraph::default_with_kmer_size(3);
    let reference = "GCACGTCA"; // CAC is unique
    let read1 = "GCACACGTCA"; // makes CAC non unique because it has a duplication
    let read2 = "CACGTCA"; // shouldn't be allowed to match CAC as start
    let mut pending = LinkedHashMap::new();
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read1.as_bytes(),
        0,
        read1.len(),
        1,
        false,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read2.as_bytes(),
        0,
        read2.len(),
        1,
        false,
    );

    assert_non_uniques(&mut assembler, vec![Kmer::new(b"CAC")].into_iter().collect(), &mut pending);
}

#[test]
fn test_counting_of_start_edges() {
    let mut assembler = ReadThreadingGraph::default_with_kmer_size(3);
    let reference = "NNNGTCAAA"; // ref has some bases before start
    let read1 = "GTCAAA"; // starts at first non N base
    let mut pending = LinkedHashMap::new();
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read1.as_bytes(),
        0,
        read1.len(),
        1,
        false,
    );
    assembler.build_graph_if_necessary(&mut pending);

    for edge in assembler.get_base_graph().graph.edge_indices() {
        let source = assembler.get_base_graph().get_edge_source(edge);
        let source_weight = assembler
            .get_base_graph()
            .graph
            .node_weight(source)
            .unwrap();
        let target = assembler.get_base_graph().get_edge_target(edge);
        let target_weight = assembler
            .get_base_graph()
            .graph
            .node_weight(target)
            .unwrap();

        let header_vertex =
            source_weight.get_suffix() == b'N' || target_weight.get_suffix() == b'N';
        if header_vertex {
            assert_eq!(
                assembler
                    .get_base_graph()
                    .graph
                    .edge_weight(edge)
                    .unwrap()
                    .get_multiplicity(),
                1,
                "Bases in the unique reference header should have multiplicity of 1"
            );
        } else {
            assert_eq!(assembler
                           .get_base_graph()
                           .graph.edge_weight(edge)
                           .unwrap().get_multiplicity(), 2,
                       "Should have multiplicity of 2 for any edge outside the ref header but got {:?} {:?} -> {:?}", assembler
                           .get_base_graph()
                           .graph.edge_weight(edge), source, target);
        }
    }
}

#[test]
fn test_counting_of_start_edges_with_multiple_prefixes() {
    let mut assembler = ReadThreadingGraph::default_with_kmer_size(3);
    assembler.set_increase_counts_through_branches(true);
    let reference = "NNNGTCAXX"; // ref has some bases before start
    let alt1 = "NNNCTCAXX"; // alt1 has SNP right after N
    let read1 = "TCAXX"; // starts right after SNP, but merges right before branch
    let mut pending = LinkedHashMap::new();
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        alt1.as_bytes(),
        0,
        alt1.len(),
        1,
        false,
    );
    assembler.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        read1.as_bytes(),
        0,
        read1.len(),
        1,
        false,
    );
    assembler.build_graph_if_necessary(&mut pending);
    // assembler.get_base_graph().print_graph("test.dot", true, 0);

    let one_count_vertices = vec!["NNN", "NNG", "NNC", "NGT", "NCT"];
    let three_count_vertices = vec!["CAX", "AXX"];

    for edge in assembler.get_base_graph().graph.edge_indices() {
        let source = assembler.get_base_graph().get_edge_source(edge);
        let _source_weight = assembler
            .get_base_graph()
            .graph
            .node_weight(source)
            .unwrap();
        let target = assembler.get_base_graph().get_edge_target(edge);
        let target_weight = assembler
            .get_base_graph()
            .graph
            .node_weight(target)
            .unwrap();

        let expected = if one_count_vertices.contains(&target_weight.get_sequence_string()) {
            1
        } else if three_count_vertices.contains(&target_weight.get_sequence_string()) {
            3
        } else {
            2
        };

        assert_eq!(
            assembler
                .get_base_graph()
                .graph.edge_weight(edge)
                .unwrap().get_multiplicity(),
            expected,
            "Should have multiplicity of {} for any edge outside the ref header but got {:?} {:?} -> {:?}",
            expected,
            assembler
                .get_base_graph()
                .graph.edge_weight(edge), source, target
        );
    }
}

#[test]
fn test_cycles_in_graph() {
    // b37 20:12655200-12655850
    let reference = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

    // SNP at 20:12655528 creates a cycle for small kmers
    let alternate = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

    let mut reads = Vec::new();
    let quals = vec![30; 100];
    for index in (0..alternate.len() - 100).step_by(20) {
        reads.push(ArtificialReadUtils::create_artificial_read(
            &alternate.as_bytes()[index..index + 100],
            &quals,
            CigarString::try_from("100M").unwrap(),
        ));
    }

    // test that there are cycles detected for small kmer
    let mut pending = LinkedHashMap::new();
    let mut rtgraph25 = ReadThreadingGraph::default_with_kmer_size(25);
    rtgraph25.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    let samples = vec!["anonymous_sample".to_string()];
    let mut count = 0;
    for read in reads.iter() {
        rtgraph25.add_read(read, &samples, &mut count, &mut pending)
    }

    rtgraph25.build_graph_if_necessary(&mut pending);
    assert!(rtgraph25.has_cycles());

    // test that there are no cycles detected for large kmer
    let mut rtgraph75 = ReadThreadingGraph::default_with_kmer_size(75);
    let mut pending = LinkedHashMap::new();
    rtgraph75.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    let samples = vec!["anonymous_sample".to_string()];
    let mut count = 0;
    for read in reads.iter() {
        rtgraph75.add_read(read, &samples, &mut count, &mut pending)
    }

    rtgraph75.build_graph_if_necessary(&mut pending);
    assert!(!rtgraph75.has_cycles());
}

// Test showing that if a read gets completely clipped in the ReadThreadingAssembler, that the assembly will not crash
#[test]
fn test_empty_read_being_added_to_graph() {
    // b37 20:12655200-12655850
    let reference = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

    // SNP at 20:12655528 creates a cycle for small kmers
    let alternate = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

    let mut reads = Vec::new();
    let quals = vec![30; 100];
    for index in (0..alternate.len() - 100).step_by(20) {
        reads.push(ArtificialReadUtils::create_artificial_read(
            &alternate.as_bytes()[index..index + 100],
            &quals,
            CigarString::try_from("100M").unwrap(),
        ));
    }
    reads.push(ReadUtils::empty_read(&reads[0]));

    // test that there are cycles detected for small kmer
    let mut rtgraph25 = ReadThreadingGraph::default_with_kmer_size(25);
    let mut pending = LinkedHashMap::new();
    rtgraph25.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    let samples = vec!["anonymous_sample".to_string()];
    let mut count = 0;
    for read in reads.iter() {
        rtgraph25.add_read(read, &samples, &mut count, &mut pending)
    }

    rtgraph25.build_graph_if_necessary(&mut pending);
    assert!(rtgraph25.has_cycles());
}

#[test]
fn test_Ns_in_reads_are_not_used_for_graph() {
    let length = 100;
    let reference = vec![b'A'; length];

    let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(25);
    let mut reads = Vec::new();
    let mut pending = LinkedHashMap::new();
    rtgraph.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_slice(),
        0,
        reference.len(),
        1,
        true,
    );

    let samples = vec!["anonymous_sample".to_string()];
    let quals = vec![30; length];
    let mut count = 0;
    // add reads with Ns at any position
    for i in 0..length {
        let mut bases = reference.clone();
        bases[i] = b'N';
        let read = ArtificialReadUtils::create_artificial_read(
            &bases,
            &quals,
            CigarString::try_from("100M").unwrap(),
        );
        reads.push(read);
    }
    for read in reads.iter() {
        rtgraph.add_read(read, &samples, &mut count, &mut pending);
    }
    rtgraph.build_graph_if_necessary(&mut pending);

    let mut graph = rtgraph.to_sequence_graph();
    let source = graph.base_graph.get_reference_source_vertex().unwrap();
    let sink = graph.base_graph.get_reference_sink_vertex().unwrap();
    let paths =
        GraphBasedKBestHaplotypeFinder::new_from_singletons(&mut graph.base_graph, source, sink)
            .find_best_haplotypes(std::usize::MAX, &graph.base_graph);
    assert_eq!(paths.len(), 1);
}

fn test_dangling_tails(
    ref_end: &str,
    alt_end: &str,
    cigar: &str,
    cigar_is_good: bool,
    merge_point_distance_from_sink: i32,
    num_leading_matches_allowed: i32,
) {
    println!(
        "params {} {} {} {} {} {}",
        ref_end,
        alt_end,
        cigar,
        cigar_is_good,
        merge_point_distance_from_sink,
        num_leading_matches_allowed
    );
    let kmer_size = 15;
    // let mut rng = ThreadRng::default();

    // construct the haplotypes
    let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let reference = common_prefix.to_string() + ref_end;
    let alternate = common_prefix.to_string() + alt_end;
    let quals = vec![30; alternate.len()];
    let read = ArtificialReadUtils::create_artificial_read(
        alternate.as_bytes(),
        &quals,
        CigarString::try_from(format!("{}M", alternate.len()).as_str()).unwrap(),
    );

    // create the graph and populate it
    let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    let mut pending = LinkedHashMap::new();
    rtgraph.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    // let reads = (0..reference.len())
    //     .map(|start| {
    //         (0..10)
    //             .map(|n| {
    //
    //                 generate_read_with_errors(
    //                     if rng.gen_range(0.0, 1.0) < 0.1 {
    //                         alternate.as_bytes()
    //                     } else {
    //                         reference.as_bytes()
    //                     },
    //                     start,
    //                     None,
    //                     0.001,
    //                     &mut rng,
    //                 )
    //             })
    //             .collect::<Vec<Vec<u8>>>()
    //     })
    //     .flat_map(|s| s)
    //     .collect::<Vec<Vec<u8>>>();
    //
    // reads.into_iter().for_each(|r| {
    //     rtgraph.add_sequence("anonymous".to_string(), 0, r.to_vec(), 0, r.len(), 1, false)
    // });

    
    let samples = vec!["anonymous".to_string()];
    let mut count = 0;
    rtgraph.add_read(&read, &samples, &mut count, &mut pending);
    rtgraph.set_min_matching_bases_to_dangling_end_recovery(num_leading_matches_allowed);
    rtgraph.build_graph_if_necessary(&mut pending);

    // confirm that we have just a single dangling tail
    let mut alt_sink = None;
    println!("reference path {:?}", &rtgraph.reference_path);
    for v in rtgraph.get_base_graph().graph.node_indices() {
        println!(
            "node {:?} node weight {:?} edges incoming {:?} outgoing {:?} is ref {}",
            v,
            rtgraph
                .get_base_graph()
                .graph
                .node_weight(v)
                .unwrap()
                .get_sequence_string(),
            rtgraph
                .get_base_graph()
                .edges_directed(v, Direction::Incoming),
            rtgraph
                .get_base_graph()
                .edges_directed(v, Direction::Outgoing),
            rtgraph.get_base_graph().is_reference_node(v)
        );
        if rtgraph.get_base_graph().is_sink(v) && !rtgraph.get_base_graph().is_reference_node(v) {
            assert!(
                alt_sink.is_none(),
                "We found more than one non-reference sink first alt sink {:?} other {:?}",
                rtgraph
                    .get_base_graph()
                    .graph
                    .node_weight(alt_sink.unwrap())
                    .unwrap()
                    .get_sequence_string(),
                rtgraph
                    .get_base_graph()
                    .graph
                    .node_weight(v)
                    .unwrap()
                    .get_sequence_string()
            );
            alt_sink = Some(v);
        }
    }

    assert!(alt_sink.is_some(), "We did not find a non-reference sink");

    // confirm that the SW alignment agrees with our expectations
    let result = rtgraph.generate_cigar_against_downwards_reference_path(
        alt_sink.unwrap(),
        0,
        4,
        false,
        &DANGLING_END_SW_PARAMETERS,
    );

    match result {
        None => {
            assert!(!cigar_is_good)
        }
        Some(result) => {
            assert_eq!(
                result.cigar.to_string(),
                cigar.to_string(),
                "SW generated cigar = {}",
                result.cigar
            );

            // confirm that the goodness of the cigar agrees with our expectations
            assert_eq!(
                ReadThreadingGraph::cigar_is_okay_to_merge(&result.cigar, false, true),
                cigar_is_good
            );

            // confirm that the tail merging works as expected
            if cigar_is_good {
                let merge_result = rtgraph.merge_dangling_tail(result);
                assert!(merge_result == 1 || merge_point_distance_from_sink == -1);

                // confirm that we created the appropriate edge
                if merge_point_distance_from_sink >= 0 {
                    let mut v = alt_sink.unwrap();
                    for _i in 0..merge_point_distance_from_sink {
                        if rtgraph.get_base_graph().in_degree_of(v) != 1 {
                            assert!(false, "Encountered vertex with multiple edges");
                        }
                        v = rtgraph
                            .get_base_graph()
                            .get_edge_source(rtgraph.get_base_graph().incoming_edge_of(v).unwrap());
                    }
                    assert!(rtgraph.get_base_graph().out_degree_of(v) > 1);
                }
            }
        }
    }
}

#[test]
fn make_dangling_tails_data() {
    // add 1M to the expected CIGAR because it includes the previous (common) base too
    test_dangling_tails("AAAAAAAAAA", "CAAA", "5M", true, 3, -1); // incomplete haplotype
    test_dangling_tails("AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10, -1); // insertion
    test_dangling_tails("CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10, -1); // deletion
    test_dangling_tails("AAAAAAAA", "CAAAAAAA", "9M", true, 7, -1); // 1 snp
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, 2, -1); // several snps
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, 2, 0); // several snps
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, 2, 1); // several snps
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, 2, 2); // several snps
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, -1, 3); // several snps (not enough bases to match)
    test_dangling_tails("AAAAAAAA", "CAAGATAA", "9M", true, -1, 4); // several snps (not enough bases to match)
    test_dangling_tails("AAAAA", "C", "1M4D1M", false, -1, -1); // funky SW alignment
    test_dangling_tails("AAAAA", "CA", "1M3D2M", false, 1, -1); // very little data
    test_dangling_tails("AAAAAAA", "CAAAAAC", "8M", true, -1, -1); // ends in mismatch
    test_dangling_tails("AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0, -1); // alignment is too complex
    test_dangling_tails("AAAAA", "YYYYY", "1M5I", false, -1, -1); // insertion
}

#[test]
fn test_forked_dangling_ends_with_suffix_code() {
    let kmer_size = 15;

    // construct the haplotypes
    let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let ref_end = "GCTAGCTAATCGTTAAGCTTTAAC";
    let alt_end1 = "GCTAGCTAAGGCG"; // two mismatches compared to the reference
    let alt_end2 = "GCTAGCTAAGCCGATGGCT";
    let reference = common_prefix.to_string() + ref_end;
    let alt1 = common_prefix.to_string() + alt_end1;
    let alt2 = common_prefix.to_string() + alt_end2;
    let quals1 = vec![30; alt1.len()];
    let quals2 = vec![30; alt2.len()];
    let read1 = ArtificialReadUtils::create_artificial_read(
        alt1.as_bytes(),
        &quals1,
        CigarString::try_from(format!("{}M", alt1.len()).as_str()).unwrap(),
    );
    let read2 = ArtificialReadUtils::create_artificial_read(
        alt2.as_bytes(),
        &quals2,
        CigarString::try_from(format!("{}M", alt2.len()).as_str()).unwrap(),
    );

    // create the graph and populate it
    let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    let mut pending = LinkedHashMap::new();
    rtgraph.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    
    let samples = vec!["anonymous".to_string()];

    rtgraph.set_min_matching_bases_to_dangling_end_recovery(1);
    let mut count = 0;
    rtgraph.add_read(&read2, &samples, &mut count, &mut pending);
    rtgraph.add_read(&read1, &samples, &mut count, &mut pending);

    rtgraph.build_graph_if_necessary(&mut pending);
    assert_eq!(rtgraph.get_base_graph().get_sinks().len(), 3);

    // Testing a degenerate case where the wrong "reference" path is selected and assuring that when
    for alt_sink in rtgraph.get_base_graph().get_sinks().into_iter() {
        if rtgraph.get_base_graph().is_reference_node(alt_sink)
            || rtgraph
                .get_base_graph()
                .graph
                .node_weight(alt_sink)
                .unwrap()
                .get_sequence_string()
                != "GCTAAGCCGATGGCT"
        {
            continue;
        } else {
            // confirm that the SW alignment agrees with our expectations
            let result = rtgraph.generate_cigar_against_downwards_reference_path(
                alt_sink,
                0,
                2,
                false,
                &DANGLING_END_SW_PARAMETERS,
            );
            assert!(result.is_some());
            let result = result.unwrap();
            assert!(ReadThreadingGraph::cigar_is_okay_to_merge(
                &result.cigar,
                false,
                true
            ));

            assert_eq!(
                min(
                    ReadThreadingGraph::longest_suffix_match(
                        result.reference_path_string.as_bytes(),
                        result.dangling_path_string.as_bytes(),
                        CigarUtils::get_reference_length(&result.cigar) as i64 - 1
                    ) as u32,
                    result.cigar.0[0].len()
                ),
                0
            );

            // confirm that the tail merging works as expected
            let merge_result = rtgraph.merge_dangling_tail(result);
            assert_eq!(merge_result, 0);
        }
    }
}

#[test]
fn test_forked_dangling_ends() {
    let kmer_size = 15;

    // construct the haplotypes
    let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let ref_end = "GCTAGCTAATCG";
    let alt_end1 = "ACTAGCTAATCG";
    let alt_end2 = "ACTAGATAATCG";
    let reference = common_prefix.to_string() + ref_end;
    let alt1 = common_prefix.to_string() + alt_end1;
    let alt2 = common_prefix.to_string() + alt_end2;
    let quals1 = vec![30; alt1.len()];
    let quals2 = vec![30; alt2.len()];
    let read1 = ArtificialReadUtils::create_artificial_read(
        alt1.as_bytes(),
        &quals1,
        CigarString::try_from(format!("{}M", alt1.len()).as_str()).unwrap(),
    );
    let read2 = ArtificialReadUtils::create_artificial_read(
        alt2.as_bytes(),
        &quals2,
        CigarString::try_from(format!("{}M", alt2.len()).as_str()).unwrap(),
    );

    // create the graph and populate it
    let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    let mut pending = LinkedHashMap::new();
    rtgraph.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    
    let samples = vec!["anonymous".to_string()];

    rtgraph.set_min_matching_bases_to_dangling_end_recovery(1);
    let mut count = 0;
    rtgraph.add_read(&read1, &samples, &mut count, &mut pending);
    rtgraph.add_read(&read2, &samples, &mut count, &mut pending);

    rtgraph.build_graph_if_necessary(&mut pending);
    assert_eq!(rtgraph.get_base_graph().get_sinks().len(), 3);

    for alt_sink in rtgraph.get_base_graph().get_sinks().into_iter() {
        if rtgraph.get_base_graph().is_reference_node(alt_sink) {
            continue;
        }

        // confirm that the SW alignment agrees with our expectations
        let result = rtgraph.generate_cigar_against_downwards_reference_path(
            alt_sink,
            0,
            4,
            true,
            &DANGLING_END_SW_PARAMETERS,
        );
        assert!(result.is_some());
        let result = result.unwrap();
        assert!(ReadThreadingGraph::cigar_is_okay_to_merge(
            &result.cigar,
            false,
            true
        ));

        // confirm that the tail merging works as expected
        let merge_result = rtgraph.merge_dangling_tail(result);
        assert_eq!(merge_result, 1);
    }

    let mut seq_graph = rtgraph.to_sequence_graph();
    seq_graph.simplify_graph("anon");
    let sources = seq_graph.base_graph.get_sources_generic().collect();
    let sinks = seq_graph.base_graph.get_sinks_generic().collect();
    let paths = GraphBasedKBestHaplotypeFinder::new(&mut seq_graph.base_graph, sources, sinks)
        .find_best_haplotypes(std::usize::MAX, &seq_graph.base_graph)
        .into_iter()
        .map(|k_best_haplotype| k_best_haplotype.path.get_bases(&seq_graph.base_graph))
        .unique()
        .sorted()
        .collect::<Vec<Vec<u8>>>();

    assert_eq!(paths.len(), 3);
    let mut expected = vec![
        reference.as_bytes().to_vec(),
        alt1.as_bytes().to_vec(),
        alt2.as_bytes().to_vec(),
    ];
    expected.sort();
    assert_eq!(paths, expected);
}

// #[test]
// fn test_whole_tail_is_insertion() {
//     let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(10);
//     let result = DanglingChainMergeHelper::new(
//
//     )
// }

#[test]
fn get_bases_for_path() {
    let kmer_size = 4;
    let test_string = "AATGGGGCAATACTA";

    let mut graph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    let mut pending = LinkedHashMap::new();
    graph.add_sequence(
        &mut pending,
        "anonymous".to_string(),
        0,
        test_string.as_bytes(),
        0,
        test_string.len(),
        1,
        true,
    );

    graph.build_graph_if_necessary(&mut pending);

    let mut vertices = Vec::new();
    let mut v = graph.get_reference_source_vertex();
    while v.is_some() {
        vertices.push(v.unwrap());
        v = graph
            .get_base_graph()
            .get_next_reference_vertex(v, false, None);
    }

    let result_for_tails = graph.get_bases_for_path(&vertices, false);
    assert_eq!(
        result_for_tails,
        std::str::from_utf8(&test_string.as_bytes()[kmer_size - 1..]).unwrap()
    );
    let result_for_heads = graph.get_bases_for_path(&vertices, true);
    assert_eq!(result_for_heads, "GTAAGGGCAATACTA".to_string()) // because the source node will be reversed
}

fn test_dangling_heads(
    reference: &str,
    alternate: &str,
    cigar: &str,
    should_be_merged: bool,
    num_leading_matches_allowed: i32,
) {
    println!(
        "params {} {} {} {} {}",
        reference, alternate, cigar, should_be_merged, num_leading_matches_allowed
    );
    let kmer_size = 5;
    let quals = vec![30; alternate.len()];
    let read = ArtificialReadUtils::create_artificial_read(
        alternate.as_bytes(),
        &quals,
        CigarString::try_from(format!("{}M", alternate.len()).as_str()).unwrap(),
    );

    // create the graph and populate it
    let mut rtgraph = ReadThreadingGraph::default_with_kmer_size(kmer_size);
    let mut pending = LinkedHashMap::new();
    rtgraph.add_sequence(
        &mut pending,
        "ref".to_string(),
        0,
        reference.as_bytes(),
        0,
        reference.len(),
        1,
        true,
    );

    
    let samples = vec!["anonymous".to_string()];
    let mut count = 0;
    rtgraph.add_read(&read, &samples, &mut count, &mut pending);
    rtgraph.set_min_matching_bases_to_dangling_end_recovery(num_leading_matches_allowed);
    rtgraph.build_graph_if_necessary(&mut pending);

    // confirm that we have just a single dangling tail
    let mut alt_source = None;
    println!("reference path {:?}", &rtgraph.reference_path);
    for v in rtgraph.get_base_graph().graph.node_indices() {
        println!(
            "node {:?} node weight {:?} edges incoming {:?} outgoing {:?} is ref {}",
            v,
            rtgraph
                .get_base_graph()
                .graph
                .node_weight(v)
                .unwrap()
                .get_sequence_string(),
            rtgraph
                .get_base_graph()
                .edges_directed(v, Direction::Incoming),
            rtgraph
                .get_base_graph()
                .edges_directed(v, Direction::Outgoing),
            rtgraph.get_base_graph().is_reference_node(v)
        );
        if rtgraph.get_base_graph().is_source(v) && !rtgraph.get_base_graph().is_reference_node(v) {
            assert!(
                alt_source.is_none(),
                "We found more than one non-reference sink first alt sink {:?} other {:?}",
                rtgraph
                    .get_base_graph()
                    .graph
                    .node_weight(alt_source.unwrap())
                    .unwrap()
                    .get_sequence_string(),
                rtgraph
                    .get_base_graph()
                    .graph
                    .node_weight(v)
                    .unwrap()
                    .get_sequence_string()
            );
            alt_source = Some(v);
        }
    }

    assert!(alt_source.is_some(), "We did not find a non-reference sink");

    // confirm that the SW alignment agrees with our expectations
    let result = rtgraph.generate_cigar_against_upwards_reference_path(
        alt_source.unwrap(),
        0,
        1,
        false,
        &DANGLING_END_SW_PARAMETERS,
    );

    match result {
        None => {
            assert!(!should_be_merged);
        }
        Some(result) => {
            assert_eq!(
                result.cigar.to_string(),
                cigar.to_string(),
                "SW generated cigar = {}",
                result.cigar
            );

            // confirm that the tail merging works as expected
            let merge_result = if num_leading_matches_allowed >= 0 {
                rtgraph.merge_dangling_head(result)
            } else {
                rtgraph.merge_dangling_head_legacy(result)
            };
            assert!(merge_result > 0 || !should_be_merged);

            // confirm that we created the appropriate bubble in the graph only if expected
            rtgraph.get_base_graph_mut().clean_non_ref_paths();
            let mut seq_graph = rtgraph.to_sequence_graph();
            let source = seq_graph.base_graph.get_reference_source_vertex().unwrap();
            let sink = seq_graph.base_graph.get_reference_sink_vertex().unwrap();
            let paths = GraphBasedKBestHaplotypeFinder::new_from_singletons(
                &mut seq_graph.base_graph,
                source,
                sink,
            )
            .find_best_haplotypes(std::usize::MAX, &seq_graph.base_graph);
            assert_eq!(paths.len(), if should_be_merged { 2 } else { 1 });

            // Asserting the contents of the path are reasonable
            let paths_converted = paths
                .into_iter()
                .map(|p| {
                    std::str::from_utf8(&p.path.get_bases(&seq_graph.base_graph))
                        .unwrap()
                        .to_string()
                })
                .collect::<Vec<String>>();
            assert!(paths_converted.contains(&reference.to_string()));
            if should_be_merged {
                assert!(paths_converted.iter().any(|p| p.contains(alternate)));
            }
        }
    }
}

#[test]
fn make_dangling_heads_data() {
    // add 1M to the expected CIGAR because it includes the last (common) base too
    test_dangling_heads("XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, -1); // 1 snp
    test_dangling_heads("XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 0); // 1 snp
    test_dangling_heads("XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 1); // 1 snp
    test_dangling_heads("XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 2); // 1 snp
    test_dangling_heads("XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", false, 3); // 1 snp

    // One SNP failing in legacy behavior (failing due to being too close to the reference start
    test_dangling_heads("YYYAACCGGTTACGT", "YAAACCGGTTACGT", "7M", false, -1); // 1 snp
    test_dangling_heads("YYYAACCGGTTACGT", "YAAACCGGTTACGT", "7M", false, 0); // 1 snp

    // Testing that indels now work
    test_dangling_heads("BBBBBBBAACCGGTTACGT", "BAACGGTTACGT", "4M1D4M", false, -1); // deletion
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "YAACGGTTACGT", "4M1D4M", true, 1); // deletion
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "YAACYCGGTTACGT", "5M1I4M", false, -1); // insertion
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "YAACYCGGTTACGT", "5M1I4M", true, 1); // insertion

    // Testing old/new behavior with multiple SNPs
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "AYYCGGTTACGT", "8M", false, -1); // 2 snps (only allow 1 SNP by default)
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "AYYCGGTTACGT", "8M", true, 1); // 2 snps
    test_dangling_heads("YYYYYYYAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", false, -1); // 2 snps (only allow 1 SNP by default)
    test_dangling_heads("YYYYYYYAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", true, 1); // 2 snps

    test_dangling_heads("YYYYYYYAACCGGTTACGT", "AYCGGTTACGT", "7M", true, -1); // very little data
    test_dangling_heads("YYYYYYYAACCGGTTACGT", "YCCGGTTACGT", "6M", true, -1); // begins in mismatch

    //Testing new behavior with indels long enough to confuse the code
    test_dangling_heads(
        "YYYYYYYAACBBBBBBCGGTTACGT",
        "AACCGGTTACGT",
        "5M6D3M",
        true,
        1,
    ); // long deletion
    test_dangling_heads(
        "YYYYYYYAACCBBBBBBGGTTACGT",
        "ATCCGGTTACGT",
        "5M6D4M",
        true,
        1,
    ); // long deletion plus close snp
    test_dangling_heads(
        "YYYYYYYAACCGGTTACGT",
        "YAACYYYYYYYCGGTTACGT",
        "5M7I4M",
        true,
        1,
    ); // 7 base insertion
    test_dangling_heads(
        "BBBBBBBAACCGGTTACGT",
        "BTACYYYYYYYCGGTTACGT",
        "5M7I4M",
        true,
        1,
    ); // 7 base with snp (NOTE: This triggers extendDanglingPathAgainstReference but should fail unless the cigar is used to predict the correct reference vertex to use for merging)
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

    result
}
