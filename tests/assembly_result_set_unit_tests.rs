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
extern crate rand;
extern crate term;

use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::assembly::assembly_result::{AssemblyResult, Status};
use lorikeet_genome::assembly::assembly_result_set::AssemblyResultSet;
use lorikeet_genome::graphs::base_edge::BaseEdgeStruct;
use lorikeet_genome::graphs::seq_graph::SeqGraph;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::byte_array_allele::Allele;
use lorikeet_genome::read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use lorikeet_genome::read_threading::read_threading_graph::ReadThreadingGraph;
use lorikeet_genome::reads::cigar_builder::CigarBuilder;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::test_utils::random_dna::RandomDNA;
use lorikeet_genome::utils::errors::BirdToolError::CigarBuilderError;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rand::rngs::ThreadRng;
use rust_htslib::bam::record::Cigar;
use std::collections::HashMap;

// #[test]
// fn test_add_reference_haplotype() {
//     let mut reference = Haplotype::new("ACGT".as_bytes(), true);
//     // let mut active_region =
//     //     AssemblyRegion::new(SimpleInterval::new(0, 0, reference.len()), true, 25, 1000, 0, 0);
//     reference.set_genome_location(SimpleInterval::new(0, 0, reference.len()));
//     let mut subject = AssemblyResultSet::<ReadThreadingGraph>::default();
//     assert!(subject.add_haplotype(reference.clone()));
//     assert!(!subject.add_haplotype(reference.clone()));
//
//     assert_eq!(subject.get_ref_haplotype(), &reference);
//     assert_eq!(subject.get_haplotypes().len(), 1);
// }


fn test_trim_to(
    haplotype_and_result_sets: HashMap<
        Haplotype<SimpleInterval>,
        AssemblyResult<SimpleInterval, ReadThreadingGraph>,
    >,
    original: AssemblyRegion,
    ref_hap: Haplotype<SimpleInterval>,
    reference: Vec<u8>,
) {
    let n = 100;
    let mut subject = AssemblyResultSet::<ReadThreadingGraph>::new(
        original.clone(),
        reference,
        SimpleInterval::new(0, 0, n + 1),
        ref_hap,
    );
    for (hap, result) in haplotype_and_result_sets.clone().into_iter() {
        let result_index = subject.add_assembly_result(result);
        subject.add_haplotype_and_assembly_result(hap, result_index);
    }
    let original_location = original.get_padded_span();
    let length = original_location.size();
    let new_location = SimpleInterval::new(
        original_location.get_contig(),
        original_location.get_start() + (length) / 2,
        original_location.get_end() - (length) / 2,
    );
    let mut new_region = original.trim_with_padded_span(new_location.clone(), new_location.clone());

    let mut original_haplotypes_by_trimmed = HashMap::new();
    for h in haplotype_and_result_sets.keys().cloned() {
        original_haplotypes_by_trimmed
            .insert(h.trim(new_region.get_padded_span()).unwrap().unwrap(), h);
    }

    let trimmed = subject.trim_to(new_region);
    for h in trimmed.get_haplotype_list() {
        assert_eq!(h.get_genome_location().unwrap(), &new_location);
        assert_eq!(h.get_bases().len(), new_location.size());
    }
}

#[test]
fn trimming_data() {
    let mut rng = ThreadRng::default();
    let mut active_region =
        AssemblyRegion::new(SimpleInterval::new(0, 1000, 1100), true, 25, 1000000, 0, 0);
    let length = active_region.get_padded_span().size();
    let mut rnd = RandomDNA::new(rng);
    let reference = rnd.next_bases(length);
    let mut reference_bases = reference.as_bytes().to_vec();
    let mut ref_hap = Haplotype::new(reference.as_bytes(), true);
    ref_hap.set_genome_location(active_region.get_padded_span());
    ref_hap.set_cigar(vec![Cigar::Match(reference_bases.len() as u32)]);

    let mut haplotype_list = Vec::new();
    haplotype_list.push(ref_hap.clone());

    for test_hap in 0..=10 {
        let mut new_bases = reference_bases.clone();
        // if test_hap % 100 == 0 {
        //     new_bases.insert(test_hap, b'A');
        //     let mut hap = Haplotype::new(new_bases.as_slice(), false);
        //     hap.set_genome_location(active_region.get_padded_span());
        //     haplotype_list.push(hap);
        // } else if test_hap % 50 == 0 {
        //     // add deletion
        //     new_bases.remove(test_hap);
        //     let mut hap = Haplotype::new(new_bases.as_slice(), false);
        //     hap.set_genome_location(active_region.get_padded_span());
        //     haplotype_list.push(hap);
        // } else {
        if new_bases[test_hap] == b'A' {
            new_bases[test_hap] = b'C'
        } else {
            new_bases[test_hap] = b'A'
        };
        let mut hap = Haplotype::new(new_bases.as_slice(), false);
        hap.set_genome_location(active_region.get_padded_span());
        hap.set_cigar(vec![Cigar::Match(new_bases.len() as u32)]);
        haplotype_list.push(hap);
        // }
    }

    let seq_graph = SeqGraph::<BaseEdgeStruct>::new(10);
    let ar = AssemblyResult::<SimpleInterval, ReadThreadingGraph>::new(
        Status::AssembledSomeVariation,
        Some(seq_graph),
        None,
    );
    let mut result = HashMap::new();
    for h in haplotype_list {
        result.insert(h, ar.clone());
    }

    test_trim_to(result, active_region, ref_hap, reference_bases)
}
