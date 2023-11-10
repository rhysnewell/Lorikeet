#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::assembly::assembly_result_set::AssemblyResultSet;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::{
    AVXMode, PCRErrorModel, PairHMMLikelihoodCalculationEngine,
};
use lorikeet_genome::processing::lorikeet_engine::ReadType;
use lorikeet_genome::read_threading::read_threading_graph::ReadThreadingGraph;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::quality_utils::QualityUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use std::collections::HashMap;

#[test]
fn test_compute_likelihoods() {
    let mut lce = PairHMMLikelihoodCalculationEngine::new(
        93,
        MathUtils::log_to_log10(QualityUtils::qual_to_error_prob_log10(45)),
        PCRErrorModel::Conservative,
        16,
        false,
        1.0,
        0.02,
        true,
        false,
        true,
        AVXMode::None,
    );

    let mut per_sample_read_list = HashMap::new();
    let n = 10;
    let mut read1 = BirdToolRead::new(
        ArtificialReadUtils::create_artificial_read_default("test", 0, 0, 10, false),
        0,
        ReadType::Short,
    );
    read1.read.set_mapq(60);
    let sample1 = 0;
    per_sample_read_list.insert(0, vec![read1.clone()]);

    let sample = vec![sample1];
    let ref_bases = vec![b'A'; n + 1];
    let mut hap1 = Haplotype::new(ref_bases.as_slice(), true);
    hap1.set_genome_location(SimpleInterval::new(
        read1.get_contig(),
        read1.get_start(),
        read1.get_end(),
    ));
    let mut assembly_result_set = AssemblyResultSet::<ReadThreadingGraph>::new(
        AssemblyRegion::new(SimpleInterval::new(0, 0, n + 1), true, 0, 100, 0, 0, 0.0),
        vec![b'A'; n + 1],
        SimpleInterval::new(0, 0, n + 1),
        hap1.clone(),
    );
    assembly_result_set.add_haplotype(hap1.clone());

    let mut bases_modified = ref_bases;
    bases_modified[5] = b'C';
    let mut hap2 = Haplotype::new(bases_modified.as_slice(), false);
    hap2.set_genome_location(SimpleInterval::new(
        read1.get_contig(),
        read1.get_start(),
        read1.get_end(),
    ));
    assembly_result_set.add_haplotype(hap2.clone());
    let mut likes =
        lce.compute_read_likelihoods(&mut assembly_result_set, sample, per_sample_read_list);

    assert_eq!(likes.alleles().len(), 2);
    assert_eq!(likes.evidence_count(), 1);

    let v1 = likes.sample_matrix(0)[[0, 0]];
    let v2 = likes.sample_matrix(0)[[1, 0]];

    assert!(
        v1 > v2,
        "Matching hap should have a higher likelihood {} -> {}",
        v1,
        v2
    );
}
