#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;


use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::haplotype::homogenous_ploidy_model::HeterogeneousPloidyModel;
use lorikeet_genome::haplotype::independent_samples_genotype_model::IndependentSamplesGenotypesModel;
use lorikeet_genome::model::allele_list::AlleleList;
use lorikeet_genome::model::byte_array_allele::Allele;

use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use rand::rngs::ThreadRng;
use rand::Rng;

fn test_calculate_likelihoods(
    ploidies: &[usize],
    allele_count: usize,
    discard_allele_count: usize,
    read_counts: &[usize],
) {
    println!(
        "Ploidies {:?}, allele count {}, discard count {}, read_counts {:?}",
        ploidies, allele_count, discard_allele_count, read_counts
    );
    let likelihoods = ReadLikelihoodsUnitTester::read_likelihoods(allele_count, read_counts);
    let genotyping_allele_list = if discard_allele_count == 0 {
        likelihoods.get_allele_list()
    } else {
        discard_alleles_at_random(likelihoods.get_allele_list(), discard_allele_count)
    };
    let sample_list = vec!["sample".to_string(); ploidies.len()];
    let ploidy_model = HeterogeneousPloidyModel::new(sample_list.clone(), ploidies.to_vec());
    let mut model = IndependentSamplesGenotypesModel::default();
    let g_likelihoods = model.calculate_likelihoods(
        &genotyping_allele_list,
        likelihoods.get_allele_list_byte_array(),
        &likelihoods,
        &ploidy_model,
        b"",
        0,
    );
    // assert!(!g_likelihoods.is_empty());
    // AlleleListUnitTester::assert_allele_list(g_likelihoods)
    let sample_count = sample_list.len();
    for i in 0..sample_count {
        let sample_likelihoods = &g_likelihoods[i];
        println!(
            "sample count {} sample_likelihoods {:?}",
            i, &sample_likelihoods
        );

        let values = sample_likelihoods.clone().get_as_vector();
        assert!(!values.is_empty());
        assert_eq!(
            values.len(),
            GenotypeLikelihoodCalculators::get_instance(
                ploidies[i],
                genotyping_allele_list.number_of_alleles()
            )
            .genotype_count as usize
        );
        for j in 0..values.len() {
            assert!(values[j] <= 0.0);
        }
    }
}

fn discard_alleles_at_random<A: Allele>(
    likelihoods: AlleleList<A>,
    discard_allele_count: usize,
) -> AlleleList<A> {
    let mut rnd = ThreadRng::default();
    let mut subset = likelihoods.list.into_iter().collect::<Vec<A>>();
    for _i in 0..discard_allele_count {
        subset.remove(rnd.gen_range(0, subset.len() - 1));
    }

    AlleleList::new(&subset)
}

lazy_static! {
    static ref ALLELE_COUNTS: Vec<Vec<usize>> = vec![
        vec![1, 0],
        vec![2, 1],
        vec![5, 2],
        vec![10, 4],
        vec![1, 0],
        vec![2, 1],
        vec![10, 7],
    ];
    static ref PLOIDIES: Vec<Vec<usize>> = vec![
        vec![1, 1, 1, 1],
        vec![1, 2, 3, 4],
        vec![2, 2, 2, 2],
        vec![2, 1, 2, 1],
        vec![1],
        vec![2],
        Vec::new(),
    ];
    static ref READ_COUNTS: Vec<Vec<usize>> = vec![
        vec![10, 100, 50, 20],
        vec![0, 100, 10, 1],
        vec![1, 2, 3, 4],
        vec![10, 20, 50, 40],
        vec![10],
        vec![20],
        Vec::new()
    ];
}

#[test]
fn ploidy_and_maximum_allele_and_read_counts_data() {
    for i in 0..PLOIDIES.len() {
        test_calculate_likelihoods(&PLOIDIES[i], ALLELE_COUNTS[i][0], 0, &READ_COUNTS[i]);
        let discard_allele_count = ALLELE_COUNTS[i][1];
        if discard_allele_count == 0 {
            continue;
        } else {
            test_calculate_likelihoods(
                &PLOIDIES[i],
                ALLELE_COUNTS[i][0],
                ALLELE_COUNTS[i][1],
                &READ_COUNTS[i],
            );
        }
    }
}
