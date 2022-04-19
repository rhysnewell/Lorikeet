#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
extern crate rust_htslib;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;
extern crate hashlink;
extern crate ndarray;
extern crate ordered_float;
extern crate rand;

use hashlink::LinkedHashMap;
use lorikeet_genome::genotype::genotype_allele_counts::GenotypeAlleleCounts;
use lorikeet_genome::genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::haplotype::haplotype_caller_genotyping_engine::HaplotypeCallerGenotypingEngine;
use lorikeet_genome::model::allele_likelihood_matrix_mapper::AlleleLikelihoodMatrixMapper;
use lorikeet_genome::model::allele_likelihoods::{AlleleLikelihoods, LOG_10_INFORMATIVE_THRESHOLD};
use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reads::read_clipper::ReadClipper;
use lorikeet_genome::test_utils::read_clipper_test_utils::ReadClipperTestUtils;
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use ndarray::Array2;
use ordered_float::OrderedFloat;
use rand::distributions::Distribution;
use rand::distributions::Normal;
use rand::{rngs::ThreadRng, Error, Rng};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::ops::Deref;
use std::sync::Mutex;
use lorikeet_genome::model::allele_list::AlleleList;

lazy_static! {
    static ref READ_COUNTS: Vec<Vec<usize>> = vec![
        Vec::new(),
        vec![100],
        vec![0],
        vec![0, 0, 0],
        vec![1, 0, 1],
        vec![100, 10, 100],
        vec![1000, 10, 100, 20, 23],
    ];
    static ref SAMPLE_SETS: Vec<Vec<String>> = vec![
        vec![format!("A"), format!("B"), format!("C"),],
        vec![format!("A"),],
        vec![
            format!("C"),
            format!("A"),
            format!("D"),
            format!("E"),
            format!("Pepper"),
            format!("Pizza")
        ],
    ];
    static ref ALLELE_SETS: Vec<Vec<ByteArrayAllele>> = vec![
        vec![
            ByteArrayAllele::new("A".as_bytes(), true),
            ByteArrayAllele::new("T".as_bytes(), false),
            ByteArrayAllele::new("C".as_bytes(), false)
        ],
        vec![ByteArrayAllele::new("A".as_bytes(), true)],
        vec![
            ByteArrayAllele::new("ATTTA".as_bytes(), false),
            ByteArrayAllele::new("A".as_bytes(), true)
        ],
        vec![
            ByteArrayAllele::new("A".as_bytes(), false),
            ByteArrayAllele::new("AT".as_bytes(), true)
        ],
        vec![
            ByteArrayAllele::new("A".as_bytes(), false),
            ByteArrayAllele::new("AT".as_bytes(), false)
        ],
    ];
}

const EPSILON: f64 = 1e-6;
const ODD_READ_START: usize = 101;
const EVEN_READ_START: usize = 1;

fn test_instantiation_and_query(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut result = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());

    assert_eq!(result.number_of_samples(), samples.len());
    assert_eq!(result.number_of_alleles(), alleles.len());
    assert_eq!(result.samples().len(), samples.len());
    assert_eq!(result.alleles().len(), alleles.len());

    test_sample_queries(&samples, &reads, &result);
    test_allele_queries(&alleles, &result);
    test_likelihood_matrix_queries(&samples, &mut result, None)
}

fn test_likelihood_filling_and_query(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut result = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    let likelihoods = fill_with_random_likelihoods(&samples, &alleles, &mut result);
    test_likelihood_matrix_queries(&samples, &mut result, Some(&likelihoods));
}

fn fill_with_random_likelihoods(
    samples: &Vec<String>,
    alleles: &Vec<ByteArrayAllele>,
    result: &mut AlleleLikelihoods<ByteArrayAllele>,
) -> Vec<Array2<f64>> {
    let mut rnd = ThreadRng::default();
    let normal = Normal::new(0.0, 1.0);
    let mut likelihoods = Vec::new();
    for s in 0..samples.len() {
        let ncols = result.sample_evidence_count(s);
        let nrows = alleles.len();
        let mut likelihood = Array2::zeros((nrows, ncols));
        for a in 0..alleles.len() {
            for r in 0..ncols {
                let val = normal.sample(&mut rnd);
                result.set(s, a, r, val);
                likelihood[[a, r]] = val;
            }
        }
        likelihoods.push(likelihood);
    }

    return likelihoods;
}

fn fill_two_with_random_likelihoods(
    samples: &Vec<String>,
    alleles: &Vec<ByteArrayAllele>,
    result_1: &mut AlleleLikelihoods<ByteArrayAllele>,
    result_2: &mut AlleleLikelihoods<ByteArrayAllele>,
) -> Vec<Array2<f64>> {
    let mut rnd = ThreadRng::default();
    let normal = Normal::new(0.0, 1.0);
    let mut likelihoods = Vec::new();
    for s in 0..samples.len() {
        let ncols = result_1.sample_evidence_count(s);
        let nrows = alleles.len();
        let mut likelihood = Array2::zeros((nrows, ncols));
        for a in 0..alleles.len() {
            for r in 0..ncols {
                let val = normal.sample(&mut rnd);
                result_1.set(s, a, r, val);
                result_2.set(s, a, r, val);
                likelihood[[a, r]] = val;
            }
        }
        likelihoods.push(likelihood);
    }

    return likelihoods;
}

fn test_likelihood_matrix_queries(
    samples: &Vec<String>,
    result: &mut AlleleLikelihoods<ByteArrayAllele>,
    likelihoods: Option<&Vec<Array2<f64>>>,
) {
    let dummy = Array2::zeros((0, 0));
    for sample in samples.iter() {
        let index_of_sample = result.index_of_sample(sample);
        let sample_read_count = result.sample_evidence_count(index_of_sample);
        let number_of_alleles = result.number_of_alleles();
        assert_eq!(result.number_of_alleles(), number_of_alleles);
        for a in 0..number_of_alleles {
            assert_eq!(
                result.sample_evidence_count(index_of_sample),
                sample_read_count
            );
            for r in 0..sample_read_count {
                if result.sample_matrix(index_of_sample)[[a, r]].is_nan() {
                    match likelihoods {
                        None => assert!(false),
                        Some(likelihoods) => {
                            assert!(likelihoods[index_of_sample][[a, r]].is_nan());
                        }
                    }
                } else {
                    assert!(
                        relative_eq!(
                            result.sample_matrix(index_of_sample)[[a, r]],
                            match likelihoods {
                                None => 0.0,
                                Some(likelihoods) => {
                                    likelihoods[index_of_sample][[a, r]]
                                }
                            },
                            epsilon = EPSILON
                        ),
                        "result {:?}, likelihoods {:?}",
                        result.sample_matrix(index_of_sample),
                        match likelihoods {
                            Some(likelihoods) => &likelihoods[index_of_sample],
                            _ => &dummy,
                        }
                    );
                }
            }
        }
    }
}

fn test_allele_queries(
    alleles: &Vec<ByteArrayAllele>,
    result: &AlleleLikelihoods<ByteArrayAllele>,
) {
    let mut allele_indices = HashSet::new();
    for allele in alleles.iter() {
        let index_of_allele = result.index_of_allele(allele).unwrap();
        assert!(index_of_allele >= 0);
        assert!(!allele_indices.contains(&index_of_allele));
        allele_indices.insert(index_of_allele);
        assert_eq!(allele, &alleles[index_of_allele]);
    }
}

fn test_sample_queries(
    samples: &Vec<String>,
    reads: &HashMap<usize, Vec<BirdToolRead>>,
    result: &AlleleLikelihoods<ByteArrayAllele>,
) {
    let mut sample_ids = HashSet::new();
    for (idx, sample) in samples.iter().enumerate() {
        let index_of_sample = result.index_of_sample(sample);
        assert!(index_of_sample >= 0);
        assert!(!sample_ids.contains(&index_of_sample));
        sample_ids.insert(index_of_sample);

        let sample_reads = result.sample_evidence(index_of_sample).unwrap();
        let sample_reads_set = sample_reads
            .iter()
            .cloned()
            .collect::<HashSet<BirdToolRead>>();
        let expected_sample_read_array = reads.get(&idx).unwrap();
        let expected_sample_reads_set = expected_sample_read_array
            .iter()
            .cloned()
            .collect::<HashSet<BirdToolRead>>();
        assert_eq!(&sample_reads_set, &expected_sample_reads_set);

        let sample_read_count = sample_reads.len();
        for r in 0..sample_read_count {
            assert_eq!(&sample_reads[r], &expected_sample_read_array[r]);
            let read_index = result
                .evidence_index(index_of_sample, &sample_reads[r])
                .unwrap();
            assert_eq!(read_index, r);
        }
    }
}

fn test_best_alleles<A: Allele>(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    fill_with_random_likelihoods(&samples, &alleles, &mut original);
    let number_of_alleles = alleles.len();
    let ref_index = original.alleles().index_of_reference();
    let ref_allele = if ref_index.is_some() {
        Some(original.alleles().get_allele(ref_index.unwrap()).unwrap().clone())
    } else {
        None
    };

    let original_alleles = original.alleles().clone();

    for s in 0..samples.len() {
        let sample_read_count = original.sample_evidence_count(s);
        let best_alleles =
            original.best_alleles_breaking_ties_main(Box::new(|allele: &ByteArrayAllele| -> i32 {
                if allele.is_reference() {
                    1
                } else {
                    0
                }
            }));

        let mut sample_matrix = original.sample_matrix(s);
        let mut best_lk_array = vec![0.0; sample_read_count];
        let mut best_index_array = vec![0; sample_read_count];
        let mut confidence_array = vec![0.0; sample_read_count];
        for r in 0..sample_read_count {
            let mut best_index_of_allele = None;
            let mut best_allele_lk = std::f64::NEG_INFINITY;
            let mut second_best_allele_lk = std::f64::NEG_INFINITY;
            for a in 0..number_of_alleles {
                let lk = sample_matrix[[a, r]];
                if lk > best_allele_lk {
                    second_best_allele_lk = best_allele_lk;
                    best_allele_lk = lk;
                    best_index_of_allele = Some(a);
                } else if lk > second_best_allele_lk {
                    second_best_allele_lk = lk;
                };
            }

            best_lk_array[r] = best_allele_lk;
            confidence_array[r] = best_allele_lk - second_best_allele_lk;
            best_index_array[r] = best_index_of_allele.unwrap();
        }

        for best_allele in best_alleles {
            if s == best_allele.sample_index {
                let read_index = best_allele.evidence_index;
                if best_allele.allele_index.is_none() {
                    continue;
                } else {
                    let ref_likelihood = if ref_index.is_some() {
                        sample_matrix[[ref_index.unwrap(), read_index]]
                    } else {
                        std::f64::NEG_INFINITY
                    };

                    let ref_override = if ref_index.is_some() {
                        if ref_index.unwrap() != best_index_array[read_index]
                            && best_lk_array[read_index] - ref_likelihood
                                < *LOG_10_INFORMATIVE_THRESHOLD
                        {
                            true
                        } else {
                            false
                        }
                    } else {
                        false
                    };

                    assert!(
                        relative_eq!(
                        if ref_override {
                            ref_likelihood
                        } else {
                            best_lk_array[read_index]
                        }, best_allele.likelihood,
                        epsilon = EPSILON
                    ),
                        "ref override {} Ref likelihood {} best_lk_array {} best_allele likelihood {}",
                        ref_override, ref_likelihood, best_lk_array[read_index], best_allele.likelihood
                    );

                    let new_ref_allele = if ref_override {
                        ref_allele.clone().unwrap()
                    } else {
                        alleles[best_index_array[read_index]].clone()
                    };
                    assert_eq!(
                        original_alleles.get_allele(best_allele.allele_index.unwrap()).unwrap(),
                        &new_ref_allele
                    );

                    assert!(relative_eq!(
                        best_allele.confidence,
                        if ref_override {
                            ref_likelihood - best_lk_array[read_index]
                        } else {
                            confidence_array[read_index]
                        },
                        epsilon = EPSILON
                    ));
                }
            }
        }
    }
}

// Not implemented does not need to be tested
// fn test_best_allele_map(samples: Vec<String>, alleles: Vec<ByteArrayAllele>, reads: HashMap<usize, Vec<BirdToolRead>>) {
//     let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
//     fill_with_random_likelihoods(&samples, &alleles, &mut original);
//     let mut expected = HashMap::new();
//
//     let number_of_alleles = alleles.len();
//     for s in 0..samples.len() {
//         let sample_read_count = original.sample_evidence_count(s);
//         let sample_matrix = original.sample_matrix(s);
//         for r in 0..sample_read_count {
//             let mut best_index_of_allele = None;
//             let mut best_allele_lk = std::f64::NEG_INFINITY;
//             let mut second_best_allele_lk = std::f64::NEG_INFINITY;
//             for a in 0..number_of_alleles {
//                 let lk = sample_matrix[[a, r]];
//                 if lk > best_allele_lk {
//                     second_best_allele_lk = best_allele_lk;
//                     best_allele_lk = lk;
//                     best_index_of_allele = Some(a);
//                 } else if lk > second_best_allele_lk {
//                     second_best_allele_lk = lk;
//                 };
//             };
//             if (best_allele_lk - second_best_allele_lk) > LOG_10_INFORMATIVE_THRESHOLD {
//                 let allele_list = expected.entry(alleles[best_index_of_allele.unwrap()].clone()).or_insert(Vec::new());
//
//             }
//         }
//     }
// }

fn test_filer_poorly_modeled_reads(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut original = make_good_and_bad_likelihoods(
        &samples,
        &alleles,
        &reads,
        Box::new(|r: usize| (r & 1) == 0),
    );
    let mut result = make_good_and_bad_likelihoods(
        &samples,
        &alleles,
        &reads,
        Box::new(|r: usize| (r & 1) == 0),
    );
    // fill_two_with_random_likelihoods(&samples, &alleles, &mut original, &mut result);

    check_evidence_to_index_map_is_correct(&result);
    result.filter_poorly_modeled_evidence(Box::new(|read| -100.0));
    check_evidence_to_index_map_is_correct(&result);

    for s in 0..samples.len() {
        let old_sample_read_count = original.sample_evidence_count(s);
        let new_sample_read_count = result.sample_evidence_count(s);
        assert_eq!(new_sample_read_count, (old_sample_read_count + 1) / 2);
        assert_eq!(
            new_sample_read_count
                + result
                    .filtered_evidence_by_sample_index
                    .get(&s)
                    .unwrap()
                    .len(),
            old_sample_read_count
        ); // assert the correct number of reads were saved int the filtered pool

        for r in 0..new_sample_read_count {
            assert_eq!(
                original.evidence_index(s, &result.sample_evidence(s).unwrap()[r]),
                Some(r * 2)
            );
            let new_sample_matrix = result.sample_matrix(s);
            let old_sample_matrix = original.sample_matrix(s);
            for a in 0..alleles.len() {
                assert_eq!(new_sample_matrix[[a, r]], old_sample_matrix[[a, r * 2]]);
            }
        }
    }
}

fn test_filter_reads_to_overlap(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    let mut result = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());

    fill_two_with_random_likelihoods(&samples, &alleles, &mut original, &mut result);

    let even_read_overlap = SimpleInterval::new(0, EVEN_READ_START, EVEN_READ_START);
    let read_qualifies_for_genotyping_predicate =
        HaplotypeCallerGenotypingEngine::compose_read_qualifies_for_genotyping_predicate();
    check_evidence_to_index_map_is_correct(&result);
    result.retain_evidence(&read_qualifies_for_genotyping_predicate, &even_read_overlap);
    check_evidence_to_index_map_is_correct(&result);

    let mut new_likelihoods = Vec::with_capacity(samples.len());
    for s in 0..samples.len() {
        let mut new_likelihood =
            Array2::zeros((alleles.len(), (original.sample_evidence_count(s) + 1) / 2));
        for a in 0..alleles.len() {
            for r in 0..new_likelihood.ncols() {
                assert_eq!(
                    result
                        .evidence_index(s, &original.sample_evidence(s).unwrap()[r << 1])
                        .unwrap(),
                    r
                );
                let mut sample_matrix = original.sample_matrix(s);
                new_likelihood[[a, r]] = sample_matrix[[a, r << 1]];
            }
        }
        new_likelihoods.push(new_likelihood);
    }

    test_likelihood_matrix_queries(&samples, &mut result, Some(&new_likelihoods));
}

fn test_filter_poorly_modeled_reads_to_overlap(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut original = make_good_and_bad_likelihoods(
        &samples,
        &alleles,
        &reads,
        Box::new(|r: usize| (r & 2) == 0),
    );
    let mut result = make_good_and_bad_likelihoods(
        &samples,
        &alleles,
        &reads,
        Box::new(|r: usize| (r & 2) == 0),
    );
    // fill_two_with_random_likelihoods(&samples, &alleles, &mut original, &mut result);
    let even_read_overlap = SimpleInterval::new(0, EVEN_READ_START, EVEN_READ_START);
    let read_qualifies_for_genotyping_predicate =
        HaplotypeCallerGenotypingEngine::compose_read_qualifies_for_genotyping_predicate();

    result.filter_poorly_modeled_evidence(Box::new(|_| -100.0));

    check_evidence_to_index_map_is_correct(&result);
    result.retain_evidence(&read_qualifies_for_genotyping_predicate, &even_read_overlap);
    check_evidence_to_index_map_is_correct(&result);

    let mut new_likelihoods = Vec::with_capacity(samples.len());

    for s in 0..samples.len() {
        let mut new_likelihood = Array2::zeros((
            alleles.len(),
            (((original.sample_evidence_count(s) + 1) / 2) + 1) / 2,
        ));

        assert_eq!(
            result.sample_evidence_count(s)
                + result
                    .filtered_evidence_by_sample_index
                    .get(&s)
                    .unwrap()
                    .len(),
            (original.sample_evidence_count(s) + 1) / 2
        ); // assert the correct number of reads were saved int the filtered pool

        for a in 0..alleles.len() {
            for r in 0..new_likelihood.ncols() {
                if r % 2 == 0 {
                    assert_eq!(
                        result
                            .evidence_index(s, &original.sample_evidence(s).unwrap()[r << 2])
                            .unwrap(),
                        r
                    );
                    let sample_matrix = original.sample_matrix(s);
                    new_likelihood[[a, r]] = sample_matrix[[a, r << 2]];
                } else {
                    assert_eq!(
                        result.evidence_index(s, &original.sample_evidence(s).unwrap()[r << 1]),
                        None
                    );
                }
            }
        }
        new_likelihoods.push(new_likelihood);
    }

    test_likelihood_matrix_queries(&samples, &mut result, Some(&new_likelihoods));
}

fn make_good_and_bad_likelihoods<A: Allele>(
    samples: &Vec<String>,
    alleles: &Vec<A>,
    reads: &HashMap<usize, Vec<BirdToolRead>>,
    reads_to_skip: Box<Fn(usize) -> bool>,
) -> AlleleLikelihoods<A> {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    for s in 0..samples.len() {
        let sample_read_count = original.sample_evidence_count(s);
        for r in 0..sample_read_count {
            if reads_to_skip(r) {
                continue;
            } else {
                for a in 0..alleles.len() {
                    original.sample_matrix(s)[[a, r]] = -10000.0;
                }
            };
        }
    }

    return original;
}

// Make sure that operations that add or remove evidence result in a correct evidence-to-index map
// Some such operations update the map; others invalidate it and leave it to regenerated lazily the next time it is queried.
// Either way, we want the queries after the mutating operations to be correct
fn check_evidence_to_index_map_is_correct<A: Allele>(subject: &AlleleLikelihoods<A>) {
    // test that evidence-to-index cache is valid after adding reads
    for s in 0..subject.number_of_samples() {
        for r in 0..subject.sample_evidence_count(s) {
            assert_eq!(
                subject
                    .evidence_index(s, &subject.sample_evidence(s).unwrap()[r])
                    .unwrap(),
                r
            );
        }
    }
}

fn test_marginalization_with_overlap(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    to_alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
    new_to_old_allele_mapping: LinkedHashMap<usize, Vec<&ByteArrayAllele>>,
) {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    let even_read_overlap = SimpleInterval::new(0, EVEN_READ_START, EVEN_READ_START);
    fill_with_random_likelihoods(&samples, &alleles, &mut original);

    let mut marginalized = original.marginalize(&new_to_old_allele_mapping);
    let read_qualifies_for_genotyping_predicate =
        HaplotypeCallerGenotypingEngine::compose_read_qualifies_for_genotyping_predicate();

    marginalized.retain_evidence(&read_qualifies_for_genotyping_predicate, &even_read_overlap);
    assert_eq!(
        new_to_old_allele_mapping.len(),
        marginalized.number_of_alleles()
    );
    let marginalized_alleles = marginalized.alleles().clone();

    for allele in marginalized_alleles.list.iter() {
        let a = marginalized_alleles.index_of_allele(allele).unwrap();
        let old_alleles = new_to_old_allele_mapping.get(&a).unwrap();
        for s in 0..samples.len() {
            let old_sample_likelihoods = original.sample_matrix(s).clone();
            let sample_likelihoods = marginalized.sample_matrix(s).clone();

            let old_sample_read_count = old_sample_likelihoods.ncols();
            let sample_read_count = marginalized.sample_evidence_count(s);
            assert_eq!(sample_read_count, (old_sample_read_count + 1) / 2);
            let sample_likelihoods = marginalized.sample_matrix(s);
            for r in 0..sample_read_count {
                let mut old_best_lk = std::f64::NEG_INFINITY;
                for old_allele in old_alleles {
                    old_best_lk = max(
                        OrderedFloat(
                            old_sample_likelihoods[[
                                original.alleles().index_of_allele(*old_allele).unwrap(),
                                r << 1,
                            ]],
                        ),
                        OrderedFloat(old_best_lk),
                    )
                    .into();
                }

                assert_eq!(
                    sample_likelihoods[[a, r]],
                    old_best_lk,
                    "Read likelihoods {:?} Allele index {}, read index {} old likelihoods {:?}",
                    sample_likelihoods.row(a),
                    a,
                    r,
                    old_sample_likelihoods.column(r)
                );
            }
        }
    }
}

fn test_marginalization(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    to_alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
    new_to_old_allele_mapping: LinkedHashMap<usize, Vec<&ByteArrayAllele>>,
) {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    let even_read_overlap = SimpleInterval::new(0, EVEN_READ_START, EVEN_READ_START);
    fill_with_random_likelihoods(&samples, &alleles, &mut original);

    println!(
        "onew to old allele mapping {:?}",
        &new_to_old_allele_mapping
    );
    let mut marginalized = original.marginalize(&new_to_old_allele_mapping, AlleleList::new(&to_alleles));

    assert_eq!(
        new_to_old_allele_mapping.len(),
        marginalized.number_of_alleles()
    );
    let marginalized_alleles = marginalized.alleles().clone();

    for allele in marginalized_alleles.list.iter() {
        let a = marginalized_alleles.index_of_allele(allele).unwrap();
        let old_alleles = new_to_old_allele_mapping.get(&a).unwrap();
        for s in 0..samples.len() {
            let old_sample_likelihoods = original.sample_matrix(s).clone();
            let sample_likelihoods = marginalized.sample_matrix(s).clone();

            let old_sample_read_count = old_sample_likelihoods.ncols();
            let sample_read_count = marginalized.sample_evidence_count(s);
            assert_eq!(sample_read_count, old_sample_read_count);
            let sample_likelihoods = marginalized.sample_matrix(s);
            for r in 0..sample_read_count {
                let mut old_best_lk = std::f64::NEG_INFINITY;
                for old_allele in old_alleles {
                    old_best_lk = max(
                        OrderedFloat(
                            old_sample_likelihoods
                                [[original.alleles().index_of_allele(*old_allele).unwrap(), r]],
                        ),
                        OrderedFloat(old_best_lk),
                    )
                    .into();
                }

                assert_eq!(
                    sample_likelihoods[[a, r]],
                    old_best_lk,
                    "Allele index {}, read index {} old likelihoods {:?} new likelihoods {:?}",
                    a,
                    r,
                    old_sample_likelihoods.column(r),
                    sample_likelihoods.column(r)
                );
            }
        }
    }
}

fn test_normalize_cap_worst_lk(
    samples: Vec<String>,
    alleles: Vec<ByteArrayAllele>,
    reads: HashMap<usize, Vec<BirdToolRead>>,
) {
    let mut original = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());
    let mut result = AlleleLikelihoods::new(alleles.clone(), samples.clone(), reads.clone());

    let original_likelihoods =
        fill_two_with_random_likelihoods(&samples, &alleles, &mut original, &mut result);
    result.normalize_likelihoods(-0.001, true);
    test_allele_queries(&alleles, &result);

    let number_of_alleles = alleles.len();
    let mut new_likelihoods = Vec::new();
    for s in 0..samples.len() {
        let sample_read_count = original.sample_evidence_count(s);
        let mut new_likelihood = Array2::zeros((number_of_alleles, sample_read_count));

        for r in 0..sample_read_count {
            let mut best_allele_lk = std::f64::NEG_INFINITY;
            //best likelihood can be alt OR reference
            for a in 0..number_of_alleles {
                best_allele_lk = max(
                    OrderedFloat(best_allele_lk),
                    OrderedFloat(original_likelihoods[s][[a, r]]),
                )
                .into();
            }
            if best_allele_lk == std::f64::NEG_INFINITY {
                for a in 0..number_of_alleles {
                    new_likelihood[[a, r]] = original_likelihoods[s][[a, r]];
                }
            } else {
                for a in 0..number_of_alleles {
                    new_likelihood[[a, r]] = max(
                        OrderedFloat(best_allele_lk - 0.001),
                        OrderedFloat(original_likelihoods[s][[a, r]]),
                    )
                    .into();
                }
            }
        }
        new_likelihoods.push(new_likelihood);
    }
    test_likelihood_matrix_queries(&samples, &mut result, Some(&new_likelihoods));
}

#[test]
fn data_for_test_instantiation_and_query() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_instantiation_and_query(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_likelihood_filling_and_query() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_likelihood_filling_and_query(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_best_alleles() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_best_alleles::<ByteArrayAllele>(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_filter_poorly_modeled_evidence() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_filer_poorly_modeled_reads(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_filter_reads_to_overlap() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_filter_reads_to_overlap(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_filter_poorly_modeled_reads_to_overlap() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_filter_poorly_modeled_reads_to_overlap(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_normalize_cap_worst_lk() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set in ALLELE_SETS.iter() {
            test_normalize_cap_worst_lk(
                sample_set.clone(),
                allele_set.clone(),
                data_set_reads(sample_set, &mut rnd),
            )
        }
    }
}

#[test]
fn data_for_test_marginalization_with_overlap() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set_1 in ALLELE_SETS.iter() {
            for allele_set_2 in ALLELE_SETS.iter() {
                if allele_set_2.len() < allele_set_1.len() {
                    test_marginalization_with_overlap(
                        sample_set.clone(),
                        allele_set_1.clone(),
                        allele_set_2.clone(),
                        data_set_reads(sample_set, &mut rnd),
                        random_allele_map(&allele_set_1, &allele_set_2),
                    );
                }
            }
        }
    }
}

#[test]
fn data_for_test_marginalization() {
    let mut rnd = ThreadRng::default();
    for sample_set in SAMPLE_SETS.iter() {
        for allele_set_1 in ALLELE_SETS.iter() {
            for allele_set_2 in ALLELE_SETS.iter() {
                if allele_set_2.len() < allele_set_1.len() {
                    println!(
                        "samples {:?} alleles1 {:?} alleles2 {:?}",
                        sample_set, allele_set_1, allele_set_2
                    );
                    test_marginalization(
                        sample_set.clone(),
                        allele_set_1.clone(),
                        allele_set_2.clone(),
                        data_set_reads(sample_set, &mut rnd),
                        random_allele_map(&allele_set_1, &allele_set_2),
                    );
                }
            }
        }
    }
}

fn data_set_reads(samples: &Vec<String>, rnd: &mut ThreadRng) -> HashMap<usize, Vec<BirdToolRead>> {
    let mut result = HashMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        let read_count = rnd.gen_range(0, 100);
        let mut reads = Vec::new();
        for r in 0..read_count {
            let alignment_start = if (r & 1) == 0 {
                EVEN_READ_START
            } else {
                ODD_READ_START
            };
            reads.push(
                ArtificialReadUtils::create_artificial_read_with_name_and_pos(
                    format!("RRR{}00{}", sample, r),
                    0,
                    alignment_start as i64,
                    "AAAAA".as_bytes(),
                    vec![30, 30, 30, 30, 30].as_slice(),
                    "5M",
                    idx,
                ),
            )
        }
        result.insert(idx, reads);
    }

    return result;
}

fn random_allele_map<'a, A: Allele>(
    from_alleles: &'a Vec<A>,
    to_alleles: &Vec<A>,
) -> LinkedHashMap<usize, Vec<&'a A>> {
    let mut result = LinkedHashMap::new();
    let mut remaining = from_alleles.iter().collect::<Vec<&'a A>>();
    let mut rnd = ThreadRng::default();
    let mut next_to_index = 0;
    for (_, _) in from_alleles.iter().enumerate() {
        let from_index_of_allele = rnd.gen_range(0, remaining.len());
        let result_entry = result.entry(next_to_index).or_insert(Vec::new());
        result_entry.push(remaining.remove(from_index_of_allele));
        next_to_index = (next_to_index + 1) % to_alleles.len();
    }

    return result;
}
