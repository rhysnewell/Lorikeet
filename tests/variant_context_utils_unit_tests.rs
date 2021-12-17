#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate lorikeet_genome;
#[macro_use]
extern crate lazy_static;
extern crate hashlink;

use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context;
use lorikeet_genome::model::variant_context::{VariantContext, VariantType};
use lorikeet_genome::model::variant_context_utils;
use lorikeet_genome::model::variant_context_utils::{VariantContextUtils, FilteredRecordMergeType, GenotypeMergeType};
use lorikeet_genome::utils::simple_interval::Locatable;
use std::collections::vec_deque::VecDeque;
use lorikeet_genome::genotype::genotype_builder::{Genotype, GenotypesContext};
use std::collections::HashSet;
use hashlink::LinkedHashSet;

#[test]
fn test_find_number_of_repetitions() {
    /*
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     */
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "GATAT".as_bytes(), false),
        2
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "GATAT".as_bytes(), true),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("A".as_bytes(), "ATATG".as_bytes(), true),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "ATATG".as_bytes(), true),
        2
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "CCC".as_bytes(),
            "CCCCCCCC".as_bytes(),
            true
        ),
        2
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "CCC".as_bytes(),
            "CCCCCCCC".as_bytes(),
            false
        ),
        2
    );

    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "ATG".as_bytes(),
            "ATGATGATGATG".as_bytes(),
            true
        ),
        4
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "G".as_bytes(),
            "ATGATGATGATG".as_bytes(),
            true
        ),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("T".as_bytes(), "T".as_bytes(), true),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "AT".as_bytes(),
            "ATGATGATCATG".as_bytes(),
            true
        ),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "CCCCCCCC".as_bytes(),
            "CCC".as_bytes(),
            true
        ),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "AT".as_bytes(), true),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "".as_bytes(), true),
        0
    ); //empty test string

    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "ATG".as_bytes(),
            "ATGATGATGATG".as_bytes(),
            false
        ),
        4
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "G".as_bytes(),
            "ATGATGATGATG".as_bytes(),
            false
        ),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("T".as_bytes(), "T".as_bytes(), false),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "AT".as_bytes(),
            "ATGATGATCATG".as_bytes(),
            false
        ),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions(
            "CCCCCCCC".as_bytes(),
            "CCC".as_bytes(),
            false
        ),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "AT".as_bytes(), false),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions("AT".as_bytes(), "".as_bytes(), false),
        0
    ); //empty test string
}

#[test]
fn test_find_number_of_repetitions_full_array() {
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "XXXATG".as_bytes(),
            3,
            3,
            "ATGATGATGATGYYY".as_bytes(),
            0,
            12,
            true
        ),
        4
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "GGGG".as_bytes(),
            0,
            1,
            "GGGGATGATGATGATG".as_bytes(),
            4,
            12,
            true
        ),
        0
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "T".as_bytes(),
            0,
            1,
            "TTTTT".as_bytes(),
            0,
            1,
            true
        ),
        1
    );

    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "AT".as_bytes(),
            0,
            2,
            "AT".as_bytes(),
            0,
            0,
            true
        ),
        0
    ); //empty test string
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "AT".as_bytes(),
            0,
            2,
            "AT".as_bytes(),
            1,
            0,
            true
        ),
        0
    ); //empty test string
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "AT".as_bytes(),
            0,
            2,
            "".as_bytes(),
            0,
            0,
            true
        ),
        0
    ); //empty test string

    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "XXXAT".as_bytes(),
            3,
            2,
            "XXXGATAT".as_bytes(),
            4,
            4,
            false
        ),
        2
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "AT".as_bytes(),
            0,
            2,
            "GATAT".as_bytes(),
            0,
            5,
            false
        ),
        2
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "ATG".as_bytes(),
            0,
            3,
            "ATGATGATGATG".as_bytes(),
            0,
            12,
            false
        ),
        4
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "ATG".as_bytes(),
            0,
            3,
            "ATGATGATGATGATG".as_bytes(),
            3,
            12,
            false
        ),
        4
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "G".as_bytes(),
            0,
            1,
            "ATGATGATGATG".as_bytes(),
            0,
            12,
            false
        ),
        1
    );
    assert_eq!(
        VariantContextUtils::find_number_of_repetitions_main(
            "G".as_bytes(),
            0,
            1,
            "ATGATGATGATGATG".as_bytes(),
            0,
            12,
            false
        ),
        1
    );
}

struct MergeAllelesTest {
    inputs: Vec<Vec<ByteArrayAllele>>,
    expected: Vec<ByteArrayAllele>,
}

impl MergeAllelesTest {
    fn new(mut arg: Vec<Vec<ByteArrayAllele>>) -> Self {
        let last = arg.remove(arg.len() - 1);
        Self {
            expected: last,
            inputs: arg
        }
    }
}

#[test]
fn merge_alleles_data() {
    let Aref = ByteArrayAllele::new(b"A", true);
    let Cref = ByteArrayAllele::new(b"C", true);
    let Gref = ByteArrayAllele::new(b"G", true);
    let Tref = ByteArrayAllele::new(b"T", true);
    let A = ByteArrayAllele::new(b"A", false);
    let T = ByteArrayAllele::new(b"T", false);
    let C = ByteArrayAllele::new(b"C", false);
    let G = ByteArrayAllele::new(b"G", false);
    let ATC = ByteArrayAllele::new(b"ATC", false);
    let ATCref = ByteArrayAllele::new(b"ATC", true);
    let ATCATC = ByteArrayAllele::new(b"ATCATC", false);
    let ATCATCT = ByteArrayAllele::new(b"ATCATCT", false);
    let ATref = ByteArrayAllele::new(b"AT", true);
    let Anoref = ByteArrayAllele::new(b"A", false);
    let GT = ByteArrayAllele::new(b"GT", false);

    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone()], vec![Aref.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone()], vec![Aref.clone()], vec![Aref.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone()], vec![Aref.clone(), T.clone()], vec![Aref.clone(), T.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone(), C.clone()], vec![Aref.clone(), T.clone()], vec![Aref.clone(), C.clone(), T.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone(), T.clone()], vec![Aref.clone(), C.clone()], vec![Aref.clone(), T.clone(), C.clone()]]));

    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone(), C.clone(), T.clone()], vec![Aref.clone(), C.clone(), T.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone(), T.clone(), C.clone()], vec![Aref.clone(), T.clone(), C.clone()]]));

    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone()], vec![Aref.clone(), ATC.clone()], vec![Aref.clone(), ATC.clone()]]));
    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone()], vec![Aref.clone(), ATC.clone(), ATCATC.clone()], vec![Aref.clone(), ATC.clone(), ATCATC.clone()]]));

    test_merge_alleles(MergeAllelesTest::new(vec![vec![Aref.clone(), ATCATC.clone()], vec![Aref.clone(), ATC.clone(), ATCATC.clone()], vec![Aref.clone(), ATCATC.clone(), ATC.clone()]]));
    test_merge_alleles(
        MergeAllelesTest::new(
            vec![
                vec![ATref.clone(), ATC.clone(), Anoref.clone(), G.clone()],
                vec![Aref.clone(), ATCATC.clone(), G.clone()],
                vec![ATref, ATC.clone(), Anoref.clone(), G.clone(), ATCATCT.clone(), GT.clone()]
            ]
        )
    );
}

fn test_merge_alleles(cfg: MergeAllelesTest) {
    let mut inputs = Vec::new();

    for (i, alleles) in cfg.inputs.iter().enumerate() {
        let name = format!("vcf{}", i);
        inputs.push(makeVC(name, alleles.clone(), None, None));
    }

    let priorities = vcs2priority(&inputs);
    let original_size = priorities.len();
    let mut merged = VariantContextUtils::simple_merge(
        inputs, Some(priorities),
        original_size,
        FilteredRecordMergeType::KeepIfAnyUnfiltered,
        GenotypeMergeType::Prioritize,
        false,
    ).unwrap();

    assert_eq!(merged.get_alleles().len(), cfg.expected.len());
    assert_eq!(merged.alleles.iter().collect::<LinkedHashSet<&ByteArrayAllele>>(), cfg.expected.iter().collect::<LinkedHashSet<&ByteArrayAllele>>());
}

fn makeVC<S: Into<String>>(source: S, alleles: Vec<ByteArrayAllele>, genotypes: Option<Vec<Genotype>>, filters: Option<HashSet<String>>) -> VariantContext {
    let start = 10;
    let stop = start + alleles[0].len() - 1;
    let mut vc =  VariantContext::build(0, start, stop, alleles);
    vc.source = source.into();

    match genotypes {
        Some(genotypes) => {
            vc.genotypes = GenotypesContext::new(genotypes)
        },
        None => {
            // pass
        }
    }

    // match filters {
    //     Some(filters) => {
    //         vc.filters = filters
    //     },
    //     None => {
    //         // pass
    //     }
    // }

    return vc
}

fn makeG<S: Into<String>>(sample: S, a1: ByteArrayAllele, a2: ByteArrayAllele, log10_p_error: f64, pls: Option<Vec<i64>>) -> Genotype {
    let mut g = Genotype::build_from_alleles(vec![a1, a2], sample.into());
    g.log10_p_error(log10_p_error);
    match pls {
        Some(pls) => {
            g.pl = pls;
        },
        None => {
            // pass
        }
    }
    return g
}

fn vcs2priority(vcs: &Vec<VariantContext>) -> Vec<String> {
    vcs.into_iter().map(|vc| vc.source.clone()).collect()
}

struct MergeGenotypesTest {
    inputs: Vec<VariantContext>,
    expected: VariantContext,
    priority: Vec<String>
}

impl MergeGenotypesTest {
    fn new(name: &str, priority: &str, mut arg: Vec<VariantContext>) -> Self {
        let last = arg.remove(arg.len() - 1);
        let priorities = priority.split(",").map(|s| s.to_string()).collect::<Vec<String>>();
        Self {
            expected: last,
            inputs: arg,
            priority: priorities
        }
    }
}

#[test]
fn merge_genotypes_data() {
    let Aref = ByteArrayAllele::new(b"A", true);
    let Cref = ByteArrayAllele::new(b"C", true);
    let Gref = ByteArrayAllele::new(b"G", true);
    let Tref = ByteArrayAllele::new(b"T", true);
    let A = ByteArrayAllele::new(b"A", false);
    let T = ByteArrayAllele::new(b"T", false);
    let C = ByteArrayAllele::new(b"C", false);
    let G = ByteArrayAllele::new(b"G", false);
    let ATC = ByteArrayAllele::new(b"ATC", false);
    let ATCref = ByteArrayAllele::new(b"ATC", true);
    let ATCATC = ByteArrayAllele::new(b"ATCATC", false);
    let ATCATCT = ByteArrayAllele::new(b"ATCATCT", false);
    let ATref = ByteArrayAllele::new(b"AT", true);
    let Anoref = ByteArrayAllele::new(b"A", false);
    let GT = ByteArrayAllele::new(b"GT", false);


    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypeByPriority-1,2", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypeByPriority-2,1", "2,1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "NonOverlappingGenotypes", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s2", Aref.clone(), T.clone(), -2.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None), makeG("s2", Aref.clone(), T.clone(), -2.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "PreserveAlleles", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), C.clone()], Some(vec![makeG("s2", Aref.clone(), C.clone(), -2.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None), makeG("s2", Aref.clone(), C.clone(), -2.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypePartialOverlap-1,2", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None), makeG("s3", Aref.clone(), T.clone(), -3.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None), makeG("s3", Aref.clone(), T.clone(), -3.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypePartialOverlap-2,1", "2,1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, None)]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None), makeG("s3", Aref.clone(), T.clone(), -3.0, None)]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, None), makeG("s3", Aref.clone(), T.clone(), -3.0, None)]), None)
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "OrderedPLs", "1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3]))]), None),
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3]))]), None),
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "OrderedPLs-3Alleles", "1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
                makeVC("1", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "OrderedPLs-3Alleles-2", "1",
            vec![
                makeVC("1", vec![Aref.clone(), C.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
                makeVC("1", vec![Aref.clone(), C.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "OrderedPLs-3Alleles-3", "1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
                makeVC("1", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s2", Aref.clone(), C.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
                makeVC("1", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6])), makeG("s2", Aref.clone(), C.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
            ]
        )
    );

    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypePartialOverlapWithPLs-2,1", "2,1",
            vec![
                makeVC("1", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -1.0, Some(vec![5, 0, 3]))]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, Some(vec![4, 0, 2])), makeG("s3", Aref.clone(), T.clone(), -3.0, Some(vec![3, 0, 2]))]), None),
                makeVC("3", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, Some(vec![4, 0, 2])), makeG("s3", Aref.clone(), T.clone(), -3.0, Some(vec![3, 0, 2]))]), None),
            ]
        )
    );

    println!("PartialOverlap");
    test_merge_genotypes(
        MergeGenotypesTest::new(
            "TakeGenotypePartialOverlapWithPLs-1,2", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), ATC.clone()], Some(vec![makeG("s1", Aref.clone(), ATC.clone(), -1.0, Some(vec![5, 0, 3]))]), None),
                makeVC("2", vec![Aref.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), T.clone(), -2.0, Some(vec![4, 0, 2])), makeG("s3", Aref.clone(), T.clone(), -3.0, Some(vec![3, 0, 2]))]), None),
                makeVC("3", vec![Aref.clone(), ATC.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), ATC.clone(), -1.0, None), makeG("s3", Aref.clone(), T.clone(), -3.0, None)]), None),
            ]
        )
    );

    println!("MultipleSamples");
    test_merge_genotypes(
        MergeGenotypesTest::new(
            "MultipleSamplePLsDifferentOrder", "1,2",
            vec![
                makeVC("1", vec![Aref.clone(), C.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), C.clone(), -1.0, Some(vec![1, 2, 3, 4, 5, 6]))]), None),
                makeVC("2", vec![Aref.clone(), T.clone(), C.clone()], Some(vec![makeG("s2", Aref.clone(), T.clone(), -2.0, Some(vec![6, 5, 4, 3, 2, 1]))]), None),
                makeVC("3", vec![Aref.clone(), C.clone(), T.clone()], Some(vec![makeG("s1", Aref.clone(), C.clone(), -1.0, None), makeG("s2", Aref.clone(), T.clone(), -2.0, None)]), None),
            ]
        )
    );
    
}

fn test_merge_genotypes(cfg: MergeGenotypesTest) {
    let merged = VariantContextUtils::simple_merge(
        cfg.inputs.clone(), Some(cfg.priority.clone()),
        cfg.priority.len(),
        FilteredRecordMergeType::KeepIfAnyUnfiltered,
        GenotypeMergeType::Prioritize,
        false,
    ).unwrap();

    assert_eq!(merged.get_alleles(), cfg.expected.get_alleles());
    assert_genotypes_are_mostly_equal(merged.get_genotypes(), cfg.expected.get_genotypes(), &cfg.priority);

}

fn assert_genotypes_are_mostly_equal(actual: &GenotypesContext, expected: &GenotypesContext, priorities: &Vec<String>) {
    if actual == expected {
        assert!(true, "genotypes are equal")
    }

    assert_eq!(actual.size(), expected.size(), "Expected and actual sizes are different: {:?} -> {:?}", expected, actual);

    for value in actual.genotypes() {
        let index = expected.genotypes().into_iter().position(|v| v.sample_name == value.sample_name).unwrap();
        let expected_value = expected.get(index);
        assert_eq!(value.gq, expected_value.gq);
        println!("value has likelihoods: {}", value.has_likelihoods());
        assert_eq!(value.has_likelihoods(), expected_value.has_likelihoods());
        if value.has_likelihoods() {
            assert_eq!(value.get_likelihoods().as_pls(), expected_value.get_likelihoods().as_pls());
        };
    };
}