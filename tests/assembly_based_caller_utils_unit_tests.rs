#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;

use gkl::smithwaterman::Parameters;
use hashlink::LinkedHashMap;
use lorikeet_genome::assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::genotype::genotype_builder::{Genotype, GenotypesContext};
use lorikeet_genome::haplotype::event_map::EventMap;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::model::variant_context::{VariantContext, VariantType};
use lorikeet_genome::processing::lorikeet_engine::ReadType;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::NEW_SW_PARAMETERS;
use lorikeet_genome::test_utils::variant_context_test_utils::VariantContextTestUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rust_htslib::{bam, bam::Read};
lazy_static! {
    static ref HAPLOTYPE_TO_REFERENCE_SW_PARAMETERS: Parameters = *NEW_SW_PARAMETERS;
}

// In finalizeRegion(), the base qualities of overlapped read clips pairs are adjusted.
// Most of read clips are clipped/copy from original reads, and the base qualities of original reads are not affected.
// However, GATK 4.0.0.x has a bug : in some cases, the clipping procedure directly returns the original reads.
// So the base qualities of original reads are changed.
// This test is added to make sure the unexpected behavior is fixed.
#[test]
fn test_finalize_region() {
    let header = SimpleInterval::new(0, 1, 100000000);
    // NOTE: These reads are mates that overlap one-another without agreement which means they should have modified base qualities after calling finalize()
    //       Read2 has a clean cigar, and thus will not be copied by the clipping code before being fed to the overlapping bases code. This test asserts that its still copied.
    let bam =
        bam::Reader::from_path("tests/resources/small/test_finalize_region_reads.sam").unwrap();
    let bam_header = bam::Header::from_template(bam.header());
    // let bam_header = bam::Header::new();

    let mut active_region = AssemblyRegion::new(
        SimpleInterval::new(0, 42596728, 42598843),
        true,
        100,
        header.size(),
        0,
        0,
    );
    assert!(active_region.is_active());
    assert!(!active_region.is_finalized());

    let org_read1 = bam::Record::from_sam(&mut bam::HeaderView::from_header(&bam_header), "HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t83\t1\t42596803\t39\t1S95M5S\t=\t42596891\t-7\tGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATAGAAT\tDAAAEDCFCCGEEDDBEDDDGCCDEDECDDFDCEECCFEECDCEDBCDBDBCC>DCECC>DBCDDBCBDDBCDDEBCCECC>DBCDBDBGC?FCCBDB>>?".as_bytes()).unwrap();
    let org_read2 = bam::Record::from_sam(&mut bam::HeaderView::from_header(&bam_header), "HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t163\t1\t42596891\t39\t101M\t=\t42596803\t7\tCTCGAATGGAATCATTTTCTACTGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTAACATGGAATGGAATCGAATGTAATCATCATCAAATGGA\t>@>:ABCDECCCEDCBBBDDBDDEBCCBEBBCBEBCBCDDCD>DECBGCDCF>CCCFCDDCBABDEDFCDCDFFDDDG?DDEGDDFDHFEGDDGECB@BAA".as_bytes()).unwrap();

    let reads = vec![
        BirdToolRead::new(org_read1.clone(), 0, ReadType::Short),
        BirdToolRead::new(org_read2.clone(), 0, ReadType::Short),
    ];

    active_region.add_all(reads.clone());
    let min_bq = 9;

    // NOTE: this test MUST be run with correctOverlappingBaseQualities enabled otherwise this test can succeed even with unsafe code
    AssemblyBasedCallerUtils::finalize_regions(
        &mut active_region,
        false,
        false,
        min_bq,
        true,
        false,
    );
}

fn test_get_variant_contexts_from_given_alleles(
    loc: usize,
    active_alleles_to_genotype: Vec<VariantContext>,
    expected_vcs_at_this_location: Vec<VariantContext>,
) {
    let vcs_at_this_position = AssemblyBasedCallerUtils::get_variant_contexts_from_given_alleles(
        loc,
        &active_alleles_to_genotype,
        true,
    );
    assert_eq!(
        vcs_at_this_position.len(),
        expected_vcs_at_this_location.len()
    );
    for i in 0..expected_vcs_at_this_location.len() {
        VariantContextTestUtils::assert_variant_contexts_are_equal(
            &vcs_at_this_position[i],
            &expected_vcs_at_this_location[i],
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(
            &vcs_at_this_position[0].source,
            &expected_vcs_at_this_location[0].source
        );
    }
}

#[test]
fn get_vcs_at_this_location_from_given_alleles_data() {
    test_get_variant_contexts_from_given_alleles(1000, Vec::new(), Vec::new());

    let mut snp_haplotype: Haplotype<SimpleInterval> =
        Haplotype::new(b"ACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let snp_alleles = vec![
        ByteArrayAllele::new(b"A", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let snp_vc_builder = VariantContext::build(20, 1000, 1000, snp_alleles);
    snp_haplotype.set_event_map(EventMap::state_for_testing(vec![snp_vc_builder.clone()]));

    // this one matches the snp haplotype above (to test duplicate removal)
    let mut snp_haplotype_duplicate: Haplotype<SimpleInterval> =
        Haplotype::new(b"ACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let snp_alleles2 = vec![
        ByteArrayAllele::new(b"A", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let snp_vc_builder2 = VariantContext::build(20, 1000, 1000, snp_alleles2);
    let snp_alleles3 = vec![
        ByteArrayAllele::new(b"T", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let snp_vc_builder3 = VariantContext::build(20, 1020, 1020, snp_alleles3);
    snp_haplotype_duplicate.set_event_map(EventMap::state_for_testing(vec![
        snp_vc_builder2,
        snp_vc_builder3,
    ]));

    let mut deletion_haplotype: Haplotype<SimpleInterval> =
        Haplotype::new(b"ACTGGTCAGGTCAACTGGTCA", false);
    let deletion_alleles = vec![
        ByteArrayAllele::new(b"ACTGGTCAACT", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let deletion_builder = VariantContext::build(20, 995, 1005, deletion_alleles.clone());
    deletion_haplotype.set_event_map(EventMap::state_for_testing(vec![deletion_builder.clone()]));

    // matches the deletion alleles above but at a different position (to catch an edge case in duplicate removal)
    let mut deletion_haplotype_false_duplicate: Haplotype<SimpleInterval> =
        Haplotype::new(b"ACTGGTCAGGTCAACTGGTCA", false);
    let deletion_alleles_false_duplicate = vec![
        ByteArrayAllele::new(b"ACTGGTCAACT", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let deletion_builder_false_duplicate =
        VariantContext::build(20, 998, 1008, deletion_alleles_false_duplicate);
    deletion_haplotype_false_duplicate.set_event_map(EventMap::state_for_testing(vec![
        deletion_builder_false_duplicate.clone(),
    ]));

    // doesn't overlap 1000
    let mut deletion_haplotype_no_span: Haplotype<SimpleInterval> =
        Haplotype::new(b"CAACTGGTCAACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let deletion_alleles_no_span = vec![
        ByteArrayAllele::new(b"GTCAA", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let deletion_vc_no_span = VariantContext::build(20, 990, 994, deletion_alleles_no_span);
    deletion_haplotype_no_span.set_event_map(EventMap::state_for_testing(vec![
        deletion_vc_no_span.clone(),
    ]));

    let mut same_loc_del_hap1: Haplotype<SimpleInterval> = Haplotype::new(b"AAAAAAAGAAA", false);
    let same_loc_del_alleles1 = vec![
        ByteArrayAllele::new(b"GTT", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let same_loc_del_vc1 = VariantContext::build(20, 10093568, 10093570, same_loc_del_alleles1);
    same_loc_del_hap1.set_event_map(EventMap::state_for_testing(vec![same_loc_del_vc1]));

    let mut same_loc_del_hap2: Haplotype<SimpleInterval> = Haplotype::new(b"AAAAAAAGTAAA", false);
    let same_loc_del_alleles2 = vec![
        ByteArrayAllele::new(b"GT", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let same_loc_del_vc2 = VariantContext::build(20, 10093568, 10093569, same_loc_del_alleles2);
    same_loc_del_hap2.set_event_map(EventMap::state_for_testing(vec![same_loc_del_vc2]));

    let mut same_loc_ins_hap1: Haplotype<SimpleInterval> = Haplotype::new(b"AAAAAAAGTTTAAA", false);
    let same_loc_ins_alleles1 = vec![
        ByteArrayAllele::new(b"G", true),
        ByteArrayAllele::new(b"GT", false),
    ];
    let same_loc_ins_vc1 = VariantContext::build(20, 10093568, 10093568, same_loc_ins_alleles1);
    same_loc_ins_hap1.set_event_map(EventMap::state_for_testing(vec![same_loc_ins_vc1]));

    let mut deletion_vc_builder_with_gts =
        VariantContext::build(20, 995, 1005, deletion_alleles.clone());
    deletion_vc_builder_with_gts.genotypes =
        GenotypesContext::new(vec![Genotype::build_from_alleles(
            deletion_alleles,
            "TEST".to_string(),
        )]);

    let mut snp_vc_expected = snp_vc_builder.clone();
    snp_vc_expected.source = "Comp0Allele0".to_string();

    test_get_variant_contexts_from_given_alleles(
        1000,
        vec![snp_vc_builder.clone()],
        vec![snp_vc_expected.clone()],
    );

    let mut del_vc_expected = deletion_builder.clone();
    del_vc_expected.source = "Comp0Allele0".to_string();
    test_get_variant_contexts_from_given_alleles(
        995,
        vec![deletion_builder.clone()],
        vec![del_vc_expected.clone()],
    );
    test_get_variant_contexts_from_given_alleles(
        1000,
        vec![deletion_builder.clone()],
        vec![del_vc_expected.clone()],
    );

    snp_vc_expected.source = "Comp1Allele0".to_string();
    test_get_variant_contexts_from_given_alleles(
        1000,
        vec![deletion_builder.clone(), snp_vc_builder.clone()],
        vec![del_vc_expected.clone(), snp_vc_expected.clone()],
    );

    // let mut del_vc_no_span_expected = deletion_vc_no_span.clone();
    // del_vc_no_span_expected.source =
    test_get_variant_contexts_from_given_alleles(
        1000,
        vec![deletion_builder.clone(), deletion_vc_no_span.clone()],
        vec![del_vc_expected.clone()],
    );
    let mut del_vc_false_duplicate_expected = deletion_builder_false_duplicate.clone();
    del_vc_false_duplicate_expected.source = "Comp1Allele0".to_string();
    // del_vc_no_span_expected.source =
    test_get_variant_contexts_from_given_alleles(
        1000,
        vec![
            deletion_builder.clone(),
            deletion_builder_false_duplicate.clone(),
            deletion_vc_no_span.clone(),
        ],
        vec![
            del_vc_expected.clone(),
            del_vc_false_duplicate_expected.clone(),
        ],
    );
}

fn test_get_event_mapper(
    merged_vc: VariantContext,
    loc: usize,
    haplotypes: Vec<Haplotype<SimpleInterval>>,
    expected_event_map: LinkedHashMap<usize, Vec<&Haplotype<SimpleInterval>>>,
) {
    let actual_event_map =
        AssemblyBasedCallerUtils::create_allele_mapper(&merged_vc, loc, &haplotypes, true);
    assert_eq!(actual_event_map.len(), expected_event_map.len());
    for key in actual_event_map.keys() {
        assert!(expected_event_map.contains_key(key));
        assert_eq!(actual_event_map.get(key), expected_event_map.get(key));
    }

    for key in expected_event_map.keys() {
        assert!(actual_event_map.contains_key(key))
    }
}

#[test]
fn get_event_mapper_data() {
    let mut ref_haplotype = Haplotype::new(b"ACTGGTCAACTAGTCAACTGGTCAACTGGTCA", true);
    ref_haplotype.set_event_map(EventMap::state_for_testing(Vec::new()));

    let mut snp_haplotype = Haplotype::new(b"ACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let ref_allele = ByteArrayAllele::new(b"A", true);
    let snp_allele = ByteArrayAllele::new(b"G", false);
    let snp_alleles = vec![ref_allele.clone(), snp_allele.clone()];
    let mut snp_vc = VariantContext::build(20, 1000, 1000, snp_alleles.clone());
    snp_vc.set_type(VariantType::Snp);
    snp_haplotype.set_event_map(EventMap::state_for_testing(vec![snp_vc.clone()]));

    let mut test1_expected_map = LinkedHashMap::new();
    for (i, a) in snp_alleles.iter().enumerate() {
        if &snp_alleles[1] == a {
            test1_expected_map.insert(i, vec![&snp_haplotype]);
        } else {
            test1_expected_map.insert(i, vec![&ref_haplotype]);
        }
    }
    test_get_event_mapper(
        snp_vc.clone(),
        snp_vc.loc.get_start(),
        vec![snp_haplotype.clone(), ref_haplotype.clone()],
        test1_expected_map,
    );
}

fn test_get_variants_contexts_from_active_haplotypes(
    haplotypes: Vec<Haplotype<SimpleInterval>>,
    loc: usize,
    expected_vcs_at_this_loc: Vec<&VariantContext>,
) {
    let vcs_at_this_position =
        AssemblyBasedCallerUtils::get_variant_contexts_from_active_haplotypes(
            loc,
            &haplotypes,
            true,
        );

    assert_eq!(
        vcs_at_this_position.len(),
        expected_vcs_at_this_loc.len(),
        "Incorrect number of vcs returned"
    );
    for i in 0..expected_vcs_at_this_loc.len() {
        VariantContextTestUtils::assert_variant_contexts_are_equal(
            vcs_at_this_position[i],
            expected_vcs_at_this_loc[i],
            Vec::new(),
            Vec::new(),
        );

        assert_eq!(
            vcs_at_this_position[i].source, expected_vcs_at_this_loc[i].source,
            "Sources differ {i}"
        );
    }
}

#[test]
fn get_variant_contexts_from_active_haplotypes_data() {
    let dummy = Vec::new();
    let dummy2 = Vec::new();
    test_get_variants_contexts_from_active_haplotypes(dummy, 1000, dummy2);

    let mut snp_haplotype = Haplotype::new(b"ACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let snp_alleles = vec![
        ByteArrayAllele::new(b"A", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let mut snp_vc = VariantContext::build(20, 1000, 1000, snp_alleles);
    snp_vc.source = "a".to_string();
    snp_vc.variant_type = Some(VariantType::Snp);
    snp_haplotype.set_event_map(EventMap::state_for_testing(vec![snp_vc.clone()]));

    // this one matches above to test dup removal
    let mut snp_haplotype_duplicate = Haplotype::new(b"ACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let snp_alleles = vec![
        ByteArrayAllele::new(b"A", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let mut snp_vc2 = VariantContext::build(20, 1000, 1000, snp_alleles);
    snp_vc2.source = "a".to_string();
    snp_vc2.variant_type = Some(VariantType::Snp);

    let snp_alleles = vec![
        ByteArrayAllele::new(b"T", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let mut snp_vc3 = VariantContext::build(20, 1020, 1020, snp_alleles);
    snp_vc3.source = "a".to_string();
    snp_vc3.variant_type = Some(VariantType::Snp);
    snp_haplotype_duplicate.set_event_map(EventMap::state_for_testing(vec![snp_vc2, snp_vc3]));

    let mut deletion_haplotype = Haplotype::new(b"ACTGGTCAGGTCAACTGGTCA", false);
    let deletion_alleles = vec![
        ByteArrayAllele::new(b"ACTGGTCAACT", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let mut deletion_vc = VariantContext::build(20, 995, 1005, deletion_alleles);
    deletion_vc.source = "a".to_string();
    deletion_vc.variant_type = Some(VariantType::Indel);
    deletion_haplotype.set_event_map(EventMap::state_for_testing(vec![deletion_vc.clone()]));

    let mut deletion_false_haplotype = Haplotype::new(b"ACTGGTCAGGTCAACTGGTCA", false);
    let deletion_false_alleles = vec![
        ByteArrayAllele::new(b"ACTGGTCAACT", true),
        ByteArrayAllele::new(b"A", false),
    ];
    let mut deletion_false_vc = VariantContext::build(20, 998, 1008, deletion_false_alleles);
    deletion_false_vc.source = "a".to_string();
    deletion_false_vc.variant_type = Some(VariantType::Indel);
    deletion_false_haplotype
        .set_event_map(EventMap::state_for_testing(vec![deletion_false_vc.clone()]));

    // doesn't overlap 1000
    let mut deletion_haplotype_no_span =
        Haplotype::new(b"CAACTGGTCAACTGGTCAACTGGTCAACTGGTCAACTGGTCA", false);
    let deletion_alleles_no_span = vec![
        ByteArrayAllele::new(b"GTCAA", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let mut deletion_vc_no_span = VariantContext::build(20, 990, 994, deletion_alleles_no_span);
    deletion_vc_no_span.source = "a".to_string();
    deletion_vc_no_span.variant_type = Some(VariantType::Indel);
    deletion_haplotype_no_span
        .set_event_map(EventMap::state_for_testing(vec![deletion_vc_no_span]));

    test_get_variants_contexts_from_active_haplotypes(
        vec![snp_haplotype.clone()],
        1000,
        vec![&snp_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![snp_haplotype.clone(), snp_haplotype_duplicate.clone()],
        1000,
        vec![&snp_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![deletion_haplotype.clone()],
        995,
        vec![&deletion_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![deletion_haplotype.clone()],
        1000,
        vec![&deletion_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![
            deletion_haplotype.clone(),
            deletion_haplotype_no_span.clone(),
        ],
        1000,
        vec![&deletion_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![
            deletion_haplotype.clone(),
            deletion_false_haplotype.clone(),
            deletion_haplotype_no_span.clone(),
        ],
        1000,
        vec![&deletion_vc, &deletion_false_vc],
    );
    test_get_variants_contexts_from_active_haplotypes(
        vec![deletion_haplotype.clone(), snp_haplotype.clone()],
        1000,
        vec![&deletion_vc, &snp_vc],
    );

    let mut same_loc_del1_haplotype = Haplotype::new(b"AAAAAAAGAAA", false);
    let same_loc_del1_alleles = vec![
        ByteArrayAllele::new(b"GTT", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let mut same_loc_del1_vc = VariantContext::build(20, 10093568, 10093570, same_loc_del1_alleles);
    same_loc_del1_vc.source = "a".to_string();
    same_loc_del1_vc.variant_type = Some(VariantType::Indel);
    same_loc_del1_haplotype
        .set_event_map(EventMap::state_for_testing(vec![same_loc_del1_vc.clone()]));

    let mut same_loc_del2_haplotype = Haplotype::new(b"AAAAAAAGTAAA", false);
    let same_loc_del2_alleles = vec![
        ByteArrayAllele::new(b"GT", true),
        ByteArrayAllele::new(b"G", false),
    ];
    let mut same_loc_del2_vc = VariantContext::build(20, 10093568, 10093569, same_loc_del2_alleles);
    same_loc_del2_vc.source = "a".to_string();
    same_loc_del2_vc.variant_type = Some(VariantType::Indel);
    same_loc_del2_haplotype
        .set_event_map(EventMap::state_for_testing(vec![same_loc_del2_vc.clone()]));

    let mut same_loc_ins1_haplotype = Haplotype::new(b"AAAAAAAGTTTAAA", false);
    let same_loc_ins1_alleles = vec![
        ByteArrayAllele::new(b"G", true),
        ByteArrayAllele::new(b"GT", false),
    ];
    let mut same_loc_ins1_vc = VariantContext::build(20, 10093568, 10093568, same_loc_ins1_alleles);
    same_loc_ins1_vc.source = "a".to_string();
    same_loc_ins1_vc.variant_type = Some(VariantType::Indel);
    same_loc_ins1_haplotype
        .set_event_map(EventMap::state_for_testing(vec![same_loc_ins1_vc.clone()]));

    test_get_variants_contexts_from_active_haplotypes(
        vec![
            same_loc_del1_haplotype.clone(),
            same_loc_del2_haplotype.clone(),
            same_loc_ins1_haplotype.clone(),
        ],
        10093568,
        vec![&same_loc_del1_vc, &same_loc_del2_vc, &same_loc_ins1_vc],
    )
}
