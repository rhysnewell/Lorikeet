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

use lorikeet_genome::assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::estimation::lorikeet_engine::ReadType;
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::smith_waterman::bindings::SWParameters;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::utils::simple_interval::SimpleInterval;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::{bam, bam::Read};
use std::convert::TryFrom;

lazy_static! {
    static ref HAPLOTYPE_TO_REFERENCE_SW_PARAMETERS: SWParameters = *NEW_SW_PARAMETERS;
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
    let mut bam =
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

    let mut reads = vec![
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
        // Vari
    }
}
