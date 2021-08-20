extern crate lorikeet_genome;
#[macro_use]
extern crate lazy_static;

use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context;
use lorikeet_genome::model::variant_context::{VariantContext, VariantType};
use lorikeet_genome::model::variant_context_utils;

static snp_loc: &str = "chr1";
const snp_loc_start: usize = 10;
const snp_loc_stop: usize = 10;

static del_loc: &str = "chr1";
const del_loc_start: usize = 20;
const del_loc_stop: usize = 22;

static ins_loc: &str = "chr1";
const ins_loc_start: usize = 20;
const ins_loc_stop: usize = 20;

lazy_static! {
    // static ref del: ByteArrayAllele = ByteArrayAllele::new("A".as_bytes(), false);
    // static ref del_ref: ByteArrayAllele = ByteArrayAllele::new("A".as_bytes(), true);
    //
    // static ref A: ByteArrayAllele = ByteArrayAllele::new("A".as_bytes(), false);
    // static ref C: ByteArrayAllele = ByteArrayAllele::new("C".as_bytes(), false);
    // static ref A_ref: ByteArrayAllele = ByteArrayAllele::new("A".as_bytes(), true);
    // static ref T: ByteArrayAllele = ByteArrayAllele::new("T".as_bytes(), false);
    // static ref T_ref: ByteArrayAllele = ByteArrayAllele::new("T".as_bytes(), true);
    //
    // static ref ATC: ByteArrayAllele = ByteArrayAllele::new("ATC".as_bytes(), false);
    // static ref ATC_ref: ByteArrayAllele = ByteArrayAllele::new("ATC".as_bytes(), true);
}

struct VariantContextUnitTest<'a> {
    del: ByteArrayAllele,
    del_ref: ByteArrayAllele,
    A: ByteArrayAllele,
    C: ByteArrayAllele,
    A_ref: ByteArrayAllele,
    T: ByteArrayAllele,
    T_ref: ByteArrayAllele,
    ATC: ByteArrayAllele,
    ATC_ref: ByteArrayAllele,
    basic_builder: VariantContext<'a>,
    snp_builder: VariantContext<'a>,
    ins_builder: VariantContext<'a>,
}

impl<'a> VariantContextUnitTest<'a> {
    pub fn new() -> Self {
        let basic_builder = VariantContext::build(
            0,
            10,
            10,
            vec![
                ByteArrayAllele::new("A".as_bytes(), true),
                ByteArrayAllele::new("T".as_bytes(), false),
            ],
        );
        let snp_builder = VariantContext::build(
            0,
            10,
            10,
            vec![
                ByteArrayAllele::new("A".as_bytes(), true),
                ByteArrayAllele::new("T".as_bytes(), false),
            ],
        );
        let ins_builder = VariantContext::build(
            0,
            20,
            22,
            vec![
                ByteArrayAllele::new("A".as_bytes(), true),
                ByteArrayAllele::new("ATC".as_bytes(), false),
            ],
        );

        Self {
            del: ByteArrayAllele::new("A".as_bytes(), false),
            del_ref: ByteArrayAllele::new("A".as_bytes(), true),
            A: ByteArrayAllele::new("A".as_bytes(), false),
            C: ByteArrayAllele::new("C".as_bytes(), false),
            A_ref: ByteArrayAllele::new("A".as_bytes(), true),
            T: ByteArrayAllele::new("T".as_bytes(), false),
            T_ref: ByteArrayAllele::new("T".as_bytes(), true),
            ATC: ByteArrayAllele::new("ATC".as_bytes(), false),
            ATC_ref: ByteArrayAllele::new("ATC".as_bytes(), true),
            basic_builder,
            snp_builder,
            ins_builder,
        }
    }
}

#[test]
fn test_determine_types() {
    let AC_ref: ByteArrayAllele = ByteArrayAllele::new("AC".as_bytes(), true);
    let AC: ByteArrayAllele = ByteArrayAllele::new("AC".as_bytes(), false);
    let AT: ByteArrayAllele = ByteArrayAllele::new("AC".as_bytes(), false);

    let C: ByteArrayAllele = ByteArrayAllele::new("C".as_bytes(), false);
    let CAT: ByteArrayAllele = ByteArrayAllele::new("CAT".as_bytes(), false);

    let TA: ByteArrayAllele = ByteArrayAllele::new("TA".as_bytes(), false);
    let TA_ref: ByteArrayAllele = ByteArrayAllele::new("TA".as_bytes(), true);
    let TC: ByteArrayAllele = ByteArrayAllele::new("TC".as_bytes(), false);

    let symbolic: ByteArrayAllele = ByteArrayAllele::new("<FOO>".as_bytes(), false);

    let mut vc_unit_test = VariantContextUnitTest::new();

    //test REF
    let alleles = vec![vc_unit_test.T_ref.clone()];
    let mut vc = VariantContext::build(0, snp_loc_stop, snp_loc_stop, alleles);
    assert_eq!(vc.get_type(), &VariantType::NoVariation);

    //test snp
    let alleles = vec![vc_unit_test.T_ref.clone(), vc_unit_test.A.clone()];
    let mut vc = VariantContext::build(0, snp_loc_stop, snp_loc_stop, alleles);
    assert_eq!(vc.get_type(), &VariantType::Snp);

    let alleles = vec![
        vc_unit_test.T_ref.clone(),
        vc_unit_test.A.clone(),
        vc_unit_test.C.clone(),
    ];
    let mut vc = VariantContext::build(0, snp_loc_stop, snp_loc_stop, alleles);
    assert_eq!(vc.get_type(), &VariantType::Snp);

    //test mnp
    let alleles = vec![AC_ref.clone(), TA.clone()];
    let mut vc = VariantContext::build(0, snp_loc_stop, snp_loc_stop + 1, alleles);
    assert_eq!(vc.get_type(), &VariantType::Mnp);

    let alleles = vec![
        vc_unit_test.ATC_ref,
        CAT.clone(),
        ByteArrayAllele::new("GGG".as_bytes(), false),
    ];
    let mut vc = VariantContext::build(0, snp_loc_stop, snp_loc_stop + 2, alleles);
    assert_eq!(vc.get_type(), &VariantType::Mnp);

    // test indels
}
