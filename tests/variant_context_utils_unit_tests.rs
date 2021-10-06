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

use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context;
use lorikeet_genome::model::variant_context::{VariantContext, VariantType};
use lorikeet_genome::model::variant_context_utils;
use lorikeet_genome::model::variant_context_utils::VariantContextUtils;
use lorikeet_genome::utils::simple_interval::Locatable;

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
