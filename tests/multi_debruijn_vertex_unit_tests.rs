#![allow(
    non_upper_case_globals,
    non_snake_case
)]

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate approx;

use lorikeet_genome::graphs::base_vertex::BaseVertex;
use lorikeet_genome::read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;

#[test]
fn merging_identicals() {
    let v1 = MultiDeBruijnVertex::new(b"fred".to_vec(), true);
    let v2 = MultiDeBruijnVertex::new(b"fred".to_vec(), true);
    assert_eq!(v1, v2);
}

#[test]
fn not_identical() {
    let mut v1 = MultiDeBruijnVertex::new(b"fred".to_vec(), false);
    let v2 = MultiDeBruijnVertex::new(b"fred".to_vec(), false);
    assert_ne!(&v1, &v2);

    assert_eq!(v1.get_kmer_size(), 4);
    assert_eq!(v1.get_suffix(), b'd');
    assert_eq!(&v1.sequence, b"fred");
    assert_eq!(v1.get_suffix_string(), "d".to_string());
    assert_eq!(v1.get_additional_sequence(true), b"fred");
    assert_eq!(v1.get_additional_sequence(false), b"d");

    assert_eq!(v1.get_kmer_size(), v2.get_kmer_size());
    assert_eq!(v1.get_additional_info(), v2.get_additional_info());
    assert_eq!(v1.has_ambiguous_sequence(), v2.has_ambiguous_sequence());
    v1.add_read("fred".to_string());
    assert_eq!(v1.get_additional_info(), "".to_string());
    v1.set_additional_info("some info".to_string());
    assert_eq!(v1.get_additional_info(), "some info".to_string());
}
