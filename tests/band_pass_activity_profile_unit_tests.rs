#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]

extern crate itertools;
extern crate lorikeet_genome;
extern crate rust_htslib;
#[macro_use]
extern crate approx;

use itertools::Itertools;
use lorikeet_genome::activity_profile::activity_profile::{ActivityProfile, Profile};
use lorikeet_genome::activity_profile::activity_profile_state::{ActivityProfileState, Type};
use lorikeet_genome::activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::haplotype::event_map::EventMap;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::byte_array_allele::ByteArrayAllele;
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::variant_context_utils::VariantContextUtils;
use lorikeet_genome::reads::cigar_utils::CigarUtils;
use lorikeet_genome::reference::reference_reader::ReferenceReader;
use lorikeet_genome::reference::reference_reader_utils::ReferenceReaderUtils;
use lorikeet_genome::utils::math_utils::MathUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use rust_htslib::bam::record::{Cigar, CigarString};
use std::cmp::{max, min};
use std::convert::TryFrom;

const MAX_PROB_PROPAGATION_DISTANCE: usize = 50;
const ACTIVE_PROB_THRESHOLD: f64 = 0.002;
static b37_reference_20_21: &str = "tests/resources/large/human_g1k_v37.20.21.fasta";

fn test_band_pass(
    start: usize,
    preceding_is_active: bool,
    n_preceding_sites: usize,
    band_pass_size: usize,
    sigma: f64,
    ref_idx: usize,
    tid: usize,
    contig_len: usize,
) {
    let mut profile = BandPassActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD,
        band_pass_size,
        sigma,
        false,
        ref_idx,
        tid,
        contig_len,
    );

    let expected_band_size = band_pass_size * 2 + 1;
    assert_eq!(
        profile.get_filtered_size(),
        band_pass_size,
        "Wrong filter size"
    );
    assert_eq!(profile.get_sigma(), sigma, "Wrong sigma");
    assert_eq!(
        profile.get_band_size(),
        expected_band_size,
        "Wrong band size"
    );

    let preceding_prob = if preceding_is_active { 1.0 } else { 0.0 };
    for i in 0..n_preceding_sites {
        let loc = SimpleInterval::new(0, i + start, i + start);
        let state = ActivityProfileState::new(loc, preceding_prob, Type::None);
        profile.add(state);
    }

    let next_loc = SimpleInterval::new(0, n_preceding_sites + start, n_preceding_sites + start);
    profile.add(ActivityProfileState::new(next_loc, 1.0, Type::None));

    if !preceding_is_active && n_preceding_sites >= band_pass_size && band_pass_size < start {
        // we have enough space that all probs fall on the genome
        let probs = profile.get_probabilities_as_array();
        assert!(
            relative_eq!(
                probs.iter().sum::<f64>(),
                1.0 * (n_preceding_sites as f64 * preceding_prob + 1.0),
                epsilon = 1e-3
            ),
            "Activity profile doesn't sum to number of non-zero prob states"
        );
    }
}

#[test]
fn make_band_pass_test() {
    let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);

    let contig_len = ref_reader.index.sequences()[0].len as usize;
    for start in vec![1, 10, 100, 1000] {
        for preceding_is_active in vec![true, false] {
            for preceding_sites in vec![0, 1, 10, 100] {
                for band_pass_size in vec![0, 1, 10, 100] {
                    for sigma in vec![1.0, 2.0, BandPassActivityProfile::DEFAULT_SIGMA] {
                        test_band_pass(
                            start,
                            preceding_is_active,
                            preceding_sites,
                            band_pass_size,
                            sigma,
                            0,
                            0,
                            contig_len,
                        );
                    }
                }
            }
        }
    }
}

fn band_pass_in_one_pass(
    profile: &BandPassActivityProfile,
    active_prob_array: &Vec<f64>,
) -> Vec<f64> {
    let mut band_pass_prob_array = vec![0.0; active_prob_array.len()];

    // apply the band pass filter for activeProbArray into filteredProbArray
    let mut gaussian_kernel = &profile.gaussian_kernel;
    for iii in 0..active_prob_array.len() {
        let kernel = &gaussian_kernel[max((profile.get_filtered_size() as i64 - iii as i64), 0)
            as usize
            ..min(
                gaussian_kernel.len(),
                profile.get_filtered_size() + active_prob_array.len() - iii,
            )];

        let active_prob_sub_array =
            &active_prob_array[max((iii as i64 - profile.get_filtered_size() as i64), 0) as usize
                ..min(
                    active_prob_array.len(),
                    iii + profile.get_filtered_size() + 1,
                )];

        band_pass_prob_array[iii] = MathUtils::dot_product(active_prob_sub_array, kernel);
    }

    return band_pass_prob_array;
}

fn test_band_pass_composition(
    band_pass_size: usize,
    integration_length: usize,
    ref_idx: usize,
    tid: usize,
    contig_len: usize,
) {
    let start = 1;
    let mut profile = BandPassActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD,
        band_pass_size,
        BandPassActivityProfile::DEFAULT_SIGMA,
        true,
        ref_idx,
        tid,
        contig_len,
    );

    let mut raw_active_probs = vec![0.0; integration_length + band_pass_size * 2];
    let raw_active_prob_len = raw_active_probs.len();

    // add a buffer so that we can get all of the band pass values
    let mut pos = start;
    let mut raw_prob_offset = 0;
    for i in 0..band_pass_size {
        let loc = SimpleInterval::new(tid, pos, pos);
        pos += 1;
        let state = ActivityProfileState::new(loc, 0.0, Type::None);
        profile.add(state);
        raw_active_probs[raw_prob_offset] = 0.0;
        raw_prob_offset += 1;
        raw_active_probs[raw_active_prob_len - raw_prob_offset] = 0.0;
    }

    for i in 0..integration_length {
        let next_loc = SimpleInterval::new(tid, pos, pos);
        pos += 1;
        profile.add(ActivityProfileState::new(next_loc, 1.0, Type::None));
        raw_active_probs[raw_prob_offset] = 1.0;
        raw_prob_offset += 1;

        for j in 0..profile.size() {
            assert!(
                profile.get_state_list()[j].is_active_prob() >= 0.0,
                "State probability < 0 at {}",
                j
            );
            assert!(
                profile.get_state_list()[j].is_active_prob() <= 1.0 + 1e-3,
                "State probability > 1 at {}",
                j
            );
        }
    }

    let expected_probs = band_pass_in_one_pass(&profile, &raw_active_probs);
    for j in 0..profile.size() {
        assert!(
            relative_eq!(
                profile.get_state_list()[j].is_active_prob(),
                expected_probs[j],
                epsilon = 1e-3
            ),
            "State probability not expected at {} got {} expected {} for profile {:?}",
            j,
            profile.get_state_list()[j].is_active_prob(),
            expected_probs[j],
            &profile
        )
    }
}

#[test]
fn make_band_pass_composition() {
    let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);

    let contig_len = ref_reader.index.sequences()[0].len as usize;
    for band_pass_size in vec![0, 1, 10, 100, BandPassActivityProfile::MAX_FILTER_SIZE] {
        for integration_length in vec![1, 10, 100, 1000] {
            test_band_pass_composition(band_pass_size, integration_length, 0, 0, contig_len);
        }
    }
}

fn test_kernel_creation(sigma: f64, max_size: usize, contig_len: usize, expected_kernel: Vec<f64>) {
    let profile = BandPassActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD,
        max_size,
        sigma,
        true,
        0,
        0,
        contig_len,
    );
    let kernel = profile.get_kernel();
    assert_eq!(kernel.len(), expected_kernel.len());
    for i in 0..kernel.len() {
        assert!(
            relative_eq!(kernel[i], expected_kernel[i], epsilon = 1e-3),
            "Kernels not equal at {}",
            i
        );
    }
}

#[test]
fn make_kernel_creation() {
    let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);

    let contig_len = ref_reader.index.sequences()[0].len as usize;
    test_kernel_creation(0.01, 1000, contig_len, vec![1.0]);
    test_kernel_creation(
        1.0,
        1000,
        contig_len,
        vec![
            0.0001338302,
            0.004431848,
            0.053990966,
            0.241970723,
            0.398942278,
            0.241970723,
            0.053990966,
            0.004431848,
            0.0001338302,
        ],
    );
    test_kernel_creation(1.0, 0, contig_len, vec![1.0]);
    test_kernel_creation(1.0, 1, contig_len, vec![0.2740686, 0.4518628, 0.2740686]);
    test_kernel_creation(
        1.0,
        2,
        contig_len,
        vec![0.05448868, 0.24420134, 0.40261995, 0.24420134, 0.05448868],
    );
    test_kernel_creation(
        5.0,
        1000,
        contig_len,
        vec![
            1.1788613551308e-05,
            2.67660451529771e-05,
            5.83893851582921e-05,
            0.000122380386022754,
            0.000246443833694604,
            0.000476817640292968,
            0.000886369682387602,
            0.00158309031659599,
            0.00271659384673712,
            0.00447890605896858,
            0.00709491856924629,
            0.0107981933026376,
            0.0157900316601788,
            0.0221841669358911,
            0.029945493127149,
            0.0388372109966426,
            0.0483941449038287,
            0.0579383105522965,
            0.0666449205783599,
            0.0736540280606647,
            0.0782085387950912,
            0.0797884560802865,
            0.0782085387950912,
            0.0736540280606647,
            0.0666449205783599,
            0.0579383105522965,
            0.0483941449038287,
            0.0388372109966426,
            0.029945493127149,
            0.0221841669358911,
            0.0157900316601788,
            0.0107981933026376,
            0.00709491856924629,
            0.00447890605896858,
            0.00271659384673712,
            0.00158309031659599,
            0.000886369682387602,
            0.000476817640292968,
            0.000246443833694604,
            0.000122380386022754,
            5.83893851582921e-05,
            2.67660451529771e-05,
            1.1788613551308e-05,
        ],
    );
    test_kernel_creation(
        17.0,
        1000,
        contig_len,
        vec![
            1.25162575710745e-05,
            1.57001772728555e-05,
            1.96260034693739e-05,
            2.44487374842009e-05,
            3.03513668801384e-05,
            3.75489089511911e-05,
            4.62928204154855e-05,
            5.68757597480354e-05,
            6.96366758708924e-05,
            8.49661819944029e-05,
            0.000103312156275406,
            0.000125185491708561,
            0.000151165896477646,
            0.000181907623161359,
            0.000218144981137171,
            0.000260697461819069,
            0.000310474281706066,
            0.000368478124457557,
            0.000435807841336874,
            0.00051365985048857,
            0.000603327960854364,
            0.000706201337376934,
            0.000823760321812988,
            0.000957569829285965,
            0.00110927005589186,
            0.00128056425833231,
            0.00147320340358764,
            0.00168896753568649,
            0.00192964376796036,
            0.00219700088266432,
            0.00249276060490197,
            0.00281856571330067,
            0.00317594525418154,
            0.00356627723683793,
            0.00399074930220799,
            0.00445031797242299,
            0.00494566720070898,
            0.00547716704583487,
            0.00604483338842317,
            0.00664828968356621,
            0.00728673180099395,
            0.00795889703644795,
            0.00866303838230695,
            0.00939690511889675,
            0.0101577307281371,
            0.010942229037054,
            0.0117465993701676,
            0.0125665413280325,
            0.0133972796167302,
            0.0142335991336574,
            0.0150698902735454,
            0.0159002041614507,
            0.0167183172536454,
            0.0175178044808441,
            0.0182921198494897,
            0.0190346831745763,
            0.0197389714002676,
            0.020398612780527,
            0.0210074820484496,
            0.0215597946062309,
            0.0220501977225941,
            0.022473856734247,
            0.0228265343139947,
            0.0231046609899767,
            0.0233053952756892,
            0.0234266719946158,
            0.0234672376502799,
            0.0234266719946158,
            0.0233053952756892,
            0.0231046609899767,
            0.0228265343139947,
            0.022473856734247,
            0.0220501977225941,
            0.0215597946062309,
            0.0210074820484496,
            0.020398612780527,
            0.0197389714002676,
            0.0190346831745763,
            0.0182921198494897,
            0.0175178044808441,
            0.0167183172536454,
            0.0159002041614507,
            0.0150698902735454,
            0.0142335991336574,
            0.0133972796167302,
            0.0125665413280325,
            0.0117465993701676,
            0.010942229037054,
            0.0101577307281371,
            0.00939690511889675,
            0.00866303838230695,
            0.00795889703644795,
            0.00728673180099395,
            0.00664828968356621,
            0.00604483338842317,
            0.00547716704583487,
            0.00494566720070898,
            0.00445031797242299,
            0.00399074930220799,
            0.00356627723683793,
            0.00317594525418154,
            0.00281856571330067,
            0.00249276060490197,
            0.00219700088266432,
            0.00192964376796036,
            0.00168896753568649,
            0.00147320340358764,
            0.00128056425833231,
            0.00110927005589186,
            0.000957569829285965,
            0.000823760321812988,
            0.000706201337376934,
            0.000603327960854364,
            0.00051365985048857,
            0.000435807841336874,
            0.000368478124457557,
            0.000310474281706066,
            0.000260697461819069,
            0.000218144981137171,
            0.000181907623161359,
            0.000151165896477646,
            0.000125185491708561,
            0.000103312156275406,
            8.49661819944029e-05,
            6.96366758708924e-05,
            5.68757597480354e-05,
            4.62928204154855e-05,
            3.75489089511911e-05,
            3.03513668801384e-05,
            2.44487374842009e-05,
            1.96260034693739e-05,
            1.57001772728555e-05,
            1.25162575710745e-05,
        ],
    );
}
