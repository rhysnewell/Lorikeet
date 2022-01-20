#![allow(
    non_upper_case_globals,
    unused_parens,
    unused_mut,
    unused_imports,
    non_snake_case
)]
extern crate lorikeet_genome;
extern crate rust_htslib;

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
const ACTIVE_PROB_THRESHOLD: f32 = 0.002;
static b37_reference_20_21: &str = "tests/resources/large/human_g1k_v37.20.21.fasta";

#[derive(Debug, Clone, Copy)]
enum ProfileType {
    BandPass,
    Basic,
}

#[derive(Debug)]
enum ActivityProfileType {
    BandPass(BandPassActivityProfile),
    Basic(ActivityProfile),
}

#[derive(Debug, Clone)]
struct BasicActivityProfileTestProvider {
    probs: Vec<f32>,
    expected_regions: Vec<AssemblyRegion>,
    extension: usize,
    region_start: SimpleInterval,
    profile_type: ProfileType,
    ref_idx: usize,
    contig_len: usize,
}

impl BasicActivityProfileTestProvider {
    pub fn new(
        profile_type: ProfileType,
        probs: Vec<f32>,
        start_active: bool,
        region_start: SimpleInterval,
        contig_len: usize,
        ref_idx: usize,
        starts_and_stops: Vec<usize>,
    ) -> BasicActivityProfileTestProvider {
        let mut result = Self {
            probs,
            profile_type,
            expected_regions: Vec::new(),
            region_start,
            extension: 0,
            ref_idx,
            contig_len,
        };

        result.to_regions(start_active, starts_and_stops);

        return result;
    }

    pub fn make_profile(&self) -> ActivityProfileType {
        match self.profile_type {
            ProfileType::Basic => {
                return ActivityProfileType::Basic(ActivityProfile::new(
                    MAX_PROB_PROPAGATION_DISTANCE,
                    ACTIVE_PROB_THRESHOLD,
                    self.ref_idx,
                    self.region_start.get_contig(),
                    self.contig_len,
                ))
            }
            ProfileType::BandPass => {
                return ActivityProfileType::BandPass(BandPassActivityProfile::new(
                    MAX_PROB_PROPAGATION_DISTANCE,
                    ACTIVE_PROB_THRESHOLD,
                    0,
                    0.01,
                    false,
                    self.ref_idx,
                    self.region_start.get_contig(),
                    self.contig_len,
                ))
            }
        }
    }

    fn to_regions(&mut self, mut is_active: bool, starts_and_stops: Vec<usize>) {
        let mut l = Vec::new();
        for i in 0..starts_and_stops.len() - 1 {
            let start = self.region_start.get_start() + starts_and_stops[i];
            let end = self.region_start.get_start() + starts_and_stops[i + 1] - 1;
            let active_loc = SimpleInterval::new(self.region_start.get_contig(), start, end);
            let r = AssemblyRegion::new(
                active_loc,
                is_active,
                self.extension,
                self.contig_len,
                self.region_start.get_contig(),
                self.ref_idx,
            );
            l.push(r);
            is_active = !is_active;
        }
        self.expected_regions = l;
    }
}

fn test_activity_profile(cfg: BasicActivityProfileTestProvider) {
    let mut profile = cfg.make_profile();

    match profile {
        ActivityProfileType::Basic(mut profile) => {
            assert!(profile.is_empty());

            for i in 0..cfg.probs.len() {
                let p = cfg.probs[i];
                let loc = SimpleInterval::new(
                    cfg.region_start.get_contig(),
                    cfg.region_start.get_start() + i,
                    cfg.region_start.get_start() + i,
                );
                profile.add(ActivityProfileState::new(loc, p, Type::None));
                assert!(
                    !profile.is_empty(),
                    "Profile should not be empty after adding state"
                );
            }

            assert_eq!(
                profile.region_start_loc.as_ref().unwrap(),
                &SimpleInterval::new(
                    cfg.region_start.get_contig(),
                    cfg.region_start.get_start(),
                    cfg.region_start.get_start()
                ),
                "Start loc should be the start of the region"
            );

            assert_eq!(
                profile.size(),
                cfg.probs.len(),
                "Should have exactly the number of states we expected to add"
            );

            assert_probs_are_equal(profile.state_list, cfg.probs);
            // TODO -- reanble tests
            //assertRegionsAreEqual(profile.createActiveRegions(0, 100), cfg.expectedRegions);
        }
        ActivityProfileType::BandPass(mut profile) => {
            assert!(profile.activity_profile.is_empty());

            for i in 0..cfg.probs.len() {
                let p = cfg.probs[i];
                let loc = SimpleInterval::new(
                    cfg.region_start.get_contig(),
                    cfg.region_start.get_start() + i,
                    cfg.region_start.get_start() + i,
                );
                profile
                    .activity_profile
                    .add(ActivityProfileState::new(loc, p, Type::None));
                assert!(
                    !profile.activity_profile.is_empty(),
                    "Profile should not be empty after adding state"
                );
            }

            assert_eq!(
                profile.activity_profile.region_start_loc.as_ref().unwrap(),
                &SimpleInterval::new(
                    cfg.region_start.get_contig(),
                    cfg.region_start.get_start(),
                    cfg.region_start.get_start()
                ),
                "Start loc should be the start of the region"
            );

            assert_eq!(
                profile.activity_profile.size(),
                cfg.probs.len(),
                "Should have exactly the number of states we expected to add"
            );

            assert_probs_are_equal(profile.activity_profile.state_list, cfg.probs);
        }
    }
}

fn assert_probs_are_equal(actual: Vec<ActivityProfileState>, expected: Vec<f32>) {
    assert_eq!(actual.len(), expected.len());

    for i in 0..actual.len() {
        assert_eq!(actual[i].is_active_prob(), expected[i])
    }
}

#[test]
fn qual_interval_tests() {
    for profile_type in vec![ProfileType::Basic, ProfileType::BandPass] {
        test_activity_profile(BasicActivityProfileTestProvider::new(
            profile_type,
            vec![1.0],
            true,
            SimpleInterval::new(0, 1, 100),
            100,
            0,
            vec![0, 1],
        ));
        test_activity_profile(BasicActivityProfileTestProvider::new(
            profile_type,
            vec![1.0, 0.0],
            true,
            SimpleInterval::new(0, 1, 100),
            100,
            0,
            vec![0, 1, 2],
        ));
        test_activity_profile(BasicActivityProfileTestProvider::new(
            profile_type,
            vec![0.0, 1.0],
            false,
            SimpleInterval::new(0, 1, 100),
            100,
            0,
            vec![0, 1, 2],
        ));
        test_activity_profile(BasicActivityProfileTestProvider::new(
            profile_type,
            vec![1.0, 0.0, 1.0],
            true,
            SimpleInterval::new(0, 1, 100),
            100,
            0,
            vec![0, 1, 2, 3],
        ));
        test_activity_profile(BasicActivityProfileTestProvider::new(
            profile_type,
            vec![1.0, 1.0, 1.0],
            true,
            SimpleInterval::new(0, 1, 100),
            100,
            0,
            vec![0, 3],
        ));
    }
}

// -------------------------------------------------------------------------------------
//
// Hardcore tests for adding to the profile and constructing active regions
//
// -------------------------------------------------------------------------------------
fn make_regions(total_region_size: usize, start_with_active: bool, n_parts: usize) -> Vec<bool> {
    let mut regions = Vec::new();

    let mut is_active = start_with_active;
    let active_region_size = max(total_region_size / n_parts, 1);
    for i in (0..total_region_size).step_by(active_region_size) {
        for j in 0..active_region_size {
            if j + i < total_region_size {
                break;
            };
            regions.push(is_active)
        }
        is_active = !is_active;
    }

    return regions;
}

fn test_region_creation(
    start: usize,
    mut probs: Vec<bool>,
    max_region_size: usize,
    n_parts: usize,
    force_conversion: bool,
    wait_until_end: bool,
    contig_len: usize,
) {
    let mut profile = ActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD,
        0,
        0,
        contig_len,
    );
    // let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);
    // let contig = ref_reader.index.sequences()[0].name.clone();
    let mut seen_sites = vec![false; probs.len()];

    let mut last_region = None;
    for i in 0..probs.len() {
        let is_active = probs[i];
        let loc = SimpleInterval::new(0, i + start, i + start);
        let mut state =
            ActivityProfileState::new(loc, if is_active { 1.0 } else { 0.0 }, Type::None);
        profile.add(state);

        if !wait_until_end {
            let regions = profile.clone().pop_ready_assembly_regions(0, 1, max_region_size, false);
            assert_good_regions(
                start,
                regions,
                max_region_size,
                &mut last_region,
                &mut probs,
                &mut seen_sites,
            );
        }
    }

    if wait_until_end || force_conversion {
        let regions = profile.clone().pop_ready_assembly_regions(0, 1, max_region_size, force_conversion);
        assert_good_regions(
            start,
            regions,
            max_region_size,
            &mut last_region,
            &mut probs,
            &mut seen_sites,
        );
    }

    for i in 0..probs.len() {
        if force_conversion
            || (i + max_region_size + profile.get_max_prob_propagation_distance() < probs.len())
        {
            // only require a site to be seen if we are forcing conversion or the site is more than maxRegionSize from the end
            assert!(seen_sites[i], "missed site {}", i);
        }
    }
}

fn assert_good_regions(
    start: usize,
    regions: Vec<AssemblyRegion>,
    max_region_size: usize,
    last_region: &mut Option<AssemblyRegion>,
    probs: &mut Vec<bool>,
    seen_sites: &mut Vec<bool>,
) {
    for region in regions {
        assert!(
            region.get_span().size() > 0,
            "Region {:?} has a bad size",
            region.get_span()
        );
        assert!(
            region.get_span().size() <= max_region_size,
            "Region {:?} has a bad size: It's bigger than the max size {}",
            region.get_span(),
            max_region_size
        );

        match last_region {
            Some(last_region) => {
                assert!(
                    region.get_span().get_start() == last_region.get_span().get_end() + 1,
                    "Region {:?} doesn't start immediately after the previous region {}",
                    &region,
                    last_region.get_span().get_end()
                );
            }
            None => {
                // pass
            }
        };

        // check that all active bases are actually active
        let region_offset = region.get_span().get_start() - start;
        assert!(
            region_offset >= 0 && region_offset < probs.len(),
            "Region {:?} has a bad offset w.r.t start",
            &region
        );
        for j in 0..region.get_span().size() {
            let site_offset = j + region_offset;
            assert_eq!(region.is_active(), probs[site_offset]);
            assert!(
                !seen_sites[site_offset],
                "Site {} in {:?} was seen already",
                j, &region
            );
            seen_sites[site_offset] = true;
        }

        *last_region = Some(region);
    }
}

#[test]
fn region_creation_tests() {
    let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);
    // ref_reader.fetch_all_by_rid(0);
    // ref_reader.index.sequences()[0].name
    // let mut seq = Vec::new();
    // ref_reader.read(&mut seq);
    let contig_len = ref_reader.index.sequences()[0].len as usize;
    for start in vec![1, 10, 100, contig_len - 100, contig_len - 10] {
        for region_size in vec![1, 10, 100, 1000, 10000] {
            for max_region_size in vec![10, 50, 200] {
                for wait_until_end in vec![false, true] {
                    for force_conversion in vec![false, true] {
                        // what do I really want to test here?  I'd like to test a few cases:
                        // -- region is all active (1.0)
                        // -- region is all inactive (0.0)
                        // -- cut the interval into 1, 2, 3, 4, 5 ... 10 regions, each with alternating activity values
                        for start_with_active in vec![true, false] {
                            for n_parts in vec![1, 2, 3, 4, 5, 7, 10, 11, 13] {
                                let region_size = min(region_size, contig_len - start);
                                let regions = make_regions(region_size, start_with_active, n_parts);
                                test_region_creation(
                                    start,
                                    regions,
                                    max_region_size,
                                    n_parts,
                                    force_conversion,
                                    wait_until_end,
                                    contig_len,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------------
//
// Hardcore tests for adding to the profile and constructing active regions
//
// -------------------------------------------------------------------------------------

fn test_soft_clips(
    start: usize,
    n_preceding_sites: usize,
    soft_clip_size: usize,
    contig_len: usize,
) {
    let mut profile = ActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD,
        0,
        0,
        contig_len,
    );
    for i in 0..n_preceding_sites {
        let loc = SimpleInterval::new(0, i + start, i + start);
        let state = ActivityProfileState::new(loc, 0.0, Type::None);
        profile.add(state);
    }

    let soft_clip_loc =
        SimpleInterval::new(0, n_preceding_sites + start, n_preceding_sites + start);
    profile.add(ActivityProfileState::new(
        soft_clip_loc.clone(),
        1.0,
        Type::HighQualitySoftClips(soft_clip_size as f32),
    ));

    let actual_num_of_soft_clips = min(soft_clip_size, profile.get_max_prob_propagation_distance());

    if n_preceding_sites == 0 {
        let profile_size = min(start + actual_num_of_soft_clips, contig_len) - start + 1;
        assert_eq!(
            profile.size(),
            profile_size,
            "Wrong number of states in the profile"
        );
    };

    for i in 0..profile.size() {
        let state = &profile.get_state_list()[i];
        let within_sc_range = state.get_loc().distance(&soft_clip_loc) <= actual_num_of_soft_clips;
        if within_sc_range {
            assert!(
                state.is_active_prob() > 0.0,
                "active prob should be changed within soft clip size"
            );
        } else {
            assert_eq!(
                state.is_active_prob(),
                0.0,
                "active prob should be changed within soft clip size"
            );
        }
    }
}

#[test]
fn run_test_soft_clips() {
    let mut ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);
    // ref_reader.fetch_all_by_rid(0);
    // ref_reader.index.sequences()[0].name
    // let mut seq = Vec::new();
    // ref_reader.read(&mut seq);
    let contig_len = ref_reader.index.sequences()[0].len as usize;
    for start in vec![
        1,
        10,
        100,
        contig_len - 100,
        contig_len - 10,
        contig_len - 1,
    ] {
        for preceding_sites in vec![0, 1, 10] {
            if preceding_sites + start < contig_len {
                for soft_clip_size in vec![1, 2, 10, 100] {
                    test_soft_clips(start, preceding_sites, soft_clip_size, contig_len)
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------------
//
// Tests to ensure we cut large active regions in the right place
//
// -------------------------------------------------------------------------------------

fn add_prob(l: &mut Vec<f32>, v: f32) {
    l.push(v)
}

fn make_gaussian(mean: f64, range: usize, sigma: f64) -> Vec<f64> {
    let mut gauss = vec![0.0; range];

    for iii in 0..range {
        gauss[iii] =
            MathUtils::normal_distribution(mean, sigma, iii as f64) + ACTIVE_PROB_THRESHOLD as f64;
    }

    return gauss;
}

fn find_cut_site_for_two_max_peaks(probs: &Vec<f32>, min_region_size: usize) -> Option<usize> {
    for i in ((min_region_size + 1)..=(probs.len() - 2))
        .into_iter()
        .rev()
    {
        let prev = probs[i - 1];
        let next = probs[i + 1];
        let cur = probs[i];
        if cur < next && cur < prev {
            return Some(i + 1);
        }
    }

    return None;
}

fn test_active_region_cuts(
    min_region_size: usize,
    max_region_size: usize,
    expected_region_size: usize,
    probs: Vec<f32>,
    contig_len: usize,
) {
    // let profile = ActivityProfileState::new()
    let mut profile = ActivityProfile::new(
        MAX_PROB_PROPAGATION_DISTANCE,
        ACTIVE_PROB_THRESHOLD as f32,
        0,
        0,
        contig_len,
    );
    for i in 0..=(max_region_size + profile.get_max_prob_propagation_distance()) {
        let loc = SimpleInterval::new(0, i + 1, i + 1);
        let prob = if i < probs.len() { probs[i] } else { 0.0 };
        let state = ActivityProfileState::new(loc, prob, Type::None);
        // println!("State: {:?}", &state);
        profile.add(state);
    }

    let regions = profile.pop_ready_assembly_regions(0, min_region_size, max_region_size, false);
    assert!(
        regions.len() >= 1,
        "Should only be one regions for this test"
    );
    let region = &regions[0];
    // println!("min size: {} max size {} expected size {} region {:?}", min_region_size, max_region_size, expected_region_size, &region);
    assert_eq!(region.get_span().get_start(), 1, "Region should start at 1");
    assert_eq!(
        region.get_span().size(),
        expected_region_size,
        "Incorrect region size; cut must have been incorrect"
    );
}

#[test]
fn make_active_region_cut_tests() {
    let ref_reader = ReferenceReaderUtils::generate_faidx(b37_reference_20_21);
    // ref_reader.fetch_all_by_rid(0);
    // ref_reader.index.sequences()[0].name
    // let mut seq = Vec::new();
    // ref_reader.read(&mut seq);
    let contig_len = ref_reader.index.sequences()[0].len as usize;

    for active_region_size in vec![10, 12, 20, 30, 40] {
        for min_region_size in vec![1, 5, 10] {
            let max_region_size = (active_region_size * 2) / 3;
            if min_region_size >= max_region_size {
                continue;
            };
            {
                // test flat activity profile
                let probs = vec![1.0; active_region_size];
                test_active_region_cuts(
                    min_region_size,
                    max_region_size,
                    max_region_size,
                    probs,
                    contig_len,
                );
            }

            {
                // test point profile is properly handled
                for end in 1..active_region_size {
                    let probs = vec![1.0; end];
                    test_active_region_cuts(
                        min_region_size,
                        max_region_size,
                        min(end, max_region_size),
                        probs,
                        contig_len,
                    );
                }
            }

            {
                // test increasing activity profile
                let mut probs = Vec::new();
                for i in 0..active_region_size {
                    add_prob(
                        &mut probs,
                        (1.0 * (i as f32 + 1.0)) / active_region_size as f32,
                    );
                }
                test_active_region_cuts(
                    min_region_size,
                    max_region_size,
                    max_region_size,
                    probs,
                    contig_len,
                );
            }

            {
                // test decreasing activity profile
                let mut probs = Vec::new();
                for i in 0..active_region_size {
                    add_prob(
                        &mut probs,
                        1.0 - (1.0 * (i as f32 + 1.0)) / active_region_size as f32,
                    );
                }
                test_active_region_cuts(
                    min_region_size,
                    max_region_size,
                    max_region_size,
                    probs,
                    contig_len,
                );
            }

            {
                // test two peaks
                for root_sigma in vec![1.0, 2.0, 3.0] {
                    for max_peak1 in 0..(active_region_size / 2) {
                        for max_peak2 in ((active_region_size / 2) + 1)..active_region_size {
                            let gauss1 =
                                make_gaussian(max_peak1 as f64, active_region_size, root_sigma);
                            let gauss2 = make_gaussian(
                                max_peak2 as f64,
                                active_region_size,
                                root_sigma + 1.0,
                            );
                            let mut probs = Vec::new();
                            for i in 0..active_region_size {
                                add_prob(&mut probs, gauss1[i] as f32 + gauss2[i] as f32);
                            }
                            let cut_site = find_cut_site_for_two_max_peaks(&probs, min_region_size);
                            match cut_site {
                                Some(cut_site) => {
                                    if cut_site < max_region_size {
                                        test_active_region_cuts(
                                            min_region_size,
                                            max_region_size,
                                            max(cut_site, min_region_size),
                                            probs,
                                            contig_len,
                                        )
                                    }
                                }
                                None => {
                                    // pass
                                }
                            }
                        }
                    }
                }
            }

            {
                // test that the lowest of two minima is taken
                // looks like a bunch of 1s, 0.5, some 1.0s, 0.75, some more 1s
                for first_min in 1..active_region_size {
                    for second_min in (first_min + 1)..active_region_size {
                        let mut probs = vec![1.0; active_region_size];
                        probs[first_min] = 0.5;
                        probs[second_min] = 0.75;
                        let mut expected_cut;
                        if first_min + 1 < min_region_size {
                            if first_min == second_min - 1 {
                                // edge case for non-min at minRegionSize
                                expected_cut = max_region_size;
                            } else {
                                expected_cut = if second_min + 1 > max_region_size {
                                    max_region_size
                                } else {
                                    if second_min + 1 < min_region_size {
                                        max_region_size
                                    } else {
                                        second_min + 1
                                    }
                                };
                            };
                        } else if first_min + 1 > max_region_size {
                            expected_cut = max_region_size;
                        } else {
                            expected_cut = first_min + 1;
                        };

                        min(first_min + 1, max_region_size);
                        test_active_region_cuts(
                            min_region_size,
                            max_region_size,
                            expected_cut,
                            probs,
                            contig_len,
                        );
                    }
                }
            }
        }
    }
}
