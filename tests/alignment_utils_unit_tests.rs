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
extern crate gkl;
extern crate itertools;
extern crate rand;
extern crate term;

use gkl::smithwaterman::{OverhangStrategy, Parameters};
use itertools::Itertools;
use lorikeet_genome::genotype::genotype_builder::Genotype;
use lorikeet_genome::genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use lorikeet_genome::genotype::genotype_likelihoods::GenotypeLikelihoods;
use lorikeet_genome::haplotype::haplotype::Haplotype;
use lorikeet_genome::model::allele_frequency_calculator::AlleleFrequencyCalculator;
use lorikeet_genome::model::allele_likelihoods::AlleleLikelihoods;
use lorikeet_genome::model::byte_array_allele::{Allele, ByteArrayAllele};
use lorikeet_genome::model::variant_context::VariantContext;
use lorikeet_genome::model::{allele_list::AlleleList, variants::SPAN_DEL_ALLELE};
use lorikeet_genome::pair_hmm::pair_hmm::PairHMM;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::{
    AVXMode, PairHMMInputScoreImputator,
};
use lorikeet_genome::processing::lorikeet_engine::ReadType;
use lorikeet_genome::reads::alignment_utils::AlignmentUtils;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reads::cigar_builder::CigarBuilder;
use lorikeet_genome::reads::cigar_utils::{
    CigarUtils, ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS,
};
use lorikeet_genome::smith_waterman::smith_waterman_aligner::{
    SmithWatermanAligner, SmithWatermanAlignmentResult, ORIGINAL_DEFAULT, STANDARD_NGS,
};
use lorikeet_genome::test_utils::read_likelihoods_unit_tester::ReadLikelihoodsUnitTester;
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::base_utils::BaseUtils;
use lorikeet_genome::utils::math_utils::{MathUtils, LOG10_ONE_HALF};
use lorikeet_genome::utils::quality_utils::QualityUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use lorikeet_genome::utils::utils::make_permutations;
use lorikeet_genome::GenomeExclusionTypes::GenomesAndContigsType;
use rand::rngs::ThreadRng;
use rand::seq::index::sample;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView, Seq};
use std::cmp::{max, min, Ordering};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::ops::Deref;
use std::sync::Mutex;

static DEBUG: bool = true;

struct AlignmentUtilsUnitTests {
    read_mapped: BirdToolRead,
    read_unmapped_flag: BirdToolRead,
    read_unknown_contig: BirdToolRead,
}

impl AlignmentUtilsUnitTests {
    fn init() -> Self {
        let read_mapped = create_mapped_read("mapped", 0);
        let mut read_unmapped = create_mapped_read("unmapped_flag", 1);
        read_unmapped.read.set_unmapped();

        let mut read_unknown_contig = create_mapped_read("unknown_contig", 2);
        read_unknown_contig.read.set_pos(3);
        read_unknown_contig.read.set_tid(-1);

        Self {
            read_mapped,
            read_unmapped_flag: read_unmapped,
            read_unknown_contig,
        }
    }
}

fn create_mapped_read(name: &str, start: i64) -> BirdToolRead {
    BirdToolRead::new(
        ArtificialReadUtils::create_artificial_read_default(
            name,
            0,
            start,
            ArtificialReadUtils::DEFAULT_READ_LENGTH,
            false,
        ),
        0,
        ReadType::Short,
    )
}

fn make_cigar_element_combinations() -> Vec<Vec<Cigar>> {
    // this functionality can be adapted to provide input data for whatever you might want in your data
    let mut cigar_elements = Vec::new();
    for size in vec![0, 10] {
        cigar_elements.push(Cigar::Match(size));
        cigar_elements.push(Cigar::Ins(size));
        cigar_elements.push(Cigar::Del(size));
        cigar_elements.push(Cigar::RefSkip(size));
        cigar_elements.push(Cigar::SoftClip(size));
        cigar_elements.push(Cigar::HardClip(size));
        cigar_elements.push(Cigar::Pad(size));
        cigar_elements.push(Cigar::Equal(size));
        cigar_elements.push(Cigar::Diff(size));
    }

    let mut combinations = Vec::new();
    for n_elements in vec![1, 2, 3] {
        combinations.extend(make_permutations(&cigar_elements, n_elements, true));
    }

    return combinations;
}

/******************************************************
 * Tests for AlignmentUtils.createReadAlignedToRef()
 ******************************************************/
fn make_haplotype_for_aligned_to_ref_test(bases: &str, cigar: &str) -> Haplotype<SimpleInterval> {
    let mut hap = Haplotype::new(bases.as_bytes(), false);
    hap.set_cigar(CigarString::try_from(cigar).unwrap().0);
    return hap;
}

fn make_read_for_aligned_to_ref_test(base_string: &str) -> BirdToolRead {
    let mut read = ArtificialReadUtils::create_artificial_read_default("my_read", 0, 0, 10, false);
    let bases = base_string.as_bytes();
    let cigar = CigarString(read.cigar().0.clone());
    read.set(
        "my_read".as_bytes(),
        Some(&cigar),
        bases,
        vec![30; bases.len()].as_slice(),
    );
    return BirdToolRead::new(read, 0, ReadType::Short);
}

fn test_read_aligned_to_ref(
    read: &BirdToolRead,
    haplotype: &Haplotype<SimpleInterval>,
    ref_haplotype: &Haplotype<SimpleInterval>,
    ref_start: usize,
    expected_read_start: i64,
    expected_read_cigar: Option<&str>,
) {
    let original_read_copy = read.clone();

    match expected_read_cigar {
        None => {
            // assert!(AlignmentUtils::create_read_aligned_to_ref(read, haplotype, ref_haplotype, ref_start, true))
        }
        Some(expected_read_cigar) => {
            let expected_cigar = CigarString::try_from(expected_read_cigar).unwrap();
            let aligned_read = AlignmentUtils::create_read_aligned_to_ref(
                read,
                haplotype,
                ref_haplotype,
                ref_start,
                true,
                AVXMode::detect_mode(),
            );
            println!(
                "aligned read {} expected cigar {}",
                aligned_read.read.cigar().to_string(),
                expected_cigar.to_string()
            );
            assert_eq!(aligned_read.read.qname(), original_read_copy.read.qname());
            assert_eq!(aligned_read.read.pos(), expected_read_start);
            assert_eq!(
                aligned_read.read.seq().as_bytes(),
                original_read_copy.read.seq().as_bytes()
            );
            assert_eq!(aligned_read.read.qual(), original_read_copy.read.qual());
            assert_eq!(
                aligned_read.read.cigar(),
                expected_cigar.into_view(expected_read_start)
            );
            assert!(aligned_read.read.aux("HC".as_bytes()).is_ok());
        }
    }
}

#[test]
fn make_read_aligned_to_ref_data() {
    let hap_bases = "ACTGAAGGTTCC";
    let all_M = make_haplotype_for_aligned_to_ref_test(hap_bases, &format!("{}M", hap_bases.len()));

    // make sure we get back a cigar of the right length
    for i in -1..hap_bases.len() as i64 {
        let mut read = make_read_for_aligned_to_ref_test(hap_bases);
        let mut bases = read.read.seq().as_bytes();
        if i != -1 {
            bases[i as usize] = b'A';
        };
        let qual = read.read.qual().to_vec();
        let cigar = read.read.cigar().take();
        let qname = read.read.qname().to_vec();
        read.update(&qname, Some(&cigar), bases, &qual);
        test_read_aligned_to_ref(
            &read,
            &all_M,
            &all_M,
            10,
            10,
            Some(all_M.get_cigar().to_string().as_str()),
        );
    }

    // make sure insertions at the front are correctly handled
    for pad_front in 1..10 {
        println!("Padding {}", pad_front);
        let read =
            make_read_for_aligned_to_ref_test(&format!("{}{}", "N".repeat(pad_front), hap_bases));
        test_read_aligned_to_ref(
            &read,
            &all_M,
            &all_M,
            10,
            10,
            Some(&format!(
                "{}I{}",
                pad_front,
                all_M.get_cigar().to_string().as_str()
            )),
        );
    }

    // make sure insertions at the back are correctly handled
    for pad_back in 1..10 {
        let read =
            make_read_for_aligned_to_ref_test(&format!("{}{}", hap_bases, "N".repeat(pad_back)));
        test_read_aligned_to_ref(
            &read,
            &all_M,
            &all_M,
            10,
            10,
            Some(&format!(
                "{}{}I",
                all_M.get_cigar().to_string().as_str(),
                pad_back
            )),
        );
    }

    // make sure refStart and hapStart are respected
    for ref_start in 1..10 {
        for hap_start in ref_start..(10 + ref_start) {
            let mut hap = Haplotype::new(all_M.get_bases(), false);
            hap.set_cigar(all_M.get_cigar().0.clone());
            hap.set_alignment_start_hap_wrt_ref(hap_start);

            let read = make_read_for_aligned_to_ref_test(hap_bases);
            test_read_aligned_to_ref(
                &read,
                &hap,
                &all_M,
                ref_start,
                (ref_start + hap_start) as i64,
                Some(all_M.get_cigar().to_string().as_str()),
            );
        }
    }

    // example case of bad alignment because SW doesn't necessarily left-align indels
    {
        let hap = "ACTGTGGGTTCCTCTTATTTTATTTCTACATCAATGTTCATATTTAACTTATTATTTTATCTTATTTTTAAATTTCTTTTATGTTGAGCCTTGATGAAAGCCATAGGTTCTCTCATATAATTGTATGTGTATGTATGTATATGTACATAATATATACATATATGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATATATACATATATG";
        let hap_cigar = format!("{}M", hap.len());
        let read_bases = "ATGTACATAATATATACATATATGTATATGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATAT";
        let read = make_read_for_aligned_to_ref_test(read_bases);
        let ref_start = 10130100;
        let hap_start = 500;
        let bad_cigar = "31M6D211M";
        let good_cigar = "28M6D214M";
        let mut bad_hap = Haplotype::new(hap.as_bytes(), false);
        bad_hap.set_cigar(CigarString::try_from(hap_cigar.as_str()).unwrap().0);
        bad_hap.set_alignment_start_hap_wrt_ref(hap_start);

        let expected_pos = 10130740;
        test_read_aligned_to_ref(
            &read,
            &bad_hap,
            &bad_hap,
            ref_start,
            expected_pos,
            Some(good_cigar),
        );
    }

    // example where left-align generates a leading deletion
    {
        let reference = "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAAAAGTGAGTGTGAATCTGGTGGCTGTGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
        // ref-to-hap:   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        let haplotype = "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAAAAGTGAGTGTGAATCTGGTGGCTGCGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
        // hap-to-read:                                                                                                                                                            ||||.|||||||||||||||||
        // aligned read:                                                                                                                                                           GCTGCTTTTGGTGTGGCTCTTT
        let read_bases = "GCTGCTTTTGGTGTGGCTCTTT";
        let hap_cigar = format!("{}M", haplotype.len());
        let read = make_read_for_aligned_to_ref_test(read_bases);
        let ref_start = 215239171;
        let hap_start = 575;
        let alignment_offset = 154;
        let good_cigar = "22M";
        let mut bad_hap = Haplotype::new(haplotype.as_bytes(), false);
        bad_hap.set_cigar(CigarString::try_from(hap_cigar.as_str()).unwrap().0);
        bad_hap.set_alignment_start_hap_wrt_ref(hap_start);
        let ref_hap = make_haplotype_for_aligned_to_ref_test(
            reference,
            format!("{}M", reference.len()).as_str(),
        );

        let expected_pos = ref_start + hap_start + alignment_offset;
        test_read_aligned_to_ref(
            &read,
            &bad_hap,
            &ref_hap,
            ref_start,
            expected_pos as i64,
            Some(good_cigar),
        );
    }

    // example where the haplotype has an indel relative to reference
    {
        let reference = "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACAGGGGAATCCCGAAGAAATGGTGGGTCCTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        // ref-to-hap:   |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||^^||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        let haplotype = "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACGGGAATCCCGAAGAAATGGTGGGTCCTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        // hap-to-read:                                                                                                                           .|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        // aligned read:                                                                                                                          CCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC
        let read_bases = "CCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        let hap_cigar = "93M2D92M";
        let read = make_read_for_aligned_to_ref_test(read_bases);
        let ref_start = 13011;
        let hap_start = 553;
        let alignment_offset = 123;
        let good_cigar = "64M";
        let mut bad_hap = Haplotype::new(haplotype.as_bytes(), false);
        bad_hap.set_cigar(CigarString::try_from(hap_cigar).unwrap().0);
        bad_hap.set_alignment_start_hap_wrt_ref(hap_start);

        let reference_hap = make_haplotype_for_aligned_to_ref_test(
            reference,
            format!("{}M", reference.len()).as_str(),
        );
        let expected_pos = ref_start + hap_start + alignment_offset;
        test_read_aligned_to_ref(
            &read,
            &bad_hap,
            &reference_hap,
            ref_start,
            expected_pos as i64,
            Some(good_cigar),
        );
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct Mutation {
    pos: usize,
    len: usize,
    operator: Cigar,
}

impl Mutation {
    fn new(pos: usize, len: usize, operator: Cigar) -> Self {
        Self { pos, len, operator }
    }

    fn get_n_mismatches(&self) -> usize {
        self.len
    }

    fn apply(&self, seq: &str) -> String {
        match self.operator {
            Cigar::Match(_) => {
                let mut bases = seq.as_bytes().to_vec();
                if self.pos < seq.len() {
                    bases[self.pos] = if bases[self.pos] == b'A' { b'C' } else { b'A' };
                };

                return std::str::from_utf8(bases.as_slice()).unwrap().to_string();
            }
            Cigar::Ins(_) => {
                let prefix = &seq[0..self.pos];
                let postfix = &seq[self.pos..seq.len()];
                return format!("{}{}{}", prefix, &"GTCAGTTA"[0..self.len], postfix);
            }
            Cigar::Del(_) => {
                let prefix = &seq[0..self.pos];
                let postfix = &seq[self.pos + self.len..seq.len()];
                return format!("{}{}", prefix, postfix);
            }
            _ => {
                panic!("Unexpected operator")
            }
        }
    }
}

impl Ord for Mutation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos.cmp(&other.pos)
    }
}

impl PartialOrd for Mutation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct MutatedSequence {
    num_mismatches: usize,
    seq: String,
}

impl MutatedSequence {
    fn new(num_mismatches: usize, seq: String) -> Self {
        Self {
            num_mismatches,
            seq,
        }
    }
}

fn mutate_sequence(hap_in: &str, mut mutations: Vec<Mutation>) -> MutatedSequence {
    mutations.sort();
    let mut mismatches = 0;
    let mut hap = hap_in.to_string();
    for mutation in mutations {
        hap = mutation.apply(hap.as_str());
        mismatches += mutation.get_n_mismatches();
    }

    return MutatedSequence::new(mismatches, hap);
}

fn test_read_aligned_to_ref_complex_alignment(
    test_index: usize,
    read: &BirdToolRead,
    reference: &str,
    haplotype: &Haplotype<SimpleInterval>,
    expected_max_mismatches: usize,
) {
    let ref_haplotype = Haplotype::new(reference.as_bytes(), true);
    let aligned_read = AlignmentUtils::create_read_aligned_to_ref(
        read,
        haplotype,
        &ref_haplotype,
        0,
        true,
        AVXMode::detect_mode(),
    );
    let mismatches = AlignmentUtils::get_mismatch_count(
        &aligned_read,
        reference.as_bytes(),
        aligned_read.get_start(),
        0,
        read.len(),
    );
    assert!(
        mismatches.num_mismatches <= expected_max_mismatches,
        "Alignment of read to ref looks broken. Expected max {} got {} for read bases {} with cigar {} reference {} haplotype with cigar {} aligned read cigar {} @ {}",
        expected_max_mismatches, mismatches.num_mismatches,
        std::str::from_utf8(read.read.seq().as_bytes().as_slice()).unwrap(), read.read.cigar().to_string(), reference, haplotype.get_cigar().to_string(), aligned_read.read.cigar().to_string(), read.get_start()
    );
}

#[test]
fn make_complex_read_aligned_to_ref() {
    let all_mutations = vec![
        Mutation::new(1, 1, Cigar::Match(0)),
        Mutation::new(2, 1, Cigar::Match(0)),
        Mutation::new(3, 1, Cigar::Ins(0)),
        Mutation::new(7, 1, Cigar::Del(0)),
    ];

    let mut i = 0;
    let reference_bases = "ACTGACTGACTG";
    let padded_reference = format!("NNNN{}NNNN", reference_bases);
    for mutations in make_permutations(&all_mutations, 3, false) {
        let hap = mutate_sequence(reference_bases, mutations);
        let mut haplotype = Haplotype::new(hap.seq.as_bytes(), false);
        println!("Padded ref {}", &padded_reference);
        let align = SmithWatermanAligner::align(
            padded_reference.as_bytes(),
            hap.seq.as_bytes(),
            &*ORIGINAL_DEFAULT,
            OverhangStrategy::SoftClip,
            AVXMode::detect_mode(),
        );
        haplotype.set_alignment_start_hap_wrt_ref(align.get_alignment_offset() as usize);
        haplotype.set_cigar(align.get_cigar().0);

        for read_mutations in make_permutations(&all_mutations, 3, false) {
            let read_bases = mutate_sequence(hap.seq.as_str(), read_mutations);
            let read = make_read_for_aligned_to_ref_test(read_bases.seq.as_str());
            test_read_aligned_to_ref_complex_alignment(
                i,
                &read,
                padded_reference.as_str(),
                &haplotype,
                hap.num_mismatches + read_bases.num_mismatches,
            );
            i += 1;
        }
    }
}

/**********************************************************
 * End of Tests for AlignmentUtils::create_read_aligned_to_ref()
 **********************************************************/

//////////////////////////////////////////
// Test AlignmentUtils.leftAlignIndel() //
//////////////////////////////////////////

fn test_with_clips_and_reference_context(
    ref_string: &str,
    read_string: &str,
    original_cigar: &str,
    expected_cigar: &str,
) {
    for leading_hard_clips in vec![0, 5] {
        for trailing_hard_clips in vec![0, 5] {
            for leading_soft_clips in vec![0, 5] {
                for trailing_soft_clips in vec![0, 5] {
                    for extra_ref_in_front in vec![0, 10] {
                        for extra_ref_in_back in vec![0, 10] {
                            let mut read_bases =
                                vec![
                                    0;
                                    read_string.len() + leading_soft_clips + trailing_soft_clips
                                ];
                            let mut ref_bases =
                                vec![0; ref_string.len() + extra_ref_in_front + extra_ref_in_back];
                            let read_len = read_bases.len();
                            let ref_len = ref_bases.len();
                            BaseUtils::fill_with_random_bases(
                                &mut read_bases,
                                0,
                                leading_soft_clips,
                            );
                            BaseUtils::fill_with_random_bases(
                                &mut read_bases,
                                leading_soft_clips + read_string.len(),
                                read_len,
                            );
                            read_bases[leading_soft_clips..leading_soft_clips + read_string.len()]
                                .clone_from_slice(read_string.as_bytes());

                            BaseUtils::fill_with_random_bases(
                                &mut ref_bases,
                                0,
                                extra_ref_in_front,
                            );
                            BaseUtils::fill_with_random_bases(
                                &mut ref_bases,
                                extra_ref_in_front + ref_string.len(),
                                ref_len,
                            );
                            ref_bases[extra_ref_in_front..extra_ref_in_front + ref_string.len()]
                                .clone_from_slice(ref_string.as_bytes());

                            let original_cigar_with_clips = format!(
                                "{}{}{}{}{}",
                                if leading_hard_clips > 0 {
                                    format!("{}H", leading_hard_clips)
                                } else {
                                    "".to_string()
                                },
                                if leading_soft_clips > 0 {
                                    format!("{}S", leading_soft_clips)
                                } else {
                                    "".to_string()
                                },
                                original_cigar,
                                if trailing_soft_clips > 0 {
                                    format!("{}S", trailing_soft_clips)
                                } else {
                                    "".to_string()
                                },
                                if trailing_hard_clips > 0 {
                                    format!("{}H", trailing_hard_clips)
                                } else {
                                    "".to_string()
                                },
                            );

                            let expected_cigar_with_clips = format!(
                                "{}{}{}{}{}",
                                if leading_hard_clips > 0 {
                                    format!("{}H", leading_hard_clips)
                                } else {
                                    "".to_string()
                                },
                                if leading_soft_clips > 0 {
                                    format!("{}S", leading_soft_clips)
                                } else {
                                    "".to_string()
                                },
                                expected_cigar,
                                if trailing_soft_clips > 0 {
                                    format!("{}S", trailing_soft_clips)
                                } else {
                                    "".to_string()
                                },
                                if trailing_hard_clips > 0 {
                                    format!("{}H", trailing_hard_clips)
                                } else {
                                    "".to_string()
                                },
                            );

                            let result = AlignmentUtils::left_align_indels(
                                CigarString::try_from(original_cigar_with_clips.as_str()).unwrap(),
                                ref_bases.as_slice(),
                                read_bases.as_slice(),
                                extra_ref_in_front as u32,
                            )
                            .cigar;

                            assert_eq!(result.to_string(), expected_cigar_with_clips);
                        }
                    }
                }
            }
        }
    }
}

fn test_left_align_indel(
    reference: &str,
    read: &str,
    original_cigar: &str,
    left_aligned_cigar: &str,
) {
    test_with_clips_and_reference_context(reference, read, original_cigar, left_aligned_cigar)
}

#[test]
fn make_left_align_indel_data() {
    // nothing happens when there is no indel
    test_left_align_indel("ACGT", "ACGT", "4M", "4M");
    test_left_align_indel("ACCT", "ACGT", "4M", "4M");
    test_left_align_indel("ACGT", "ACAT", "2M1X1M", "2M1X1M");

    // one insertion already left-aligned
    test_left_align_indel("AAATTT", "AAACCCTTT", "3M3I3M", "3M3I3M");
    test_left_align_indel("CCCTTT", "AAACCCTTT", "3I6M", "3I6M");
    test_left_align_indel("AAACCC", "AAACCCTTT", "6M3I", "6M3I");
    test_left_align_indel("AAACCC", "AAACCGTTT", "6M3I", "6M3I");

    // one deletion already left-aligned
    test_left_align_indel("AAACCCTTT", "AAATTT", "3M3D3M", "3M3D3M");

    // one insertion not left-aligned in homopolymer and STR with greater unit length
    test_left_align_indel("AAACCCTTT", "AAACCCCCCTTT", "5M3I4M", "3M3I6M");
    test_left_align_indel("AAACCCTTT", "AAACCCCCCTTT", "6M3I3M", "3M3I6M");
    test_left_align_indel("AAACCCTTT", "AAGCCCCCCTGT", "6M3I3M", "3M3I6M");
    test_left_align_indel("AAACGCGCGCGTTT", "AAACGCGCGCGCGCGTTT", "7M4I7M", "3M4I11M");
    test_left_align_indel("CCGCCG", "CCGCCGCCG", "6M3I", "3I6M");
    test_left_align_indel("ACCGCCG", "TCCGCCGCCG", "7M3I", "1M3I6M");

    // one deletion not left-aligned in homopolymer and STR with greater unit length
    test_left_align_indel("AAACCCCCCTTT", "AAACCCTTT", "5M3D4M", "3M3D6M");
    test_left_align_indel("AAACCCCCCTTT", "AAACCCTTT", "6M3D3M", "3M3D6M");
    test_left_align_indel("AAACGCGCGCGCGCGTTT", "AAACGCGCGCGTTT", "7M4D7M", "3M4D11M");

    //multiple separated indels
    test_left_align_indel(
        "AAACCCTTTGGGAAA",
        "AAACCCCCCTTTGGGGGGAAA",
        "6M3I6M3I3M",
        "3M3I6M3I6M",
    );
    test_left_align_indel(
        "AAACCCTTTGGGGGGAAA",
        "AAACCCCCCTTTGGGAAA",
        "6M3I6M3D3M",
        "3M3I6M3D6M",
    );

    // multiple indels in the same STR that combine or cancel
    test_left_align_indel("AAACCCCCTTT", "AAACCCCCTTT", "4M3I3D4M", "11M");
    test_left_align_indel("AAACCCCCTTT", "AAACCCCCTTT", "4M3D3I4M", "11M");
    test_left_align_indel("AAACCCCCTTT", "AAACCCCCTTT", "3M3I2M3D3M", "11M");
    test_left_align_indel("AACGCGCGCGTT", "AACGCGCGCGCGCGTT", "2M2I8M2I2M", "2M4I10M");
    test_left_align_indel("AACGCGCGCGCGCGTT", "AACGCGCGCGTT", "2M2D8M2D2M", "2M4D10M");
}

//////////////////////////////////////////
// Test AlignmentUtils.trimCigarByReference() //
//////////////////////////////////////////
fn test_trim_cigar(cigar_string: &str, start: usize, length: usize, expected_cigar_string: &str) {
    let cigar = CigarString::try_from(cigar_string).unwrap();
    let expected_cigar_raw = CigarString::try_from(expected_cigar_string).unwrap();

    // trimming throws error if all but deletion elements are trimmed
    if expected_cigar_raw.len() == 1
        && CigarUtils::cigar_elements_are_same_type(&expected_cigar_raw.0[0], &Some(Cigar::Del(0)))
    {
        // pass
    } else {
        let mut expected_cigar = CigarBuilder::new(true);
        expected_cigar.add_all(expected_cigar_raw.0);
        let expected_cigar = expected_cigar.make(false).unwrap();
        let actual_cigar =
            AlignmentUtils::trim_cigar_by_reference(cigar, start as u32, length as u32).cigar;
        assert_eq!(actual_cigar, expected_cigar);
    }
}

#[test]
fn make_trim_cigar_data() {
    for op in vec!["D", "=", "X", "M"] {
        for my_length in 1..6 {
            for start in 0..my_length - 1 {
                for end in start..my_length {
                    let length = end - start + 1;

                    for pad_op in vec!["D", "M"] {
                        for left_pad in 0..2 {
                            for right_pad in 0..2 {
                                test_trim_cigar(
                                    format!(
                                        "{}{}{}{}",
                                        if left_pad > 0 {
                                            format!("{}{}", left_pad, pad_op)
                                        } else {
                                            "".to_string()
                                        },
                                        my_length,
                                        op,
                                        if right_pad > 0 {
                                            format!("{}{}", right_pad, pad_op)
                                        } else {
                                            "".to_string()
                                        },
                                    )
                                    .as_str(),
                                    start + left_pad,
                                    end + left_pad,
                                    format!("{}{}", length, op).as_str(),
                                )
                            }
                        }
                    }
                }
            }
        }
    }

    for left_pad in vec![0, 1, 2, 5] {
        for right_pad in vec![0, 1, 2, 5] {
            let length = left_pad + right_pad;
            if length > 0 {
                for ins_size in vec![1, 10] {
                    for start in 0..=left_pad {
                        for stop in left_pad..length {
                            let left_pad_remaining = left_pad - start;
                            let right_pad_remaining = stop - left_pad + 1;
                            let ins_c = format!("{}I", ins_size);
                            test_trim_cigar(
                                format!("{}M{}{}M", left_pad, ins_c, right_pad).as_str(),
                                start,
                                stop,
                                format!(
                                    "{}{}{}",
                                    if left_pad_remaining > 0 {
                                        format!("{}M", left_pad_remaining)
                                    } else {
                                        "".to_string()
                                    },
                                    ins_c,
                                    if right_pad_remaining > 0 {
                                        format!("{}M", right_pad_remaining)
                                    } else {
                                        "".to_string()
                                    }
                                )
                                .as_str(),
                            )
                        }
                    }
                }
            }
        }
    }

    test_trim_cigar("3M2D4M", 0, 8, "3M2D4M");
    test_trim_cigar("3M2D4M", 2, 8, "1M2D4M");
    test_trim_cigar("3M2D4M", 2, 6, "1M2D2M");
    test_trim_cigar("3M2D4M", 3, 6, "2D2M");
    test_trim_cigar("3M2D4M", 4, 6, "1D2M");
    test_trim_cigar("3M2D4M", 5, 6, "2M");
    test_trim_cigar("3M2D4M", 6, 6, "1M");

    test_trim_cigar("2M3I4M", 0, 5, "2M3I4M");
    test_trim_cigar("2M3I4M", 1, 5, "1M3I4M");
    test_trim_cigar("2M3I4M", 1, 4, "1M3I3M");
    test_trim_cigar("2M3I4M", 2, 4, "3I3M");
    test_trim_cigar("2M3I4M", 2, 3, "3I2M");
    test_trim_cigar("2M3I4M", 2, 2, "3I1M");
    test_trim_cigar("2M3I4M", 3, 4, "2M");
    test_trim_cigar("2M3I4M", 3, 3, "1M");
    test_trim_cigar("2M3I4M", 4, 4, "1M");

    // this doesn't work -- but I'm not sure it should
    // test_trim_cigar("2M3I4M", 2, 1, "3I");
}

fn test_trim_cigar_by_bases(
    cigar_string: &str,
    start: usize,
    length: usize,
    expected_cigar_string: &str,
) {
    let cigar = CigarString::try_from(cigar_string).unwrap();
    let expected_cigar = CigarString::try_from(expected_cigar_string).unwrap();
    let actual_cigar =
        AlignmentUtils::trim_cigar_by_bases(cigar, start as u32, length as u32).cigar;
    assert_eq!(actual_cigar, expected_cigar);
}

#[test]
fn make_trim_cigar_by_bases_data() {
    test_trim_cigar_by_bases("2M3I4M", 0, 8, "2M3I4M");
    test_trim_cigar_by_bases("2M3I4M", 1, 8, "1M3I4M");
    test_trim_cigar_by_bases("2M3I4M", 2, 8, "3I4M");
    test_trim_cigar_by_bases("2M3I4M", 3, 8, "2I4M");
    test_trim_cigar_by_bases("2M3I4M", 4, 8, "1I4M");
    test_trim_cigar_by_bases("2M3I4M", 4, 7, "1I3M");
    test_trim_cigar_by_bases("2M3I4M", 4, 6, "1I2M");
    test_trim_cigar_by_bases("2M3I4M", 4, 5, "1I1M");
    test_trim_cigar_by_bases("2M3I4M", 4, 4, "1I");
    test_trim_cigar_by_bases("2M3I4M", 5, 5, "1M");
    test_trim_cigar_by_bases("2M2D2I", 0, 3, "2M2I");
    test_trim_cigar_by_bases("2M2D2I", 1, 3, "1M2I");
    test_trim_cigar_by_bases("2M2D2I", 2, 3, "2I");
    test_trim_cigar_by_bases("2M2D2I", 3, 3, "1I");
    test_trim_cigar_by_bases("2M2D2I", 2, 2, "1I");
    test_trim_cigar_by_bases("2M2D2I", 1, 2, "1M1I");
    test_trim_cigar_by_bases("2M2D2I", 0, 1, "2M");
    test_trim_cigar_by_bases("2M2D2I", 1, 1, "1M");
}

fn test_apply_cigar_to_cigar(
    first_to_second_string: &str,
    second_to_third_string: &str,
    expected_cigar_string: &str,
) {
    let first_to_second = CigarString::try_from(first_to_second_string).unwrap();
    let second_to_third = CigarString::try_from(second_to_third_string).unwrap();
    let expected_cigar = CigarString::try_from(expected_cigar_string).unwrap();
    let actual_cigar = AlignmentUtils::apply_cigar_to_cigar(&first_to_second, &second_to_third);
    assert_eq!(actual_cigar, expected_cigar);
}

#[test]
fn make_apply_cigar_to_cigar_data() {
    for i in 1..5 {
        test_apply_cigar_to_cigar(
            format!("{}M", i).as_str(),
            format!("{}M", i).as_str(),
            format!("{}M", i).as_str(),
        );
    }

    //        * ref   : ACGTAC
    //        * hap   : AC---C  - 2M3D1M
    //        * read  : AC---C  - 3M
    //        * result: AG---C => 2M3D
    test_apply_cigar_to_cigar("3M", "2M3D1M", "2M3D1M");

    //        * ref   : ACxG-TA
    //        * hap   : AC-G-TA  - 2M1D3M
    //        * read  : AC-GxTA  - 3M1I2M
    //        * result: AC-GxTA => 2M1D1M1I2M
    test_apply_cigar_to_cigar("3M1I2M", "2M1D3M", "2M1D1M1I2M");

    //        * ref   : ACGTA
    //        * hap   : A-GTA  - 1M1D3M
    //        * read  : A--TA  - 1M1D2M
    //        * result: A--TA => 1M2D2M
    test_apply_cigar_to_cigar("1M1D2M", "1M1D3M", "1M2D2M");

    //        * ref   : ACG-TA
    //        * hap   : A-GxTA  - 1M1D1M1I2M
    //        * read  : A---TA  - 1M2D2M
    //        * result: A---TA => 1M2D2M
    test_apply_cigar_to_cigar("1M2D2M", "1M1D1M1I2M", "1M2D2M");

    //        * ref   : A-CGTA
    //        * hap   : A-CGTA  - 5M
    //        * read  : AxCGTA  - 1M1I4M
    //        * result: AxCGTA => 1M1I4M
    test_apply_cigar_to_cigar("1M1I4M", "5M", "1M1I4M");

    //        * ref   : ACGTA
    //        * hap   : ACGTA  - 5M
    //        * read  : A--TA  - 1M2D2M
    //        * result: A--TA => 1M2D2M
    test_apply_cigar_to_cigar("1M2D2M", "5M", "1M2D2M");

    //        * ref   : AC-GTA
    //        * hap   : ACxGTA  - 2M1I3M
    //        * read  : A--GTA  - 1M2D3M
    //        * result: A--GTA => 1M1D3M
    test_apply_cigar_to_cigar("108M14D24M2M18I29M92M1000M", "2M1I3M", "2M1I3M");
}

fn test_append_clipped_elements_from_original_cigar(
    original_cigar_string: &str,
    shifted_cigar_string: &str,
    expected_string: &str,
) {
    let original_cigar = CigarString::try_from(original_cigar_string).unwrap();
    let shifted_cigar = CigarString::try_from(shifted_cigar_string).unwrap();
    let expected_cigar = CigarString::try_from(expected_string).unwrap();
    let actual_cigar = AlignmentUtils::append_clipped_elements_from_cigar_to_cigar(
        shifted_cigar,
        original_cigar.into_view(0),
    );

    assert_eq!(actual_cigar, expected_cigar);
}

#[test]
fn make_test_append_clipped_elements_from_original_cigar() {
    test_append_clipped_elements_from_original_cigar("30M", "30M", "30M");
    test_append_clipped_elements_from_original_cigar("30M", "15M6I15M", "15M6I15M");
    test_append_clipped_elements_from_original_cigar("5S30M", "30M", "5S30M");
    test_append_clipped_elements_from_original_cigar("5H30M", "30M", "5H30M");
    test_append_clipped_elements_from_original_cigar("5H5S30M", "30M", "5H5S30M");
    test_append_clipped_elements_from_original_cigar("30M5H", "30M", "30M5H");
    test_append_clipped_elements_from_original_cigar("30M5S", "30M", "30M5S");
    test_append_clipped_elements_from_original_cigar("10H30M5S5H", "30M", "10H30M5S5H");
    test_append_clipped_elements_from_original_cigar(
        "10H10M6D6M6D50M5S5H",
        "10M6I50M",
        "10H10M6I50M5S5H",
    );
}

fn test_get_bases_covering_ref_interval(
    bases_string: &str,
    ref_start: usize,
    ref_end: usize,
    cigar_string: &str,
    expected: Option<&str>,
) {
    let actual_bytes = AlignmentUtils::get_bases_covering_ref_interval(
        ref_start,
        ref_end,
        bases_string.as_bytes(),
        0,
        &CigarString::try_from(cigar_string).unwrap(),
    );

    match expected {
        None => {
            assert!(actual_bytes.is_none())
        }
        Some(expected) => {
            assert_eq!(
                std::str::from_utf8(actual_bytes.unwrap()).unwrap(),
                expected
            )
        }
    }
}

#[test]
fn make_get_bases_covering_ref_interval_data() {
    // matches
    // 0123
    // ACGT
    test_get_bases_covering_ref_interval("ACGT", 0, 3, "4M", Some("ACGT"));
    test_get_bases_covering_ref_interval("ACGT", 1, 3, "4M", Some("CGT"));
    test_get_bases_covering_ref_interval("ACGT", 1, 2, "4M", Some("CG"));
    test_get_bases_covering_ref_interval("ACGT", 1, 1, "4M", Some("C"));

    // deletions
    // 012345
    // AC--GT
    test_get_bases_covering_ref_interval("ACGT", 0, 5, "2M2D2M", Some("ACGT"));
    test_get_bases_covering_ref_interval("ACGT", 1, 5, "2M2D2M", Some("CGT"));
    test_get_bases_covering_ref_interval("ACGT", 2, 5, "2M2D2M", None);
    test_get_bases_covering_ref_interval("ACGT", 3, 5, "2M2D2M", None);
    test_get_bases_covering_ref_interval("ACGT", 4, 5, "2M2D2M", Some("GT"));
    test_get_bases_covering_ref_interval("ACGT", 5, 5, "2M2D2M", Some("T"));
    test_get_bases_covering_ref_interval("ACGT", 0, 4, "2M2D2M", Some("ACG"));
    test_get_bases_covering_ref_interval("ACGT", 0, 3, "2M2D2M", None);
    test_get_bases_covering_ref_interval("ACGT", 0, 2, "2M2D2M", None);
    test_get_bases_covering_ref_interval("ACGT", 0, 1, "2M2D2M", Some("AC"));
    test_get_bases_covering_ref_interval("ACGT", 0, 0, "2M2D2M", Some("A"));

    // insertions
    // 01--23
    // ACTTGT
    test_get_bases_covering_ref_interval("ACTTGT", 0, 3, "2M2I2M", Some("ACTTGT"));
    test_get_bases_covering_ref_interval("ACTTGT", 1, 3, "2M2I2M", Some("CTTGT"));
    test_get_bases_covering_ref_interval("ACTTGT", 2, 3, "2M2I2M", Some("GT"));
    test_get_bases_covering_ref_interval("ACTTGT", 3, 3, "2M2I2M", Some("T"));
    test_get_bases_covering_ref_interval("ACTTGT", 0, 2, "2M2I2M", Some("ACTTG"));
    test_get_bases_covering_ref_interval("ACTTGT", 0, 1, "2M2I2M", Some("AC"));
    test_get_bases_covering_ref_interval("ACTTGT", 1, 2, "2M2I2M", Some("CTTG"));
    test_get_bases_covering_ref_interval("ACTTGT", 2, 2, "2M2I2M", Some("G"));
    test_get_bases_covering_ref_interval("ACTTGT", 1, 1, "2M2I2M", Some("C"));

    // leading and terminal insertions - test that they are excluded
    test_get_bases_covering_ref_interval("ACTTGT", 0, 3, "2I4M", Some("TTGT"));
    test_get_bases_covering_ref_interval("ACTTGT", 0, 3, "4M2I", Some("ACTT"));
    test_get_bases_covering_ref_interval("ACGT", 0, 1, "2M2I", Some("AC"));
    test_get_bases_covering_ref_interval("ACGT", 1, 1, "2M2I", Some("C"));
    test_get_bases_covering_ref_interval("ACGT", 0, 0, "2M2I", Some("A"));

    // Weird edge case breaking function
    test_get_bases_covering_ref_interval(
        "GATGAAAATGACCTGCCCCCCCGTATCAGAAAGAACTATTGGAACTCCAAGGAACGAAGAGTGGACAGGCCTATAATCATCTCAGAAAACGCCGTTGACAGAATCTACGAGATTCCGCACGCATGGCGTTTGATCAAAACGTTGCCATTTGTAGTATTCGAGGAATTTGGCGCCCGGATCAACCTTACGGTACTCGACCTGGCTGCCAAGTGGTTTGCCACGCAAGACTCACTATCTCGGATCAACCAAAATCCTGCCTTGGCTTTCTACTACTCAAGAAATGACTCGCTGGACGTCGATTACGAACAGGTGCGTCGGCTGAACAAGCCTAAAGAAATGGACAGGGAGTGGATATTTTCTTTAATAGACGAGATCGCCGGCGAGGGAGAAGAGCAGAAGACATGCTCAAGCTCGTCGCAATCTCTGCGCCTGAAGACGTTCGTGAGAGTCTCGACAAACT",
        93,
        460,
        "460M",
        Some("GTTGACAGAATCTACGAGATTCCGCACGCATGGCGTTTGATCAAAACGTTGCCATTTGTAGTATTCGAGGAATTTGGCGCCCGGATCAACCTTACGGTACTCGACCTGGCTGCCAAGTGGTTTGCCACGCAAGACTCACTATCTCGGATCAACCAAAATCCTGCCTTGGCTTTCTACTACTCAAGAAATGACTCGCTGGACGTCGATTACGAACAGGTGCGTCGGCTGAACAAGCCTAAAGAAATGGACAGGGAGTGGATATTTTCTTTAATAGACGAGATCGCCGGCGAGGGAGAAGAGCAGAAGACATGCTCAAGCTCGTCGCAATCTCTGCGCCTGAAGACGTTCGTGAGAGTCTCGACAAACT")
    );
}

fn test_read_start_on_reference_haplotype(
    cigar: &str,
    read_start_on_haplotype: u32,
    expected_offset_in_ref: u32,
) {
    let haplotype_vs_ref_cigar = CigarString::try_from(cigar).unwrap();
    let offset_in_ref = AlignmentUtils::read_start_on_reference_haplotype(
        &haplotype_vs_ref_cigar,
        read_start_on_haplotype,
    );
    assert_eq!(offset_in_ref, expected_offset_in_ref);
}

#[test]
fn make_read_start_on_reference_haplotype_data() {
    test_read_start_on_reference_haplotype("30M5D20M", 50, 55);
    test_read_start_on_reference_haplotype("30M5I20M", 50, 45);
    test_read_start_on_reference_haplotype("55M", 50, 50);
    test_read_start_on_reference_haplotype("30M5D30M5D30M", 80, 90);
    test_read_start_on_reference_haplotype("30M5D30M5I30M", 80, 80);
    test_read_start_on_reference_haplotype("30M5D30M5I30M", 80, 80);
}
