// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Deref, Range};

use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};

//use crate::grammar;

//pub mod model.evidence;
//pub mod likelihood;
//pub mod modes;
//pub mod sample;
//
//#[derive(Debug, Clone)]
//pub struct Contamination {
//    pub by: usize,
//    pub fraction: f64,
//}
//
//#[derive(Ord, Eq, PartialOrd, PartialEq, Clone, Debug)]
//pub struct Event {
//    pub name: String,
//    pub vafs: grammar::VAFTree,
//    pub strand_bias: StrandBias,
//}
//
//impl Event {
//    pub fn is_artifact(&self) -> bool {
//        self.strand_bias != StrandBias::None
//    }
//}

pub type AlleleFreq = NotNan<f64>;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord)]
pub enum StrandBias {
    None,
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None
    }
}

impl StrandBias {
    pub fn is_some(&self) -> bool {
        if let StrandBias::None = self {
            false
        } else {
            true
        }
    }

    pub fn forward_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_one(),
            StrandBias::Reverse => LogProb::ln_zero(),
        }
    }

    pub fn reverse_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_zero(),
            StrandBias::Reverse => LogProb::ln_one(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumString, EnumIter, IntoStaticStr)]
pub enum VariantType {
    #[strum(serialize = "INS")]
    Insertion(Option<Range<u32>>),
    #[strum(serialize = "DEL")]
    Deletion(Option<Range<u32>>),
    #[strum(serialize = "SNV")]
    SNV,
    #[strum(serialize = "MNV")]
    MNV,
    #[strum(serialize = "REF")]
    None, // site with no suggested alternative allele
}

impl From<&str> for VariantType {
    fn from(string: &str) -> VariantType {
        match string {
            "INS" => VariantType::Insertion(None),
            "DEL" => VariantType::Deletion(None),
            "SNV" => VariantType::SNV,
            "REF" => VariantType::None,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum SVType {
    Dup
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub struct SV {
    sv: SVType,
    len: u32,
    start: u32,
    end: u32,
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    SNV(u8),
    MNV(Vec<u8>),
    SV(SV),
    None,
}

impl Variant {

    pub fn has_fragment_evidence(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SV(_) => true,
            &Variant::SNV(_) => false,
            &Variant::MNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_single_base(&self) -> bool {
        match self {
            &Variant::SNV(_) | &Variant::None => true,
            _ => false,
        }
    }

    pub fn is_snv(&self) -> bool {
        match self {
            &Variant::SNV(_) => true,
            _ => false,
        }
    }

    pub fn is_indel(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SV(_) => false,
            &Variant::SNV(_) => false,
            &Variant::MNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            }
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            (&Variant::MNV(_), &VariantType::MNV) => true,
            (&Variant::None, &VariantType::None) => true,
            _ => false,
        }
    }

    pub fn end(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length,
            &Variant::Insertion(_) => start + 1, // end of insertion is the next regular base
            &Variant::SV(SV) => SV.end,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32,
        }
    }

    pub fn centerpoint(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length / 2,
            &Variant::Insertion(_) => start, // end of insertion is the next regular base
            &Variant::SV(SV) => (SV.start + SV.len) / 2,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32 / 2,
        }
    }

    pub fn len(&self) -> u32 {
        match self {
            &Variant::Deletion(l) => l,
            &Variant::Insertion(ref s) => s.len() as u32,
            &Variant::SV(SV) => SV.len,
            &Variant::SNV(_) => 1,
            &Variant::MNV(ref alt) => alt.len() as u32,
            &Variant::None => 1,
        }
    }
}

#[cfg(test)]
mod tests {
//    use crate::model::model.evidence::{observation::ObservationBuilder, Observation};
//    use crate::utils;

//    use bio::stats::LogProb;
//
//    pub fn observation(prob_mapping: LogProb, prob_alt: LogProb, prob_ref: LogProb) -> Observation {
//        ObservationBuilder::default()
//            .prob_mapping_mismapping(prob_mapping)
//            .prob_alt(prob_alt)
//            .prob_ref(prob_ref)
//            .prob_missed_allele(utils::max_prob(prob_ref, prob_alt))
//            .prob_sample_alt(LogProb::ln_one())
//            .prob_overlap(LogProb::ln_one())
//            .prob_any_strand(LogProb::ln_one())
//            .forward_strand(true)
//            .reverse_strand(true)
//            .build()
//            .unwrap()
//    }

    // fn setup_pairwise_test<'a>(
    // ) -> PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::TumorNormalModel> {
    //     let insert_size = InsertSize {
    //         mean: 250.0,
    //         sd: 50.0,
    //     };
    //     let prior_model = priors::TumorNormalModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, Prob(0.001));
    //     let case_sample = Sample::new(
    //         bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
    //         true,
    //         AlignmentProperties::default(insert_size),
    //         LatentVariableModel::new(1.0),
    //         constants::PROB_ILLUMINA_INS,
    //         constants::PROB_ILLUMINA_DEL,
    //         Prob(0.0),
    //         Prob(0.0),
    //         10,
    //         500,
    //         &[],
    //     );
    //     let control_sample = Sample::new(
    //         bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
    //         true,
    //         AlignmentProperties::default(insert_size),
    //         LatentVariableModel::new(1.0),
    //         constants::PROB_ILLUMINA_INS,
    //         constants::PROB_ILLUMINA_DEL,
    //         Prob(0.0),
    //         Prob(0.0),
    //         10,
    //         500,
    //         &[],
    //     );
    //
    //     let model = PairCaller::new(case_sample, control_sample, prior_model);
    //
    //     model
    // }
    //
    //
    // /// scenario 1: same pileup -> germline call
    // #[test]
    // fn test_same_pileup() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         observations.clone(),
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     // germline
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     assert_relative_eq!(p_germline.exp(), 1.0);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.0);
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert!(*pileup.map_allele_freqs().0 >= 0.97);
    // }
    //
    // /// scenario 2: empty control pileup -> somatic call
    // #[test]
    // fn test_empty_control_pileup() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         vec![],
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     let (af_case, af_control) = pileup.map_allele_freqs();
    //     // we have no model.evidence for germline, but an allele frequency of 1 is most likely with a germline variant!
    //     assert!(p_germline > p_somatic);
    //     // germline < somatic
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert!(*af_case >= 0.97);
    //     assert!(*af_control >= 0.97);
    // }
    //
    // /// scenario 3: subclonal variant
    // #[test]
    // fn test_subclonal() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //     for _ in 0..50 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //             LogProb::ln_one(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         vec![],
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.9985, epsilon = 0.01);
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert_relative_eq!(*pileup.map_allele_freqs().0, 0.09, epsilon = 0.03);
    // }
    //
    // /// scenario 4: absent variant
    // #[test]
    // fn test_absent() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let tumor_ref = ContinuousAlleleFreqs::inclusive(0.0..0.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..10 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //             LogProb::ln_one(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         observations.clone(),
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);
    //
    //     // germline
    //     assert_relative_eq!(p_germline.exp(), 0.0, epsilon = 0.02);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.0, epsilon = 0.02);
    //     // absent
    //     assert_relative_eq!(p_absent.exp(), 1.0, epsilon = 0.02);
    // }
    //
    // #[test]
    // fn test_allele_freq_observable() {
    //     let afs = ContinuousAlleleFreqs::left_exclusive(0.0..0.5);
    //     assert_relative_eq!(*afs.observable_min(15), 0.0666666666666667);
    //     assert_relative_eq!(*afs.observable_max(15), 0.4666666666666667);
    //
    //     let afs = ContinuousAlleleFreqs::right_exclusive(0.0..0.5);
    //     assert_relative_eq!(*afs.observable_min(12), 0.0);
    //     assert_relative_eq!(*afs.observable_max(12), 0.4166666666666667);
    // }
}
