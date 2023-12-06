use hashlink::{LinkedHashMap, LinkedHashSet};
use rand::distributions::{Distribution, Normal};
use rand::rngs::ThreadRng;
use std::cmp::Ordering;

use crate::genotype::genotype_builder::{AttributeObject, Genotype, GenotypesContext};
use crate::model::allele_likelihoods::AlleleLikelihoods;
use crate::model::byte_array_allele::Allele;
use crate::model::variant_context::VariantContext;
use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::read_utils::ReadUtils;
use crate::utils::math_utils::MathUtils;

/// Determine whether the annotation appears in the info or format field of the VCF
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum AnnotationType {
    Info,
    Format,
}

impl Ord for AnnotationType {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other {
            Ordering::Equal
        } else if self == &Self::Format {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

impl PartialOrd for AnnotationType {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Possible Variant Annotation types, each enum branch takes the tag for each field
/// i.e. VariantAnnotation::Depth("DP")
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum VariantAnnotations {
    Depth,
    AlleleFraction,
    AlleleCount,
    MappingQuality,
    BaseQuality,
    DepthPerAlleleBySample,
    QualByDepth,
    MLEAC,
    MLEAF,
    PhredLikelihoods,
    GenotypeQuality,
    Genotype,
    VariantGroup,
    Strain,
    Qualified,
}

/// The actual annotation struct, Holds all information about an annotation
#[derive(Debug, Clone)]
pub struct Annotation {
    annotation_type: AnnotationType,
    annotation: VariantAnnotations,
    value: AttributeObject,
}

impl VariantAnnotations {
    const MAX_QD_BEFORE_FIXING: f64 = 45.0; // This is 35.0 in GATK
    const IDEAL_HIGH_QD: f64 = 45.0; // this is 35.0 in GATK
    const JITTER_SIGMA: f64 = 3.0;

    pub fn to_key(&self) -> &str {
        match self {
            Self::Depth => "DP",
            Self::AlleleFraction => "AF",
            Self::AlleleCount => "AC",
            Self::MappingQuality => "MQ",
            Self::BaseQuality => "BQ",
            Self::DepthPerAlleleBySample => "AD",
            Self::QualByDepth => "QD",
            Self::MLEAC => "MLEAC",
            Self::MLEAF => "MLEAF",
            Self::PhredLikelihoods => "PL",
            Self::GenotypeQuality => "GQ",
            Self::Genotype => "GT",
            Self::VariantGroup => "VG",
            Self::Strain => "ST",
            Self::Qualified => "QF",
        }
    }

    pub fn annotate<A: Allele>(
        &self,
        vc: &mut VariantContext,
        genotype: Option<&mut Genotype>,
        likelihoods: &mut AlleleLikelihoods<A>,
        annotation_type: AnnotationType,
    ) -> AttributeObject {
        match self {
            Self::Depth => {
                // debug!("Depth");
                match annotation_type {
                    AnnotationType::Format => {
                        let genotype = genotype.unwrap();
                        if !genotype.has_ad() {
                            // if there is no AD value calculate it now using likelihoods
                            Self::DepthPerAlleleBySample.annotate(
                                vc,
                                Some(genotype),
                                likelihoods,
                                annotation_type,
                            );
                        };
                        let total_ad: i32 = genotype.ad.iter().sum();
                        genotype.dp = total_ad;
                        return AttributeObject::None;
                    }
                    AnnotationType::Info => {
                        return AttributeObject::UnsizedInteger(vc.genotypes.get_dp() as usize);
                    }
                }
            }
            Self::AlleleFraction => {
                let genotype = genotype.unwrap();
                // debug!("Allele Fraction");
                if genotype.has_ad() {
                    let allele_fractions = MathUtils::normalize_sum_to_one(
                        genotype
                            .get_ad()
                            .into_iter()
                            .map(|a| *a as f64)
                            .collect::<Vec<f64>>(),
                    );
                    genotype.attribute(
                        self.to_key().to_string(),
                        AttributeObject::Vecf64(allele_fractions),
                    );
                    return AttributeObject::None;
                } else {
                    // if there is no AD value calculate it now using likelihoods
                    Self::DepthPerAlleleBySample.annotate(
                        vc,
                        Some(genotype),
                        likelihoods,
                        annotation_type,
                    );
                    let allele_fractions = MathUtils::normalize_sum_to_one(
                        genotype
                            .get_ad()
                            .into_iter()
                            .map(|a| *a as f64)
                            .collect::<Vec<f64>>(),
                    );
                    genotype.attribute(
                        self.to_key().to_string(),
                        AttributeObject::Vecf64(allele_fractions),
                    );
                    return AttributeObject::None;
                }
            }
            Self::AlleleCount => {
                let genotype = genotype.unwrap();
                // debug!("Allele Count");
                if genotype.has_ad() {
                    let allele_counts = genotype.get_ad().into_iter().filter(|ad| **ad > 0).count();
                    genotype.attribute(
                        self.to_key().to_string(),
                        AttributeObject::UnsizedInteger(allele_counts),
                    );
                    return AttributeObject::None;
                } else {
                    // if there is no AD value calculate it now using likelihoods
                    Self::DepthPerAlleleBySample.annotate(
                        vc,
                        Some(genotype),
                        likelihoods,
                        annotation_type,
                    );
                    let allele_counts = genotype.get_ad().into_iter().filter(|ad| **ad > 0).count();
                    genotype.attribute(
                        self.to_key().to_string(),
                        AttributeObject::UnsizedInteger(allele_counts),
                    );
                    return AttributeObject::None;
                }
            }
            Self::MappingQuality | Self::BaseQuality => {
                let mut values: LinkedHashMap<usize, Vec<u8>> = LinkedHashMap::new();

                likelihoods
                    .best_alleles_breaking_ties_main(Box::new(|allele: &A| {
                        if allele.is_reference() {
                            1
                        } else {
                            0
                        }
                    }))
                    .into_iter()
                    .filter(|ba| {
                        ba.is_informative()
                            && Self::is_usable_read(
                                &likelihoods
                                    .evidence_by_sample_index
                                    .get(&ba.sample_index)
                                    .unwrap()[ba.evidence_index],
                            )
                    })
                    .for_each(|ba| {
                        let value = values.entry(ba.allele_index.unwrap()).or_insert(Vec::new());
                        match self.get_value_u8(
                            &likelihoods
                                .evidence_by_sample_index
                                .get(&ba.sample_index)
                                .unwrap()[ba.evidence_index],
                            vc,
                        ) {
                            None => {} // pass
                            Some(val) => value.push(val),
                        }
                    });

                let statistics = (0..vc.alleles.len())
                    .into_iter()
                    .map(|index| {
                        let mut vals = values.entry(index).or_insert(Vec::new());
                        if vals.len() > 0 {
                            MathUtils::median(&mut vals)
                        } else {
                            30
                        }
                    })
                    .collect::<Vec<u8>>();

                return AttributeObject::VecU8(statistics);
            }
            Self::DepthPerAlleleBySample => {
                let genotype = genotype.unwrap();
                let alleles = vc.alleles.clone().into_iter().collect::<LinkedHashSet<_>>();
                // debug!("Depth per allele alleles {:?}", &alleles);
                // alleles.iter().for_each(|a| {
                //     // type difference mean we can't check if this allele is in the array at this point
                //     assert!(
                //         likelihoods.alleles.contains_allele(a),
                //         "Likelihoods {:?} does not contain {:?}",
                //         likelihoods.alleles,
                //         a
                //     )
                // });
                if likelihoods.number_of_alleles() <= 1 {
                    return AttributeObject::None;
                }
                let mut allele_counts = LinkedHashMap::new();
                // let mut subset = LinkedHashMap::new();
                for (allele_index, _allele) in alleles.iter().enumerate() {
                    allele_counts.insert(allele_index, 0);
                    // subset.insert(allele_index, vec![allele]);
                }

                let sample_index = likelihoods
                    .samples
                    .iter()
                    .position(|s| s == &genotype.sample_name)
                    .unwrap_or(0);
                likelihoods
                    .best_alleles_breaking_ties_for_sample(sample_index)
                    .into_iter()
                    .filter(|ba| ba.is_informative())
                    .for_each(|ba| {
                        let count = allele_counts.entry(ba.allele_index.unwrap()).or_insert(0);
                        *count += 1;
                    });
                let mut counts = vec![0; vc.alleles.len()];
                counts[0] = *allele_counts.get(&vc.get_reference_and_index().0).unwrap();
                // debug!("Allele counts {:?}", &allele_counts);
                for (vec_index, (allele_index, _)) in vc
                    .get_alternate_alleles_with_index()
                    .into_iter()
                    .enumerate()
                {
                    // debug!(
                    //     "{} {} \n {:?}",
                    //     vec_index,
                    //     allele_index,
                    //     allele_counts.get(&allele_index)
                    // );
                    counts[vec_index + 1] = *allele_counts.get(&allele_index).unwrap();
                }

                // debug!("{:?}", &counts);
                genotype.ad = counts;

                return AttributeObject::None;
            }
            Self::QualByDepth => {
                // debug!("Qual by depth");
                // debug!(
                //     "vc log10_p_error {} {}",
                //     vc.log10_p_error,
                //     vc.has_log10_p_error()
                // );
                if !vc.has_log10_p_error() {
                    return AttributeObject::None;
                }

                let genotypes = vc.get_genotypes_mut();
                // debug!("genotypes empty {}", genotypes.is_empty());
                if genotypes.is_empty() {
                    return AttributeObject::None;
                }

                let depth = Self::get_depth(genotypes, &likelihoods);
                if depth == 0 {
                    return AttributeObject::None;
                }

                let qual = -10.0 * vc.log10_p_error;
                let mut QD = qual / (depth as f64);
                // debug!(
                //     "Log10 p error {} depth {} QD {}",
                //     vc.log10_p_error, depth, QD
                // );

                QD = Self::fix_too_high_qd(QD);
                // debug!("Updated QD {}", QD);

                return AttributeObject::f64(QD);
            }
            Self::MLEAF
            | Self::MLEAC
            | Self::PhredLikelihoods
            | Self::Genotype
            | Self::GenotypeQuality
            | Self::Strain
            | Self::VariantGroup
            | Self::Qualified => {
                // These are returned in genotype contexts already
                // Or calculated elsewhere i.e. Strain & Qualified
                AttributeObject::None
            }
        }
    }

    fn get_value_u8(&self, read: &BirdToolRead, vc: &VariantContext) -> Option<u8> {
        let return_val = match self {
            Self::MappingQuality => Some(read.read.mapq()),
            Self::BaseQuality => {
                ReadUtils::get_read_base_quality_at_reference_coordinate(read, vc.loc.start)
            }
            _ => panic!("u8 read value not appropriate for {:?}", &self),
        };

        return return_val;
    }

    fn is_usable_read(read: &BirdToolRead) -> bool {
        read.read.mapq() != 0
    }

    pub fn get_depth<A: Allele>(
        genotypes: &mut GenotypesContext,
        likelihoods: &AlleleLikelihoods<A>,
    ) -> i32 {
        let mut depth = 0;
        let mut AD_restrict_depth = 0;

        for genotype in genotypes.genotypes_mut() {
            // we care only about variant calls with likelihoods
            if !genotype.is_het() && !genotype.is_hom_var() && !genotype.is_hom_ref() {
                // debug!(
                //     "Skipping: {} {} {:?}",
                //     genotype.is_het(),
                //     genotype.is_hom_var(),
                //     genotype.genotype_type
                // );
                continue;
            }

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if genotype.has_ad() {
                let total_ad: i32 = genotype.ad.iter().sum();
                genotype.dp = total_ad;
                if total_ad != 0 {
                    if total_ad - genotype.ad[0] > 0 {
                        AD_restrict_depth += total_ad;
                    }
                    depth += total_ad;
                    continue;
                }
            }

            // if there is no AD value or it is a dummy value, we want to look to other means to get the depth
            if let Some(sample_index) = likelihoods.index_of_sample(genotype.sample_name) {
                depth += likelihoods.sample_evidence_count(sample_index) as i32;
            } else if genotype.has_dp() {
                depth += genotype.dp;
            }
        }

        if AD_restrict_depth > 0 {
            depth = AD_restrict_depth;
        }

        return depth;
    }

    /**
     * The haplotype caller generates very high quality scores when multiple events are on the
     * same haplotype.  This causes some very good variants to have unusually high QD values,
     * and VQSR will filter these out.  This code looks at the QD value, and if it is above
     * threshold we map it down to the mean high QD value, with some jittering
     *
     * @param QD the raw QD score
     * @return a QD value
     */
    pub fn fix_too_high_qd(qd: f64) -> f64 {
        if qd < Self::MAX_QD_BEFORE_FIXING {
            return qd;
        } else {
            let mut rng = ThreadRng::default();
            let normal = Normal::new(0.0, 1.0);
            return Self::IDEAL_HIGH_QD + normal.sample(&mut rng) * Self::JITTER_SIGMA;
        }
    }

    /// Generates a string holding the information for the VariantAnnotation to be inserted
    /// into the VCF header
    pub fn header_record(&self, annotation_type: &AnnotationType) -> String {
        match self {
            VariantAnnotations::Depth => match annotation_type {
                AnnotationType::Info => {
                    format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">", self.to_key())
                }
                AnnotationType::Format => {
                    format!("##FORMAT=<ID={},Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">", self.to_key())
                }
            },
            VariantAnnotations::BaseQuality => {
                format!("##INFO=<ID={},Number=R,Type=Integer,Description=\"Median PHRED-scaled Base Quality of the variant\">", self.to_key())
            }
            VariantAnnotations::MappingQuality => {
                format!(
                    "##INFO=<ID={},Number=R,Type=Float,Description=\"RMS Mapping Quality\">",
                    self.to_key()
                )
            }
            VariantAnnotations::AlleleFraction => {
                format!("##INFO=<ID={},Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">", self.to_key())
            }
            VariantAnnotations::AlleleCount => {
                format!("##INFO=<ID={},Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">", self.to_key())
            }
            VariantAnnotations::QualByDepth => {
                format!("##INFO=<ID={},Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">", self.to_key())
            }
            VariantAnnotations::Qualified => {
                format!("##INFO=<ID={},Number=1,Type=String,Description=\"Whether the variant passed quality checks to be included in ANI  analyses\">", self.to_key())
            }
            VariantAnnotations::DepthPerAlleleBySample => {
                format!("##FORMAT=<ID={},Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", self.to_key())
            }
            VariantAnnotations::MLEAF => {
                format!("##INFO=<ID={},Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">", self.to_key())
            }
            VariantAnnotations::MLEAC => {
                format!("##INFO=<ID={},Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">", self.to_key())
            }
            VariantAnnotations::GenotypeQuality => {
                format!(
                    "##FORMAT=<ID={},Number=1,Type=Integer,Description=\"Genotype Quality\">",
                    self.to_key()
                )
            }
            VariantAnnotations::Genotype => {
                format!(
                    "##FORMAT=<ID={},Number=1,Type=String,Description=\"Genotype\">",
                    self.to_key()
                )
            }
            VariantAnnotations::PhredLikelihoods => {
                format!("##FORMAT=<ID={},Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">", self.to_key())
            }
            VariantAnnotations::VariantGroup => {
                format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"The base variant group assigned to this variant after clustering\">", self.to_key())
            }
            VariantAnnotations::Strain => {
                format!("##INFO=<ID={},Number=N,Type=Integer,Description=\"A list of potential strain ids associated with this variant location\">", self.to_key())
            }
        }
    }
}

impl Annotation {
    pub fn new(annotation: VariantAnnotations, annotation_type: AnnotationType) -> Self {
        Self {
            annotation_type,
            annotation,
            value: AttributeObject::None,
        }
    }

    pub fn new_with_value(
        annotation: VariantAnnotations,
        annotation_type: AnnotationType,
        value: AttributeObject,
    ) -> Self {
        Self {
            annotation_type,
            annotation,
            value,
        }
    }

    pub fn annotate<A: Allele>(
        mut self,
        vc: &mut VariantContext,
        genotype: Option<&mut Genotype>,
        likelihoods: &mut AlleleLikelihoods<A>,
    ) -> Self {
        self.value = self
            .annotation
            .annotate(vc, genotype, likelihoods, self.annotation_type);
        self
    }

    pub fn get_key(&self) -> &str {
        self.annotation.to_key()
    }

    pub fn get_object(self) -> AttributeObject {
        self.value
    }

    pub fn annotation_type(&self) -> &AnnotationType {
        &self.annotation_type
    }

    pub fn generate_header_record(&self) -> String {
        self.annotation.header_record(&self.annotation_type)
    }
}

impl Ord for Annotation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.annotation_type.cmp(&other.annotation_type)
    }
}

impl PartialOrd for Annotation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Annotation {
    fn eq(&self, other: &Self) -> bool {
        self.annotation_type == other.annotation_type && self.annotation == other.annotation
    }
}

impl Eq for Annotation {}
