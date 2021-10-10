use annotator::variant_annotation::{Annotation, AnnotationType, VariantAnnotations};
use genotype::genotype_builder::{AttributeObject, Genotype, GenotypesContext};
use haplotype::haplotype::Haplotype;
use hashlink::LinkedHashMap;
use model::allele_likelihoods::AlleleLikelihoods;
use model::byte_array_allele::Allele;
use model::variant_context::VariantContext;
use rust_htslib::bcf::Header;
use std::collections::HashMap;
use utils::simple_interval::SimpleInterval;

/**
 * The class responsible for computing annotations for variants.
 * Annotations are auto-discovered - ie, any class that extends {@link VariantAnnotation} and
 * lives in this package is treated as an annotation and the engine will attempt to create instances of it
 * by calling the non-arg constructor (loading will fail if there is no no-arg constructor).
 */
#[derive(Debug)]
pub struct VariantAnnotationEngine {}

impl VariantAnnotationEngine {
    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     * @param vc the variant context to annotate
     * @param features context containing the features that overlap the given variant
     * @param ref the reference context of the variant to annotate or null if there is none
     * @param readLikelihoods readLikelihoods indexed by sample, allele, and read within sample. May be null
     * @param addAnnot function that indicates if the given annotation type should be added to the variant
     *
     */
    pub fn annotate_context(
        vc: &VariantContext,
        read_likelihoods: &mut AlleleLikelihoods<Haplotype<SimpleInterval>>,
        add_annotation: Box<dyn Fn(&Annotation) -> bool>,
    ) -> VariantContext {
        // annotate genotypes, creating another new VC in the process
        let mut builder = VariantContext::build_from_vc(vc);
        // genotype context annotation here
        builder.genotypes = Self::add_genotype_annotations(&mut builder, read_likelihoods);
        debug!(
            "genotypes {:?} empty {}",
            &builder.genotypes,
            builder.genotypes.is_empty()
        );
        let info_annot_map =
            Self::add_info_annotations(&mut builder, read_likelihoods, add_annotation);

        builder.attributes(info_annot_map);

        return builder;
    }

    fn add_genotype_annotations(
        vc: &mut VariantContext,
        likelihoods: &mut AlleleLikelihoods<Haplotype<SimpleInterval>>,
    ) -> GenotypesContext {
        let mut genotypes = GenotypesContext::create(vc.get_n_samples());

        for g_index in 0..vc.genotypes.genotypes().len() {
            let mut gb = vc.genotypes.genotypes()[g_index].clone();
            for genotype_annotation in Self::genotype_annotations() {
                genotype_annotation.annotate(vc, Some(&mut gb), likelihoods);
            }

            genotypes.add(gb);
        }

        return genotypes;
    }

    fn add_info_annotations(
        vc: &mut VariantContext,
        likelihoods: &mut AlleleLikelihoods<Haplotype<SimpleInterval>>,
        add_annotation: Box<dyn Fn(&Annotation) -> bool>,
    ) -> LinkedHashMap<String, AttributeObject> {
        let mut info_annot_map = LinkedHashMap::new();
        for mut annotation in Self::vc_annotations() {
            if add_annotation(&annotation) {
                let annotation_result = annotation.annotate(vc, None, likelihoods);
                info_annot_map.insert(
                    annotation_result.get_key().to_string(),
                    annotation_result.get_object(),
                );
            }
        }

        return info_annot_map;
    }

    /// Annotations added to the VariantContext
    pub fn vc_annotations() -> Vec<Annotation> {
        vec![
            Annotation::new(VariantAnnotations::Depth, AnnotationType::Info),
            Annotation::new(VariantAnnotations::QualByDepth, AnnotationType::Info),
            Annotation::new(VariantAnnotations::MappingQuality, AnnotationType::Info),
            Annotation::new(VariantAnnotations::BaseQuality, AnnotationType::Info),
        ]
    }

    /// Annotations added to the Genotype of VariantContexts
    pub fn genotype_annotations() -> Vec<Annotation> {
        vec![
            Annotation::new(
                VariantAnnotations::DepthPerAlleleBySample,
                AnnotationType::Format,
            ),
            Annotation::new(VariantAnnotations::AlleleFraction, AnnotationType::Info),
            Annotation::new(VariantAnnotations::AlleleCount, AnnotationType::Info),
        ]
    }

    /// Annotations that are precalculated or calculated through other annotations
    pub fn precalculated_annotations() -> Vec<Annotation> {
        vec![
            Annotation::new(VariantAnnotations::Depth, AnnotationType::Format),
            Annotation::new(VariantAnnotations::Genotype, AnnotationType::Format),
            Annotation::new(VariantAnnotations::GenotypeQuality, AnnotationType::Format),
            Annotation::new(VariantAnnotations::PhredLikelihoods, AnnotationType::Format),
            Annotation::new(VariantAnnotations::MLEAC, AnnotationType::Info),
            Annotation::new(VariantAnnotations::MLEAF, AnnotationType::Info),
        ]
    }

    /// Sorted list of annotations. Format field annotations appear first
    fn all_annotations() -> Vec<Annotation> {
        let mut annotations = Self::precalculated_annotations();
        annotations.extend(Self::genotype_annotations());
        annotations.extend(Self::vc_annotations());
        annotations.sort();
        return annotations;
    }

    /// Populates a given VCF header with all possible annotation fields and info
    pub fn populate_vcf_header(header: &mut Header) {
        for annotation in Self::all_annotations() {
            header.push_record(annotation.generate_header_record().as_bytes());
        }
    }
}
