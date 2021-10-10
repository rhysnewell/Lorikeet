use annotator::tandem_repeat::TandemRepeat;
use assembly::assembly_region::AssemblyRegion;
use model::variant_context::VariantContext;
use rayon::prelude::*;
use reference::reference_reader::ReferenceReader;
use std::cmp::{max, min};
use std::collections::BTreeSet;
use utils::simple_interval::{Locatable, SimpleInterval};

#[derive(Debug, Clone)]
pub struct AssemblyRegionTrimmer {
    assembly_region_padding: usize,
    indel_padding_for_genotyping: usize,
    snp_padding_for_genotyping: usize,
    str_padding_for_genotyping: usize,
    max_extension_into_region_padding: usize,
}

/**
 * Helper component to manage active region trimming
 *
 * <p/>
 * It receives the user arguments that controls trimming and also performs the trimming region calculation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
impl AssemblyRegionTrimmer {
    pub fn new(
        assembly_region_padding: usize,
        indel_padding_for_genotyping: usize,
        snp_padding_for_genotyping: usize,
        str_padding_for_genotyping: usize,
        max_extension_into_region_padding: usize,
    ) -> AssemblyRegionTrimmer {
        Self {
            assembly_region_padding,
            indel_padding_for_genotyping,
            snp_padding_for_genotyping,
            str_padding_for_genotyping,
            max_extension_into_region_padding,
        }
    }

    fn no_variation(target_region: AssemblyRegion) -> AssemblyRegionTrimmerResult {
        return AssemblyRegionTrimmerResult::new(target_region, None, None);
    }

    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered later.
     *
     * @param region the genome location range to trim.
     * @param variants list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code region} are simply ignored.
     * @param referenceContext
     * @return never {@code null}.
     */
    pub fn trim(
        &self,
        ref_idx: usize,
        region: AssemblyRegion,
        variants: BTreeSet<VariantContext>,
        legacy_trimming: bool,
        reference_reader: &ReferenceReader,
        reference_bases: &Vec<u8>,
    ) -> AssemblyRegionTrimmerResult {
        if legacy_trimming {
            return self.trim_legacy(region, variants, reference_reader);
        };

        let variants_in_region = variants
            .into_iter()
            .filter(|variant| region.get_span().overlaps(&variant.loc))
            .collect::<Vec<VariantContext>>();

        if variants_in_region.is_empty() {
            return Self::no_variation(region);
        };

        let mut min_start = variants_in_region
            .iter()
            .map(|vc| vc.loc.get_start())
            .min()
            .unwrap_or(0);
        let mut max_end = variants_in_region
            .iter()
            .map(|vc| vc.loc.get_end())
            .max()
            .unwrap_or(0);
        let variant_span = SimpleInterval::new(region.get_contig(), min_start, max_end)
            .intersect(&region.active_span);

        for mut vc in variants_in_region {
            let mut padding = self.snp_padding_for_genotyping;
            if vc.is_indel() {
                padding = self.indel_padding_for_genotyping;
                let num_repeats_and_unit =
                    TandemRepeat::get_num_tandem_repeat_units(&region, &mut vc, reference_bases);
                match num_repeats_and_unit {
                    Some(num_repeats_and_unit) => {
                        let repeat_length = num_repeats_and_unit.1.len();
                        let most_repeats = num_repeats_and_unit.0.iter().max().unwrap_or(&0);
                        let longest_str = most_repeats * repeat_length;
                        padding = self.str_padding_for_genotyping + longest_str;
                    }
                    None => {
                        // pass
                    }
                }
            }

            min_start = min(
                min_start,
                vc.loc.get_start().checked_sub(padding).unwrap_or(0),
            );
            max_end = max(max_end, vc.loc.get_end() + padding);
        }

        let padded_variant_span = SimpleInterval::new(region.get_contig(), min_start, max_end)
            .intersect(&region.get_padded_span());
        debug!(
            "Padded and trimmed the region to this span: {:?}",
            &padded_variant_span
        );

        return AssemblyRegionTrimmerResult::new(
            region,
            Some(variant_span),
            Some(padded_variant_span),
        );
    }

    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered latter.
     *
     * @param originalRegion the genome location range to trim.
     * @param allVariantsWithinExtendedRegion list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code originalRegion} are simply ignored.
     * @return never {@code null}.
     */
    pub fn trim_legacy(
        &self,
        original_region: AssemblyRegion,
        all_variants_within_extended_region: BTreeSet<VariantContext>,
        reference_reader: &ReferenceReader,
    ) -> AssemblyRegionTrimmerResult {
        if all_variants_within_extended_region.is_empty() {
            return Self::no_variation(original_region);
        };

        let mut within_active_region = Vec::new();
        // let original_region_range = original_region.get_span();
        let mut found_non_snp = false;
        let mut variant_span = None;

        for mut vc in all_variants_within_extended_region {
            if original_region.get_span().overlaps(&vc.loc) {
                found_non_snp = found_non_snp || !vc.is_snp();
                variant_span = if variant_span.is_none() {
                    Some(vc.loc.clone())
                } else {
                    Some(variant_span.as_mut().unwrap().span_with(&vc.loc))
                };
                within_active_region.push(vc);
            }
        }

        let padding = if found_non_snp {
            self.indel_padding_for_genotyping
        } else {
            self.snp_padding_for_genotyping
        };

        // we don't actually have anything in the region after skipping out variants that don't overlap
        // the region's full location
        if variant_span.is_none() {
            return Self::no_variation(original_region);
        };

        let target_length =
            reference_reader.get_contig_length(original_region.get_contig()) as usize;
        let maximum_span = original_region
            .get_span()
            .expand_within_contig(self.max_extension_into_region_padding, target_length);
        let ideal_span = variant_span
            .as_ref()
            .unwrap()
            .expand_within_contig(self.max_extension_into_region_padding, target_length);
        let final_span = maximum_span
            .intersect(&ideal_span)
            .merge_with_contiguous(variant_span.as_ref().unwrap())
            .unwrap();

        //      TODO disable this code with ERC
        //        // Make double sure that, if we are emitting GVCF we won't call non-variable positions beyond the target active region span.
        //        // In regular call we don't do so so we don't care and we want to maintain behavior, so the conditional.
        //        final SimpleInterval callableSpan = emitReferenceConfidence ? variantSpan.intersect(originalRegionRange) : variantSpan;

        return AssemblyRegionTrimmerResult::new(original_region, variant_span, Some(final_span));
    }
}

#[derive(Debug)]
pub struct AssemblyRegionTrimmerResult {
    // Holds the input active region.
    pub(crate) original_region: AssemblyRegion,
    // Holds the smaller range that contain all relevant callable variants in the
    // input active region (not considering the extension).
    pub(crate) variant_span: Option<SimpleInterval>,
    // The trimmed variant region span including the extension.
    pub(crate) padded_span: Option<SimpleInterval>,
}

impl AssemblyRegionTrimmerResult {
    /**
     * Creates a trimming result given all its properties.
     * @param originalRegion the original active region.
     * @param variantSpan variant containing span without padding.
     * @param paddedSpan final trimmed variant span including the extension.
     */
    pub fn new(
        original_region: AssemblyRegion,
        variant_span: Option<SimpleInterval>,
        padded_span: Option<SimpleInterval>,
    ) -> AssemblyRegionTrimmerResult {
        if variant_span.is_some() && padded_span.is_some() {
            assert!(
                padded_span
                    .as_ref()
                    .unwrap()
                    .contains(variant_span.as_ref().unwrap()),
                "the padded span must include the variant span"
            )
        };

        AssemblyRegionTrimmerResult {
            original_region,
            variant_span,
            padded_span,
        }
    }

    /**
     * Checks whether there is any variation present in the target region.
     *
     * @return {@code true} if there is any variant, {@code false} otherwise.
     */
    pub fn is_variation_present(&self) -> bool {
        self.variant_span.is_some()
    }

    /**
     * Returns the trimmed variant containing region
     *
     * @throws IllegalStateException if there is no variation detected.
     *
     * @return never {@code null}.
     */
    pub fn get_variant_region(self) -> AssemblyRegion {
        assert!(
            self.variant_span.is_some(),
            "There is no variation present {:?}",
            &self
        );
        return self
            .original_region
            .trim_with_padded_span(self.variant_span.unwrap(), self.padded_span.unwrap());
    }

    /**
     *  Returns the trimmed out left non-variant region.
     *  <p/>
     *  Notice that in case of no variation, the whole original region is considered the left flanking region.
     *
     *  @throws IllegalStateException if there is not such as left flanking region.
     */
    pub fn non_variant_left_flank_region(
        self,
        assembly_region_padding: usize,
    ) -> Option<AssemblyRegion> {
        if !self.is_variation_present() {
            return Some(self.original_region);
        } else if self.original_region.get_start() < self.variant_span.as_ref().unwrap().get_start()
        {
            let left_flank = SimpleInterval::new(
                self.original_region.get_contig(),
                self.original_region.get_start(),
                self.variant_span.as_ref().unwrap().get_start() - 1,
            );
            return Some(
                self.original_region
                    .trim_with_padding(left_flank, assembly_region_padding),
            );
        } else {
            None
        }
    }

    /**
     *  Returns the trimmed out right non-variant region.
     */
    pub fn non_variant_right_flank_region(
        self,
        assembly_region_padding: usize,
    ) -> Option<AssemblyRegion> {
        if !self.is_variation_present() {
            return Some(self.original_region);
        } else if self.variant_span.as_ref().unwrap().get_end() < self.original_region.get_end() {
            let right_flank = SimpleInterval::new(
                self.original_region.get_contig(),
                self.variant_span.as_ref().unwrap().get_end() + 1,
                self.original_region.get_end(),
            );
            return Some(
                self.original_region
                    .trim_with_padding(right_flank, assembly_region_padding),
            );
        } else {
            return None;
        }
    }
}
