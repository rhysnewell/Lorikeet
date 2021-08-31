use haplotype::haplotype::Haplotype;
use model::byte_array_allele::Allele;
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_builder::{CigarBuilder, CigarBuilderResult};
use reads::cigar_utils::CigarUtils;
use reads::cigar_utils::ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS;
use reads::read_clipper::ReadClipper;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use smith_waterman::bindings::SWOverhangStrategy;
use smith_waterman::smith_waterman_aligner::SmithWatermanAligner;
use std::ops::Range;
use utils::simple_interval::SimpleInterval;

lazy_static! {
    pub static ref HAPLOTYPE_TAG: String = format!("HC");
}

pub struct AlignmentUtils {}

impl AlignmentUtils {
    /**
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     * @param isInformative true if the read is differentially informative for one of the haplotypes
     *
     * @param aligner
     * @throws IllegalArgumentException if {@code originalRead} is {@code null} or {@code haplotype} is {@code null} or it
     *   does not have a Cigar or the {@code referenceStart} is invalid (less than 1).
     *
     * @return a GATKRead aligned to reference. Never {@code null}.
     */
    pub fn create_read_aligned_to_ref(
        original_read: &BirdToolRead,
        haplotype: &Haplotype<SimpleInterval>,
        ref_haplotype: &Haplotype<SimpleInterval>,
        reference_start: usize,
        is_informative: bool,
    ) -> BirdToolRead {
        assert!(
            reference_start >= 1,
            "Reference start must be >= 1, but got {}",
            reference_start
        );

        // compute the smith-waterman alignment of read -> haplotype //TODO use more efficient than the read clipper here
        let read_minus_soft_clips = ReadClipper::new(original_read).hard_clip_soft_clipped_bases();
        let soft_clipped_bases = original_read.len() - read_minus_soft_clips.len();
        let read_to_haplotype_sw_alignment = SmithWatermanAligner::align(
            haplotype.get_bases(),
            read_minus_soft_clips.read.seq().encoded,
            *ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
            SWOverhangStrategy::SoftClip,
        );

        if read_to_haplotype_sw_alignment.alignment_offset == -1 {
            // sw can fail (reasons not clear) so if it happens just don't realign the read
            return original_read.clone();
        }

        let mut sw_cigar_builder = CigarBuilder::new(true);
        sw_cigar_builder
            .add_all(read_to_haplotype_sw_alignment.cigar.0)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder"));
        let sw_cigar = sw_cigar_builder
            .make(false)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder"));

        // since we're modifying the read we need to clone it
        let mut copied_read = original_read.clone();

        // TODO: Figure out best way to add HashCode to read
        // if is_informative {
        //     copied_read.read.push_aux(b"HC", Aux::)
        // }

        // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
        let right_padded_haplotype_vs_ref_cigar = haplotype
            .get_consolidated_padded_ciagr(1000)
            .unwrap_or_else(|_| panic!("Unhandled error in haplotype"));

        // this computes the number of reference bases before the read starts, based on the haplotype vs ref cigar
        // This cigar corresponds exactly to the readToRefCigarRaw, below.  One might wonder whether readToRefCigarRaw and
        // readToRefCigarClean ever imply different starts, which could occur if if the former has a leading deletion.  However,
        // according to the logic of applyCigarToCigar, this can only happen if the read has a leading deletion wrt its best haplotype,
        // which our SW aligner won't do, or if the read starts on a haplotype base that is in a deletion wrt to reference, which is nonsensical
        // since a base that exists is not a deletion.  Thus, there is nothing to worry about, in contrast to below where we do check
        // whether left-alignment shifted the start position.
        // TODO: Check that this is non-negative
        let read_start_on_reference_haplotype = Self::read_start_on_reference_haplotype(
            &right_padded_haplotype_vs_ref_cigar,
            read_to_haplotype_sw_alignment.alignment_offset as u32,
        );

        let read_start_on_reference = reference_start
            + haplotype.alignment_start_hap_wrt_ref
            + read_start_on_reference_haplotype as usize;

        // compute the read -> ref alignment by mapping read -> hap -> ref from the
        // SW of read -> hap mapped through the given by hap -> ref
        // this is the sub-cigar of the haplotype-to-ref alignment,
        // with cigar elements before the read start removed.
        // Elements after the read end are kept.
        let read_length = CigarUtils::get_read_length(&right_padded_haplotype_vs_ref_cigar);
        let haplotype_to_ref = Self::trim_cigar_by_bases(
            right_padded_haplotype_vs_ref_cigar,
            read_to_haplotype_sw_alignment.alignment_offset as u32,
            read_length - 1,
        )
        .cigar;

        let read_to_ref_cigar = Self::apply_cigar_to_cigar(&sw_cigar, &haplotype_to_ref);
        let left_aligned_read_to_ref_cigar_result = Self::left_align_indels(
            read_to_ref_cigar,
            ref_haplotype.get_bases(),
            read_minus_soft_clips.read.seq().encoded,
            read_start_on_reference_haplotype,
        );
        let left_aligned_read_to_ref_cigar = &left_aligned_read_to_ref_cigar_result.cigar;
        // it's possible that left-alignment shifted a deletion to the beginning of a read and
        // removed it, shifting the first aligned base to the right
        copied_read.read.set_pos(
            (read_start_on_reference
                + left_aligned_read_to_ref_cigar_result.leading_deletion_bases_removed as usize)
                as i64,
        );

        // the SW Cigar does not contain the hard clips of the original read
        // Here we reconcile the aligned read (that has had any softclips removed)
        // with its softclipped bases
        let original_cigar = original_read.read.cigar();
        let new_cigar = Self::append_clipped_elements_from_cigar_to_cigar(
            left_aligned_read_to_ref_cigar.clone(),
            original_cigar,
        );
        copied_read.read.set(
            original_read.read.qname(),
            Some(&new_cigar),
            original_read.read.seq().encoded,
            original_read.read.qual(),
        );

        if CigarUtils::get_read_length(left_aligned_read_to_ref_cigar) + soft_clipped_bases as u32
            != copied_read.len() as u32
        {
            panic!(
                "Cigar {:?} with read length {} != read length {} for hap {:?}",
                left_aligned_read_to_ref_cigar,
                CigarUtils::get_read_length(left_aligned_read_to_ref_cigar),
                copied_read.len(),
                ref_haplotype
            );
        };

        return copied_read;
    }

    /**
     * Helper method that handles the work of re-appending clipped bases from the original cigar to the new one
     *
     * @param cigarToHaveClippedElementsAdded cigar to attach softclips to
     * @param originalClippedCigar cigar to check for clipped bases
     * @return a new cigar that has had the clipped elements from the original appended to either end
     */
    fn append_clipped_elements_from_cigar_to_cigar(
        cigar_to_have_clipped_elements_added: CigarString,
        original_clipped_cigar: CigarStringView,
    ) -> CigarString {
        let mut first_index = 0;
        let mut last_index = original_clipped_cigar.len() - 1;
        let mut first_element = original_clipped_cigar.0[0].clone();
        let mut last_element = original_clipped_cigar.0[last_index].clone();
        let mut read_to_ref_cigar_elements_with_hard_clips = Vec::new();

        while CigarUtils::is_clipping(&first_element) && first_index != last_index {
            read_to_ref_cigar_elements_with_hard_clips.push(first_element);
            first_index += 1;
            first_element = original_clipped_cigar.0[first_index].clone();
        }

        read_to_ref_cigar_elements_with_hard_clips
            .par_extend(cigar_to_have_clipped_elements_added.0.into_par_iter());

        let mut end_cigar_elements_to_reverse = Vec::new();
        while CigarUtils::is_clipping(&last_element) && first_index != last_index {
            end_cigar_elements_to_reverse.push(last_element);
            last_index -= 1;
            last_element = original_clipped_cigar.0[last_index].clone();
        }

        // Reverse the order to preserve the original soft/hardclip ordering in
        // mixed clipping cases where softclips precede hardclips
        end_cigar_elements_to_reverse.reverse();
        read_to_ref_cigar_elements_with_hard_clips.par_extend(end_cigar_elements_to_reverse);

        return CigarString::from(read_to_ref_cigar_elements_with_hard_clips);
    }

    /**
     * Generate a new Cigar that maps the operations of the first cigar through those in a second
     *
     * For example, if first is 5M and the second is 2M1I2M then the result is 2M1I2M.
     * However, if first is 1M2D3M and second is 2M1I3M this results in a cigar X
     *
     * ref   : AC-GTA
     * hap   : ACxGTA  - 2M1I3M
     * read  : A--GTA  - 1M2D3M
     * result: A--GTA => 1M1D3M
     *
     * ref   : ACxG-TA
     * hap   : AC-G-TA  - 2M1D3M
     * read  : AC-GxTA  - 3M1I2M
     * result: AC-GxTA => 2M1D1M1I2M
     *
     * ref   : ACGTA
     * hap   : ACGTA  - 5M
     * read  : A-GTA  - 1M1I3M
     * result: A-GTA => 1M1I3M
     *
     * ref   : ACGTAC
     * hap   : AC---C  - 2M3D1M
     * read  : AC---C  - 3M
     * result: AG---C => 2M3D
     *
     * The constraint here is that both cigars should imply that the result have the same number of
     * reference bases (i.e.g, cigar.getReferenceLength() are equals).
     *
     * @param firstToSecond the cigar mapping hap1 -> hap2
     * @param secondToThird the cigar mapping hap2 -> hap3
     * @return A cigar mapping hap1 -> hap3
     */
    pub fn apply_cigar_to_cigar(
        first_to_second: &CigarString,
        second_to_third: &CigarString,
    ) -> CigarString {
        let mut new_elements = CigarBuilder::new(true);
        let n_elements_12 = first_to_second.0.len();
        let n_elements_23 = second_to_third.0.len();

        let mut cigar_12_I = 0;
        let mut cigar_23_I = 0;
        let mut elt_12_I = 0;
        let mut elt_23_I = 0;

        while cigar_12_I < n_elements_12 && cigar_23_I < n_elements_23 {
            let elt_12 = &first_to_second.0[cigar_12_I];
            let elt_23 = &second_to_third.0[cigar_23_I];

            let transform = CigarPairTransform::new(elt_12, elt_23);

            elt_12_I += transform.advance_12;
            elt_23_I += transform.advance_23;

            if transform.op_13.is_some() {
                new_elements.add(transform.op_13.unwrap());
            };

            // if have exhausted our current element, advance to the next one
            if elt_12_I == elt_12.len() {
                cigar_12_I += 1;
                elt_12_I = 0;
            };

            if elt_23_I == elt_23.len() {
                cigar_23_I += 1;
                elt_23_I = 0;
            };
        }

        return new_elements
            .make(false)
            .unwrap_or_else(|_| panic!("Unhandled error in cigar builder"));
    }

    fn read_start_on_reference_haplotype(
        haplotype_vs_ref_cigar: &CigarString,
        read_start_on_haplotype: u32,
    ) -> u32 {
        if read_start_on_haplotype == 0 {
            return 0;
        } else {
            // move forward in the haplotype vs ref cigar until we have consumed enough haplotype bases to reach the read start
            // the number of reference bases consumed during this traversal gives us the reference start
            let mut ref_bases_consumed_before_read_start = 0;
            let mut haplotype_bases_consumed = 0;
            for element in haplotype_vs_ref_cigar.0.iter() {
                ref_bases_consumed_before_read_start += Self::length_on_reference(element);
                haplotype_bases_consumed += Self::length_on_read(element);

                if haplotype_bases_consumed >= read_start_on_haplotype {
                    let excess = if CigarUtils::cigar_consumes_reference_bases(element) {
                        haplotype_bases_consumed - read_start_on_haplotype
                    } else {
                        0
                    };

                    return ref_bases_consumed_before_read_start - excess;
                };
            }

            panic!("Cigar doesn't reach the read start");
        }
    }

    /**
     * Trim cigar down to one that starts at start base in the cigar and extends to (inclusive) end base
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
     * @return a new Cigar containing == start - end + 1 reads
     */
    pub fn trim_cigar_by_bases(cigar: CigarString, start: u32, end: u32) -> CigarBuilderResult {
        return Self::trim_cigar(cigar, start, end, false);
    }

    /**
     * Workhorse for trimCigarByBases and trimCigarByReference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
     * @param byReference should start and end be interpreted as position in the reference or the read to trim to/from?
     * @return a non-null cigar
     */
    fn trim_cigar(
        cigar: CigarString,
        start: u32,
        end: u32,
        by_reference: bool,
    ) -> CigarBuilderResult {
        assert!(end >= start, "End position cannot be before start position");

        let mut new_elements = CigarBuilder::new(true);

        // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read (otherwise) coordinates
        let mut element_start; // inclusive
        let mut element_end = 0; // exclusive -- start of next element

        for elt in cigar.iter() {
            element_start = element_end;
            element_end = element_start
                + (if by_reference {
                    Self::length_on_reference(elt)
                } else {
                    Self::length_on_read(elt)
                });

            // we are careful to include zero-length elements at both ends, that is,
            // elements with elementStart == elementEnd == start and elementStart == elementEnd == end + 1
            if element_end < start || (element_end == start && element_start < start) {
                continue;
            } else if element_start > end && element_end > end + 1 {
                break;
            }

            let overlap_length = if element_end == element_start {
                elt.len()
            } else {
                std::cmp::min(end + 1, element_end) - std::cmp::max(start, element_start)
            };
            new_elements.add(CigarUtils::new_cigar_from_operator_and_length(
                elt,
                overlap_length,
            ));
        }

        assert!(
            element_end > end,
            "Cigar elements don't reach end position (inclusive)"
        );

        return new_elements.make_and_record_deletions_removed_result();
    }

    fn length_on_read(element: &Cigar) -> u32 {
        return if CigarUtils::cigar_consumes_read_bases(element) {
            element.len()
        } else {
            0
        };
    }

    fn length_on_reference(element: &Cigar) -> u32 {
        return if CigarUtils::cigar_consumes_reference_bases(element) {
            element.len()
        } else {
            0
        };
    }

    /**
     * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>readStart</code> on the <code>ref</code> and specified by its <code>cigar</code>.
     * <p/>
     * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
     * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
     * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
     * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
     *
     * Soft-clipped bases in the cigar are presumed to correspond to bases in the byte[] of read sequence.  That is, this method
     * assumes that inputs are precise about the distinction between hard clips (removed from the read sequence) and soft clips
     * (kept in the read sequence but not aligned).  For example, with the inputs {cigar: 2S2M2I, read sequence: TTAAAA, ref sequence: GGAA, read start: 2}
     * the method lines up the AAAA (2M2I) of the read with the AA of the ref and left-aligns the indel to yield a cigar of
     * 2S2I2M.
     *
     * @param cigar     structure of the original alignment
     * @param ref    reference sequence the read is aligned to
     * @param read   read sequence
     * @param readStart  0-based position on ref of the first aligned base in the read sequence
     * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    pub fn left_align_indels(
        cigar: CigarString,
        ref_seq: &[u8],
        read: &[u8],
        read_start: u32,
    ) -> CigarBuilderResult {
        if cigar.0.par_iter().all(|elem| !CigarUtils::is_indel(elem)) {
            return CigarBuilderResult::new(cigar, 0, 0);
        }

        // we need reference bases from the start of the read to the rightmost indel
        let last_indel = (0..cigar.0.len())
            .into_par_iter()
            .filter(|n| CigarUtils::is_indel(&cigar.0[*n]))
            .max()
            .unwrap_or(0);
        let necessary_ref_length = read_start
            + cigar
                .0
                .iter()
                .take(last_indel + 1)
                .map(|e| Self::length_on_reference(e))
                .sum::<u32>();
        assert!(
            necessary_ref_length <= ref_seq.len() as u32,
            "Read goes past end of reference"
        );

        // at this point, we are one base past the end of the read.  Now we traverse the cigar from right to left
        let mut result_right_to_left = Vec::new();
        let ref_length = CigarUtils::get_reference_length(&cigar);
        let mut ref_indel_range =
            (read_start + ref_length) as i32..(read_start + ref_length) as i32;
        let mut read_indel_range = read.len() as i32..read.len() as i32;
        for n in (0..cigar.0.len()).into_iter().rev() {
            let element = cigar.0[n];
            // if it's an indel, just accumulate the read and ref bases consumed.  We won't shift the indel until we hit an alignment
            // block or the read start.
            if CigarUtils::is_indel(&element) {
                ref_indel_range.start -= Self::length_on_reference(&element) as i32;
                read_indel_range.start -= Self::length_on_read(&element) as i32;
            } else if ref_indel_range.len() == 0 && read_indel_range.len() == 0 {
                ref_indel_range.start -= Self::length_on_reference(&element) as i32;
                read_indel_range.start -= Self::length_on_read(&element) as i32;
                result_right_to_left.push(element);
            } else {
                let max_shift = if CigarUtils::is_alignment(&element) {
                    element.len()
                } else {
                    0
                };
                let shifts = Self::normalize_alleles(
                    &vec![ref_seq, read],
                    vec![&mut ref_indel_range, &mut read_indel_range],
                    max_shift,
                    true,
                );

                // account for new match alignments on the right due to left-alignment
                result_right_to_left.push(Cigar::Match(shifts.1 as u32));

                // emit if we didn't go all the way to the start of an alignment block
                // OR we have reached clips OR we have reached the start of the cigar
                let emit_indel =
                    n == 0 || shifts.0 < max_shift as i32 || !CigarUtils::is_alignment(&element);
                let new_match_left_due_to_trimming = if shifts.0 < 0 { -shifts.0 } else { 0 };
                let remaining_bases_on_left = if shifts.0 < 0 {
                    element.len() as i32
                } else {
                    element.len() as i32 - shifts.0
                };

                if emit_indel {
                    // some of this alignment block remains after left-alignment -- emit the indel
                    result_right_to_left.push(Cigar::Del(ref_indel_range.len() as u32));
                    result_right_to_left.push(Cigar::Ins(read_indel_range.len() as u32));
                    ref_indel_range.end -= ref_indel_range.len() as i32; // ref indel range is now empty and points to start of left-aligned indel
                    read_indel_range.end -= read_indel_range.len() as i32; // read indel range is now empty and points to start of left-aligned indel

                    ref_indel_range.start -= new_match_left_due_to_trimming
                        + if CigarUtils::cigar_consumes_reference_bases(&element) {
                            remaining_bases_on_left
                        } else {
                            0
                        };
                    ref_indel_range.end -= new_match_left_due_to_trimming
                        + if CigarUtils::cigar_consumes_reference_bases(&element) {
                            remaining_bases_on_left
                        } else {
                            0
                        };

                    read_indel_range.start -= new_match_left_due_to_trimming
                        + if CigarUtils::cigar_consumes_read_bases(&element) {
                            remaining_bases_on_left
                        } else {
                            0
                        };
                    read_indel_range.end -= new_match_left_due_to_trimming
                        + if CigarUtils::cigar_consumes_read_bases(&element) {
                            remaining_bases_on_left
                        } else {
                            0
                        };
                }

                result_right_to_left.push(Cigar::Match(new_match_left_due_to_trimming as u32));
                result_right_to_left.push(CigarUtils::new_cigar_from_operator_and_length(
                    &element,
                    remaining_bases_on_left as u32,
                ));
            }
        }

        result_right_to_left.push(Cigar::Del(ref_indel_range.len() as u32));
        result_right_to_left.push(Cigar::Ins(read_indel_range.len() as u32));
        assert!(
            read_indel_range.start == 0,
            "Given cigar does not account for all bases of the read"
        );
        let mut cigar_builder = CigarBuilder::new(true);
        result_right_to_left.reverse();
        cigar_builder.add_all(result_right_to_left);

        return cigar_builder.make_and_record_deletions_removed_result();
    }

    /**
     *  Example usage:  reference = GAAT, read = GAAAT (insertion of one A) and we initially consider the insertion of the A to occur before
     *  the T.  Thus the reference range of this allele is [3,3) (no bases) and the read range is [3,4).  This will be left-aligned so that
     *  the insertion occurs after the G, so that the ranges become [1,1) and [1,2) and the returned shifts are 2 bases for both the start and end
     *  of the range.
     *
     *  If the given allele ranges are not parsimonious, for example [3,4) and [3,5) in the above example to include the common T in both alleles,
     *  the resulting ranges will be shifted by different amounts.  In this case, the shifts are 2 bases in the front and 3 at the end.
     *
     *  Note that we use the convention that the ref allele in the case of an alt insertion, or the alt allele in case of a deletion, is represented
     *  by [n,n) where n is the last aligned coordinate before the indel.  This makes sense when you think in terms of alignment CIGARs:
     *  suppose for example we have a 5M5I5M read with start 100.  Then the match bases are [100,105) on the ref and [0,5) on the read and the inserted bases are
     *  [5,10) on the read and [5,5) on the reference.
     *
     * @param sequences bases of sequences containing different alleles -- could be reference, a haplotype, a read, or subsequences thereof
     * @param bounds    initial ranges (inclusive start, exclusive end) of alleles in same order as {@code sequences}
     *                  Note that this method adjusts these ranges as a side effect
     * @param maxShift  maximum allowable shift left.  This may be less than the amount demanded by the array bounds.  For example, when
     *                  left-aligning a read with multiple indels, we don't want to realign one indel past another (if they "collide" we merge
     *                  them into a single indel and continue -- see {@link AlignmentUtils::leftAlignIndels}
     * @param trim      If true, remove redundant shared bases at the start and end of alleles
     * @return          The number of bases the alleles were shifted left such that they still represented the same event.
     */
    pub fn normalize_alleles(
        sequences: &Vec<&[u8]>,
        mut bounds: Vec<&mut Range<i32>>,
        max_shift: u32,
        trim: bool,
    ) -> (i32, i32) {
        assert!(!sequences.is_empty());
        assert!(
            sequences.len() == bounds.len(),
            "Must have one initial allele range per sequence"
        );
        bounds.par_iter().for_each(|bound| {
            assert!(
                max_shift <= bound.start as u32,
                "maxShift goes past the start of a sequence"
            )
        });

        let mut start_shift = 0;
        let mut end_shift = 0;

        // consume any redundant shared bases at the end of the alleles
        let mut min_size = bounds.par_iter().map(|bound| bound.len()).min().unwrap() as i32;
        while trim && min_size > 0 && Self::last_base_on_right_is_same(sequences, &bounds) {
            bounds.par_iter_mut().for_each(|bound| {
                bound.end -= 1;
            });
            min_size -= 1;
            end_shift += 1;
        }

        while trim && min_size > 0 && Self::first_base_on_left_is_same(sequences, &bounds) {
            bounds.par_iter_mut().for_each(|bound| {
                bound.start += 1;
            });
            min_size -= 1;
            start_shift -= 1;
        }

        // we shift left as long as the last bases on the right are equal among all sequences and the next bases on the left are all equal.
        // if a sequence is empty (eg the reference relative to an insertion alt allele) the last base on the right is the next base on the left
        while trim
            && min_size > 0
            && Self::next_base_on_left_is_same(sequences, &bounds)
            && Self::last_base_on_right_is_same(sequences, &bounds)
        {
            bounds.par_iter_mut().for_each(|bound| {
                bound.start -= 1;
                bound.end -= 1;
            });
            start_shift += 1;
            end_shift += 1;
        }

        return (start_shift, end_shift);
    }

    // do all sequences share a common base at the end of the given index range
    fn last_base_on_right_is_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut Range<i32>>) -> bool {
        let last_base_on_right = sequences[0][(bounds[0].end - 1) as usize];
        return (0..sequences.len()).into_par_iter().all(|n| {
            if sequences[n][(bounds[n].end - 1) as usize] != last_base_on_right {
                false
            } else {
                true
            }
        });
    }

    // do all sequences share a common first base
    fn first_base_on_left_is_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut Range<i32>>) -> bool {
        let first_base_on_left = sequences[0][bounds[0].start as usize];
        return (0..sequences.len()).into_par_iter().all(|n| {
            if sequences[n][bounds[n].start as usize] != first_base_on_left {
                false
            } else {
                true
            }
        });
    }

    // do all sequences share a common base just before the given index ranges
    fn next_base_on_left_is_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut Range<i32>>) -> bool {
        let next_base_on_left = sequences[0][(bounds[0].start - 1) as usize];
        return (0..sequences.len()).into_par_iter().all(|n| {
            if sequences[n][(bounds[n].start - 1) as usize] != next_base_on_left {
                false
            } else {
                true
            }
        });
    }

    /**
     * Returns the number of bases in a read minus the number of softclipped bases.
     */
    pub fn unclipped_read_length(read: &BirdToolRead) -> usize {
        let mut soft_clipped_bases = 0;
        for element in read.read.cigar().iter() {
            match element {
                Cigar::SoftClip(len) => {
                    soft_clipped_bases += len;
                }
                _ => {
                    // pass
                }
            };
        }

        return read.len() - soft_clipped_bases as usize;
    }

    /**
     * Removing a trailing deletion from the incoming cigar if present
     *
     * @param c the cigar we want to update
     * @return a non-null Cigar
     */
    pub fn remove_trailing_deletions(c: CigarString) -> CigarString {
        match c.0[c.len() - 1] {
            Cigar::Del(_) => CigarString::from(c.0[0..(c.len() - 1)].to_vec()),
            _ => c,
        }
    }

    /**
     * Find the last occurrence of the query sequence in the reference sequence
     *
     * Returns the index of the last occurrence or -1 if the query sequence is not found
     *
     * @param reference the reference sequence
     * @param query the query sequence
     */
    pub fn last_index_of(reference: &[u8], query: &[u8]) -> Option<usize> {
        let query_length = query.len();
        // start search from the last possible matching position and search to the left
        for r in (0..(reference.len().checked_sub(query_length).unwrap_or(0) + 1))
            .into_iter()
            .rev()
        {
            let mut q = 0;
            while q < query_length && reference[r + q] == query[q] {
                q += 1;
            }

            if q == query_length {
                return Some(r);
            }
        }

        return None;
    }

    /**
     * Get the byte[] from bases that cover the reference interval refStart -> refEnd given the
     * alignment of bases to the reference (basesToRefCigar) and the start offset of the bases on the reference
     *
     * refStart and refEnd are 0 based offsets that we want to obtain.  In the client code, if the reference
     * bases start at position X and you want Y -> Z, refStart should be Y - X and refEnd should be Z - X.
     *
     * If refStart or refEnd would start or end the new bases within a deletion, this function will return null
     *
     * If the cigar starts with an insertion, the inserted bases are considered as coming before the start position and
     * are therefore excluded from the result.  That is getBasesCoveringRefInterval(0, 3, "ACTTGT", 0, 2I4M) should yield "TTGT".
     *
     * @param bases
     * @param refStart
     * @param refEnd
     * @param basesStartOnRef where does the bases array start w.r.t. the reference start?  For example, bases[0] of
     *                        could be at refStart == 0 if basesStartOnRef == 0, but it could just as easily be at
     *                        10 (meaning bases doesn't fully span the reference), which would be indicated by basesStartOnRef == 10.
     *                        It's not trivial to eliminate this parameter because it's tied up with the cigar
     * @param basesToRefCigar the cigar that maps the bases to the reference genome
     * @return a byte[] containing the bases covering this interval, or null if we would start or end within a deletion
     */
    pub fn get_bases_covering_ref_interval<'a>(
        ref_start: usize,
        ref_end: usize,
        bases: &'a [u8],
        bases_start_on_ref: usize,
        bases_to_ref_cigar: &CigarString,
    ) -> Option<&'a [u8]> {
        assert!(
            ref_end >= ref_start,
            "Bad start or stop: start {} stop {}",
            ref_start,
            ref_end
        );
        if bases.len() != CigarUtils::get_read_length(bases_to_ref_cigar) as usize {
            panic!("Mismatch between length in reference bases and cigar length");
        };

        let mut ref_pos = bases_start_on_ref;
        let mut bases_pos = 0;
        let mut bases_start = None;
        let mut bases_stop = None;
        let mut done = false;

        let mut iii = 0;
        while iii < bases_to_ref_cigar.0.len() && !done {
            let ce = &bases_to_ref_cigar.0[iii];
            match ce {
                &Cigar::Ins(len) => {
                    bases_pos += len;
                }
                &Cigar::Match(len) | &Cigar::Diff(len) | &Cigar::Equal(len) => {
                    for i in 0..len {
                        if ref_pos == ref_start {
                            bases_start = Some(bases_pos);
                        };
                        if ref_pos == ref_end {
                            bases_stop = Some(bases_pos);
                            done = true;
                            break;
                        }
                        ref_pos += 1;
                        bases_pos += 1;
                    }
                }
                &Cigar::Del(len) => {
                    for i in 0..len {
                        if ref_pos == ref_end || ref_pos == ref_start {
                            // if we ever reach a ref position that is either a start or an end, we fail
                            return None;
                        }
                        ref_pos += 1;
                    }
                }
                _ => {
                    panic!("Unsupported operator {:?}", ce)
                }
            }
            iii += 1;
        }

        if bases_start.is_none() || bases_stop.is_none() {
            panic!(
                "Never found start {:?} or stop {:?} given cigar {:?}",
                bases_start, bases_stop, bases_to_ref_cigar
            );
        };
        return Some(&bases[bases_start.unwrap() as usize..bases_stop.unwrap() as usize + 1]);
    }

    /**
     * Trim cigar down to one that starts at start reference on the left and extends to end on the reference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases on the reference?  The first position is 0
     * @param end Where should we stop keeping bases on the reference?  The maximum value is cigar.getReferenceLength()
     * @return a new Cigar with reference length == start - end + 1
     */
    pub fn trim_cigar_by_reference(cigar: CigarString, start: u32, end: u32) -> CigarBuilderResult {
        return Self::trim_cigar(cigar, start, end, true);
    }
}

/**
 * transformations that project one alignment state through another
 *
 * Think about this as a state machine, where we have:
 *
 * bases3 : xxx A zzz
 * bases2 : xxx B zzz
 * bases1 : xxx C zzz
 *
 * where A, B and C are alignment states of a three way alignment.  We want to capture
 * the transition from operation mapping 1 -> 2 and an operation mapping 2 -> 3 and its
 * associated mapping from 1 -> 3 and the advancement of the cigar states of 1->2 and 2->3.
 *
 * Imagine that A, B, and C are all equivalent (so that op12 = M and op23 = M).  This implies
 * a mapping of 1->3 of M, and in this case the next states to consider in the 3 way alignment
 * are the subsequent states in 1 and 2 (so that advance12 and advance23 are both 1).
 *
 * Obviously not all of the states and their associated transitions are so simple.  Suppose instead
 * that op12 = I, and op23 = M.  What does this look like:
 *
 * bases3 : xxx - A zzz
 * bases2 : xxx - B zzz
 * bases1 : xxx I C zzz
 *
 * It means that op13 must be an insertion (as we have an extra base in 1 thats not present in 2 and
 * so not present in 3).  We advance the cigar in 1 by 1 (as we've consumed one base in 1 for the I)
 * but we haven't yet found the base corresponding to the M of op23.  So we don't advance23.
 */
struct CigarPairTransform {
    op_13: Option<Cigar>,
    advance_12: u32,
    advance_23: u32,
}

impl CigarPairTransform {
    pub fn new(op_12: &Cigar, op_23: &Cigar) -> CigarPairTransform {
        match op_12 {
            //
            // op12 is a match
            //
            &Cigar::Match(_) | &Cigar::Equal(_) | &Cigar::Diff(_) => {
                match op_23 {
                    &Cigar::Match(_) | &Cigar::Equal(_) | &Cigar::Diff(_) => {
                        // 3: xxx B yyy
                        // ^^^^^^^^^^^^
                        // 2: xxx M yyy
                        // 1: xxx M yyy
                        return Self::new_inner(Some(Cigar::Match(1)), 1, 1);
                    }
                    &Cigar::Ins(_) | &Cigar::SoftClip(_) => {
                        // 3: xxx I yyy
                        // ^^^^^^^^^^^^
                        // 2: xxx I yyy
                        // 1: xxx M yyy
                        return Self::new_inner(Some(Cigar::Ins(1)), 1, 1);
                    }
                    &Cigar::Del(_) => {
                        // 3: xxx D yyy
                        // ^^^^^^^^^^^^
                        // 2: xxx D yyy
                        // 1: xxx M yyy
                        return Self::new_inner(Some(Cigar::Del(1)), 0, 1);
                    }
                    _ => panic!("Unexpected stat {:?}", op_12),
                };
            }
            //
            // op12 is a insertion
            //
            &Cigar::Ins(_) | &Cigar::SoftClip(_) => match op_23 {
                &Cigar::Match(_) | &Cigar::Equal(_) | &Cigar::Diff(_) => {
                    return CigarPairTransform::new_inner(Some(Cigar::Ins(1)), 1, 0)
                }
                &Cigar::Ins(_) | &Cigar::SoftClip(_) => {
                    return CigarPairTransform::new_inner(Some(Cigar::Ins(1)), 1, 0)
                }
                &Cigar::Del(_) => return CigarPairTransform::new_inner(Some(Cigar::Ins(1)), 1, 0),
                _ => panic!("Unexpected stat {:?}", op_12),
            },
            //
            // op12 is a deletion
            //
            &Cigar::Del(_) => {
                match op_23 {
                    &Cigar::Match(_) | &Cigar::Equal(_) | &Cigar::Diff(_) => {
                        // 3: xxx D M yyy
                        // ^^^^^^^^^^^^
                        // 2: xxx M yyy
                        // 1: xxx D yyy
                        return Self::new_inner(Some(Cigar::Del(1)), 1, 1);
                    }
                    &Cigar::Ins(_) | &Cigar::SoftClip(_) => {
                        // 3: xxx X yyy => no-op, we skip emitting anything here
                        // ^^^^^^^^^^^^
                        // 2: xxx I yyy
                        // 1: xxx D yyy
                        return Self::new_inner(None, 1, 1);
                    }
                    &Cigar::Del(_) => {
                        // 3: xxx D2 D1 yyy
                        // ^^^^^^^^^^^^
                        // 2: xxx D2 yyy
                        // 1: xxx D1 yyy
                        return Self::new_inner(Some(Cigar::Del(1)), 0, 1);
                    }
                    _ => panic!("Unexpected stat {:?}", op_12),
                };
            }
            _ => panic!("Unexpected stat {:?}", op_12),
        };
    }

    fn new_inner(op_13: Option<Cigar>, advance_12: u32, advance_23: u32) -> CigarPairTransform {
        Self {
            op_13,
            advance_12,
            advance_23,
        }
    }
}
