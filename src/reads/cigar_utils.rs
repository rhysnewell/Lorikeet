use gkl::smithwaterman::{OverhangStrategy, Parameters};
use pair_hmm::pair_hmm_likelihood_calculation_engine::AVXMode;
use rayon::prelude::*;
use reads::alignment_utils::AlignmentUtils;
use reads::cigar_builder::CigarBuilder;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use smith_waterman::smith_waterman_aligner::{SmithWatermanAligner, SmithWatermanAlignmentResult};

lazy_static! {
    static ref SW_PAD: String = format!("NNNNNNNNNN");
    // FROM GATK COMMENTS:
    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    pub static ref NEW_SW_PARAMETERS: Parameters = Parameters::new(200, -150, -260, -11);
    // FROM GATK COMMENTS:
    // In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
    // The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
    // Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
    // would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
    // Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
    // we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
    pub static ref ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS: Parameters = Parameters::new(10, -15, -30, -5);
}

pub struct CigarUtils {}

impl CigarUtils {

    pub fn is_valid(cigar: &CigarString) -> bool {
        if cigar.0.is_empty() {
            return true
        } else {
            let mut seen_real_operator = false;
            for (i, cig) in cigar.0.iter().enumerate() {
                if cig.len() == 0 {
                    // can't have cigar operator of length 0
                    return false
                }

                if Self::is_clipping(cig) {
                    if Self::cigar_is_hard_clip(cig) && (i != 0 && i != cigar.0.len() - 1){
                        return false
                    } else if i != 0 && i != cigar.0.len() - 1 {
                        if i == 1 {
                            if cigar.0.len() != 3
                                || (!Self::cigar_is_hard_clip(&cigar.0[0])
                                && !Self::cigar_is_hard_clip(&cigar.0[2])) {
                                // Soft Clip CIGAR can only be inside of hard clipping operator
                                // if not at start or end of cigar
                                return false
                            }
                        } else if i == cigar.0.len() - 2 {
                            if !Self::cigar_is_hard_clip(&cigar.0[cigar.0.len() - 1]) {
                                return false
                            }
                        } else {
                            return false
                        }
                    }
                } else if Self::is_real_operator(cig) {
                    seen_real_operator = true;
                } else if Self::is_padding_operator(cig) && i != 0 {
                    if i == cigar.0.len() - 1 {
                        // padding operator not valid at end of cigar
                        return false
                    } else if !Self::is_real_operator(&cigar.0[i - 1])
                        || !Self::is_real_operator(&cigar.0[i + 1]) {
                        // padding operator not between real operators
                        return false
                    }
                }
            }

            if !seen_real_operator {
                // no real operators (M|I|D|N|X|=)
                return false
            }
        }
        return true
    }

    pub fn is_padding_operator(cig: &Cigar) -> bool {
        match cig {
            Cigar::Pad(_) => true,
            _ => false,
        }
    }

    pub fn is_real_operator(cig: &Cigar) -> bool {
        match cig {
            Cigar::Match(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_)
            | Cigar::Ins(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_) => true,
            _ => false,
        }
    }


    pub fn cigar_consumes_read_bases(cig: &Cigar) -> bool {
        // Consumes read bases
        match cig {
            Cigar::Match(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_)
            | Cigar::Ins(_)
            | Cigar::SoftClip(_) => true,
            _ => false,
        }
    }

    pub fn cigar_consumes_reference_bases(cig: &Cigar) -> bool {
        // consumes reference bases
        match cig {
            Cigar::Match(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => true,
            _ => false,
        }
    }

    pub fn cigar_is_soft_clip(cig: &Cigar) -> bool {
        match cig {
            Cigar::SoftClip(_) => true,
            _ => false,
        }
    }

    pub fn cigar_is_hard_clip(cig: &Cigar) -> bool {
        match cig {
            Cigar::HardClip(_) => true,
            _ => false
        }
    }

    /**
     * Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin
     * @param start initial index to clip within read bases, inclusive
     * @param stop final index to clip within read bases exclusive
     * @param clippingOperator      type of clipping -- must be either hard clip or soft clip
     */
    pub fn clip_cigar(
        cigar: &CigarStringView,
        start: u32,
        stop: u32,
        clipping_operator: Cigar,
    ) -> CigarString {
        let clip_left = start == 0;

        let mut new_cigar = CigarBuilder::new(true);

        let mut element_start = 0;
        for element in cigar.iter() {
            match element {
                // copy hard clips
                Cigar::HardClip(len) => new_cigar.add(Cigar::HardClip(*len)).unwrap(),
                Cigar::SoftClip(len)
                | Cigar::Diff(len)
                | Cigar::Equal(len)
                | Cigar::RefSkip(len)
                | Cigar::Del(len)
                | Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::Pad(len) => {
                    let element_end = element_start
                        + if CigarUtils::cigar_consumes_read_bases(element) {
                            *len
                        } else {
                            0
                        };

                    // element precedes start or follows end of clip, copy it to new cigar
                    if element_end <= start || element_start >= stop {
                        // edge case: deletions at edge of clipping are meaningless and we skip them
                        if CigarUtils::cigar_consumes_read_bases(element)
                            || (element_start != start && element_start != stop)
                        {
                            new_cigar.add(element.clone()).unwrap()
                        }
                    } else {
                        // otherwise, some or all of the element is soft-clipped
                        let unclipped_length = if clip_left {
                            element_end.checked_sub(stop)
                        } else {
                            start.checked_sub(element_start)
                        };

                        match unclipped_length {
                            None => {
                                // Totally clipped
                                if CigarUtils::cigar_consumes_read_bases(element) {
                                    new_cigar
                                        .add(CigarUtils::cigar_from_element_and_length(
                                            &clipping_operator,
                                            element.len(),
                                        ))
                                        .unwrap()
                                }
                            }
                            Some(unclipped_length) => {
                                if unclipped_length == 0 {
                                    // Totally clipped
                                    if CigarUtils::cigar_consumes_read_bases(element) {
                                        new_cigar
                                            .add(CigarUtils::cigar_from_element_and_length(
                                                &clipping_operator,
                                                element.len(),
                                            ))
                                            .unwrap()
                                    }
                                } else {
                                    let clipped_length = len.checked_sub(unclipped_length).unwrap();
                                    if clip_left {
                                        new_cigar
                                            .add(CigarUtils::cigar_from_element_and_length(
                                                &clipping_operator,
                                                clipped_length,
                                            ))
                                            .unwrap();
                                        new_cigar
                                            .add(CigarUtils::cigar_from_element_and_length(
                                                element,
                                                unclipped_length,
                                            ))
                                            .unwrap();
                                    } else {
                                        new_cigar
                                            .add(CigarUtils::cigar_from_element_and_length(
                                                element,
                                                unclipped_length,
                                            ))
                                            .unwrap();
                                        new_cigar
                                            .add(CigarUtils::cigar_from_element_and_length(
                                                &clipping_operator,
                                                clipped_length,
                                            ))
                                            .unwrap();
                                    }
                                }
                            }
                        }
                    };
                    element_start = element_end
                }
            }
        }
        return new_cigar.make(false).unwrap();
    }

    /**
     * replace soft clips (S) with match (M) operators, normalizing the result by all the transformations of the {@link CigarBuilder} class:
     * merging consecutive identical operators and removing zero-length elements.  For example 10S10M -> 20M and 10S10M10I10I -> 20M20I.
     */
    pub fn revert_soft_clips(cigar: &CigarStringView) -> CigarString {
        let mut builder = CigarBuilder::new(true);
        for element in cigar.iter() {
            match element {
                Cigar::SoftClip(length) => builder
                    .add(CigarUtils::cigar_from_element_and_length(
                        &Cigar::Match(0),
                        *length,
                    ))
                    .unwrap(),
                _ => builder.add(element.clone()).unwrap(),
            }
        }
        return builder.make(false).unwrap();
    }

    /**
     * How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips
     */
    pub fn alignment_start_shift(cigar: &CigarStringView, num_clipped: i64) -> i64 {
        let mut ref_bases_clipped = 0;

        let mut element_start = 0; // this and elementEnd are indices in the read's bases
        for element in cigar.iter() {
            match element {
                // copy hard clips
                Cigar::HardClip(len) => continue,
                Cigar::SoftClip(len)
                | Cigar::Diff(len)
                | Cigar::Equal(len)
                | Cigar::RefSkip(len)
                | Cigar::Del(len)
                | Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::Pad(len) => {
                    let element_end = element_start
                        + if CigarUtils::cigar_consumes_read_bases(element) {
                            *len as i64
                        } else {
                            0
                        };

                    if element_end <= num_clipped {
                        // totally within clipped span -- this includes deletions immediately following clipping
                        ref_bases_clipped += if CigarUtils::cigar_consumes_reference_bases(element)
                        {
                            *len as i64
                        } else {
                            0
                        };
                    } else if element_start < num_clipped {
                        // clip in middle of element, which means the element necessarily consumes read bases
                        let clipped_length = num_clipped - element_start;
                        ref_bases_clipped += if CigarUtils::cigar_consumes_reference_bases(element)
                        {
                            clipped_length
                        } else {
                            0
                        };
                        break;
                    }
                    element_start = element_end;
                }
            }
        }
        return ref_bases_clipped;
    }

    pub fn cigar_from_element_and_length(cigar: &Cigar, length: u32) -> Cigar {
        match cigar {
            Cigar::Pad(_) => return Cigar::Pad(length),
            Cigar::Ins(_) => return Cigar::Ins(length),
            Cigar::Match(_) => return Cigar::Match(length),
            Cigar::Del(_) => return Cigar::Del(length),
            Cigar::RefSkip(_) => return Cigar::RefSkip(length),
            Cigar::Equal(_) => return Cigar::Equal(length),
            Cigar::Diff(_) => return Cigar::Diff(length),
            Cigar::SoftClip(_) => return Cigar::SoftClip(length),
            Cigar::HardClip(_) => return Cigar::HardClip(length),
        }
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence.
     *
     * This assumes that the reference and alt sequence are haplotypes derived from a de Bruijn graph or SeqGraph and have the same
     * ref source and ref sink vertices.  That is, the alt sequence start and end are assumed anchored to the reference start and end, which
     * occur at the ends of the padded assembly region.  Hence, unlike read alignment, there is no concept of a start or end coordinate here.
     * Furthermore, it is important to note that in the rare case that the alt cigar begins or ends with a deletion, we must keep the leading
     * or trailing deletion in order to maintain the original reference span of the alt haplotype.  This can occur, for example, when the ref
     * haplotype starts with N repeats of a long sequence and the alt haplotype starts with N-1 repeats.
     *
     * @param aligner
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    pub fn calculate_cigar(
        ref_seq: &[u8],
        alt_seq: &[u8],
        strategy: OverhangStrategy,
        sw_parameters: &Parameters,
        avx_mode: AVXMode,
    ) -> Option<CigarString> {
        if alt_seq.len() == 0 {
            // horrible edge case from the unit tests, where this path has no bases
            return Some(CigarString::from(vec![Cigar::Del(ref_seq.len() as u32)]));
        }

        //Note: this is a performance optimization.
        // If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
        // Furthermore, if their lengths are equal and their element-by-element comparison yields two or fewer mismatches
        // it's also a trivial M-only CIGAR, because in order to have equal length one would need at least one insertion and
        // one deletion, in which case two substitutions is a better alignment.
        if alt_seq.len() == ref_seq.len() {
            let mismatch_count = (0..ref_seq.len())
                .into_iter()
                .map(|n| if alt_seq[n] == ref_seq[n] { 0 } else { 1 })
                .sum::<usize>();

            if mismatch_count <= 2 {
                let matching = CigarString::from(vec![Cigar::Match(ref_seq.len() as u32)]);
                return Some(matching);
            }
        }

        let padded_ref = format!(
            "{}{}{}",
            *SW_PAD,
            std::str::from_utf8(ref_seq).unwrap(),
            *SW_PAD
        );
        let padded_path = format!(
            "{}{}{}",
            *SW_PAD,
            std::str::from_utf8(alt_seq).unwrap(),
            *SW_PAD
        );
        let alignment = SmithWatermanAligner::align(
            padded_ref.as_bytes(),
            padded_path.as_bytes(),
            sw_parameters,
            strategy,
            avx_mode,
        );


        if Self::is_s_w_failure(&alignment) {
            return None;
        }

        // cut off the padding bases
        let base_start = SW_PAD.len();
        let base_end = padded_path.len() - SW_PAD.len() - 1; // -1 becuase it is inclusive

        let mut trimmed_cigar_and_deletions_removed = AlignmentUtils::trim_cigar_by_bases(
            alignment.cigar,
            base_start as u32,
            base_end as u32,
        );

        let mut non_standard = trimmed_cigar_and_deletions_removed.cigar.0;
        if trimmed_cigar_and_deletions_removed.trailing_deletion_bases_removed > 0 {
            non_standard.push(Cigar::Del(
                trimmed_cigar_and_deletions_removed.trailing_deletion_bases_removed,
            ));
        }

        let mut left_alignment_result = AlignmentUtils::left_align_indels(
            CigarString::from(non_standard),
            ref_seq,
            alt_seq,
            trimmed_cigar_and_deletions_removed.leading_deletion_bases_removed,
        );

        // we must account for possible leading deletions removed when trimming the padding and when left-aligning
        // trailing deletions removed when trimming have already been restored for left-alignment, but left-alingment may have removed them again.
        let total_leading_deletions_removed = trimmed_cigar_and_deletions_removed
            .leading_deletion_bases_removed
            + left_alignment_result.leading_deletion_bases_removed;
        let total_trailing_deletions_removed =
            left_alignment_result.trailing_deletion_bases_removed;

        if total_leading_deletions_removed == 0 && total_trailing_deletions_removed == 0 {
            return Some(left_alignment_result.cigar);
        } else {
            let mut result_elements = Vec::new();
            if total_leading_deletions_removed > 0 {
                result_elements.push(Cigar::Del(total_leading_deletions_removed))
            }

            result_elements.extend(left_alignment_result.cigar.0);
            if total_trailing_deletions_removed > 0 {
                result_elements.push(Cigar::Del(total_trailing_deletions_removed))
            }
            return Some(CigarString::from(result_elements));
        }
    }

    pub fn is_indel(cigar: &Cigar) -> bool {
        match cigar {
            Cigar::Del(_) | Cigar::Ins(_) => true,
            _ => false,
        }
    }

    /**
     * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
     */
    fn is_s_w_failure(alignment: &SmithWatermanAlignmentResult) -> bool {
        // check that the alignment starts at the first base, which it should given the padding
        if alignment.alignment_offset > 0 {
            return true;
        }

        // check that we aren't getting any S operators (which would be very bad downstream)
        for ce in alignment.cigar.0.iter() {
            match ce {
                Cigar::SoftClip(_) => {
                    return true;
                    // soft clips at the end of the alignment are really insertions
                }
                _ => continue,
            }
        }

        return false;
    }

    /**
     * Check to see if two cigar operators are of the same type regardless of length
     */
    pub fn cigar_elements_are_same_type(this: &Cigar, other: &Option<Cigar>) -> bool {
        match other {
            None => false,
            Some(element) => match element {
                Cigar::SoftClip(_) => match this {
                    Cigar::SoftClip(_) => true,
                    _ => false,
                },
                Cigar::HardClip(_) => match this {
                    Cigar::HardClip(_) => true,
                    _ => false,
                },
                Cigar::Match(_) => match this {
                    Cigar::Match(_) => true,
                    _ => false,
                },
                Cigar::Diff(_) => match this {
                    Cigar::Diff(_) => true,
                    _ => false,
                },
                Cigar::Equal(_) => match this {
                    Cigar::Equal(_) => true,
                    _ => false,
                },
                Cigar::RefSkip(_) => match this {
                    Cigar::RefSkip(_) => true,
                    _ => false,
                },
                Cigar::Ins(_) => match this {
                    Cigar::Ins(_) => true,
                    _ => false,
                },
                Cigar::Del(_) => match this {
                    Cigar::Del(_) => true,
                    _ => false,
                },
                Cigar::Pad(_) => match this {
                    Cigar::Pad(_) => true,
                    _ => false,
                },
            },
        }
    }

    pub fn is_clipping(element: &Cigar) -> bool {
        match element {
            Cigar::SoftClip(_) | Cigar::HardClip(_) => true,
            _ => false,
        }
    }

    /**
     * Combine two cigar operators of equal type and add their lengths creating new operator
     */
    pub fn combine_cigar_operators(this: &Cigar, other: &Cigar) -> Option<Cigar> {
        match other {
            &Cigar::SoftClip(other_len) => match this {
                &Cigar::SoftClip(this_len) => Some(Cigar::SoftClip(other_len + this_len)),
                _ => None,
            },
            &Cigar::HardClip(other_len) => match this {
                &Cigar::HardClip(this_len) => Some(Cigar::HardClip(other_len + this_len)),
                _ => None,
            },
            &Cigar::Match(other_len) => match this {
                &Cigar::Match(this_len) => Some(Cigar::Match(other_len + this_len)),
                _ => None,
            },
            &Cigar::Equal(other_len) => match this {
                &Cigar::Equal(this_len) => Some(Cigar::Equal(other_len + this_len)),
                _ => None,
            },
            &Cigar::Diff(other_len) => match this {
                &Cigar::Diff(this_len) => Some(Cigar::Diff(other_len + this_len)),
                _ => None,
            },
            &Cigar::RefSkip(other_len) => match this {
                &Cigar::RefSkip(this_len) => Some(Cigar::RefSkip(other_len + this_len)),
                _ => None,
            },
            &Cigar::Ins(other_len) => match this {
                &Cigar::Ins(this_len) => Some(Cigar::Ins(other_len + this_len)),
                _ => None,
            },
            &Cigar::Del(other_len) => match this {
                &Cigar::Del(this_len) => Some(Cigar::Del(other_len + this_len)),
                _ => None,
            },
            &Cigar::Pad(other_len) => match this {
                &Cigar::Pad(this_len) => Some(Cigar::Pad(other_len + this_len)),
                _ => None,
            },
        }
    }

    /**
     * Creates a new cigar element using the same operator as the provided cigar but with the new
     * provided length
     */
    pub fn new_cigar_from_operator_and_length(old_cigar: &Cigar, new_length: u32) -> Cigar {
        match old_cigar {
            Cigar::SoftClip(_) => Cigar::SoftClip(new_length),
            Cigar::HardClip(_) => Cigar::HardClip(new_length),
            Cigar::Match(_) => Cigar::Match(new_length),
            Cigar::Equal(_) => Cigar::Equal(new_length),
            Cigar::Diff(_) => Cigar::Diff(new_length),
            Cigar::RefSkip(_) => Cigar::RefSkip(new_length),
            Cigar::Ins(_) => Cigar::Ins(new_length),
            Cigar::Del(_) => Cigar::Del(new_length),
            Cigar::Pad(_) => Cigar::Pad(new_length),
        }
    }

    /**
     * @return The number of reference bases that the read covers, excluding padding.
     */
    pub fn get_reference_length(cigar: &CigarString) -> u32 {
        let length = cigar
            .0
            .iter()
            .map(|elem| match elem {
                Cigar::Match(len)
                | Cigar::Del(len)
                | Cigar::RefSkip(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len) => *len,
                _ => 0,
            })
            .sum::<u32>();

        return length;
    }

    /**
     * @return The number of reference bases that the read covers, including padding.
     */
    pub fn get_padded_reference_length(cigar: &CigarString) -> u32 {
        let length = cigar
            .0
            .iter()
            .map(|elem| match elem {
                Cigar::Match(len)
                | Cigar::Del(len)
                | Cigar::RefSkip(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len)
                | Cigar::Pad(len) => *len,
                _ => 0,
            })
            .sum::<u32>();

        return length;
    }

    /**
     * @return The number of read bases that the read covers.
     */
    pub fn get_read_length(cigar: &CigarString) -> u32 {
        let length = cigar
            .0
            .iter()
            .filter(|elem| CigarUtils::cigar_consumes_read_bases(elem))
            .map(|elem| elem.len())
            .sum::<u32>();

        return length;
    }

    /** Returns true if the operator is a M, a X or a EQ */
    pub fn is_alignment(cigar: &Cigar) -> bool {
        match cigar {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
            _ => false,
        }
    }

    /**
     * A good Cigar object obeys the following rules:
     *  - is valid as per SAM spec {@link Cigar#isValid(String, long)}.
     *  - has no consecutive I/D elements
     *  - does not start or end with deletions (with or without preceding clips).
     */
    pub fn is_good(c: &CigarString) -> bool {
        return !(Self::has_consecutive_indels(&c.0)
            || Self::starts_or_ends_with_deletion_ignoring_clips(&c.0));
    }

    /**
     * Checks if cigar has consecutive I/D elements.
     */
    pub fn has_consecutive_indels(elems: &Vec<Cigar>) -> bool {
        let mut prev_indel = false;
        for elem in elems {
            match elem {
                &Cigar::Ins(_) | &Cigar::Del(_) => {
                    let is_indel = true;
                    if prev_indel && is_indel {
                        return true;
                    };
                    prev_indel = is_indel;
                }
                _ => {
                    prev_indel = false;
                }
            };
        }

        return false;
    }

    /**
     * Checks if cigar starts with a deletion (ignoring any clips at the beginning).
     */
    pub fn starts_or_ends_with_deletion_ignoring_clips(elems: &Vec<Cigar>) -> bool {
        for left_side in &[true, false] {
            if *left_side {
                for elem in elems.iter() {
                    match elem {
                        Cigar::Del(_) => return true,
                        _ => {
                            if !Self::is_clipping(elem) {
                                break;
                            }
                        }
                    }
                }
            } else {
                for elem in elems.iter().rev() {
                    match elem {
                        Cigar::Del(_) => return true,
                        _ => {
                            if !Self::is_clipping(elem) {
                                break;
                            }
                        }
                    }
                }
            }
        }

        return false;
    }
}
