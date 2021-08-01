use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_builder::{CigarBuilder, CigarBuilderResult};
use reads::cigar_utils::CigarUtils;
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use std::ops::Range;

pub struct AlignmentUtils {}

impl AlignmentUtils {
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
    pub fn get_bases_covering_ref_interval(
        ref_start: usize,
        ref_end: usize,
        bases: &Vec<u8>,
        bases_start_on_ref: usize,
        bases_to_ref_cigar: &CigarString,
    ) -> Option<Vec<u8>> {
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
        }

        if bases_start.is_none() || bases_stop.is_none() {
            panic!(
                "Mever found start {:?} or stop {:?} given cigar {:?}",
                bases_start, bases_stop, bases_to_ref_cigar
            );
        };
        return Some(
            bases[bases_start.unwrap() as usize..bases_stop.unwrap() as usize + 1].to_vec(),
        );
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
