use rayon::prelude::*;
use reads::cigar_builder::{CigarBuilderResult, CigarBuilder};
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};
use reads::cigar_utils::CigarUtils;

pub struct AlignmentUtils {

}

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
        return Self::trim_cigar(cigar, start, end, false)
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
    fn trim_cigar(cigar: CigarString, start: u32, end: u32, by_reference: bool) -> CigarBuilderResult {
        assert!(end >= start, "End position cannot be before start position");

        let mut new_elements = CigarBuilder::new(true);

        // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read (otherwise) coordinates
        let mut element_start; // inclusive
        let mut element_end = 0; // exclusive -- start of next element

        for elt in cigar.iter() {
            element_start = element_end;
            element_end = element_start + (if by_reference { Self::length_on_reference(elt) } else { Self::length_on_read(elt) });

            // we are careful to include zero-length elements at both ends, that is,
            // elements with elementStart == elementEnd == start and elementStart == elementEnd == end + 1
            if element_end < start || (element_end == start && element_start < start) {
                continue
            } else if element_start > end && element_end > end + 1 {
                break
            }

            let overlap_length = if element_end == element_start { elt.len() } else { std::cmp::min(end + 1, element_end) - std::cmp::max(start, element_start) };
            new_elements.add(CigarUtils::new_cigar_from_operator_and_length(elt, overlap_length));
        }

        assert!(element_end > end, "Cigar elements don't reach end position (inclusive)");

        return new_elements.make_and_record_deletions_removed_result()
    }

    fn length_on_read(element: &Cigar) -> u32 {
        return if CigarUtils::cigar_consumes_read_bases(element) { element.len() } else { 0 }
    }

    fn length_on_reference(element: &Cigar) -> u32 {
        return if CigarUtils::cigar_consumes_reference_bases(element) { element.len() } else { 0 }
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
    pub fn left_align_indels(cigar: CigarString, ref_seq: &[u8], read: &[u8], read_start: u32) -> CigarBuilderResult {
        if &cigar.0.par_iter().all(|elem| !CigarUtils::is_indel(elem)) {
            return CigarBuilderResult::new(cigar, 0, 0)
        }

        // we need reference bases from the start of the read to the rightmost indel
        let last_indel = (0..cigar.len()).into_par_iter().filter(|n| CigarUtils::is_indel(&cigar.0[n])).max();
        let necessary_ref_length = read_start + cigar.0.iter().take(last_indel + 1).map(|e| Self::length_on_reference(e)).sum();
        assert!(necessary_ref_length <= ref_seq.len(), "Read goes past end of reference");

        // at this point, we are one base past the end of the read.  Now we traverse the cigar from right to left
        let mut result_right_to_left = Vec::new();
        let ref_length = CigarUtils::get_reference_length(&cigar);
        let ref_indel_range = (read_start + ref_length..read_start + ref_length);
        let read_indel_range = (read.len()..read.len());
    }
}