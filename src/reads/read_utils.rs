use reads::bird_tool_reads::BirdToolRead;
use rust_htslib::bam::record::{CigarStringView, Cigar, CigarString};

pub struct ReadUtils {}

impl ReadUtils {
    /**
     * Find the 0-based index within a read base array corresponding to a given 1-based position in the reference, along with the cigar operator of
     * the element containing that base.  If the reference coordinate occurs within a deletion, the first index after the deletion is returned.
     * Note that this treats soft-clipped bases as if they align with the reference, which is useful for hard-clipping reads with soft clips.
     *
     * @param alignmentStart        The soft start of the read on the reference
     * @param cigar                 The read's cigar
     * @param refCoord              The target reference coordinate
     * @return                      If the reference coordinate occurs before the read start or after the read end {@code CLIPPING_GOAL_NOT_REACHED};
     *                              if the reference coordinate falls within an alignment block of the read's cigar, the corresponding read coordinate;
     *                              if the reference coordinate falls within a deletion, the first read coordinate after the deletion.  Note: if the last cigar element is
     *                              a deletion (which isn't meaningful), it returns {@code CLIPPING_GOAL_NOT_REACHED}.
     */
    pub fn get_read_index_for_reference_coordinate(
        alignment_start: usize, cigar: CigarStringView, ref_coord: usize
    ) -> (Option<usize>, Option<Cigar>) {
        if ref_coord < alignment_start {
            return (None, None)
        }

        let mut first_read_pos_of_element = 0;
        let mut first_ref_pos_of_element = alignment_start;
        let mut last_read_pos_of_element = 0;
        let mut last_ref_pos_of_element = alignment_start;

        for cig in cigar.iter() {
            first_read_pos_of_element = last_read_pos_of_element;
            first_ref_pos_of_element = last_ref_pos_of_element;
            last_read_pos_of_element += if Self::cigar_consumes_read_bases(cig) { cig.len() as usize } else { 0 };
            last_ref_pos_of_element += if Self::cigar_consumes_reference_bases(cig) || Self::cigar_is_soft_clip(cig) { cig.len() as usize } else { 0 };

            if first_ref_pos_of_element <= ref_coord && ref_coord < last_ref_pos_of_element { // refCoord falls within this cigar element
                let read_pos_at_ref_coord = first_read_pos_of_element +
                    (if Self::cigar_consumes_read_bases(cig) { ref_coord.checked_sub(first_ref_pos_of_element).unwrap_or(0) } else { 0 });
                return (Some(read_pos_at_ref_coord), Some(cig.clone()))
            }
        }

        return (None, None)
    }

    /**
     * Returns the index within the read's bases array corresponding to the requested reference coordinate -- or the read coordinate immediately preceding
     * a deletion in which the reference coordinate falls -- along with the cigar operator in which the reference coordinate occurs.
     */
    pub fn get_read_index_for_reference_coordinate_from_read(
        read: &BirdToolRead, ref_coord: usize
    ) -> (Option<usize>, Option<Cigar>) {
        Self::get_read_index_for_reference_coordinate(read.get_soft_start(), read.read.cigar(), ref_coord)
    }

    pub fn empty_read(read: &BirdToolRead) -> BirdToolRead {
        let mut empty_read = read.clone();

        empty_read.read.set_mate_unmapped();
        empty_read.read.set_unmapped();
        empty_read.read.set_mapq(0);
        empty_read.read.set(read.read.qname(), Some(&CigarString(vec![""])), &[0], &[0]);

        return empty_read
    }
}