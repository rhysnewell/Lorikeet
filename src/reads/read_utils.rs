use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_utils::CigarUtils;
use rust_htslib::bam::record::{Aux, AuxArray, Cigar, CigarString, CigarStringView};
use std::cmp::Ordering;
use utils::simple_interval::Locatable;

pub struct ReadUtils {}

impl ReadUtils {
    pub const CANNOT_COMPUTE_ADAPTOR_BOUNDARY: usize = std::usize::MIN;

    /**
     * The default quality score for an insertion or deletion, if
     * none are provided for this read.
     */
    pub const DEFAULT_INSERTION_DELETION_QUAL: u8 = 45;

    /**
     * Find the 0-based index within a read base array corresponding to a given 0-based position in the reference, along with the cigar operator of
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
        alignment_start: usize,
        cigar: CigarStringView,
        ref_coord: usize,
    ) -> (Option<usize>, Option<Cigar>) {
        if ref_coord < alignment_start {
            return (None, None);
        }

        // let mut current_read_pos = alignment_start;
        // let mut current_read_index = 0;
        //
        // for cig in cigar.iter() {
        //     if ref_coord - 1 <= current_read_pos {
        //         let read_index_at_position = current_read_index as i64 +
        //             if CigarUtils::cigar_consumes_read_bases(cig) {
        //                 ref_coord as i64 - current_read_pos as i64
        //             } else {
        //                 0
        //             };
        //
        //         println!("Found read index {} for ref coord {}", read_index_at_position, ref_coord);
        //         return (Some(read_index_at_position as usize), Some(cig.clone()))
        //     };
        //
        //     current_read_pos += if CigarUtils::cigar_consumes_reference_bases(cig)
        //         || CigarUtils::cigar_is_soft_clip(cig) {
        //         cig.len() as usize
        //     } else {
        //         0
        //     };
        //
        //     current_read_index += if CigarUtils::cigar_consumes_read_bases(cig) {
        //         cig.len() as usize
        //     } else {
        //         0
        //     }
        //
        //
        // }

        let mut first_read_pos_of_element = 0; // inclusive
        let mut first_ref_pos_of_element = alignment_start; // inclusive
        let mut last_read_pos_of_element = 0; // exclusive
        let mut last_ref_pos_of_element = alignment_start; // exclusive

        // advance forward through all the cigar elements until we bracket the reference coordinate
        for cig in cigar.iter() {
            first_read_pos_of_element = last_read_pos_of_element;
            first_ref_pos_of_element = last_ref_pos_of_element;
            last_read_pos_of_element += if CigarUtils::cigar_consumes_read_bases(cig) {
                cig.len() as usize
            } else {
                0
            };
            last_ref_pos_of_element += if CigarUtils::cigar_consumes_reference_bases(cig)
                || CigarUtils::cigar_is_soft_clip(cig)
            {
                cig.len() as usize
            } else {
                0
            };

            if first_ref_pos_of_element <= ref_coord && ref_coord < last_ref_pos_of_element {
                // refCoord falls within this cigar element
                let read_pos_at_ref_coord = first_read_pos_of_element as i64
                    + (if CigarUtils::cigar_consumes_read_bases(cig) {
                        (ref_coord) as i64 - first_ref_pos_of_element as i64
                    // ref_coord.checked_sub(first_ref_pos_of_element).unwrap_or(0)
                    } else {
                        0
                    });

                return (Some(read_pos_at_ref_coord as usize), Some(cig.clone()));
            }
        }
        return (None, None);
    }

    /**
     * Returns the index within the read's bases array corresponding to the requested reference coordinate -- or the read coordinate immediately preceding
     * a deletion in which the reference coordinate falls -- along with the cigar operator in which the reference coordinate occurs.
     */
    pub fn get_read_index_for_reference_coordinate_from_read(
        read: &BirdToolRead,
        ref_coord: usize,
    ) -> (Option<usize>, Option<Cigar>) {
        Self::get_read_index_for_reference_coordinate(
            read.get_soft_start().unwrap(),
            read.read.cigar(),
            ref_coord,
        )
    }

    pub fn empty_read(read: &BirdToolRead) -> BirdToolRead {
        let mut empty_read = read.clone();

        empty_read.read.set_mate_unmapped();
        empty_read.read.set_unmapped();
        empty_read.read.set_mapq(0);
        empty_read
            .read
            .set(read.read.qname(), Some(&CigarString(Vec::new())), &[], &[]);

        return empty_read;
    }

    pub fn empty_read_mut(empty_read: &mut BirdToolRead) {
        let qname = empty_read.read.qname().to_vec();
        empty_read.read.set_mate_unmapped();
        empty_read.read.set_unmapped();
        empty_read.read.set_mapq(0);
        empty_read
            .read
            .set(&qname[..], Some(&CigarString(Vec::new())), &[], &[]);
    }

    pub fn empty_read_with_quals_and_bases_mut(
        empty_read: &mut BirdToolRead,
        quals: &[u8],
        bases: &[u8],
    ) {
        let qname = empty_read.read.qname().to_vec();
        empty_read.read.set_mate_unmapped();
        empty_read.read.set_unmapped();
        empty_read.read.set_mapq(0);
        empty_read.read.set(
            &qname[..],
            Some(&CigarString(empty_read.read.cigar().0.clone())),
            bases,
            quals,
        );
    }

    pub fn empty_read_with_quals_and_bases(
        read: &BirdToolRead,
        quals: &[u8],
        bases: &[u8],
    ) -> BirdToolRead {
        let mut empty_read = read.clone();

        empty_read.read.set_mate_unmapped();
        empty_read.read.set_unmapped();
        empty_read.read.set_mapq(0);
        empty_read.read.set(
            read.read.qname(),
            Some(&CigarString(read.read.cigar().0.clone())),
            bases,
            quals,
        );

        return empty_read;
    }

    pub fn set_insertion_base_qualities(read: &mut BirdToolRead, base_ins_quals: &Vec<u8>) {
        read.read
            .push_aux(b"BI", Aux::ArrayU8(AuxArray::from(base_ins_quals)));
    }

    pub fn set_deletion_base_qualities(read: &mut BirdToolRead, base_del_quals: &Vec<u8>) {
        read.read
            .push_aux(b"BD", Aux::ArrayU8(AuxArray::from(base_del_quals)));
    }

    pub fn compare_coordinates(first: &BirdToolRead, second: &BirdToolRead) -> Ordering {
        let first_ref_index = first.read.tid();
        let second_ref_index = second.read.tid();

        if first_ref_index == -1 {
            if second_ref_index == -1 {
                return Ordering::Equal;
            } else {
                return Ordering::Greater;
            }
        } else if second_ref_index == -1 {
            return Ordering::Less;
        }

        if first_ref_index != second_ref_index {
            return first_ref_index.cmp(&second_ref_index);
        } else {
            return first.get_start().cmp(&second.get_start());
        }
    }

    /**
     * Can the adaptor sequence of read be reliably removed from the read based on the alignment of
     * read and its mate?
     *
     * @param read the read to check
     * @return true if it can, false otherwise
     */
    pub fn has_well_defined_fragment_size(read: &BirdToolRead) -> bool {
        if read.read.insert_size() == 0 {
            // no adaptors in reads with mates in another chromosome or unmapped pai
            return false;
        }

        if !read.read.is_paired() {
            // only reads that are paired can be adaptor trimmed
            return false;
        }

        if read.read.is_unmapped() || read.read.is_mate_unmapped() {
            // only reads when both reads are mapped can be trimmed
            return false;
        }

        if read.read.is_reverse() == read.read.is_mate_reverse() {
            // sanity check on isProperlyPaired to ensure that read1 and read2 aren't on the same strand
            return false;
        }

        if read.read.is_reverse() {
            // we're on the negative strand, so our read runs right to left
            return read.get_end() as i64 > read.read.mpos();
        } else {
            // we're on the positive strand, so our mate should be to our right (his start + insert size should be past our start)
            return read.get_start() as i64 <= read.read.mpos() + read.read.insert_size();
        }
    }

    /**
     * Finds the adaptor boundary around the read and returns the first base inside the adaptor that is closest to
     * the read boundary. If the read is in the positive strand, this is the first base after the end of the
     * fragment (Picard calls it 'insert'), if the read is in the negative strand, this is the first base before the
     * beginning of the fragment.
     *
     * There are two cases we need to treat here:
     *
     * 1) Our read is in the reverse strand :
     *
     *     <----------------------| *
     *   |--------------------->
     *
     *   in these cases, the adaptor boundary is at the mate start (minus one)
     *
     * 2) Our read is in the forward strand :
     *
     *   |---------------------->   *
     *     <----------------------|
     *
     *   in these cases the adaptor boundary is at the start of the read plus the inferred insert size (plus one)
     *
     * @param read the read being tested for the adaptor boundary
     * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read.
     * CANNOT_COMPUTE_ADAPTOR_BOUNDARY if the read is unmapped or the mate is mapped to another contig.
     */
    pub fn get_adaptor_boundary(read: &BirdToolRead) -> usize {
        if !Self::has_well_defined_fragment_size(read) {
            return Self::CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
        } else if read.read.is_reverse() {
            return read.read.mpos() as usize - 1;
        } else {
            let insert_size = read.read.insert_size().abs() as usize;
            return read.get_start() + insert_size;
        }
    }

    /**
     * Is a base inside a read?
     *
     * @param read                the read to evaluate
     * @param referenceCoordinate the reference coordinate of the base to test
     * @return true if it is inside the read, false otherwise.
     */
    pub fn is_inside_read(read: &BirdToolRead, reference_coordinate: usize) -> bool {
        return reference_coordinate >= read.get_start() && reference_coordinate <= read.get_end();
    }

    /**
     * Default utility to query the base insertion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base insertion quality array
     */
    pub fn get_base_insertion_qualities(read: &BirdToolRead) -> Vec<u8> {
        debug!("read aux tag {:?}", read.read.aux(b"BI"));
        match read.read.aux(b"BI") {
            Ok(tag) => {
                match tag {
                    Aux::ArrayU8(quals) => return quals.iter().collect::<Vec<u8>>(),
                    _ => {
                        // Some day in the future when base insertion and base deletion quals exist the samtools API will
                        // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                        vec![Self::DEFAULT_INSERTION_DELETION_QUAL; read.read.qual().len()]
                    }
                }
            }
            _ => {
                // Some day in the future when base insertion and base deletion quals exist the samtools API will
                // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                vec![Self::DEFAULT_INSERTION_DELETION_QUAL; read.read.qual().len()]
            }
        }
    }

    /**
     * Default utility to query the base Deletion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base insertion quality array
     */
    pub fn get_base_deletion_qualities(read: &BirdToolRead) -> Vec<u8> {
        debug!("read aux tag {:?}", read.read.aux(b"BD"));
        match read.read.aux(b"BD") {
            Ok(tag) => {
                match tag {
                    Aux::ArrayU8(quals) => return quals.iter().collect::<Vec<u8>>(),
                    _ => {
                        // Some day in the future when base insertion and base deletion quals exist the samtools API will
                        // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                        vec![Self::DEFAULT_INSERTION_DELETION_QUAL; read.read.qual().len()]
                    }
                }
            }
            _ => {
                // Some day in the future when base insertion and base deletion quals exist the samtools API will
                // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                vec![Self::DEFAULT_INSERTION_DELETION_QUAL; read.read.qual().len()]
            }
        }
    }
}
