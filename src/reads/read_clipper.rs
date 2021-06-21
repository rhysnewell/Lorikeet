use rust_htslib::bam::Record;
use reads::bird_tool_reads::BirdToolRead;
use bio_types::sequence::SequenceRead;
use reads::read_utils::ReadUtils;

/**
 * A comprehensive clipping tool.
 *
 * General Contract:
 *  - All clipping operations return a new read with the clipped bases requested, it never modifies the original read.
 *  - If a read is fully clipped, return an empty SAMRecord, never null.
 *  - When hard clipping, add cigar operator H for every *reference base* removed (i.e. Matches, SoftClips and Deletions, but *not* insertions). See Hard Clipping notes for details.
 *
 *
 * There are several types of clipping to use:
 *
 * Write N's:
 *   Change the bases to N's in the desired region. This can be applied anywhere in the read.
 *
 * Write Q0's:
 *   Change the quality of the bases in the desired region to Q0. This can be applied anywhere in the read.
 *
 * Write both N's and Q0's:
 *   Same as the two independent operations, put together.
 *
 * Soft Clipping:
 *   Do not change the read, just mark the reads as soft clipped in the Cigar String
 *   and adjust the alignment start and end of the read.
 *
 * Hard Clipping:
 *   Creates a new read without the hard clipped bases (and base qualities). The cigar string
 *   will be updated with the cigar operator H for every reference base removed (i.e. Matches,
 *   Soft clipped bases and deletions, but *not* insertions). This contract with the cigar
 *   is necessary to allow read.getUnclippedStart() / End() to recover the original alignment
 *   of the read (before clipping).
 *
 */
pub struct ReadClipper {}

impl ReadClipper {

    /**
     * Hard clip the read to the variable region (from refStart to refStop)
     *
     * @param read     the read to be clipped
     * @param refStart the beginning of the variant region (inclusive)
     * @param refStop  the end of the variant region (inclusive)
     * @return the read hard clipped to the variant region (Could return an empty, unmapped read)
     */
    pub fn hard_clip_to_region(
        read: BirdToolRead, ref_start: usize, ref_stop: usize
    ) -> BirdToolRead {
        let start = read.get_start();
        let end = read.get_end();
        return Self::hard_clip_to_region_with_alignment_interval()
    }

    fn hard_clip_to_region_with_alignment_interval(
        read: BirdToolRead, ref_start: usize, ref_stop: usize, alignment_start: usize, alignment_stop: usize
    ) -> BirdToolRead {
        if alignment_start <= ref_stop && alignment_stop >= ref_start {
            if alignment_start < ref_start && alignment_stop > ref_stop {
                return Self::
            }
        }
    }

    /**
     * Generic functionality to  clip a read, used internally by hardClipByReferenceCoordinatesLeftTail
     * and hardClipByReferenceCoordinatesRightTail. Should not be used directly.
     *
     * Note, it REQUIRES you to give the directionality of your hard clip (i.e. whether you're clipping the
     * left of right tail) by specifying either refStart == None or refStop == None.
     *
     * @param refStart  first base to clip (inclusive)
     * @param refStop last base to clip (inclusive)
     * @param clippingOp clipping operation to be performed
     * @return a new read, without the clipped bases (May return empty, unclipped reads)
     */
    fn clip_by_reference_coordinates(
        read: &BirdToolRead, ref_start: Option<usize>, ref_stop: Option<usize>, clipping_op: ClippingRepresentation
    ) -> BirdToolRead {
        if read.read.is_empty() {
            return read
        }
        if clipping_op == ClippingRepresentation::SOFTCLIP_BASES && read.read.is_unmapped() {
            panic!("Cannot soft-clip read {:?} by reference coordinates because it is unmapped", read)
        }

        let mut start;
        let mut stop;

        // Determine the read coordinate to start and stop hard clipping
        match ref_start {
            None => {
                match ref_stop {
                    None => {
                        panic!("Only one of ref_start or ref_stop can be None, not both.")
                    },
                    Some(ref_stop) => {
                        start = Some(0);
                        let stop_pos_and_operator = ReadUtils::get_read_index_for_reference_coordinate_from_read(read, ref_stop);
                        match stop_pos_and_operator.0 {
                            Some(pos) => {
                                stop = Some(pos - (if ReadUtils::cigar_consumes_read_bases(stop_pos_and_operator.0.unwrap()) { 0 } else { 1 }));
                            },
                            None => {
                                stop = None;
                            }
                        }
                    }
                }
            },
            Some(ref_start) => {
                match ref_stop {
                    None => {
                        start = ReadUtils::get_read_index_for_reference_coordinate_from_read(read, ref_start).0.unwrap();
                        stop = read.read.seq_len() as usize - 1;
                    },
                    Some(ref_stop) => {
                        panic!("Either ref_start or ref_stop needs to be None")
                    }
                }
            }

        }
    }
}


/**
 * How should we represent a clipped bases in a read?
 */
pub enum ClippingRepresentation {
    /** Clipped bases are changed to Ns */
    WRITE_NS,

    /** Clipped bases are changed to have Q0 quality score */
    WRITE_Q0S,

    /** Clipped bases are change to have both an N base and a Q0 quality score */
    WRITE_NS_Q0S,

    /**
     * Change the read's cigar string to soft clip (S, see sam-spec) away the bases.
     * Note that this can only be applied to cases where the clipped bases occur
     * at the start or end of a read.
     */
    SOFTCLIP_BASES,

    /**
     * WARNING: THIS OPTION IS STILL UNDER DEVELOPMENT AND IS NOT SUPPORTED.
     *
     * Change the read's cigar string to hard clip (H, see sam-spec) away the bases.
     * Hard clipping, unlike soft clipping, actually removes bases from the read,
     * reducing the resulting file's size but introducing an irrevesible (i.e.,
     * lossy) operation.  Note that this can only be applied to cases where the clipped
     * bases occur at the start or end of a read.
     */
    HARDCLIP_BASES,

    /**
     * Turn all soft-clipped bases into matches
     */
    REVERT_SOFTCLIPPED_BASES,
}