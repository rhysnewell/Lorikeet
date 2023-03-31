use num::traits::AsPrimitive;
use rust_htslib::bam::record::{Cigar, CigarString};

use crate::utils::simple_interval::Locatable;
use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::cigar_utils::CigarUtils;
use crate::reads::read_clipper::ClippingRepresentation;
use crate::reads::read_utils::ReadUtils;

#[derive(Debug)]
pub struct ClippingOp {
    pub start: usize,
    pub stop: usize,
}

impl ClippingOp {
    pub fn new(start: usize, stop: usize) -> ClippingOp {
        ClippingOp { start, stop }
    }

    pub fn get_length(&self) -> usize {
        (self.stop + 1).saturating_sub(self.start)
    }

    /**
     * Clips the bases in read according to this operation's start and stop.  Uses the clipping
     * representation used is the one provided by algorithm argument.
     *  @param algorithm    clipping algorithm to use
     * @param originalRead the read to be clipped
     */
    pub fn apply(&self, algorithm: &ClippingRepresentation, original_read: &mut BirdToolRead) {
        match algorithm {
            // This note probably doesn't apply in Rust as this function has explicit reference of the
            // the BirdToolRead. Kept for posterity
            // important note:
            //   it's not safe to call read.getBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the BirdToolRead
            ClippingRepresentation::WriteNs => {
                Self::apply_write_ns(original_read);
            }
            ClippingRepresentation::WriteQ0s => {
                Self::apply_write_q0s(original_read);
            }
            ClippingRepresentation::WriteNsQ0s => {
                Self::apply_write_ns_and_q0s(original_read);
            }
            ClippingRepresentation::HardclipBases => {
                ClippingOp::apply_hard_clip_bases(original_read, self.start, self.stop);
            }
            ClippingRepresentation::SoftclipBases => {
                self.apply_soft_clip_bases(original_read);
            }
            ClippingRepresentation::RevertSoftclippedBases => {
                self.apply_revert_soft_clipped_bases(original_read);
            }
        }
    }

    fn apply_soft_clip_bases(&self, read: &mut BirdToolRead) {
        assert!(
            !read.read.is_unmapped(),
            "ReadClipper cannot soft clip unmapped reads"
        );

        let new_length = read.len() - (self.stop - self.start + 1);
        if read.read.seq_len() <= 2 || new_length == 0 {
            // pass
        } else {
            let stop = std::cmp::min(self.stop, self.start + read.read.seq_len() - 2);

            assert!(self.start < 1 || stop == read.read.seq_len() - 1,
                    "Cannot apply soft clipping operator to the middle of a read: {:?} to be clipped at {}-{}", read.read.qname(), self.start, stop);

            let old_cigar = read.read.cigar();
            let new_cigar = CigarUtils::clip_cigar(
                &old_cigar,
                self.start as u32,
                (stop + 1) as u32,
                Cigar::SoftClip(0),
            );

            let qname = read.read.qname().to_vec();
            let bases = read.read.seq().as_bytes();
            let quals = read.read.qual().to_vec();

            read.read
                .set(&qname[..], Some(&new_cigar), &bases[..], &quals[..]);

            let alignment_start_shift = if self.start == 0 {
                CigarUtils::alignment_start_shift(&old_cigar, (self.stop + 1) as i64)
            } else {
                0
            };
            let new_start = read.get_start() as i64 + alignment_start_shift;

            read.read.set_pos(new_start);
        }
    }

    fn apply_revert_soft_clipped_bases(&self, read: &mut BirdToolRead) {
        if read.read.cigar().is_empty()
            || !(CigarUtils::is_clipping(read.read.cigar().first().unwrap())
                || CigarUtils::is_clipping(read.read.cigar().last().unwrap()))
        {
            // pass
        } else {
            let unclipped_cigar = CigarUtils::revert_soft_clips(&read.read.cigar());
            let qname = read.read.qname().to_vec();
            let bases = read.read.seq().as_bytes();
            let quals = read.read.qual().to_vec();
            let new_start = read.get_soft_start_i64();

            read.read
                .set(&qname[..], Some(&unclipped_cigar), &bases[..], &quals[..]);

            if new_start <= 0 {
                // if the start of the unclipped read occurs before the contig,
                // we must hard clip away the bases since we cannot represent reads with
                // negative or 0 alignment start values in the SAMRecord (e.g., 0 means unaligned)

                // We cannot set the read to temporarily have a negative start position, as our Read
                // interface will not allow it. Instead, since we know that the final start position will
                // be 0 after the hard clip operation, set it to 0 explicitly. We have to set it twice:
                // once before the hard clip (to reset the alignment stop / read length in read implementations
                // that cache these values, such as SAMRecord), and again after the hard clip.

                read.read.set_pos(0);
                ClippingOp::apply_hard_clip_bases(read, 0, -new_start as usize);

                // Reset the position to 0 again only if we didn't end up with an empty, unmapped read after hard clipping.
                // See https://github.com/broadinstitute/gatk/issues/3845
                if !read.read.is_unmapped() {
                    read.read.set_pos(0);
                }
            } else {
                read.read.set_pos(new_start as i64);
            }
        }
    }

    fn apply_write_q0s(original_read: &mut BirdToolRead) {
        let new_quals = vec![0; original_read.read.qual().len()];
        let qname = original_read.read.qname().to_vec();
        let bases = original_read.read.seq().as_bytes();
        original_read.read.set(
            &qname[..],
            Some(&original_read.read.cigar().take()),
            &bases[..],
            &new_quals,
        );
    }

    fn apply_write_ns(original_read: &mut BirdToolRead) {
        let new_bases = vec!['N'.as_(); original_read.len()];
        let qname = original_read.read.qname().to_vec();
        let quals = original_read.read.qual().to_vec();
        original_read.read.set(
            &qname[..],
            Some(&original_read.read.cigar().take()),
            &new_bases,
            &quals[..],
        );
    }

    fn apply_write_ns_and_q0s(original_read: &mut BirdToolRead) {
        let new_quals = vec![0; original_read.read.qual().len()];

        let new_bases = vec!['N'.as_(); original_read.len()];
        let quals = original_read.read.qual().to_vec();
        original_read.read.set(
            &quals[..],
            Some(&original_read.read.cigar().take()),
            &new_bases,
            &new_quals,
        );
    }

    // fn overwrite_from_start_to_stop(arr: &mut [u8], new_val: u8) {
    //     arr.iter_mut().for_each(|val| *val = new_val);
    // }

    /**
     * Hard clip bases from read, from start to stop in base coordinates
     * <p>
     * If start == 0, then we will clip from the front of the read, otherwise we clip
     * from the right.  If start == 0 and stop == 10, this would clip out the first
     * 10 bases of the read.
     * <p>
     * Note that this function works with reads with negative alignment starts, in order to
     * allow us to hardClip reads that have had their soft clips reverted and so might have
     * negative alignment starts
     * <p>
     * Works properly with reduced reads and insertion/deletion base qualities
     * <p>
     *
     * @param read  a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop  a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down (Could be an empty, unmapped read)
     */
    fn apply_hard_clip_bases(read: &mut BirdToolRead, start: usize, stop: usize) {
        let new_length = read.len() - (stop - start + 1);
        // If the new read is going to be empty, return an empty read now. This avoids initializing the new
        // read with invalid values below in certain cases (such as a negative alignment start).
        if new_length == 0 {
            ReadUtils::empty_read_mut(read);
        } else {
            // If the read is unmapped there is no Cigar string and neither should we create a new cigar string
            let cigar = read.read.cigar();
            let new_cigar = if read.read.is_unmapped() {
                CigarString(vec![Cigar::Match(0)])
            } else {
                CigarUtils::clip_cigar(&cigar, start as u32, stop as u32 + 1, Cigar::HardClip(0))
            };

            let copy_start = if start == 0 { stop + 1 } else { 0 };

            let name = read.read.qname().to_vec();
            let new_bases =
                read.read.seq().as_bytes()[copy_start..(new_length + copy_start)].to_vec();
            let new_quals = read.read.qual()[copy_start..(new_length + copy_start)].to_vec();

            // let mut hard_clipped_read = read.clone();

            read.update(&name[..], Some(&new_cigar), new_bases, &new_quals);

            if start == 0 && !read.read.is_unmapped() {
                read.read.set_pos(
                    read.read.pos() + CigarUtils::alignment_start_shift(&cigar, stop as i64 + 1),
                )
            }

            // BSQR shenanigans happen here but we won't worry about that
        }
    }
}
