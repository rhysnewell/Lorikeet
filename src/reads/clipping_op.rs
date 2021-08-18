use num::traits::AsPrimitive;
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_utils::CigarUtils;
use reads::read_clipper::ClippingRepresentation;
use reads::read_utils::ReadUtils;
use rust_htslib::bam::record::{Cigar, CigarString};
use utils::simple_interval::Locatable;

pub struct ClippingOp {
    pub start: usize,
    pub stop: usize,
}

impl ClippingOp {
    pub fn new(start: usize, stop: usize) -> ClippingOp {
        ClippingOp { start, stop }
    }

    pub fn get_length(&self) -> usize {
        self.stop.checked_sub(self.start).unwrap_or(0) + 1
    }

    /**
     * Clips the bases in read according to this operation's start and stop.  Uses the clipping
     * representation used is the one provided by algorithm argument.
     *  @param algorithm    clipping algorithm to use
     * @param originalRead the read to be clipped
     */
    pub fn apply(
        &self,
        algorithm: &ClippingRepresentation,
        original_read: &BirdToolRead,
    ) -> BirdToolRead {
        match algorithm {
            // This note probably doesn't apply in Rust as this function has explicit reference of the
            // the BirdToolRead. Kept for posterity
            // important note:
            //   it's not safe to call read.getBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the BirdToolRead
            &ClippingRepresentation::WriteNs => {
                let mut read_copied = original_read.clone();
                Self::apply_write_ns(&mut read_copied, &original_read);
                return read_copied;
            }
            &ClippingRepresentation::WriteQ0s => {
                let mut read_copied = original_read.clone();
                Self::apply_write_q0s(&mut read_copied, &original_read);
                return read_copied;
            }
            &ClippingRepresentation::WriteNsQ0s => {
                let mut read_copied = original_read.clone();
                Self::apply_write_ns_and_q0s(&mut read_copied, &original_read);
                return read_copied;
            }
            &ClippingRepresentation::HardclipBases => {
                return ClippingOp::apply_hard_clip_bases(original_read, self.start, self.stop)
            }
            &ClippingRepresentation::SoftclipBases => {
                return self.apply_soft_clip_bases(original_read)
            }
            &ClippingRepresentation::RevertSoftclippedBases => {
                return self.apply_revert_soft_clipped_bases(original_read)
            }
        }
    }

    fn apply_soft_clip_bases(&self, read: &BirdToolRead) -> BirdToolRead {
        assert!(
            !read.read.is_unmapped(),
            "ReadClipper cannot soft clip unmapped reads"
        );

        if read.read.seq_len() <= 2 {
            return read.clone();
        }

        let stop = std::cmp::min(self.stop, self.start + read.read.seq_len() - 2);

        assert!(self.start <= 0 || stop == read.read.seq_len() - 1,
                "Cannot apply soft clipping operator to the middle of a read: {:?} to be clipped at {}-{}", read.read.qname(), self.start, stop);

        let old_cigar = read.read.cigar();
        let mut new_cigar = CigarUtils::clip_cigar(
            &old_cigar,
            self.start as u32,
            (stop + 1) as u32,
            Cigar::SoftClip(0),
        );

        let mut read_copied = read.clone();
        read_copied.read.set(
            read.read.qname(),
            Some(&new_cigar),
            read.read.seq().encoded,
            read.read.qual(),
        );

        let alignment_start_shift = if self.start == 0 {
            CigarUtils::alignment_start_shift(&old_cigar, (stop + 1) as i64)
        } else {
            0
        };
        let new_start = read_copied.get_start() as i64 + alignment_start_shift;

        read_copied.read.set_pos(new_start);

        return read_copied;
    }

    fn apply_revert_soft_clipped_bases(&self, read: &BirdToolRead) -> BirdToolRead {
        let original_cigar = read.read.cigar();

        if original_cigar.is_empty()
            || !(((original_cigar.leading_hardclips() + original_cigar.leading_softclips()) > 0)
                || (original_cigar.trailing_hardclips() + original_cigar.trailing_softclips() > 0))
        {
            return read.clone();
        }

        let mut unclipped = read.clone();
        let mut unclipped_cigar = CigarUtils::revert_soft_clips(&original_cigar);

        unclipped.read.set(
            read.read.qname(),
            Some(&unclipped_cigar),
            read.read.seq().encoded,
            read.read.qual(),
        );

        let new_start = read.get_soft_start();

        match new_start {
            None => {
                // if the start of the unclipped read occurs before the contig,
                // we must hard clip away the bases since we cannot represent reads with
                // negative or 0 alignment start values in the SAMRecord (e.g., 0 means unaligned)

                // We cannot set the read to temporarily have a negative start position, as our Read
                // interface will not allow it. Instead, since we know that the final start position will
                // be 0 after the hard clip operation, set it to 0 explicitly. We have to set it twice:
                // once before the hard clip (to reset the alignment stop / read length in read implementations
                // that cache these values, such as SAMRecord), and again after the hard clip.
                let new_start =
                    (read.get_start() as i64) - (read.read.cigar().leading_softclips() as i64);

                let new_end = if new_start < 0 { -new_start } else { 0 };
                unclipped.read.set_pos(0);
                unclipped = ClippingOp::apply_hard_clip_bases(&unclipped, 0, new_end as usize);

                // Reset the position to 1 again only if we didn't end up with an empty, unmapped read after hard clipping.
                // See https://github.com/broadinstitute/gatk/issues/3845
                if !unclipped.read.is_unmapped() {
                    unclipped.read.set_pos(0);
                }

                return unclipped;
            }
            Some(new_start) => {
                unclipped.read.set_pos(new_start as i64);
                return unclipped;
            }
        }
    }

    fn apply_write_q0s(read_copied: &mut BirdToolRead, original_read: &BirdToolRead) {
        let new_quals = vec![0; original_read.read.qual().len()];
        read_copied.read.set(
            original_read.read.qname(),
            Some(&original_read.read.cigar().take()),
            original_read.read.seq().encoded,
            &new_quals,
        );
    }

    fn apply_write_ns(read_copied: &mut BirdToolRead, original_read: &BirdToolRead) {
        let new_bases = vec!['N'.as_(); original_read.len()];
        read_copied.read.set(
            original_read.read.qname(),
            Some(&original_read.read.cigar().take()),
            &new_bases,
            original_read.read.qual(),
        );
    }

    fn apply_write_ns_and_q0s(read_copied: &mut BirdToolRead, original_read: &BirdToolRead) {
        let new_quals = vec![0; original_read.read.qual().len()];

        let new_bases = vec!['N'.as_(); original_read.len()];

        read_copied.read.set(
            original_read.read.qname(),
            Some(&original_read.read.cigar().take()),
            &new_bases,
            &new_quals,
        );
    }

    fn overwrite_from_start_to_stop(arr: &mut [u8], new_val: u8) {
        arr.par_iter_mut().for_each(|val| *val = new_val);
    }

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
     * Note: this method does not assume that the read is directly modifiable
     * and makes a copy of it.
     *
     * @param read  a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop  a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down (Could be an empty, unmapped read)
     */
    fn apply_hard_clip_bases(read: &BirdToolRead, start: usize, stop: usize) -> BirdToolRead {
        let new_length = read.read.seq_len() - (stop - start + 1);

        // If the new read is going to be empty, return an empty read now. This avoids initializing the new
        // read with invalid values below in certain cases (such as a negative alignment start).
        // See https://github.com/broadinstitute/gatk/issues/3466
        if new_length == 0 {
            return ReadUtils::empty_read(read);
        }

        // If the read is unmapped there is no Cigar string and neither should we create a new cigar string
        let cigar = read.read.cigar();
        let new_cigar = if read.read.is_unmapped() {
            CigarString(vec![Cigar::Match(0)])
        } else {
            CigarUtils::clip_cigar(&cigar, start as u32, (stop + 1) as u32, Cigar::HardClip(0))
        };

        let copy_start = if start == 0 { stop + 1 } else { 0 };

        let mut new_bases = &read.read.seq().encoded[copy_start..new_length];
        let mut new_quals = &read.read.qual()[copy_start..new_length];

        let mut hard_clipped_read = read.clone();

        hard_clipped_read
            .read
            .set(read.read.qname(), Some(&new_cigar), new_bases, new_quals);

        if start == 0 && !read.read.is_unmapped() {
            hard_clipped_read.read.set_pos(
                read.read.pos() + CigarUtils::alignment_start_shift(&cigar, (stop + 1) as i64),
            )
        }

        // BSQR shenanigans happen here but we won't worry about that

        return hard_clipped_read;
    }
}
