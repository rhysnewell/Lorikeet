use rust_htslib::bam::record::{Cigar, CigarString, Record};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::hash::{Hash, Hasher};
use std::ops::Deref;

use crate::processing::lorikeet_engine::ReadType;
use crate::reads::cigar_utils::CigarUtils;
use crate::reads::read_utils::ReadUtils;
use crate::utils::simple_interval::Locatable;

/**
 * Unified read interface for use throughout the Bird Tools.
 *
 * Adapter classes implementing this interface exist for rust_htslib's {@link Record}
 *
 * All BirdToolRead methods that return mutable reference types make defensive copies, with the exception
 * of the conversion method {@link #convertToSAMRecord}.
 *
 * Note that {@link #getContig} and {@link #getStart} will not expose nominal positions assigned to unmapped
 * reads for sorting purposes -- for unmapped reads, these methods will always return {@code null} or 0,
 * respectively. To access positions assigned to unmapped reads for sorting purposes, use {@link #getAssignedContig}
 * and {@link #getAssignedStart}.
 */
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BirdToolRead {
    pub read: Record,
    pub sample_index: usize,
    pub read_type: ReadType,
    pub transient_attributes: HashMap<String, Vec<u8>>,
    pub bases: Vec<u8>,
}

impl BirdToolRead {
    pub fn new(read: Record, sample_index: usize, read_type: ReadType) -> BirdToolRead {
        let bases = read.seq().as_bytes();
        BirdToolRead {
            read,
            sample_index,
            read_type,
            transient_attributes: HashMap::new(),
            bases,
        }
    }

    pub fn update(
        &mut self,
        name: &[u8],
        cigar: Option<&CigarString>,
        bases: Vec<u8>,
        quals: &[u8],
    ) {
        self.read.set(name, cigar, &bases, quals);
        self.bases = bases;
    }

    pub fn seq(&self) -> &[u8] {
        self.bases.as_slice()
    }

    pub fn get_contig(&self) -> usize {
        self.read.tid() as usize
    }

    // pub fn get_end(&self) -> usize {
    //     self.read.reference_end() as usize
    // }
    /**
     * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped start of the read taking soft clips (but not hard clips) into account
     */
    pub fn get_soft_start(&self) -> Result<usize, <usize as TryFrom<i64>>::Error> {
        let mut start = self.get_start() as i64;
        // let new_start = start.checked_sub(self.read.cigar().leading_softclips() as usize);
        for cig in self.read.cigar().iter() {
            match cig {
                &Cigar::SoftClip(len) => {
                    start -= len as i64;
                }
                &Cigar::HardClip(_) => continue,
                _ => break,
            }
        }
        return start.try_into();
    }

    pub fn get_soft_start_i64(&self) -> i64 {
        let mut start = self.get_start() as i64;
        // let new_start = start.checked_sub(self.read.cigar().leading_softclips() as usize);
        for cig in self.read.cigar().iter() {
            match cig {
                &Cigar::SoftClip(len) => {
                    start -= len as i64;
                }
                &Cigar::HardClip(_) => continue,
                _ => break,
            }
        }
        return start;
    }

    /**
     * Calculates the reference coordinate for the end of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedEnd() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped end of the read taking soft clips (but not hard clips) into account
     */
    pub fn get_soft_end(&self) -> usize {
        let mut found_aligned_base = false;
        let mut soft_end = self.get_end();
        let cigs = self.read.cigar();

        for cig in cigs.iter().rev() {
            match cig {
                &Cigar::SoftClip(len) => soft_end += len as usize,
                &Cigar::HardClip(_) => continue,
                _ => {
                    found_aligned_base = true;
                    break;
                }
            }
        }

        if !found_aligned_base {
            soft_end = self.get_end();
        }

        return soft_end;
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
     * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read).
     * CANNOT_COMPUTE_ADAPTOR_BOUNDARY if the read is unmapped or the mate is mapped to another contig.
     */
    pub fn get_adaptor_boundary(&self) -> usize {
        return ReadUtils::get_adaptor_boundary(&self);
    }

    pub fn name(&self) -> &[u8] {
        self.read.qname()
    }

    pub fn len(&self) -> usize {
        self.read.seq_len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn set_transient_attribute(&mut self, tag: String, val: Vec<u8>) {
        self.transient_attributes.insert(tag, val);
    }

    /**
     * @return the alignment start (0-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     *
     * Invalid to call on an unmapped read.
     */
    pub fn get_unclipped_start(&self) -> usize {
        let mut unclipped_start = self.get_start();
        for cig in self.read.cigar().iter() {
            match cig {
                &Cigar::HardClip(len) | &Cigar::SoftClip(len) => unclipped_start -= len as usize,
                _ => break,
            }
        }

        return unclipped_start;
    }

    /**
     * @param alignmentEnd The end (1-based) of the alignment
     * @param cigar        The cigar containing the alignment information
     * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     * <p/>
     * Invalid to call on an unmapped read.
     * Invalid to call with cigar = null
     */
    pub fn get_unclipped_end(&self) -> usize {
        let mut unclipped_end = self.get_end();
        for cig in self.read.cigar().iter().rev() {
            match cig {
                &Cigar::HardClip(len) | &Cigar::SoftClip(len) => {
                    unclipped_end += len as usize;
                }
                _ => break,
            }
        }

        return unclipped_end;
    }
}

impl Hash for BirdToolRead {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.read.qname().hash(state);
        self.read.qual().hash(state);
        self.sample_index.hash(state);
        self.bases.hash(state);
    }
}

impl Locatable for BirdToolRead {
    fn tid(&self) -> i32 {
        self.read.tid()
    }

    fn get_start(&self) -> usize {
        self.read.pos() as usize // because it is 0-based whereas SimpleInterval should be 1-based?
    }

    fn get_end(&self) -> usize {
        self.get_start()
            + (CigarUtils::get_reference_length(self.read.cigar().deref()) as usize)
                .checked_sub(1)
                .unwrap_or(0)
        // because it is 0-based whereas SimpleInterval should be 1-based?
    }

    fn get_length_on_reference(&self) -> usize {
        CigarUtils::get_reference_length(self.read.cigar().deref()) as usize
    }
}

/**
 * Comparator for sorting Reads by coordinate.
 *
 * Uses the various other fields in a read to break ties for reads that share
 * the same location.
 *
 * Ordering is almost identical to the {@link htsjdk.samtools.SAMRecordCoordinateComparator},
 * modulo a few subtle differences in tie-breaking rules for reads that share the same
 * position. This comparator will produce an ordering consistent with coordinate ordering
 * in a bam file, including interleaving unmapped reads assigned the positions of their
 * mates with the mapped reads.
 */
impl Ord for BirdToolRead {
    fn cmp(&self, other: &Self) -> Ordering {
        let result = ReadUtils::compare_coordinates(&self, &other);

        if result != Ordering::Equal {
            return result;
        }

        if self.read.is_reverse() != other.read.is_reverse() {
            if self.read.is_reverse() {
                return Ordering::Greater;
            } else {
                return Ordering::Less;
            }
        }

        if !self.read.qname().is_empty() && !other.read.qname().is_empty() {
            let result = self.read.qname().cmp(&other.read.qname());
            if result != Ordering::Equal {
                return result;
            }
        }

        let result = self.read.flags().cmp(&other.read.flags());
        if result != Ordering::Equal {
            return result;
        }

        let result = self.read.mapq().cmp(&other.read.mapq());
        if result != Ordering::Equal {
            return result;
        }

        if self.read.is_paired() && other.read.is_paired() {
            let result = self.read.mtid().cmp(&other.read.mtid());
            if result != Ordering::Equal {
                return result;
            }

            let result = self.read.mpos().cmp(&other.read.mpos());
            if result != Ordering::Equal {
                return result;
            }
        }

        let result = self.read.seq_len().cmp(&other.read.seq_len());
        return result;
    }
}

impl PartialOrd for BirdToolRead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
