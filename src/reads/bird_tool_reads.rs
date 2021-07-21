use rust_htslib::bam::record::Record;
use rust_htslib::bam::ext::BamRecordExtensions;
use utils::simple_interval::Locatable;
use std::cmp::Ordering;
use reads::read_utils::ReadUtils;
use std::hash::{Hasher, Hash};


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
}

impl BirdToolRead {
    pub fn new(read: Record, sample_index: usize) -> BirdToolRead {
        BirdToolRead {
            read,
            sample_index
        }
    }

    // pub fn get_start(&self) -> usize {
    //     self.read.reference_start() as usize
    // }

    pub fn get_contig(&self) -> usize {
        self.read.tid() as usize
    }

    // pub fn get_end(&self) -> usize {
    //     self.read.reference_end() as usize
    // }

    pub fn get_soft_start(&self) -> Option<usize> {
        let mut start = self.get_start();
        let new_start = start.checked_sub(self.read.cigar().leading_softclips() as usize);
        return new_start
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
        return ReadUtils::get_adaptor_boundary(&self)
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

}

impl Hash for BirdToolRead {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.read.qname().hash(state);
        self.read.qual().hash(state);
        self.sample_index.hash(state);
        self.read.seq().encoded.hash(state);
    }
}

impl Locatable for BirdToolRead {
    fn tid(&self) -> i32 {
        self.read.tid()
    }

    fn get_start(&self) -> usize {
        self.read.reference_start() as usize
    }

    fn get_end(&self) -> usize {
        self.read.reference_end() as usize
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
            return result
        }

        if self.read.is_reverse() != other.read.is_reverse() {
            if self.read.is_reverse() { return Ordering::Greater } else { return Ordering::Less }
        }

        if !self.read.qname().is_empty() && !other.read.qname().is_empty() {
            let result = self.read.qname().cmp(&other.read.qname());
            if result != Ordering::Equal {
                return result
            }
        }

        let result = self.read.flags().cmp(&other.read.flags());
        if result != Ordering::Equal {
            return result
        }

        let result = self.read.mapq().cmp(&other.read.mapq());
        if result != Ordering::Equal {
            return result
        }

        if self.read.is_paired() && other.read.is_paired() {
            let result = self.read.mtid().cmp(&other.read.mtid());
            if result != Ordering::Equal {
                return result
            }

            let result = self.read.mpos().cmp(&other.read.mpos());
            if result != Ordering::Equal {
                return result
            }
        }

        let result = self.read.seq_len().cmp(&other.read.seq_len());
        return result
    }
}

impl PartialOrd for BirdToolRead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}