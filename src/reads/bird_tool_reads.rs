use rust_htslib::bam::record::{Record, CigarString, Cigar, CigarStringView};
use rust_htslib::bam::ext::BamRecordExtensions;
use std::fs::soft_link;


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
#[derive(Debug, Clone)]
pub struct BirdToolRead {
    pub read: Record
}

impl BirdToolRead {
    pub fn new(read: Record) -> BirdToolRead {
        BirdToolRead {
            read
        }
    }

    pub fn get_start(&self) -> usize {
        self.read.reference_start() as usize
    }

    pub fn get_contig(&self) -> usize {
        self.read.tid() as usize
    }

    pub fn get_end(&self) -> usize {
        self.read.reference_end() as usize
    }

    pub fn get_soft_start(&self) -> Option<usize> {
        let mut start = self.get_start();
        start.checked_sub(self.read.cigar().leading_softclips());
        return start
    }
}

