use utils::simple_interval::{SimpleInterval, CoordMath};

pub struct IntervalUtils {}

impl IntervalUtils {

    /**
     * Create a new interval, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @param contigLength length of the contig
     * @return a valid interval over contig, or null if a meaningful interval cannot be created
     */
    pub fn trim_interval_to_contig(tid: usize, start: usize, stop: usize, contig_length: usize) -> Option<SimpleInterval> {
        assert!(contig_length >= 1, "Contig length should be at least 1 but was {}", contig_length);
        let bounded_start = std::cmp::max(1, start);
        let bounded_stop = std::cmp::min(contig_length, stop);

        if (bounded_start > contig_length) || (bounded_stop < 1) {
            return None
        } else {
            return Some(SimpleInterval::new(tid, bounded_start, bounded_stop))
        }
    }
}