use std::cmp::{max, min, Ordering};
use std::fmt::Debug;
use std::hash::Hash;
use utils::errors::BirdToolError;
use utils::interval_utils::IntervalUtils;

/**
* Minimal immutable class representing a 0-based closed ended genomic interval
* SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
*
*@warning 0 length intervals are NOT currently allowed, but support may be added in the future
*/
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct SimpleInterval {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) tid: usize,
}

impl SimpleInterval {
    pub const CONTIG_SEPARATOR: char = ':';
    pub const START_END_SEPARATOR: char = '-';
    pub const MAG_SEPARATOR: char = '~';
    pub const END_OF_CONTIG: char = '+';

    /**
     * Create a new immutable 0-based interval of the form [start, end]
     * @param contig the name of the contig, must not be null
     * @param start  0-based inclusive start position
     * @param end  0-based inclusive end position
     */
    pub fn new(tid: usize, start: usize, end: usize) -> SimpleInterval {
        SimpleInterval { start, end, tid }
    }

    pub fn get_contig(&self) -> usize {
        self.tid
    }

    // pub fn get_start(&self) -> usize {
    //     self.start
    // }
    //
    // pub fn get_end(&self) -> usize {
    //     self.end
    // }

    /**
     * Determine if this is on the same contig as other
     * this must be equivalent to this.getContig().equals(other.getContig()) but may be implemented more efficiently
     *
     * @return true iff this.getContig().equals(other.getContig())
     */
    pub fn contigs_match<L: Locatable>(&self, other: &L) -> bool {
        self.get_contig() == other.tid() as usize
    }

    /**
     * @return number of bases covered by this interval (will always be > 0)
     */
    pub fn size(&self) -> usize {
        self.end - self.start + 1
    }

    /**
     * Determines whether this interval comes within {@code distance} of overlapping the provided locatable.
     * When distance = 0 this is equal to {@link #overlaps(Locatable)}
     *
     * @param other interval to check
     * @param distance how many bases may be between the two intervals for us to still consider them overlapping.
     * @return true if this interval overlaps other, otherwise false
     */
    pub fn within_distance_of(&self, other: &SimpleInterval, distance: usize) -> bool {
        // Avoid usize going below 0
        let other_start = if distance > other.get_start() {
            0
        } else {
            other.get_start() - distance
        };
        self.contigs_match(other)
            && CoordMath::overlaps(
                self.get_start(),
                self.get_end(),
                other_start,
                other.get_end() + distance,
            )
    }

    /**
     * Determines whether this interval contains the entire region represented by other
     * (in other words, whether it covers it).
     *
     *
     * @param other interval to check
     * @return true if this interval contains all of the base positions spanned by other, otherwise false
     */
    pub fn contains(&self, other: &SimpleInterval) -> bool {
        self.contigs_match(other)
            && CoordMath::encloses(
                self.get_start(),
                self.get_end(),
                other.get_start(),
                other.get_end(),
            )
    }

    /**
     * Returns a new SimpleInterval that represents the region between the endpoints of this and other.
     *
     * Unlike {@link #mergeWithContiguous}, the two intervals do not need to be contiguous
     *
     * @param other the other interval with which to calculate the span
     * @return a new SimpleInterval that represents the region between the endpoints of this and other.
     */
    pub fn span_with(&self, other: &Self) -> SimpleInterval {
        if !self.contigs_match(other) {
            panic!("Cannot get span for intervals on different contigs")
        } else {
            return SimpleInterval::new(
                self.tid,
                min(self.start, other.start),
                max(self.end, other.end),
            );
        }
    }

    /**
     * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
     * directions, bounded by the contig start/stop if necessary.
     *
     * @param padding amount to expand this interval
     * @param contigLength length of this interval's contig
     * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
     *         directions, bounded by the contig start/stop if necessary.
     */
    pub fn expand_within_contig(&self, padding: usize, contig_length: usize) -> SimpleInterval {
        let mut start = self.get_start();
        if start < padding {
            start = 0
        } else {
            start = start - padding
        }

        IntervalUtils::trim_interval_to_contig(self.tid, start, self.end + padding, contig_length)
            .unwrap()
    }

    /**
     * Returns the intersection of the two intervals. The intervals must overlap or IllegalArgumentException will be thrown.
     */
    pub fn intersect(&self, that: &Self) -> SimpleInterval {
        assert!(
            self.overlaps(that),
            "The two intervals need to overlap {:?} and {:?}",
            &self,
            that
        );

        return SimpleInterval::new(
            self.get_contig(),
            max(self.start, that.start),
            min(self.end, that.end),
        );
    }

    /**
     * Returns a new SimpleInterval that represents the entire span of this and that.  Requires that
     * this and that SimpleInterval are contiguous.
     */
    pub fn merge_with_contiguous(&self, that: &Self) -> Result<SimpleInterval, BirdToolError> {
        if !self.contiguous(that) {
            return Err(BirdToolError::NonContiguousIntervals(format!(
                "The two intervals need to be contiguous: {:?} {:?}",
                &self, that
            )));
        };

        return Ok(SimpleInterval::new(
            self.get_contig(),
            min(self.get_start(), that.get_start()),
            max(self.get_end(), that.get_end()),
        ));
    }

    fn contiguous(&self, that: &Self) -> bool {
        return self.get_contig() == that.get_contig()
            && self.get_start() <= that.get_end() + 1
            && that.get_start() <= self.get_end() + 1;
    }

    /**
     * Determines whether this interval comes within "margin" of overlapping the provided locatable.
     * This is the same as plain overlaps if margin=0.
     *
     * @param other interval to check
     * @param margin how many bases may be between the two interval for us to still consider them overlapping; must be non-negative
     * @return true if this interval overlaps other, otherwise false
     * @throws IllegalArgumentException if margin is negative
     */
    pub fn overlaps_with_margin<L: Locatable>(&self, other: &L, margin: usize) -> bool {
        self.contigs_match(other)
            && self.start <= other.get_end() + margin
            && other.get_start().checked_sub(margin).unwrap_or(0) <= self.end
    }
}

// The priority queue depends on `Ord`.
// Explicitly implement the trait so the queue becomes a min-heap
// instead of a max-heap.
impl Ord for SimpleInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.tid
            .cmp(&other.tid)
            .then_with(|| other.end.cmp(&self.end))
            .then_with(|| self.start.cmp(&other.start))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for SimpleInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub struct CoordMath {}

impl CoordMath {
    pub fn get_length(start: usize, end: usize) -> usize {
        return ((end + 1).checked_sub(start).unwrap_or(0));
    }

    pub fn get_start(end: usize, length: usize) -> usize {
        return (end).checked_sub(length).unwrap_or(0);
    }

    pub fn get_end(start: usize, length: usize) -> usize {
        return (start + length) - 1;
    }

    /**
     * Checks to see if the two sets of coordinates have any overlap.
     */
    pub fn overlaps(start: usize, end: usize, start2: usize, end2: usize) -> bool {
        (start2 >= start && start2 <= end)
            || (end2 >= start && end2 <= end)
            || CoordMath::encloses(start2, end2, start, end)
    }

    /** Returns true if the "inner" coords and totally enclosed by the "outer" coords. */
    pub fn encloses(
        outer_start: usize,
        outer_end: usize,
        inner_start: usize,
        inner_end: usize,
    ) -> bool {
        inner_start >= outer_start && inner_end <= outer_end
    }

    /**
     * Determines the amount of overlap between two coordinate ranges. Assumes that the two ranges
     * actually do overlap and therefore may produce strange results when they do not!
     */
    pub fn get_overlap(start: usize, end: usize, start2: usize, end2: usize) -> usize {
        CoordMath::get_length(std::cmp::max(start, start2), std::cmp::min(end, end2))
    }

    /**
     * Determines the read cycle number for the base
     *
     *  @param isNegativeStrand true if the read is negative strand
     *  @param readLength
     *  @param readBaseIndex the 0-based index of the read base in question
     */
    pub fn get_cycle(
        is_negative_strand: bool,
        read_length: usize,
        read_base_index: usize,
    ) -> usize {
        if is_negative_strand {
            read_length - read_base_index
        } else {
            read_base_index + 1
        }
    }
}

pub trait Locatable:
    Clone + Debug + Hash + Eq + PartialEq + Ord + PartialOrd + Send + Sync
{
    fn tid(&self) -> i32;

    fn get_start(&self) -> usize;

    fn get_end(&self) -> usize;

    fn overlaps<L: Locatable>(&self, other: &L) -> bool {
        (other.get_start() >= self.get_start() && other.get_start() <= self.get_end())
            || (other.get_end() >= self.get_start() && other.get_end() <= self.get_end())
            || CoordMath::encloses(
                other.get_start(),
                other.get_end(),
                self.get_start(),
                self.get_end(),
            )
    }

    fn contains<L: Locatable>(&self, other: &L) -> bool {
        return self.tid() == other.tid()
            && self.get_start() <= other.get_start()
            && self.get_end() >= other.get_end();
    }

    fn distance<L: Locatable>(&self, other: &L) -> usize {
        if self.tid() == other.tid() {
            max(self.get_start(), other.get_start()) - min(self.get_start(), other.get_start())
        } else {
            usize::MAX
        }
    }

    fn get_length_on_reference(&self) -> usize {
        return CoordMath::get_length(self.get_start(), self.get_end());
    }
}

impl Locatable for SimpleInterval {
    fn tid(&self) -> i32 {
        self.tid as i32
    }

    fn get_start(&self) -> usize {
        self.start
    }

    fn get_end(&self) -> usize {
        self.end
    }

    fn overlaps<L: Locatable>(&self, other: &L) -> bool {
        self.overlaps_with_margin(other, 0)
    }
}
