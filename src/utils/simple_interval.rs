use std::cmp::Ordering;
// #[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
// pub struct Locus {
//     pub chrom: String,
//     pub start: u32,
//     pub end: u32,
//     pub activity: f32,
// }
//
// impl Locus {
//     pub fn new(chrom: String, start: u32, end: u32) -> Locus {
//         Locus {
//             chrom: chrom,
//             start: start,
//             end: end,
//             activity: 0.0,
//         }
//     }
// }
//
// impl fmt::Display for Locus {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
//     }
// }

/**
* Minimal immutable class representing a 1-based closed ended genomic interval
* SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
*
*@warning 0 length intervals are NOT currently allowed, but support may be added in the future
*/
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SimpleInterval {
    start: usize,
    end: usize,
    tid: usize,
}

impl SimpleInterval {
    pub const CONTIG_SEPARATOR: char = ':';
    pub const START_END_SEPARATOR: char = '-';
    pub const MAG_SEPARATOR: char = '~';
    pub const END_OF_CONTIG: char = '+';

    /**
     * Create a new immutable 1-based interval of the form [start, end]
     * @param contig the name of the contig, must not be null
     * @param start  1-based inclusive start position
     * @param end  1-based inclusive end position
     */
    pub fn new(tid: usize, start: usize, end: usize) -> SimpleInterval {
        SimpleInterval {
            start,
            end,
            tid
        }
    }

    pub fn get_contig(&self) -> usize {
        self.tid
    }

    pub fn get_start(&self) -> usize {
        self.start
    }

    pub fn get_end(&self) -> usize {
        self.end
    }

    /**
     * Determine if this is on the same contig as other
     * this must be equivalent to this.getContig().equals(other.getContig()) but may be implemented more efficiently
     *
     * @return true iff this.getContig().equals(other.getContig())
     */
    pub fn contigs_match(&self, other: &SimpleInterval) -> bool {
        self.get_contig() == other.get_contig()
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
        self.contigs_match(other) && CoordMath::overlaps(self.get_start(), self.get_end(), other_start, other.get_end() + distance)
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
        self.contigs_match(other) && CoordMath::encloses(self.get_start(), self.get_end(), other.get_start(), other.get_end())
    }
}

// The priority queue depends on `Ord`.
// Explicitly implement the trait so the queue becomes a min-heap
// instead of a max-heap.
impl Ord for SimpleInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        other.end.cmp(&self.end).then_with(||self.start.cmp(&other.start))
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
        return (end - start) + 1
    }

    pub fn get_start(end: usize, length: usize) -> usize {
        return (end - length) + 1
    }

    pub fn get_end(start: usize, length: usize) -> usize {
        return (start + length) - 1
    }

    /**
     * Checks to see if the two sets of coordinates have any overlap.
     */
    pub fn overlaps(start: usize, end: usize, start2: usize, end2: usize) -> bool {
        (start2 >= start && start2 <= end) || (end2 >= start && end2 <= end)
        || CoordMath::encloses(start2, end2, start, end)
    }

    /** Returns true if the "inner" coords and totally enclosed by the "outer" coords. */
    pub fn encloses(outer_start: usize, outer_end: usize, inner_start: usize, inner_end: usize) -> bool {
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
    pub fn get_cycle(is_negative_strand: bool, read_length: usize, read_base_index: usize) -> usize {
        if is_negative_strand {
            read_length - read_base_index
        } else {
            read_base_index + 1
        }
    }
}