
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
#[derive(Clone, Debug)]
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

    /**
     * @return number of bases covered by this interval (will always be > 0)
     */
    pub fn size(&self) -> usize {
        self.end - self.start + 1
    }
}