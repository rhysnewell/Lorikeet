use clap::ArgMatches;

use crate::utils::simple_interval::SimpleInterval;

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
    pub fn trim_interval_to_contig(
        tid: usize,
        start: usize,
        stop: usize,
        contig_length: usize,
    ) -> Option<SimpleInterval> {
        assert!(
            contig_length >= 1,
            "Contig length should be at least 1 but was {}",
            contig_length
        );
        let bounded_start = start;
        let bounded_stop = std::cmp::min(contig_length, stop);

        if bounded_start > contig_length {
            return None;
        } else {
            return Some(SimpleInterval::new(tid, bounded_start, bounded_stop));
        }
    }

    pub fn parse_limiting_interval(args: &ArgMatches) -> Option<SimpleInterval> {
        if args.contains_id("limiting-interval") {
            let interval_str = args.get_one::<String>("limiting-interval").unwrap();
            let split = interval_str.split('-').collect::<Vec<&str>>();
            if split.len() == 1 {
                None
            } else {
                let mut split_iter = split.into_iter();
                let start = split_iter.next().unwrap().parse().unwrap();
                let end = split_iter.next().unwrap().parse().unwrap();
                Some(SimpleInterval::new(0, start, end))
            }
        } else {
            None
        }
    }
}
