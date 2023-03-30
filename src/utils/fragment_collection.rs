use reads::bird_tool_reads::BirdToolRead;
use std::collections::HashMap;
use utils::simple_interval::Locatable;

/**
 * Represents the results of the reads -> fragment calculation.
 *
 * Contains singleton -- objects whose underlying reads do not overlap their mate pair
 * Contains overlappingPairs -- objects whose underlying reads do overlap their mate pair
 */
pub struct FragmentCollection {
    singletons: Vec<BirdToolRead>,
    overlapping_pairs: Vec<(BirdToolRead, BirdToolRead)>,
}

impl FragmentCollection {
    pub fn new(
        singletons: Vec<BirdToolRead>,
        overlapping_pairs: Vec<(BirdToolRead, BirdToolRead)>,
    ) -> Self {
        Self {
            singletons,
            overlapping_pairs,
        }
    }

    pub fn consume(self) -> (Vec<BirdToolRead>, Vec<(BirdToolRead, BirdToolRead)>) {
        (self.singletons, self.overlapping_pairs)
    }

    pub fn create(read_containing_objects: Vec<BirdToolRead>) -> Self {
        let mut singletons = Vec::with_capacity(read_containing_objects.len());
        let mut overlapping = Vec::with_capacity(read_containing_objects.len());
        let mut name_map: HashMap<Vec<u8>, BirdToolRead> =
            HashMap::with_capacity(read_containing_objects.len());

        let mut last_start = -1;
        for read in read_containing_objects {
            if (read.get_start() as i64) < last_start {
                panic!("FragmentUtils.create assumes that the incoming objects are ordered by SAMRecord \
                alignment start, but saw a read {} with alignment start {} before the previous start {}",
                       std::str::from_utf8(read.read.qname()).unwrap(), read.get_start(), last_start);
            };
            last_start = read.get_start() as i64;

            if !read.read.is_paired()
                || read.read.is_mate_unmapped()
                || read.read.mpos() == -1
                || read.read.mpos() > read.get_end() as i64
            {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                singletons.push(read);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                let read_name = read.read.qname().to_vec();
                match name_map.remove(&read_name) {
                    Some(pe1) => {
                        // assumes we have at most 2 reads per fragment
                        overlapping.push((pe1, read));
                    }
                    None => {
                        name_map.insert(read_name, read);
                    }
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if !name_map.is_empty() {
            singletons.extend(name_map.into_values());
        }

        Self::new(singletons, overlapping)
    }
}
