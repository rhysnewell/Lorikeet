use anyhow::Result;
use rust_htslib::bam::record::{Cigar, CigarString};
use std::convert::TryFrom;

use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::cigar_builder::CigarBuilder;
use crate::reads::cigar_utils::CigarUtils;
use crate::utils::artificial_read_utils::ArtificialReadUtils;

lazy_static! {
    pub static ref BASES: Vec<u8> = vec!['A' as u8, 'C' as u8, 'T' as u8, 'G' as u8];
    pub static ref QUALS: Vec<u8> = vec![2, 15, 25, 30];
    pub static ref LEADING_CLIPS: Vec<Vec<Cigar>> = vec![
        Vec::new(),
        vec![Cigar::HardClip(1)],
        vec![Cigar::SoftClip(1)],
        vec![Cigar::HardClip(1), Cigar::SoftClip(1)]
    ];
    pub static ref TRAILING_CLIPS: Vec<Vec<Cigar>> = LEADING_CLIPS
        .clone()
        .into_iter()
        .rev()
        .collect::<Vec<Vec<Cigar>>>();
    pub static ref CORE_CIGAR_ELEMENTS: Vec<Cigar> =
        vec![Cigar::Ins(1), Cigar::Del(1), Cigar::Match(1)];
    pub static ref CORE_CIGAR_ELEMENTS_INCLUDING_SKIPS: Vec<Cigar> = vec![
        Cigar::Ins(1),
        Cigar::Del(1),
        Cigar::Match(1),
        Cigar::RefSkip(1)
    ];
}

pub struct ReadClipperTestUtils {}

impl ReadClipperTestUtils {
    //Should contain all the utils needed for tests to mass produce
    //reads, cigars, and other needed classes
    pub fn make_read_from_str(cigar_string: &str, length_change: i64) -> BirdToolRead {
        let cigar = CigarString::try_from(cigar_string).unwrap();
        let mut read_length = CigarUtils::get_read_length(&cigar) as i64;
        if read_length >= -length_change {
            read_length += length_change;
        }

        ArtificialReadUtils::create_artificial_read(
            &Self::array_from_array_with_length(&*BASES, read_length as usize),
            &Self::array_from_array_with_length(&*QUALS, read_length as usize),
            cigar,
        )
    }

    pub fn make_read_from_cigar(cigar_string: CigarString, length_change: i64) -> BirdToolRead {
        let mut read_length = CigarUtils::get_read_length(&cigar_string) as i64;
        if read_length >= -length_change {
            read_length += length_change;
        }

        ArtificialReadUtils::create_artificial_read(
            &Self::array_from_array_with_length(&*BASES, read_length as usize),
            &Self::array_from_array_with_length(&*QUALS, read_length as usize),
            cigar_string,
        )
    }

    pub fn array_from_array_with_length(array: &[u8], length: usize) -> Vec<u8> {
        let mut output = vec![0; length];
        for j in 0..length {
            output[j] = array[j % array.len()]
        }

        return output;
    }

    /**
     * This function generates every valid permutation of cigar strings (with a given set of cigarElement) with a given length.
     * <p>
     * A valid cigar object obeys the following rules:
     * - No Hard/Soft clips in the middle of the read
     * - No deletions in the beginning / end of the read
     * - No repeated adjacent element (e.g. 1M2M, this should be 3M)
     * - No consecutive I/D elements
     *
     * @param maximumCigarElements the maximum number of elements in the cigar
     * @return a list with all valid Cigar objects
     */
    pub fn generate_cigar_list(
        maximum_cigar_elements: usize,
        include_skips: bool,
    ) -> Result<Vec<CigarString>> {
        let core_elements = if include_skips {
            CORE_CIGAR_ELEMENTS_INCLUDING_SKIPS.to_vec()
        } else {
            CORE_CIGAR_ELEMENTS.to_vec()
        };

        let mut cigar_list = Vec::new();
        for leading_clip_elements in LEADING_CLIPS.iter() {
            for trailing_clip_elements in TRAILING_CLIPS.iter() {
                let maximum_core_elements = maximum_cigar_elements as i64
                    - leading_clip_elements.len() as i64
                    - trailing_clip_elements.len() as i64;
                if maximum_core_elements < 1 {
                    continue;
                };

                let num_core_elements = core_elements.len();
                let mut cigar_combination = vec![0; maximum_core_elements as usize];

                let mut current_index = 0;
                loop {
                    let core_cigar = Self::create_cigar_from_combination(
                        cigar_combination.clone(),
                        core_elements.clone(),
                    ); // create the cigar
                    if CigarUtils::is_good(&core_cigar) {
                        let mut builder = CigarBuilder::new(true);
                        // Some of the tests create CIGAR elements that cause errors here. It is intentional, but we
                        // don't want to return the error. Under normal circumstances, this should never happen and
                        // the error will be caught and lorikeet will panic.
                        match builder.add_all(leading_clip_elements.clone()) {
                            Ok(_) => {}
                            Err(_) => {}
                        };
                        match builder.add_all(core_cigar.0) {
                            Ok(_) => {}
                            Err(_) => {}
                        };
                        match builder.add_all(trailing_clip_elements.clone()) {
                            Ok(_) => {}
                            Err(_) => {}
                        };
                        match builder.make(false) {
                            Ok(cigar) => cigar_list.push(cigar),
                            Err(_) => {} // do nothing
                        }
                    }

                    let mut current_index_changed = false;
                    while current_index < maximum_core_elements as usize
                        && cigar_combination[current_index] == num_core_elements - 1
                    {
                        current_index += 1; // find the next index to increment
                        current_index_changed = true; // keep track of the fact that we have changed indices!
                    }

                    if current_index == maximum_core_elements as usize {
                        // if we hit the end of the array, we're done.
                        break;
                    };

                    cigar_combination[current_index] += 1; // otherwise advance the current index
                    if current_index_changed {
                        // if we have changed index, then...
                        for i in 0..current_index {
                            cigar_combination[i] = 0; // reset everything from 0->currentIndex
                        }
                        current_index = 0; // go back to the first index
                    }
                }
            }
        }

        return Ok(cigar_list);
    }

    fn create_cigar_from_combination(
        cigar_combination: Vec<usize>,
        cigar_elements: Vec<Cigar>,
    ) -> CigarString {
        CigarString::from(
            cigar_combination
                .into_iter()
                .map(|n| cigar_elements[n].clone())
                .collect::<Vec<Cigar>>(),
        )
    }
}
