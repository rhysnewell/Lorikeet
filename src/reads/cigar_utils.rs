use rust_htslib::bam::record::{CigarStringView, Cigar, CigarString};

pub struct CigarUtils {}

impl CigarUtils {


    pub fn cigar_consumes_read_bases(cig: &Cigar) -> bool {
        // Consumes read bases
        match cig {
            Cigar::Match(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_)
            | Cigar::Ins(_)
            | Cigar::SoftClip(_) => true,
            _ => false
        }
    }

    pub fn cigar_consumes_reference_bases(cig: &Cigar) -> bool {
        // consumes reference bases
        match cig {
            Cigar::Match(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => true,
            _ => false
        }
    }

    pub fn cigar_is_soft_clip(cig: &Cigar) -> bool {
        match cig {
            Cigar::SoftClip(_) => true,
            _ => false,
        }
    }

    /**
     * Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin
     * @param start initial index to clip within read bases, inclusive
     * @param stop final index to clip within read bases exclusive
     * @param clippingOperator      type of clipping -- must be either hard clip or soft clip
     */
    pub fn clip_cigar(cigar: &CigarStringView, start: u32, stop: u32, clipping_operator: Cigar) -> CigarString {
        let clip_left = start == 0;

        let mut new_cigar = Vec::new();

        let mut element_start = 0;
        for element in cigar.iter() {
            match element {
                // copy hard clips
                Cigar::HardClip(len) => {
                    new_cigar.push(Cigar::HardClip(*len))
                },
                Cigar::SoftClip(len)
                | Cigar::Diff(len)
                | Cigar::Equal(len)
                | Cigar::RefSkip(len)
                | Cigar::Del(len)
                | Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::Pad(len) => {
                    let element_end = element_start + if CigarUtils::cigar_consumes_read_bases(element) { *len } else { 0 };

                    // element precedes start or follows end of clip, copy it to new cigar
                    if element_end <= start || element_start >= stop {
                        // edge case: deletions at edge of clipping are meaningless and we skip them
                        if CigarUtils::cigar_consumes_read_bases(element) ||
                            (element_start != start && element_start != stop) {
                            new_cigar.push(element.clone())
                        }
                    } else { // otherwise, some or all of the element is soft-clipped
                        let unclipped_length = if clip_left { element_end.checked_sub(stop) } else { start.checked_sub(element_start) };
                        match unclipped_length {
                            None => {
                                // Totally clipped
                                if CigarUtils::cigar_consumes_read_bases(element) {
                                    new_cigar.push(element.clone())
                                }
                            },
                            Some(unclipped_length) => {
                                let clipped_length = len.checked_sub(unclipped_length).unwrap();
                                if clip_left {
                                    new_cigar.push(CigarUtils::cigar_from_element_and_length(&clipping_operator, clipped_length));
                                    new_cigar.push(CigarUtils::cigar_from_element_and_length(element, unclipped_length));
                                } else {
                                    new_cigar.push(CigarUtils::cigar_from_element_and_length(element, unclipped_length));
                                    new_cigar.push(CigarUtils::cigar_from_element_and_length(&clipping_operator, clipped_length));
                                }
                            }
                        }
                    };
                    element_start = element_end
                }
            }
        }
        return CigarString(new_cigar)
    }

    /**
     * replace soft clips (S) with match (M) operators, normalizing the result by all the transformations of the {@link CigarBuilder} class:
     * merging consecutive identical operators and removing zero-length elements.  For example 10S10M -> 20M and 10S10M10I10I -> 20M20I.
     */
    pub fn revert_soft_clips(cigar: &CigarStringView) -> CigarString {
        let mut builder = Vec::new();
        for element in cigar.iter() {
            match element {
                Cigar::SoftClip(length) => {
                    builder.push(CigarUtils::cigar_from_element_and_length(&Cigar::Match(0), *length))
                },
                _ => {
                    builder.push(element.clone())
                }
            }
        }
        CigarString::from(builder)
    }

    /**
     * How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips
     */
    pub fn alignment_start_shift(cigar: &CigarStringView, num_clipped: i64) -> i64 {
        let ref_bases_clipped = 0;

        let element_start = 0; // this and elementEnd are indices in the read's bases
        for element in cigar.iter() {
            match element {
                // copy hard clips
                Cigar::HardClip(len) => {
                    continue
                },
                Cigar::SoftClip(len)
                | Cigar::Diff(len)
                | Cigar::Equal(len)
                | Cigar::RefSkip(len)
                | Cigar::Del(len)
                | Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::Pad(len) => {
                    let element_end = element_start + if CigarUtils::cigar_consumes_read_bases(element) { *len as i64 } else { 0 };

                    if element_end <= num_clipped { // totally within clipped span -- this includes deletions immediately following clipping
                        ref_bases_clipped += if CigarUtils::cigar_consumes_reference_bases(element) { *len as i64 } else { 0 };
                    } else if element_start < num_clipped { // clip in middle of element, which means the element necessarily consumes read bases
                        let clipped_length = num_clipped - element_start;
                        ref_bases_clipped += if CigarUtils::cigar_consumes_reference_bases(element) { clipped_length } else { 0 };
                    }
                    element_start = element_end;
                }
            }
        }
        return ref_bases_clipped
    }

    pub fn cigar_from_element_and_length(cigar: &Cigar, length: u32) -> Cigar {
        match cigar {
            Cigar::Pad(_) => {
                return Cigar::Pad(length)
            },
            Cigar::Ins(_) => {
                return Cigar::Ins(length)
            },
            Cigar::Match(_) => {
                return Cigar::Match(length)
            },
            Cigar::Del(_) => {
                return Cigar::Del(length)
            },
            Cigar::RefSkip(_) => {
                return Cigar::RefSkip(length)
            },
            Cigar::Equal(_) => {
                return Cigar::Equal(length)
            },
            Cigar::Diff(_) => {
                return Cigar::Diff(length)
            },
            Cigar::SoftClip(_) => {
                return Cigar::SoftClip(length)
            },
            Cigar::HardClip(_) => {
                return Cigar::HardClip(length)
            }
        }
    }
}