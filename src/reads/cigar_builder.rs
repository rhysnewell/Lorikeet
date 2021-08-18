use reads::cigar_utils::CigarUtils;
use rust_htslib::bam::record::{Cigar, CigarString};

#[derive(Debug, Eq, PartialEq)]
enum Section {
    LeftHardClip,
    LeftSoftClip,
    Middle,
    RightSoftClip,
    RightHardClip,
}

/**
 * This class allows code that manipulates cigars to do so naively by handling complications such as merging consecutive
 * identical operators within the builder.  A CigarBuilder takes care of the following:
 *
 * 1)  Merging consecutive identical operators, eg 10M5M -> 15M
 * 2)  Eliminating leading and trailing deletions, eg 10D10M -> 10M and 10M10D -> 10M
 * 3)  Shifting deletions to the left of adjacent insertions, eg 10M1ID10D -> 10M10D10I
 * 4)  Validating the overall structure of [hard clip] [soft clip] non-clip [soft clip] [hard clip]
 *
 * Edge cases, such as removing a deletion that immediately follows a leading insertion, *are* handled correctly.  See the unit tests.
 *
 * Leading and trailing deletions may be kept by using the non-default CigarBuilder(false) constructor.
 *
 * All of this is achieved simply by invoking add() repeatedly, followed by make().
 */
pub struct CigarBuilder {
    cigar_elements: Vec<Cigar>,
    // track the last operator so we can merge consecutive elements with the same operator
    // for example, adding 3M and 4M is equivalent to adding 7M
    // also we ignore leading deletions so for example 10S + 5D = 10S
    last_operator: Option<Cigar>,
    section: Section,
    remove_deletions_at_ends: bool,
    leading_deletion_bases_removed: u32,
    trailing_deletion_bases_removed: u32,
    trailing_deletion_bases_removed_in_make: u32,
}

impl CigarBuilder {
    pub fn new(remove_deletions_at_ends: bool) -> Self {
        Self {
            remove_deletions_at_ends,
            cigar_elements: Vec::new(),
            last_operator: None,
            section: Section::LeftHardClip,
            leading_deletion_bases_removed: 0,
            trailing_deletion_bases_removed: 0,
            trailing_deletion_bases_removed_in_make: 0,
        }
    }

    pub fn add(&mut self, element: Cigar) {
        if element.len() > 0 {
            if self.remove_deletions_at_ends
                && match element {
                    Cigar::Del(_) => true,
                    _ => false,
                }
            {
                if match self.last_operator {
                    None => true,
                    Some(operator) => match operator {
                        Cigar::SoftClip(_) | Cigar::HardClip(_) => true,
                        Cigar::Ins(length) => {
                            if self.cigar_elements.len() == 1
                                || CigarUtils::is_clipping(
                                    &self.cigar_elements[self.cigar_elements.len() - 2],
                                )
                            {
                                true
                            } else {
                                false
                            }
                        }
                        _ => false,
                    },
                } {
                    self.leading_deletion_bases_removed += element.len();
                };
            }

            self.advance_section_and_validate_cigar_order(&element);

            if CigarUtils::cigar_elements_are_same_type(&element, &self.last_operator) {
                let n = self.cigar_elements.len() - 1;
                self.cigar_elements[n] =
                    CigarUtils::combine_cigar_operators(&element, &self.cigar_elements[n])
                        .unwrap_or(self.cigar_elements[n]);
            } else {
                match self.last_operator {
                    None => {
                        self.cigar_elements.push(element.clone());
                        self.last_operator = Some(element);
                    }
                    Some(_) => {
                        if CigarUtils::is_clipping(&element) {
                            // if we have just started clipping on the right and realize the last operator was a deletion, remove it
                            // if we have just started clipping on the right and the last two operators were a deletion and insertion, remove the deletion
                            let cigar_elements_len = self.cigar_elements.len();
                            if self.remove_deletions_at_ends
                                && !CigarUtils::cigar_consumes_read_bases(
                                    &self.last_operator.unwrap(),
                                )
                                && !CigarUtils::is_clipping(&self.last_operator.unwrap())
                            {
                                self.trailing_deletion_bases_removed +=
                                    self.cigar_elements[cigar_elements_len - 1].len();
                                self.cigar_elements[cigar_elements_len - 1] = element.clone();
                                self.last_operator = Some(element);
                            } else if self.remove_deletions_at_ends
                                && self.last_two_elements_were_deletion_and_insertion()
                            {
                                self.trailing_deletion_bases_removed +=
                                    self.cigar_elements[cigar_elements_len - 2].len();
                                self.cigar_elements[cigar_elements_len - 2] =
                                    self.cigar_elements[cigar_elements_len - 1].clone();
                                self.cigar_elements[cigar_elements_len - 1] = element;
                                // self.last_operator = Some(element);
                            } else {
                                self.cigar_elements.push(element.clone());
                                self.last_operator = Some(element);
                            }
                        } else {
                            match element {
                                Cigar::Del(_) => {
                                    match self.last_operator {
                                        None => {
                                            self.cigar_elements.push(element.clone());
                                            self.last_operator = Some(element);
                                        }
                                        Some(last_operator) => {
                                            match last_operator {
                                                Cigar::Ins(_) => {
                                                    // The order of deletion and insertion elements is arbitrary, so to standardize we shift deletions to the left
                                                    // that is, we place the deletion before the insertion and shift the insertion right
                                                    // if the element before the insertion is another deletion, we merge in the new deletion
                                                    // note that the last operator remains an insertion
                                                    let size = self.cigar_elements.len();
                                                    if size > 1
                                                        && match self.cigar_elements[size - 2] {
                                                            Cigar::Del(_) => true,
                                                            _ => false,
                                                        }
                                                    {
                                                        self.cigar_elements[size - 2] = Cigar::Del(
                                                            self.cigar_elements[size - 2].len()
                                                                + element.len(),
                                                        );
                                                        // self.last_operator = Some(element);
                                                    } else {
                                                        self.cigar_elements
                                                            .insert(size - 1, element);
                                                        // self.last_operator = Some(element);
                                                    }
                                                }
                                                _ => {
                                                    self.cigar_elements.push(element.clone());
                                                    self.last_operator = Some(element);
                                                }
                                            }
                                        }
                                    }
                                }
                                _ => {
                                    self.cigar_elements.push(element.clone());
                                    self.last_operator = Some(element);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn add_all(&mut self, elements: Vec<Cigar>) {
        for element in elements {
            self.add(element)
        }
    }

    fn last_two_elements_were_deletion_and_insertion(&self) -> bool {
        match self.last_operator {
            None => false,
            Some(operator) => {
                if self.cigar_elements.len() > 1 {
                    match operator {
                        Cigar::Ins(_) => match self.cigar_elements[self.cigar_elements.len() - 2] {
                            Cigar::Del(_) => true,
                            _ => false,
                        },
                        _ => false,
                    }
                } else {
                    false
                }
            }
        }
    }

    // validate that cigar structure is hard clip, soft clip, unclipped, soft clip, hard clip
    fn advance_section_and_validate_cigar_order(&mut self, operator: &Cigar) {
        match operator {
            Cigar::HardClip(_) => {
                match self.section {
                    Section::LeftSoftClip | Section::Middle | Section::RightSoftClip => {
                        self.section = Section::RightHardClip
                    }
                    _ => {
                        // Do nothing?
                    }
                }
            }
            Cigar::SoftClip(_) => {
                match self.section {
                    Section::RightHardClip => {
                        panic!("Cigar has already reached its right hard clip");
                    }
                    Section::LeftHardClip => self.section = Section::LeftSoftClip,
                    Section::Middle => self.section = Section::RightSoftClip,
                    _ => {
                        // do nothing
                    }
                }
            }
            _ => {
                match self.section {
                    Section::RightSoftClip | Section::RightHardClip => {
                        panic!("Section has alread reached right clip")
                    }
                    Section::LeftHardClip | Section::LeftSoftClip => self.section = Section::Middle,
                    _ => {
                        // do nothing
                    }
                }
            }
        }
    }

    pub fn make(mut self, allow_empty: bool) -> CigarString {
        assert!(
            !(self.section == Section::LeftSoftClip
                && match self.cigar_elements[0] {
                    Cigar::SoftClip(_) => true,
                    _ => false,
                }),
            "Cigar is completely soft clipped"
        );

        self.trailing_deletion_bases_removed_in_make = 0;
        if self.remove_deletions_at_ends
            && match self.last_operator {
                Some(element) => match element {
                    Cigar::Del(_) => true,
                    _ => false,
                },
                None => {
                    panic!("Last element cannot be None at this point");
                }
            }
        {
            self.trailing_deletion_bases_removed_in_make =
                self.cigar_elements[self.cigar_elements.len() - 1].len();
            self.cigar_elements.remove(self.cigar_elements.len() - 1);
        } else if self.remove_deletions_at_ends
            && self.last_two_elements_were_deletion_and_insertion()
        {
            self.trailing_deletion_bases_removed_in_make =
                self.cigar_elements[self.cigar_elements.len() - 2].len();
            self.cigar_elements.remove(self.cigar_elements.len() - 2);
        }

        assert!(
            allow_empty || !self.cigar_elements.is_empty(),
            "No cigar elements left after removing leading and trailing deletions."
        );

        return CigarString::from(self.cigar_elements);
    }

    pub fn make_and_record_deletions_removed_result(mut self) -> CigarBuilderResult {
        let leading_deletion_bases_removed = self.leading_deletion_bases_removed;
        let trailing_deletion_bases_removed = self.get_trailing_deletion_bases_removed();
        let cigar = self.make(false);
        return CigarBuilderResult::new(
            cigar,
            leading_deletion_bases_removed,
            trailing_deletion_bases_removed,
        );
    }

    /**
     * Count the number of leading deletion bases that have been removed by this builder and that will not show up in any call to make().
     * Note that all leading deletions are removed prior to calling make().  For example, successively adding 3S2D10I7D10M would result in
     * the 2D and 7D elements being discarded, for a total of 9 removed deletion bases.
     */
    pub fn get_leading_deletion_bases_removed(&self) -> u32 {
        self.leading_deletion_bases_removed
    }

    /**
     * Counts the number of trailing deletion bases that were removed in the last call to make().  These may be removed
     * before or during make().  For example, adding 3M and 3D does not removed the 3D because the builder does not know that 3D
     * is a terminal element.  If make() is then called, the builder will record the discarded 3D and this method will return 3.
     * Subsequently adding 3M, calling make(), and then calling this method will result in 0.
     */
    pub fn get_trailing_deletion_bases_removed(&self) -> u32 {
        self.trailing_deletion_bases_removed + self.trailing_deletion_bases_removed_in_make
    }
}

pub struct CigarBuilderResult {
    pub cigar: CigarString,
    pub leading_deletion_bases_removed: u32,
    pub trailing_deletion_bases_removed: u32,
}

impl CigarBuilderResult {
    pub fn new(
        cigar: CigarString,
        leading_deletion_bases_removed: u32,
        trailing_deletion_bases_removed: u32,
    ) -> Self {
        Self {
            cigar,
            leading_deletion_bases_removed,
            trailing_deletion_bases_removed,
        }
    }
}
