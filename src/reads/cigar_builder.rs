use rust_htslib::bam::record::{CigarString, Cigar};

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
    leading_deletion_bases_removed: i64,
    trailing_deletion_bases_removed: i64,
    trailing_deletion_bases_removed_in_make: i64,
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
            trailing_deletion_bases_removed_in_make: 0
        }
    }

    pub fn add(&mut self, element: Cigar) {
        if element.len() > 0 {
            if self.remove_deletions_at_ends && match element { Cigar::Del(_) => true, _ => false } {
                if match self.last_operator {
                    None => true,
                    Some(operator) => {
                        match operator {
                            Cigar::SoftClip(_)
                            | Cigar::HardClip(_) => {
                                true
                            },
                            _ => false,
                        }
                    },
                } {
                    self.leading_deletion_bases_removed += element.len();
                }
            }
        }
    }
}