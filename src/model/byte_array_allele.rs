use model::variants;
use rayon::prelude::*;
use utils::vcf_constants::VCFConstants;

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ByteArrayAllele {
    pub(crate) is_ref: bool,
    pub(crate) is_no_call: bool,
    pub(crate) is_symbolic: bool,
    pub(crate) bases: Vec<u8>,
}

impl ByteArrayAllele {
    const SINGLE_BREAKEND_INDICATOR: char = '.';
    const BREAKEND_EXTENDING_RIGHT: char = '[';
    const BREAKEND_EXTENDING_LEFT: char = ']';
    const SYMBOLIC_ALLELE_START: char = '<';
    const SYMBOLIC_ALLELE_END: char = '>';

    pub const NO_CALL: char = '.';
    pub const SPAN_DEL: char = '*';

    pub fn new(bases: &[u8], is_ref: bool) -> ByteArrayAllele {
        if Self::would_be_null_allele(bases) {
            panic!("Null alleles are not supported")
        }

        if Self::would_be_no_call_allele(bases) {
            if is_ref {
                panic!("Cannot tag a no call allele as the reference allele")
            } else {
                return ByteArrayAllele {
                    bases: bases.to_ascii_uppercase(),
                    is_ref: false,
                    is_no_call: true,
                    is_symbolic: false,
                };
            }
        }

        if Self::would_be_symbolic_allele(bases) {
            if is_ref {
                panic!("Cannot tag a no call allele as the reference allele")
            } else {
                return ByteArrayAllele {
                    bases: bases.to_ascii_uppercase(),
                    is_ref: false,
                    is_no_call: false,
                    is_symbolic: true,
                };
            }
        }

        if !Self::acceptable_allele_bases(bases, is_ref) {
            panic!(
                "Unexpected base in allele bases {} ",
                String::from_utf8_lossy(bases).to_string()
            )
        } else {
            return ByteArrayAllele {
                bases: bases.to_ascii_uppercase(),
                is_ref,
                is_no_call: false,
                is_symbolic: false,
            };
        }
    }

    pub fn len(&self) -> usize {
        return if self.is_symbolic {
            0
        } else {
            self.bases.len()
        };
    }

    pub fn get_bases(&self) -> &Vec<u8> {
        return if self.is_symbolic {
            &*variants::EMPTY_ALLELE_BASES
        } else {
            &self.bases
        };
    }

    pub fn fake(is_ref: bool) -> ByteArrayAllele {
        Self::new(".".as_bytes(), is_ref)
    }

    pub fn create_fake_alleles() -> Vec<ByteArrayAllele> {
        let alleles = vec![Self::fake(true), Self::fake(false)];

        return alleles;
    }

    pub fn no_call() -> ByteArrayAllele {
        Self {
            bases: vec![Self::NO_CALL as u8],
            is_ref: false,
            is_no_call: true,
            is_symbolic: false,
        }
    }

    pub fn extend(left: &ByteArrayAllele, right: &[u8]) -> ByteArrayAllele {
        if left.is_symbolic {
            panic!("Cannot extend a symbolic allele");
        };

        let mut bases = vec![0; left.len() + right.len()];
        bases[0..left.len()].clone_from_slice(&left.get_bases()[0..left.len()]);
        bases[left.len()..].clone_from_slice(right);

        return Self::new(&bases, left.is_ref);
    }

    pub fn would_be_null_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::NULL_ALLELE
            || bases.len() == 0;
    }

    pub fn would_be_no_call_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::NO_CALL_ALLELE;
    }

    pub fn would_be_star_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::SPANNING_DELETION_ALLELE;
    }

    pub fn would_be_symbolic_allele(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false;
        } else {
            return bases[0] == Self::SYMBOLIC_ALLELE_START as u8
                || bases[bases.len() - 1] == Self::SYMBOLIC_ALLELE_END as u8
                || Self::would_be_breakpoint(bases)
                || Self::would_be_single_breakend(bases);
        }
    }

    pub fn would_be_breakpoint(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false;
        }
        return bases.iter().par_bridge().any(|base| {
            *base as char == Self::BREAKEND_EXTENDING_LEFT
                || *base as char == Self::BREAKEND_EXTENDING_RIGHT
        });
    }

    pub fn would_be_single_breakend(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false;
        } else {
            return bases[0] == Self::SINGLE_BREAKEND_INDICATOR as u8
                || bases[bases.len() - 1] == Self::SINGLE_BREAKEND_INDICATOR as u8;
        }
    }

    pub fn acceptable_allele_bases(bases: &[u8], is_ref: bool) -> bool {
        if Self::would_be_null_allele(bases) {
            return false;
        } else if Self::would_be_no_call_allele(bases) || Self::would_be_symbolic_allele(bases) {
            return true;
        } else if Self::would_be_star_allele(bases) {
            return !is_ref;
        } else {
            // return true if there are any unacceptable bases, so take conjugate value
            !bases.iter().par_bridge().any(|base| {
                let base = *base as char;
                match base {
                    'A' | 'C' | 'T' | 'G' | 'a' | 'c' | 't' | 'g' | 'N' | 'n' => false,
                    _ => true,
                }
            })
        }
    }
}

pub trait Allele: Eq + PartialEq + Clone {
    fn is_reference(&self) -> bool;

    fn length(&self) -> usize;

    fn is_symbolic(&self) -> bool;

    fn is_called(&self) -> bool;

    fn is_no_call(&self) -> bool;
}

impl Allele for ByteArrayAllele {
    fn is_reference(&self) -> bool {
        self.is_ref
    }

    fn length(&self) -> usize {
        self.len()
    }

    fn is_symbolic(&self) -> bool {
        self.is_symbolic
    }

    fn is_called(&self) -> bool {
        !self.is_no_call
    }

    fn is_no_call(&self) -> bool {
        self.is_no_call
    }
}
