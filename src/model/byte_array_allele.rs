use model::variants::Allele;

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ByteArrayAllele {
    pub(crate) is_ref: bool,
    is_no_call: bool,
    is_symbolic: bool,
    pub(crate) bases: Vec<u8>,
}

impl ByteArrayAllele {
    pub fn new(bases: &[u8], is_ref: bool) -> ByteArrayAllele {
        if Allele::would_be_null_allele(bases) {
            panic!("Null alleles are not supported")
        }

        if Allele::would_be_no_call_allele(bases) {
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

        if Allele::would_be_symbolic_allele(bases) {
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

        if !Allele::acceptable_allele_bases(bases, is_ref) {
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
}
