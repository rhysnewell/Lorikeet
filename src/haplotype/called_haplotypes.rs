use crate::model::variant_context::VariantContext;


pub struct CalledHaplotypes {
    pub(crate) calls: Vec<VariantContext>,
    // pub(crate) called_haplotypes: HashSet<Haplotype<SimpleInterval>>,
}

impl CalledHaplotypes {
    pub fn new(
        calls: Vec<VariantContext>,
        // called_haplotypes: HashSet<Haplotype<SimpleInterval>>,
    ) -> CalledHaplotypes {
        Self {
            calls,
            // called_haplotypes,
        }
    }
}
