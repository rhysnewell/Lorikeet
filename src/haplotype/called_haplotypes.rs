use haplotype::haplotype::Haplotype;
use model::variant_context::VariantContext;
use std::collections::HashSet;
use utils::simple_interval::SimpleInterval;

pub struct CalledHaplotypes {
    pub(crate) calls: Vec<VariantContext>,
    pub(crate) called_haplotypes: HashSet<Haplotype<SimpleInterval>>,
}

impl CalledHaplotypes {
    pub fn new(
        calls: Vec<VariantContext>,
        called_haplotypes: HashSet<Haplotype<SimpleInterval>>,
    ) -> CalledHaplotypes {
        Self {
            calls,
            called_haplotypes,
        }
    }
}
