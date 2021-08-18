use haplotype::haplotype::Haplotype;
use model::variant_context::VariantContext;
use std::collections::HashSet;
use utils::simple_interval::SimpleInterval;

pub struct CalledHaplotypes<'a> {
    pub(crate) calls: Vec<VariantContext<'a>>,
    pub(crate) called_haplotypes: HashSet<Haplotype<'a, SimpleInterval>>,
}

impl<'a> CalledHaplotypes<'a> {
    pub fn new(
        calls: Vec<VariantContext<'a>>,
        called_haplotypes: HashSet<Haplotype<'a, SimpleInterval>>,
    ) -> CalledHaplotypes<'a> {
        Self {
            calls,
            called_haplotypes,
        }
    }
}
