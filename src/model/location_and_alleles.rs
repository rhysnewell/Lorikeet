use model::byte_array_allele::Allele;

/**
 * This class exists to allow VariantContext objects to be compared based only on their location and set of alleles,
 * providing a more liberal equals method so that VariantContext objects can be placed into a Set
 * which retains only VCs that have non-redundant location and Allele lists.
 */
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct LocationAndAlleles<A: Allele> {
    loc: usize,
    alleles: Vec<A>,
}

impl<A: Allele> LocationAndAlleles<A> {
    pub fn new(loc: usize, alleles: Vec<A>) -> LocationAndAlleles<A> {
        Self { loc, alleles }
    }

    pub fn get_loc(&self) -> usize {
        self.loc
    }

    pub fn get_alleles(&self) -> &Vec<A> {
        &self.alleles
    }
}
