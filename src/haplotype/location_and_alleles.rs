use model::byte_array_allele::Allele;

/**
 * This class exists to allow VariantContext objects to be compared based only on their location and set of alleles,
 * providing a more liberal equals method so that VariantContext objects can be placed into a Set
 * which retains only VCs that have non-redundant location and Allele lists.
 */
#[derive(Debug, Clone, Eq, Ord, PartialOrd, PartialEq, Hash)]
pub struct LocationsAndAlleles<'a, A: Allele> {
    loc: usize,
    alleles: &'a Vec<A>,
}

impl<'a, A: Allele> LocationsAndAlleles<'a, A> {
    pub fn new(loc: usize, alleles: &'a Vec<A>) -> LocationsAndAlleles<'a, A> {
        Self { loc, alleles }
    }

    pub fn get_loc(&self) -> usize {
        self.loc
    }

    pub fn get_alleles(&self) -> &Vec<A> {
        &self.alleles
    }
}
