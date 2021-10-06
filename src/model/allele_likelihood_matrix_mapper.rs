use model::allele_list::Permutation;
use model::byte_array_allele::Allele;

#[derive(Debug)]
pub struct AlleleLikelihoodMatrixMapper<A: Allele> {
    pub(crate) permutation: Permutation<A>,
}

impl<A: Allele> AlleleLikelihoodMatrixMapper<A> {
    /**
     * Constructs a new mapper given an allele-list permutation.
     * @param permutation the requested permutation.
     *
     * @throws IllegalArgumentException if {@code permutation} is {@code null}.
     */
    pub fn new(permutation: Permutation<A>) -> AlleleLikelihoodMatrixMapper<A> {
        Self { permutation }
    }
}
