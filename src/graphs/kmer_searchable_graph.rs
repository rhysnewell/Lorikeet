use crate::assembly::kmer::Kmer;
use crate::graphs::base_edge::BaseEdge;
use crate::graphs::base_vertex::BaseVertex;

/**
 * Common interface for those graphs that implement vertex by kmer look-up.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Rhys Newell &lt;rhys.newell@hdr.qut.edu.au&gt; for Rust implementation
 */
pub trait KmerSearchableGraph<V: BaseVertex, E: BaseEdge> {
    /**
     * Returns the vertex that represents or contains the last base of a given kmer.
     * @param k the query kmer.
     *
     * @throws NullPointerException if {@code k} is {@code null}.
     * @return {@code null} if there is no such a kmer in the graph or it is not unique.
     */
    fn find_kmer(&self, k: Kmer) -> V;

    /**
     * The kmer-size of indexed kmers.
     *
     * @return greater than 0.
     */
    fn get_kmer_size(&self) -> usize;
}
