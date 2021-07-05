/**
 * Collection of read assembly using several kmerSizes.
 *
 * <p>
 *     There could be a different assembly per each kmerSize. In turn, haplotypes are result of one of those
 *     assemblies.
 * </p>
 *
 * <p>
 *     Where there is more than one possible kmerSize that generates a haplotype we consider the smaller one.
 * </p>
 *
 * @original_author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 * @author Rhys Newell; rhys.newell@hdr.qut.edu.au; Rust re-implementation
 */
pub struct AssemblyResultSet {
    assembly_result_by_kmer_size: HashMap<usize, Assemb>
}