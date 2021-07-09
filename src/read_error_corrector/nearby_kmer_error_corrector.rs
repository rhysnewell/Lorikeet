/**
 * Utility class that error-corrects reads.
 * Main idea: An error in a read will appear as a bubble in a k-mer (de Bruijn) graph and such bubble will have very low multiplicity.
 * Hence, read errors will appear as "sparse" kmers with very little support.
 * Historically, the most common approach to error-correct reads before assembly has been to first compute the kmer spectrum of the reads,
 * defined as the kmer composition of a set of reads along with the multiplicity of each kmer.
 * First-generation correctors like the Euler corrector (Pevzner 2001) mapped low frequency kmers (kmers appearing say below N times)
 * into high frequency ones that lied within a certain Hamming or edit distance.
 * This is doable, but has some drawbacks:
 * - Kmers used for error correction become tied to kmers used for graph building.
 * - Hence, large kmers (desireable for graph building because they can resolve repeats better) are a hindrance for error correction,
 * because they are seen less often.
 * - After error correction, there is no guarantee that a sequence of kmers corresponds to an "actual" read.
 *
 * An error-corrected set of reads also makes a much smoother graph without the need to resolving so many bubbles.
 *
 * Idea hence is to correct reads based on their kmer content, but in a context independent from graph building.
 * In order to do this, the following steps are taken:
 * - The k-mer spectrum of a set of reads in computed. However, we are at freedom to choose the most convenient k-mer size (typicially around
 * read length /2).
 * - We partition the set of observed k-mers into "solid" kmers which have multiplicity > M, and "insolid" ones otherwise (Pevzner 2001).
 *
 * - Main idea of the algorithm is to try to substitute a sequence of bases in a read by a sequence better supported by kmers.
 * - For each "unsolid" kmer observed in reads, we try to find a "solid" kmer within a maximum Hamming distance.
 * - If such solid kmer exists, then this unsolid kmer is "correctable", otherwise, uncorrectable.
 * - For each read, then:
 * -- Walk through  read and visit all kmers.
 * -- If kmer is solid, continue to next kmer.
 * -- If not, and if it's correctable (i.e. there exists a mapping from an unsolid kmer to a solid kmer within a given Hamming distance),
 *    add the bases and offsets corresponding to differing positions between unsolid and solid kmer to correction list.
 * -- At the end, each base in read will have a list of corrections associated with it. We can then choose to correct or not.
 *    If read has only consistent corrections, then we can correct base to common base in corrections.
 *
 *    TODO:
 *    todo Q: WHAT QUALITY TO USE??
 *    todo how do we deal with mate pairs?
 *
 *
 */
pub struct NearbyKmerErrorCorrector {
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
}