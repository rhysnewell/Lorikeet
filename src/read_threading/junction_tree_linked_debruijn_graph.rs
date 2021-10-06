// TODO: Implement this

/**
 * Experimental version of the ReadThreadingGraph with support for threading reads to generate JunctionTrees for resolving
 * connectivity information at longer ranges.
 *
 * Note that many of the non-DeBruijn graph alterations that are made to ReadThreadingGraph are not made here:
 * - Non-Unique kmers are not duplicated by this graph
 * - Kmers are not Zipped together to form a SeqGraph
 * - The reference path is stored in its entirety rather than being calculated on the fly
 *
 * For ease of debugging, this graph supports the method {@link #printSimplifiedGraph(File, int)}} which generates a SequenceGraph and
 * adds the junction trees to the output .dot file.
 */