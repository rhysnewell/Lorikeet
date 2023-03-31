use petgraph::stable_graph::NodeIndex;

use crate::graphs::base_edge::BaseEdge;
use crate::graphs::base_graph::BaseGraph;
use crate::graphs::base_vertex::BaseVertex;

/**
 * Utility functions used in the graphs package
 */
pub struct GraphUtils {}

impl GraphUtils {
    /**
     * Compute the maximum shared prefix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @return the number of shared bytes common at the start of all bytes
     */
    pub fn common_maximum_prefix_length(list_of_bytes: &Vec<&[u8]>) -> usize {
        let min_length = Self::min_kmer_length(list_of_bytes);
        for i in 0..min_length {
            let b = list_of_bytes[0][i];
            for j in 1..list_of_bytes.len() {
                if b != list_of_bytes[j][i] {
                    return i;
                }
            }
        }

        return min_length;
    }

    /**
     * Compute the maximum shared suffix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @param minLength the min. length among all byte[] in listOfBytes
     * @return the number of shared bytes common at the end of all bytes
     */
    pub fn common_maximum_suffix_length(list_of_bytes: &Vec<&[u8]>, min_length: usize) -> usize {
        for suffix_len in 0..min_length {
            let b = list_of_bytes[0][list_of_bytes[0].len() - suffix_len - 1];
            for j in 1..list_of_bytes.len() {
                if b != list_of_bytes[j][list_of_bytes[j].len() - suffix_len - 1] {
                    return suffix_len;
                }
            }
        }

        return min_length;
    }

    /**
     * Get the list of kmers as byte[] from the vertices in the graph
     *
     * @param vertices a collection of vertices
     * @return a list of their kmers in order of the iterator on vertices
     */
    pub fn get_kmers<'a, V: BaseVertex, E: BaseEdge, I: IntoIterator<Item = &'a NodeIndex>>(
        vertices: I,
        graph: &'a BaseGraph<V, E>,
    ) -> Vec<&'a [u8]> {
        vertices
            .into_iter()
            .map(|v| graph.get_sequence_from_index(*v))
            .collect::<Vec<&[u8]>>()
    }

    /**
     * Get the minimum length of a collection of byte[]
     *
     * @param kmers a list of kmers whose .length min we want
     * @return the min of the kmers, if kmers is empty the result is 0
     */
    pub fn min_kmer_length(kmers: &Vec<&[u8]>) -> usize {
        return kmers.into_iter().map(|k| k.len()).min().unwrap_or(0);
    }
}
