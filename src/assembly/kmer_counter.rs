use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::assembly::kmer::Kmer;

/**
 * generic utility class that counts kmers
 *
 * Basically you add kmers to the counter, and it tells you how many occurrences of each kmer it's seen.
 */
pub struct KmerCounter {
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    counts_by_kmer: HashMap<Kmer, CountedKmer>,
    kmer_length: usize,
}

impl KmerCounter {
    /**
     * Create a new kmer counter
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     */
    pub fn new(kmer_length: usize) -> KmerCounter {
        Self {
            counts_by_kmer: HashMap::new(),
            kmer_length,
        }
    }

    /**
     * Get the count of kmer in this kmer counter
     * @param kmer a non-null counter to get
     * @return a positive integer
     */
    pub fn get_kmer_count(&self, kmer: &Kmer) -> usize {
        match self.counts_by_kmer.get(kmer) {
            Some(count) => count.count,
            None => 0,
        }
    }

    pub fn get_counted_kmers(&self) -> Vec<&CountedKmer> {
        self.counts_by_kmer.values().collect::<Vec<&CountedKmer>>()
    }

    /**
     * Get kmers that have minCount or greater in this counter
     * @param minCount only return kmers with count >= this value
     * @return a non-null collection of kmers
     */
    pub fn get_kmers_with_counts_at_least(&self, min_count: usize) -> Vec<&CountedKmer> {
        return self
            .get_counted_kmers()
            .into_par_iter()
            .filter(|counted_kmer| counted_kmer.count >= min_count)
            .collect::<Vec<&CountedKmer>>();
    }

    /**
     * Remove all current counts, resetting the counter to an empty state
     */
    pub fn clear(&mut self) {
        self.counts_by_kmer.clear()
    }

    /**
     * Add a kmer that occurred kmerCount times
     *
     * @param kmer a kmer
     * @param kmerCount the number of occurrences
     */
    pub fn add_kmer(&mut self, kmer: Kmer, kmer_count: usize) {
        assert!(
            kmer.len() == self.kmer_length,
            "Incorrect kmer length, {} expected {}",
            kmer.len(),
            self.kmer_length
        );
        let mut counts_from_map = self
            .counts_by_kmer
            .entry(kmer.clone())
            .or_insert_with(|| CountedKmer::new(kmer));
        counts_from_map.count += kmer_count;
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct CountedKmer {
    kmer: Kmer,
    count: usize,
}

impl CountedKmer {
    pub fn new(kmer: Kmer) -> CountedKmer {
        CountedKmer { kmer, count: 0 }
    }

    pub fn get_kmer(&self) -> &Kmer {
        &self.kmer
    }

    pub fn get_count(&self) -> usize {
        self.count
    }

    pub fn to_string(&self) -> String {
        return format!(
            "CountedKmer{{kmer=\'{}\', count={}}}",
            self.kmer.to_string(),
            self.count
        );
    }
}

impl Ord for CountedKmer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.count.cmp(&other.count)
    }
}

impl PartialOrd for CountedKmer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hash for CountedKmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kmer.hash(state)
    }
}
