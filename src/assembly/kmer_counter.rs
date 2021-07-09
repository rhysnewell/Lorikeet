use std::collections::HashMap;
use assembly::kmer::Kmer;
use std::cmp::Ordering;

/**
 * generic utility class that counts kmers
 *
 * Basically you add kmers to the counter, and it tells you how many occurrences of each kmer it's seen.
 */
pub struct KmerCounter<'a> {
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    counts_by_kmer: HashMap<Kmer<'a>, CountedKmer<'a>>,
    kmer_length: usize,
}

impl<'a> KmerCounter<'a> {
    /**
     * Create a new kmer counter
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     */
    pub fn new(kmer_length: usize) -> KmerCounter<'a> {
        Self {
            counts_by_kmer: HashMap::new(),
            kmer_length
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
            None => 0
        }
    }

    pub fn get_counted_kmers(&self) -> &'a Vec<CountedKmer<'a>> {
        self.counts_by_kmer.values()
    }

    /**
     * Get kmers that have minCount or greater in this counter
     * @param minCount only return kmers with count >= this value
     * @return a non-null collection of kmers
     */
    pub fn get_kmers_with_counts_at_least(&self, min_count: usize) -> Vec<&'a Kmer<'a>> {
        let mut result = Vec::new();
        return self.get_counted_kmers().into_par_iter().filter()(|counted_kmer| {
            counted_kmer.count >= min_count
        }).collect::<Vec<&Kmer>>()
    }
}

#[derive(Debug, Eq, PartialEq, Hash)]
struct CountedKmer<'a> {
    kmer: Kmer<'a>,
    count: usize,
}

impl<'a> CountedKmer<'a> {
    pub fn new(kmer: Kmer) -> CountedKmer<'a> {
        CountedKmer {
            kmer,
            count: 0
        }
    }

    pub fn get_kmer(&self) -> &Kmer<'a> {
        &self.kmer
    }

    pub fn get_count(&self) -> usize {
        self.count
    }

    pub fn to_string(&self) -> String {
            return format!("CountedKmer{{kmer=\'{}\', count={}}}", self.kmer.to_string(), self.count)
        }
    }

impl Ord for CountedKmer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.count.cmp(other.count)
    }
}

impl PartialOrd for CountedKmer {
    fn cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}