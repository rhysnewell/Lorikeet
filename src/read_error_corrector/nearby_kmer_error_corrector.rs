use assembly::kmer::Kmer;
use assembly::kmer_counter::KmerCounter;
use rayon::prelude::*;
use read_error_corrector::read_error_corrector::ReadErrorCorrector;
use reads::bird_tool_reads::BirdToolRead;
use reads::read_clipper::ReadClipper;
use std::collections::HashMap;
use std::ops::Deref;
use utils::base_utils::BaseUtils;
use utils::quality_utils::QualityUtils;

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
    pub counts_by_kmer: KmerCounter,
    kmer_correction_map: HashMap<Kmer, Kmer>,
    kmer_differing_bases: HashMap<Kmer, (Vec<usize>, Vec<u8>)>,
    kmer_length: usize,
    trim_low_quality_bases: bool,
    min_tail_quality: u8,
    max_mismatches_to_correct: usize,
    quality_of_corrected_bases: u8,
    max_observations_for_kmer_to_be_correctable: usize,
    max_homopolymer_length_in_region: usize,
    min_observations_for_kmer_to_be_solid: usize,
    do_in_place_error_correction: bool,
    read_error_correction_stats: ReadErrorCorrectionStats,
}

struct ReadErrorCorrectionStats {
    num_reads_corrected: usize,
    num_reads_uncorrected: usize,
    num_bases_corrected: usize,
    num_solid_kmers: usize,
    num_uncorrectable_kmers: usize,
    num_corrected_kmers: usize,
}

impl ReadErrorCorrectionStats {
    fn new() -> Self {
        Self {
            num_reads_corrected: 0,
            num_reads_uncorrected: 0,
            num_bases_corrected: 0,
            num_solid_kmers: 0,
            num_uncorrectable_kmers: 0,
            num_corrected_kmers: 0,
        }
    }
}

impl NearbyKmerErrorCorrector {
    pub const MAX_MISMATCHES_TO_CORRECT: usize = 2;
    pub const QUALITY_OF_CORRECTED_BASES: u8 = 30;
    pub const MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE: usize = 1;
    pub const TRIM_LOW_QUAL_TAILS: bool = false;
    pub const DONT_CORRECT_IN_LONG_HOMOPOLYMERS: bool = false;
    pub const MAX_HOMOPOLYMER_THRESHOLD: usize = 12;

    /**
     * Create a new kmer corrector
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     * @param maxMismatchesToCorrect e >= 0
     * @param qualityOfCorrectedBases  Bases to be corrected will be assigned this quality
     */
    pub fn new(
        kmer_length: usize,
        max_mismatches_to_correct: usize,
        max_observations_for_kmer_to_be_correctable: usize,
        quality_of_corrected_bases: u8,
        min_observations_for_kmer_to_be_solid: usize,
        trim_low_quality_bases: bool,
        min_tail_quality: u8,
        full_reference_with_padding: &[u8],
    ) -> NearbyKmerErrorCorrector {
        assert!(
            quality_of_corrected_bases >= 2
                && quality_of_corrected_bases <= QualityUtils::MAX_REASONABLE_Q_SCORE,
            "Quality of corrected bases must be >= 2 and <= 60 but got {}",
            quality_of_corrected_bases
        );

        let counts_by_kmer = KmerCounter::new(kmer_length);
        // when region has long homopolymers, we may want not to correct reads, since assessment is complicated,
        // so we may decide to skip error correction in these regions
        let max_homopolymer_length_in_region =
            Self::compute_max_hl_len(full_reference_with_padding);

        Self {
            kmer_length,
            counts_by_kmer,
            max_mismatches_to_correct,
            quality_of_corrected_bases,
            min_observations_for_kmer_to_be_solid,
            max_homopolymer_length_in_region,
            trim_low_quality_bases,
            min_tail_quality,
            max_observations_for_kmer_to_be_correctable,
            read_error_correction_stats: ReadErrorCorrectionStats::new(),
            do_in_place_error_correction: false,
            kmer_correction_map: HashMap::new(),
            kmer_differing_bases: HashMap::new(),
        }
    }

    /**
     * Simple constructor with sensible defaults
     * @param kmerLength            K-mer length for error correction (not necessarily the same as for assembly graph)
     * @param minTailQuality        Minimum tail quality: remaining bases with Q's below this value are hard-clipped after correction
     * @param debug                 Output debug information
     */
    pub fn default(
        kmer_length: usize,
        min_tail_quality: u8,
        min_observations_for_kmer_to_be_solid: usize,
        full_reference_with_padding: &[u8],
    ) -> NearbyKmerErrorCorrector {
        Self::new(
            kmer_length,
            Self::MAX_MISMATCHES_TO_CORRECT,
            Self::MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE,
            Self::QUALITY_OF_CORRECTED_BASES,
            min_observations_for_kmer_to_be_solid,
            Self::TRIM_LOW_QUAL_TAILS,
            min_tail_quality,
            full_reference_with_padding,
        )
    }

    /**
     * Main entry routine to add all kmers in a read to the read map counter
     * @param read                        Read to add bases
     */
    pub fn add_read_kmers(&mut self, read: &BirdToolRead) {
        if Self::DONT_CORRECT_IN_LONG_HOMOPOLYMERS
            && self.max_homopolymer_length_in_region > Self::MAX_HOMOPOLYMER_THRESHOLD
        {
            // pass
        } else {
            // TODO: Change KMER to be a refernce to a sequence to avoid cloning
            let read_bases = read.read.seq().as_bytes();
            for offset in 0..read_bases.len().saturating_sub(self.kmer_length) {
                self.counts_by_kmer.add_kmer(
                    Kmer::new_with_start_and_length(
                        read_bases.clone(),
                        offset,
                        self.kmer_length,
                    ),
                    1,
                )
            }
        }
    }

    /**
     * Do actual read correction based on k-mer map. First, loop through stored k-mers to get a list of possible corrections
     * for each position in the read. Then correct read based on all possible consistent corrections.
     * @param inputRead                               Read to correct
     * @return                                        Corrected read (can be same reference as input if doInplaceErrorCorrection is set)
     */
    fn correct_read(&mut self, input_read: BirdToolRead) -> BirdToolRead {
        // do actual correction
        let mut corrected = false;
        let mut corrected_bases = input_read.read.seq().as_bytes();
        let mut corrected_quals = input_read.read.qual().to_vec();

        // array to store list of possible corrections for read
        let correction_set = self.build_correction_map(&corrected_bases);

        for offset in 0..corrected_bases.len() {
            let b = correction_set.get_consensus_correction(offset);
            match b {
                None => continue,
                Some(b) => {
                    if b != corrected_bases[offset] {
                        corrected_bases[offset] = b;
                        corrected_quals[offset] = self.quality_of_corrected_bases;
                        corrected = true;
                        self.read_error_correction_stats.num_bases_corrected += 1;
                    }
                }
            }
        }

        if corrected {
            self.read_error_correction_stats.num_reads_corrected += 1;
            let name = input_read.read.qname().to_vec();
            let mut corrected_read = input_read;
            corrected_read.update(
                name.as_slice(),
                Some(corrected_read.read.cigar().deref()),
                corrected_bases,
                &corrected_quals,
            );

            return corrected_read;
        } else {
            self.read_error_correction_stats.num_reads_uncorrected += 1;
            return input_read;
        }
    }

    /**
     * Build correction map for each of the bases in read.
     * For each of the constituent kmers in read:
     *
     * a) See whether the kmer has been mapped to a corrected kmer.
     *
     * b) If so, get list of differing positions and corresponding bases.
     *
     * c) Add then list of new bases to index in correction list.
     *
     * Correction list is of read size, and holds a list of bases to correct.
     *
     * @param correctedBases                        Bases to attempt to correct
     * @return                                      CorrectionSet object.
     */
    fn build_correction_map(&self, corrected_bases: &[u8]) -> CorrectionSet {
        // array to store list of possible corrections for read
        let mut correction_set = CorrectionSet::new(corrected_bases.len());

        for offset in 0..corrected_bases.len().saturating_sub(self.kmer_length) {
            let kmer =
                Kmer::new_with_start_and_length(corrected_bases.to_vec(), offset, self.kmer_length);
            debug!("Kmer setup: {:?}", &kmer);
            let new_kmer = self.kmer_correction_map.get(&kmer);
            match new_kmer {
                None => continue,
                Some(new_kmer) => {
                    if new_kmer == &kmer {
                        continue;
                    } else {
                        let differing_positions = self.kmer_differing_bases.get(&kmer).unwrap();

                        for k in 0..differing_positions.0.len() {
                            // get list of differing positions for corrected kmer
                            // for each of these, add correction candidate to correction set
                            correction_set
                                .add(offset + differing_positions.0[k], differing_positions.1[k]);
                        }
                    }
                }
            }
        }

        return correction_set;
    }

    /**
     * Top-level entry point that adds a collection of reads to our kmer list.
     * For each read in list, its constituent kmers will be logged in our kmer table.
     * @param reads
     */
    pub fn add_reads_to_kmers(&mut self, reads: Vec<&BirdToolRead>) {
        for read in reads {
            self.add_read_kmers(read);
        }
    }

    /**
     * For each kmer we've seen, do the following:
     * a) If kmer count > threshold1, this kmer is good, so correction map will be to itself.
     * b) If kmer count <= threshold2, this kmer is bad.
     *    In that case, loop through all other kmers. If kmer is good, compute distance, and get minimal distance.
     *    If such distance is < some threshold, map to this kmer, and record differing positions and bases.
     *
     */
    fn compute_kmer_correction_map(&mut self) {
        for stored_kmer in self.counts_by_kmer.get_counted_kmers().into_iter() {
            if stored_kmer.get_count() >= self.min_observations_for_kmer_to_be_solid {
                // this kmer is good: map to itself
                debug!("Kmer {:?} count {:?}", stored_kmer.get_kmer(), stored_kmer.get_count());
                self.kmer_correction_map.insert(
                    stored_kmer.get_kmer().clone(),
                    stored_kmer.get_kmer().clone(),
                );
                self.kmer_differing_bases
                    .insert(stored_kmer.get_kmer().clone(), (vec![0], vec![0]));
                self.read_error_correction_stats.num_solid_kmers += 1;
            } else if stored_kmer.get_count() <= self.max_observations_for_kmer_to_be_correctable {
                // loop now thru all other kmers to find nearest neighbor
                let nearest_neighbour = self.find_nearest_neighbour(stored_kmer.get_kmer());
                match nearest_neighbour.0 {
                    None => {
                        self.read_error_correction_stats.num_uncorrectable_kmers += 1;
                    }
                    Some(neighbour) => {
                        debug!("Stored kmer {:?}", stored_kmer.get_kmer());
                        debug!("Neighbour {:?}", &neighbour);
                        self.kmer_correction_map
                            .insert(stored_kmer.get_kmer().clone(), neighbour);
                        self.kmer_differing_bases
                            .insert(stored_kmer.get_kmer().clone(), nearest_neighbour.1);
                        self.read_error_correction_stats.num_corrected_kmers += 1;
                    }
                }
            }
        }
    }

    /**
     * Finds nearest neighbor of a given k-mer, among a list of counted K-mers, up to a given distance.
     * If many k-mers share same closest distance, an arbitrary k-mer is picked
     * @param kmer                        K-mer of interest
     * @param countsByKMer                KMerCounter storing set of counted k-mers (may include kmer of interest)
     * @param maxDistance                 Maximum distance to search
     * @return                            Pair of values: closest K-mer in Hamming distance and list of differing bases.
     *                                      If no neighbor can be found up to given distance, returns null
     */
    fn find_nearest_neighbour(
        &self,
        kmer: &Kmer,
        // counts_by_kmer: &KmerCounter,
        // max_distance: usize,
    ) -> (Option<Kmer>, (Vec<usize>, Vec<u8>)) {
        // assert!(max_distance >= 1, "max_distance must be >=1, got {}", self.max_mismatches_to_correct);
        let mut minimum_distance = std::i32::MAX;
        let mut closest_kmer = None;

        let mut differing_indices = vec![0; self.max_mismatches_to_correct + 1];
        let mut differing_bases = vec![0; self.max_mismatches_to_correct + 1];

        let mut closest_differing_indices = vec![0; self.max_mismatches_to_correct + 1];
        let mut closest_differing_bases = vec![0; self.max_mismatches_to_correct + 1];

        for candidate_kmer in self.counts_by_kmer.get_counted_kmers() {
            // skip if candidate set includes test kmer
            if candidate_kmer.get_kmer() == kmer {
                continue;
            }

            let hamming_distance = kmer.get_differing_positions(
                candidate_kmer.get_kmer(),
                self.max_mismatches_to_correct,
                &mut differing_indices,
                &mut differing_bases,
            );

            debug!("Hamming distance: {}", hamming_distance);
            if hamming_distance < 0 {
                // can't compare kmer? skip
                continue;
            }

            if hamming_distance < minimum_distance {
                minimum_distance = hamming_distance;
                closest_kmer = Some(candidate_kmer.get_kmer().clone());
                closest_differing_bases = differing_bases.clone();
                closest_differing_indices = differing_indices.clone();
            }
        }

        return (
            closest_kmer,
            (closest_differing_indices, closest_differing_bases),
        );
    }

    /**
     * experimental function to compute max homopolymer length in a given reference context
     * @param fullReferenceWithPadding                Reference context of interest
     * @return                                        Max homopolymer length in region
     */
    fn compute_max_hl_len(full_reference_with_padding: &[u8]) -> usize {
        let mut left_run = 1;
        let mut max_run = 1;

        for i in 1..full_reference_with_padding.len() {
            if full_reference_with_padding[i] == full_reference_with_padding[i - 1] {
                left_run += 1;
            } else {
                left_run = 1;
            }
        }

        if left_run > max_run {
            max_run = left_run
        }

        return max_run;
    }
}

/**
 * Wrapper utility class that holds, for each position in read, a list of bytes representing candidate corrections.
 * So, a read ACAGT where the middle A has found to be errorful might look like:
 * 0: {}
 * 1: {}
 * 2: {'C','C','C'}
 * 3: {}
 * 4: {}
 *
 * It's up to the method getConsensusCorrection()  to decide how to use the correction sets for each position.
 * By default, only strict consensus is allowed right now.
 *
 */
pub struct CorrectionSet {
    size: usize,
    pub corrections: Vec<Vec<u8>>,
}

impl CorrectionSet {
    /**
     * Main class constructor.
     * @param size      Size of correction set, needs to be set equal to the read being corrected
     */
    pub fn new(size: usize) -> Self {
        let mut corrections = Vec::with_capacity(size);
        for k in 0..size {
            corrections.push(Vec::with_capacity(k));
        }

        Self { size, corrections }
    }

    /**
     * Add a base to this correction set at a particular offset, measured from the start of the read
     * @param offset                          Offset from start of read
     * @param base                            base to be added to list of corrections at this offset
     */
    pub fn add(&mut self, offset: usize, base: u8) {
        if offset >= self.size {
            panic!("Bad entry into CorrectionSet: offset > size")
        }
        if !BaseUtils::is_regular_base(base) {
            // pass since no irregular base correction
        } else {
            let mut stored_bytes = &mut self.corrections[offset];
            stored_bytes.push(base)
        }
    }

    /**
     * Get consensus correction for a particular offset. In this implementation, it just boils down to seeing if
     * byte list associated with offset has identical values. If so, return this base, otherwise return null.
     * @param offset
     * @return                                 Consensus base, or null if no consensus possible.
     */
    pub fn get_consensus_correction(&self, offset: usize) -> Option<u8> {
        if offset >= self.size {
            panic!("Bad entry into CorrectionSet: offset > size")
        }
        let stored_bytes = &self.corrections[offset];
        if stored_bytes.len() == 0 {
            return None;
        } else {
            let last_base = stored_bytes.last().unwrap();
            if stored_bytes.iter().any(|b| b != last_base) {
                return None;
            } else {
                return Some(*last_base);
            }
        }
    }
}

impl ReadErrorCorrector for NearbyKmerErrorCorrector {
    /**
     * Correct a collection of reads based on stored k-mer counts
     * @param reads
     */
    fn correct_reads(&mut self, reads: Vec<BirdToolRead>) -> Vec<BirdToolRead> {
        let mut corrected_reads = Vec::with_capacity(reads.len());

        if Self::DONT_CORRECT_IN_LONG_HOMOPOLYMERS
            && self.max_homopolymer_length_in_region > Self::MAX_HOMOPOLYMER_THRESHOLD
        {
            return reads.to_vec(); // Can't correct so just return the reads
        } else {
            self.compute_kmer_correction_map();
            for read in reads {
                let corrected_read = self.correct_read(read);
                if self.trim_low_quality_bases {
                    corrected_reads.push(
                        ReadClipper::new(&corrected_read)
                            .hard_clip_low_qual_ends(self.min_tail_quality),
                    )
                } else {
                    corrected_reads.push(corrected_read)
                }
            }
            debug!(
                "Number of corrected bases: {}",
                self.read_error_correction_stats.num_bases_corrected
            );
            debug!(
                "Number of corrected reads: {}",
                self.read_error_correction_stats.num_reads_corrected
            );
            debug!(
                "Number of skipped reads: {}",
                self.read_error_correction_stats.num_reads_uncorrected
            );
            debug!(
                "Number of solid kmers: {}",
                self.read_error_correction_stats.num_solid_kmers
            );
            debug!(
                "Number of corrected kmers: {}",
                self.read_error_correction_stats.num_corrected_kmers
            );
            debug!(
                "Number of uncorrectable kmers: {}",
                self.read_error_correction_stats.num_uncorrectable_kmers
            );

            return corrected_reads;
        }
    }
}
