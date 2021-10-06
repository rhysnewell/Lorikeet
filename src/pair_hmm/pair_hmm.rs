use model::allele_likelihoods::AlleleLikelihoods;
use model::byte_array_allele::Allele;
use ndarray::prelude::*;
use ndarray::{Array2, Zip};
use pair_hmm::pair_hmm_likelihood_calculation_engine::PairHMMInputScoreImputator;
use pair_hmm::pair_hmm_model::PairHMMModel;
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use rust_htslib::bam::record::Seq;
use utils::quality_utils::QualityUtils;

lazy_static! {
    static ref INITIAL_CONDITION: f64 = 2.0_f64.powf(1020.0);
    static ref INITIAL_CONDITION_LOG10: f64 = (*INITIAL_CONDITION).log10();
    static ref LOG10_3: f64 = 3.0_f64.log10();
}

/**
 * Class for performing the pair HMM for global alignment. Figure 4.1 in Durbin 1998 book.
 */
pub struct PairHMM {
    constants_are_initialized: bool,
    previous_haplotype_length: Option<usize>,
    hap_start_index: Option<usize>,
    max_haplotype_length: usize,
    max_read_length: usize,
    padded_max_read_length: usize,
    padded_max_haplotype_length: usize,
    padded_read_length: Option<usize>,
    padded_haplotype_length: Option<usize>,
    initialized: bool,
    do_not_use_tristate_correction: bool,
    m_log_likelihood_array: Vec<f64>,
    do_profiling: bool,
    transition: Array2<f64>,
    prior: Array2<f64>,
    match_matrix: Array2<f64>,
    insertion_matrix: Array2<f64>,
    deletion_matrix: Array2<f64>,
    model: PairHMMModel,
    do_exact_log10: bool,
    logless: bool,
}

impl PairHMM {
    const TRISTATE_CORRECTION: f64 = 3.0;
    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     *
     * Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
     *
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     * @throws IllegalArgumentException if haplotypeMaxLength is less than or equal to zero
     */
    pub fn initialize(max_read_length: usize, haplotype_max_length: usize) -> PairHMM {
        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for
        // initial conditions and + 1 to consider the final base in a non-global alignment
        let padded_max_read_length = max_read_length + 1;
        let padded_max_haplotype_length = haplotype_max_length + 1;

        let mut initialized = true;
        if max_read_length == 0 && haplotype_max_length == 0 {
            initialized = false;
        }

        PairHMM {
            max_read_length,
            max_haplotype_length: haplotype_max_length,
            padded_max_read_length,
            padded_max_haplotype_length,
            match_matrix: Array2::zeros((padded_max_read_length, padded_max_haplotype_length)),
            insertion_matrix: Array2::zeros((padded_max_read_length, padded_max_haplotype_length)),
            deletion_matrix: Array2::zeros((padded_max_read_length, padded_max_haplotype_length)),
            transition: PairHMMModel::create_transition_matrix(max_read_length),
            prior: Array2::zeros((padded_max_read_length, padded_max_haplotype_length)),
            initialized,
            do_not_use_tristate_correction: false,
            m_log_likelihood_array: Vec::new(),
            do_profiling: true,
            previous_haplotype_length: None,
            constants_are_initialized: false,
            model: PairHMMModel::new(),
            padded_haplotype_length: None,
            padded_read_length: None,
            hap_start_index: None,
            do_exact_log10: false,
            logless: true,
        }
    }

    /**
     * Set the do_exact_log10 parameter
     *
     * @param do_exact_log10 should the log calculations be exact (slow) or approximate (faster)
     */
    pub fn set_do_exact_log10(&mut self, do_exact_log10: bool) {
        if do_exact_log10 {
            self.logless = false
        }
        self.do_exact_log10 = do_exact_log10
    }

    pub fn set_logless(&mut self, logless: bool) {
        if logless {
            if self.do_exact_log10 {
                self.do_exact_log10 = false
            }
        }
        self.logless = logless
    }

    /**
     * Only used for debugging purposes
     */
    pub fn do_not_use_tristate_correction(&mut self) {
        self.do_not_use_tristate_correction = true
    }

    fn reinitialize(&mut self, max_read_length: usize, haplotype_max_length: usize) {
        let padded_max_read_length = max_read_length + 1;
        let padded_max_haplotype_length = haplotype_max_length + 1;
        self.max_read_length = max_read_length;
        self.max_haplotype_length = haplotype_max_length;
        self.padded_max_read_length = padded_max_read_length;
        self.padded_max_haplotype_length = padded_max_haplotype_length;
        // self.match_matrix = Array2::zeros((padded_max_read_length, padded_max_haplotype_length));
        // self.insertion_matrix =
        //     Array2::zeros((padded_max_read_length, padded_max_haplotype_length));
        // self.deletion_matrix = Array2::zeros((padded_max_read_length, padded_max_haplotype_length));
        // self.transition = PairHMMModel::create_transition_matrix(max_read_length);
        // self.prior = Array2::zeros((padded_max_read_length, padded_max_haplotype_length));
        self.initialized = true;
    }

    /**
     *  Given a list of reads and haplotypes, for every read compute the total probability of said read arising from
     *  each haplotype given base substitution, insertion, and deletion probabilities.
     *
     * @param processedReads reads to analyze instead of the ones present in the destination read-likelihoods.
     * @param logLikelihoods where to store the log likelihoods where position [a][r] is reserved for the log likelihood of {@code reads[r]}
     *             conditional to {@code alleles[a]}.
     */
    pub fn compute_log10_likelihoods<A: Allele>(
        &mut self,
        sample_index: usize,
        allele_likelihoods: &mut AlleleLikelihoods<A>,
        processed_reads: Vec<BirdToolRead>,
        input_score_imputator: &PairHMMInputScoreImputator,
    ) {
        if !processed_reads.is_empty() {
            // (re)initialize the pairHMM only if necessary
            let max_read_length = processed_reads
                .par_iter()
                .map(|r| r.len())
                .max()
                .unwrap_or(0);
            let max_haplotype_length = allele_likelihoods
                .alleles
                .as_list_of_alleles()
                .par_iter()
                .map(|hap| hap.length())
                .max()
                .unwrap_or(0);
            if !self.initialized
                || max_haplotype_length > self.max_haplotype_length
                || max_read_length > self.max_read_length
            {
                self.reinitialize(max_read_length, max_haplotype_length);
            };

            let read_count = processed_reads.len();
            let allele_count = allele_likelihoods.alleles.len();

            self.m_log_likelihood_array = vec![0.0; read_count * allele_count];
            let mut idx = 0;
            let mut read_index = 0;
            for read in processed_reads {
                let read_bases = read.read.seq().as_bytes();
                let read_quals = read.read.qual();
                let read_ins_quals = input_score_imputator.ins_open_penalties(&read);
                let read_del_quals = input_score_imputator.del_open_penalties(&read);
                let overall_gcp = input_score_imputator.gap_continuation_penalties(&read);
                // peek at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
                let mut is_first_haplotype = true;
                for a in 0..allele_count {
                    let allele = &allele_likelihoods.alleles.get_allele(a);
                    let allele_bases = allele.get_bases();
                    let next_allele_bases = if a == allele_count - 1 {
                        None
                    } else {
                        Some(allele_likelihoods.alleles.get_allele(a + 1).get_bases())
                    };
                    let lk = self.compute_read_likelihood_given_haplotype_log10(
                        allele_bases,
                        read_bases.as_slice(),
                        read_quals,
                        &read_ins_quals,
                        &read_del_quals,
                        &overall_gcp,
                        is_first_haplotype,
                        next_allele_bases,
                    );

                    allele_likelihoods.values_by_sample_index[sample_index][[a, read_index]] = lk;
                    self.m_log_likelihood_array[idx] = lk;
                    if is_first_haplotype {
                        is_first_haplotype = false;
                    }
                    idx += 1;
                }
                read_index += 1;
            }
        }
    }

    pub fn get_log_likelihood_array(&self) -> &Vec<f64> {
        &self.m_log_likelihood_array
    }

    /**
     * Compute the total probability of read arising from haplotypeBases given base substitution, insertion, and deletion
     * probabilities.
     *
     * Note on using hapStartIndex.  This allows you to compute the exact true likelihood of a full haplotypes
     * given a read, assuming that the previous calculation read over a full haplotype, recaching the read values,
     * starting only at the place where the new haplotype bases and the previous haplotype bases differ.  This
     * index is 0-based, and can be computed with findFirstPositionWhereHaplotypesDiffer given the two haplotypes.
     * Note that this assumes that the read and all associated quals values are the same.
     *
     * @param haplotypeBases the full sequence (in standard SAM encoding) of the haplotype, must be >= than read bases in length
     * @param readBases the bases (in standard encoding) of the read, must be <= haplotype bases in length
     * @param readQuals the phred-scaled per base substitution quality scores of read.  Must be the same length as readBases
     * @param insertionGOP the phred-scaled per base insertion quality scores of read.  Must be the same length as readBases
     * @param deletionGOP the phred-scaled per base deletion quality scores of read.  Must be the same length as readBases
     * @param overallGCP the phred-scaled gap continuation penalties scores of read.  Must be the same length as readBases
     * @param recacheReadValues if false, we don't recalculate any cached results, assuming that readBases and its associated
     *                          parameters are the same, and only the haplotype bases are changing underneath us
     * @throws IllegalStateException  if did not call initialize() beforehand
     * @throws IllegalArgumentException haplotypeBases is null or greater than maxHaplotypeLength
     * @throws IllegalArgumentException readBases is null or greater than maxReadLength
     * @throws IllegalArgumentException readBases, readQuals, insertionGOP, deletionGOP and overallGCP are not the same size
     * @return the log10 probability of read coming from the haplotype under the provided error model
     */
    pub fn compute_read_likelihood_given_haplotype_log10(
        &mut self,
        haplotype_bases: &[u8],
        read_bases: &[u8],
        read_quals: &[u8],
        insertion_gop: &[u8],
        deletion_gop: &[u8],
        overall_gcp: &[u8],
        recache_read_values: bool,
        next_haplotype_bases: Option<&[u8]>,
        // current_allele_index: usize,
    ) -> f64 {
        assert!(
            self.initialized,
            "Must call initialize before calling compute_read_likelihood_given_haplotype_log10"
        );
        assert!(
            haplotype_bases.len() <= self.max_haplotype_length,
            "Haplotype bases is too long"
        );
        assert!(
            read_quals.len() == read_bases.len(),
            "Read bases and read quals aren't the same size"
        );
        assert!(
            insertion_gop.len() == read_bases.len(),
            "Read bases and insertion gcp aren't the same size"
        );
        assert!(
            deletion_gop.len() == read_bases.len(),
            "Read bases and deletion gcp aren't the same size"
        );
        assert!(
            overall_gcp.len() == read_bases.len(),
            "Read bases and overal GCP aren't the same size"
        );

        self.padded_read_length = Some(read_bases.len() + 1);
        self.padded_haplotype_length = Some(haplotype_bases.len() + 1);
        self.hap_start_index = if recache_read_values {
            Some(0)
        } else {
            self.hap_start_index
        };

        // Pre-compute the difference between the current haplotype and the next one to be run
        // Looking ahead is necessary for the ArrayLoglessPairHMM implementation
        let next_hap_start_index = match next_haplotype_bases {
            None => 0,
            Some(next_haplotype_bases) => {
                if haplotype_bases.len() != next_haplotype_bases.len() {
                    0
                } else {
                    Self::find_first_position_where_haplotypes_differ(
                        haplotype_bases,
                        next_haplotype_bases,
                    )
                }
            }
        };

        let result = self.sub_compute_read_likelihood_given_haplotype_log10(
            haplotype_bases,
            read_bases,
            read_quals,
            insertion_gop,
            deletion_gop,
            overall_gcp,
            self.hap_start_index.unwrap(),
            recache_read_values,
            next_hap_start_index,
        );

        assert!(
            result <= 0.0,
            "PairHmm Log Probability cannot be greater than 0.0"
        );
        // Warning: This assumes no downstream modification of the haplotype bases (saves us from copying the array). It is okay for the haplotype caller.
        // self.previous_haplotype_bases = haplotype_bases; // this won't work due to reference only existing inside this function

        // For the next iteration, the hapStartIndex for the next haploytpe becomes the index for the current haplotype
        // The array implementation has to look ahead to the next haplotype to store caching info. It cannot do this if nextHapStart is before hapStart
        self.hap_start_index = match self.hap_start_index {
            None => Some(next_hap_start_index),
            Some(hap_start_index) => {
                if next_hap_start_index < hap_start_index {
                    Some(0)
                } else {
                    Some(next_hap_start_index)
                }
            }
        };

        self.previous_haplotype_length = Some(haplotype_bases.len());

        return result;
    }

    pub fn sub_compute_read_likelihood_given_haplotype_log10(
        &mut self,
        haplotype_bases: &[u8],
        read_bases: &[u8],
        read_quals: &[u8],
        insertion_gop: &[u8],
        deletion_gop: &[u8],
        overall_gcp: &[u8],
        hap_start_index: usize,
        recache_read_values: bool,
        next_hap_start_index: usize,
    ) -> f64 {
        match self.previous_haplotype_length {
            None => {
                let initial_value = *INITIAL_CONDITION / haplotype_bases.len() as f64;
                // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
                // self.deletion_matrix.row_mut(0).fill(initial_value);
                self.deletion_matrix.row_mut(0).fill(initial_value);
            }
            Some(previous_haplotype_length) => {
                if previous_haplotype_length != haplotype_bases.len() {
                    let initial_value = *INITIAL_CONDITION / haplotype_bases.len() as f64;
                    // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
                    self.deletion_matrix.row_mut(0).fill(initial_value);
                };
            }
        };

        if !self.constants_are_initialized || recache_read_values {
            self.initialize_probabilities(insertion_gop, deletion_gop, overall_gcp);
            self.constants_are_initialized = true;
        };

        self.initialize_priors(haplotype_bases, read_bases, read_quals, hap_start_index);

        // let padded_haplotype_length = self.padded_haplotype_length.unwrap();
        // // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
        // let prior = &self.prior;
        // let transition = &self.transition;
        // let mut match_matrix = &self.match_matrix;
        // let mut insertion_matrix = &self.insertion_matrix;
        // let mut deletion_matrix = &self.deletion_matrix;
        // Zip::indexed(&mut self.match_matrix)
        //     .and(&mut self.insertion_matrix)
        //     .and(&mut self.deletion_matrix)
        //     .apply(|(i, j), match_val, insertion_val, deletion_val| {
        //         if i > 0 && (j > hap_start_index && j < padded_haplotype_length) {
        //             *match_val = prior[[i, j]]
        //                 * (match_matrix[[i - 1, j - 1]]
        //                 * transition[[i, PairHMMModel::match_to_match]]
        //                 + insertion_matrix[[i - 1, j - 1]]
        //                 * transition[[i, PairHMMModel::indel_to_match]]
        //                 + deletion_matrix[[i - 1, j - 1]]
        //                 * transition[[i, PairHMMModel::indel_to_match]]);
        //
        //             *insertion_val = match_matrix[[i - 1, j]]
        //                 * transition[[i, PairHMMModel::match_to_insertion]]
        //                 + insertion_matrix[[i - 1, j]]
        //                 * transition[[i, PairHMMModel::insertion_to_insertion]];
        //
        //             *deletion_val = match_matrix[[i, j - 1]]
        //                 * transition[[i, PairHMMModel::match_to_deletion]]
        //                 + deletion_matrix[[i, j - 1]]
        //                 * transition[[i, PairHMMModel::deletion_to_deletion]];
        //         }
        //     });
        // TODO: This part is slow for large reads and haplotypes. It cannot be parallelized due to
        //       requiring previous cells be populated with values before the next cell can be calculated
        //       Unsure what to do to speed this up?
        // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
        for i in 1..self.padded_read_length.unwrap() {
            for j in (hap_start_index + 1)..self.padded_haplotype_length.unwrap() {
                self.match_matrix[[i, j]] = self.prior[[i, j]]
                    * (self.match_matrix[[i - 1, j - 1]]
                        * self.transition[[i, PairHMMModel::match_to_match]]
                        + self.insertion_matrix[[i - 1, j - 1]]
                            * self.transition[[i, PairHMMModel::indel_to_match]]
                        + self.deletion_matrix[[i - 1, j - 1]]
                            * self.transition[[i, PairHMMModel::indel_to_match]]);

                self.insertion_matrix[[i, j]] = self.match_matrix[[i - 1, j]]
                    * self.transition[[i, PairHMMModel::match_to_insertion]]
                    + self.insertion_matrix[[i - 1, j]]
                        * self.transition[[i, PairHMMModel::insertion_to_insertion]];

                self.deletion_matrix[[i, j]] = self.match_matrix[[i, j - 1]]
                    * self.transition[[i, PairHMMModel::match_to_deletion]]
                    + self.deletion_matrix[[i, j - 1]]
                        * self.transition[[i, PairHMMModel::deletion_to_deletion]];
            }
        }

        // final log probability is the log10 sum of the last element in the Match and Insertion state arrays
        // this way we ignore all paths that ended in deletions! (huge)
        // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
        let end_i = self.padded_read_length.unwrap() - 1;
        // let mut final_sum_probabilities = 0.0;
        // for j in 1..self.padded_haplotype_length.unwrap() {
        //     final_sum_probabilities += self.match_matrix[[end_i, j]] + self.insertion_matrix[[end_i, j]];
        // };
        let match_matrix = &self.match_matrix;
        let insertion_matrix = &self.insertion_matrix;
        // potential parallel implementation
        let final_sum_probabilities: f64 = (1..self.padded_haplotype_length.unwrap())
            .into_par_iter()
            .fold_with(0.0_f64, |a: f64, j: usize| {
                a + (match_matrix[[end_i, j]] + insertion_matrix[[end_i, j]])
            })
            .sum::<f64>();

        return final_sum_probabilities.log10() - *INITIAL_CONDITION_LOG10;
    }

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    fn initialize_priors(
        &mut self,
        haplotype_bases: &[u8],
        read_bases: &[u8],
        read_quals: &[u8],
        start_index: usize,
    ) {
        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.
        let do_not_use_tristate_correction = self.do_not_use_tristate_correction;

        // Potential parallel implementation
        Zip::indexed(&mut self.prior).par_for_each(|(i, j), value| {
            if (i > 0 && i <= read_bases.len())
                && (j >= start_index + 1 && j <= haplotype_bases.len())
            {
                let x = read_bases[i - 1];
                let qual = read_quals[i - 1];
                let y = haplotype_bases[j - 1];
                *value = if x == y || x == b'N' || y == b'N' {
                    QualityUtils::qual_to_prob(qual)
                } else {
                    QualityUtils::qual_to_error_prob(qual)
                        / (if do_not_use_tristate_correction {
                            1.0
                        } else {
                            Self::TRISTATE_CORRECTION
                        })
                }
            }
        });

        // Non parallel implementation
        // for i in 0..read_bases.len() {
        //     let x = read_bases[i];
        //     let qual = read_quals[i];
        //     for j in start_index..haplotype_bases.len() {
        //         let y = haplotype_bases[j];
        //         self.prior[[i + 1, j + 1]] = if x == y || x == b'N' || y == b'N' {
        //             QualityUtils::qual_to_prob(qual)
        //         } else {
        //             QualityUtils::qual_to_error_prob(qual) / (if self.do_not_use_tristate_correction {
        //                 1.0
        //             } else {
        //                 Self::TRISTATE_CORRECTION
        //             })
        //         };
        //     }
        // }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    fn initialize_probabilities(
        &mut self,
        insertion_gop: &[u8],
        deletion_gop: &[u8],
        overall_gcp: &[u8],
    ) {
        self.model.qual_to_trans_probs_with_array(
            &mut self.transition,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        );
    }

    /**
     * Compute the first position at which two haplotypes differ
     *
     * If the haplotypes are exact copies of each other, returns the min length of the two haplotypes.
     *
     * @param haplotype1 the first haplotype1
     * @param haplotype2 the second haplotype1
     * @throws IllegalArgumentException if haplotype1 or haplotype2 are null or zero length
     * @return the index of the first position in haplotype1 and haplotype2 where the byte isn't the same
     */
    pub fn find_first_position_where_haplotypes_differ(
        haplotype1: &[u8],
        haplotype2: &[u8],
    ) -> usize {
        match (0..std::cmp::min(haplotype1.len(), haplotype2.len()))
            .into_par_iter()
            .position_first(|i| haplotype1[i] != haplotype2[i])
        {
            Some(i) => return i,
            None => return std::cmp::min(haplotype1.len(), haplotype2.len()),
        }
    }
}
