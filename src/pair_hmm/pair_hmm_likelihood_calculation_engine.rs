use assembly::assembly_result_set::AssemblyResultSet;
use haplotype::haplotype::Haplotype;
use model::allele_likelihoods::AlleleLikelihoods;
use model::variant_context_utils::VariantContextUtils;
use ordered_float::OrderedFloat;
use pair_hmm::pair_hmm::PairHMM;
use rayon::prelude::*;
use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use reads::bird_tool_reads::BirdToolRead;
use reads::read_clipper::ReadClipper;
use reads::read_utils::ReadUtils;
use std::cmp::{max, min};
use std::collections::HashMap;
use utils::quality_utils::QualityUtils;
use utils::simple_interval::SimpleInterval;

lazy_static! {
    // table used for disqualifying reads for genotyping
    // Format for each row of table: baseQ, mean, variance
    // Actual threshold is calculated over the length of the read as:
    // sum(means) + K * sqrt(sum(variances))
    static ref dynamic_read_qual_thresh_lookup_table: Vec<f64> = vec![
                    //baseQ,mean,variance
                    1.0,  5.996842844, 0.196616587, 2.0,  5.870018422, 1.388545569, 3.0,  5.401558531, 5.641990128,
                    4.0,  4.818940919, 10.33176216, 5.0,  4.218758304, 14.25799688, 6.0,  3.646319832, 17.02880749,
                    7.0,  3.122346753, 18.64537883, 8.0,  2.654731979, 19.27521677, 9.0,  2.244479156, 19.13584613,
                    10.0, 1.88893867,  18.43922003, 11.0, 1.583645342, 17.36842261, 12.0, 1.3233807, 16.07088712,
                    13.0, 1.102785365, 14.65952563, 14.0, 0.916703025, 13.21718577, 15.0, 0.760361881, 11.80207947,
                    16.0, 0.629457387, 10.45304833, 17.0, 0.520175654, 9.194183767, 18.0, 0.42918208,  8.038657241,
                    19.0, 0.353590663, 6.991779595, 20.0, 0.290923699, 6.053379213, 21.0, 0.23906788,  5.219610436,
                    22.0, 0.196230431, 4.484302033, 23.0, 0.160897421, 3.839943445, 24.0, 0.131795374, 3.27839108,
                    25.0, 0.1078567,   2.791361596, 26.0, 0.088189063, 2.370765375, 27.0, 0.072048567, 2.008921719,
                    28.0, 0.058816518, 1.698687797, 29.0, 0.047979438, 1.433525748, 30.0, 0.039111985, 1.207526336,
                    31.0, 0.031862437, 1.015402928, 32.0, 0.025940415, 0.852465956, 33.0, 0.021106532, 0.714585285,
                    34.0, 0.017163711, 0.598145851, 35.0, 0.013949904, 0.500000349, 36.0, 0.011332027, 0.41742159,
                    37.0, 0.009200898, 0.348056286, 38.0, 0.007467036, 0.289881373, 39.0, 0.006057179, 0.241163527,
                    40.0, 0.004911394, 0.200422214
                ];
}

pub struct PairHMMLikelihoodCalculationEngine {
    constant_gcp: u8,
    log10_global_read_mismapping_rate: f64,
    pair_hmm: PairHMM,
    dynamic_disqualification: bool,
    read_disqualification_scale: f64,
    expected_error_rate_per_base: f64,
    disable_cap_read_qualities_to_mapq: bool,
    symmetrically_normalize_alleles_to_reference: bool,
    modify_soft_clipped_bases: bool,
    pcr_error_model: PCRErrorModel,
    base_quality_score_threshold: u8,
    pcr_indel_error_model_cache: Vec<u8>,
    input_score_imputator: PairHMMInputScoreImputator,
}

#[derive(Debug, Copy, Clone)]
pub enum PCRErrorModel {
    /** no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used */
    None = 0,
    /** a most aggressive model will be applied that sacrifices true positives in order to remove more false positives */
    Hostile = 1,
    /** a more aggressive model will be applied that sacrifices true positives in order to remove more false positives */
    Aggresive = 2,
    /** a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives */
    Conservative = 3,
}

impl PCRErrorModel {
    /**
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer).
     */
    pub fn new(args: &clap::ArgMatches) -> PCRErrorModel {
        let pcr_error_model_arg = args
            .value_of("pcr-indel-model")
            .unwrap()
            .to_ascii_lowercase();
        let pcr_error_model = match pcr_error_model_arg.as_str() {
            "none" => PCRErrorModel::None,
            "hostile" => PCRErrorModel::Hostile,
            "aggressive" => PCRErrorModel::Aggresive,
            "conservative" => PCRErrorModel::Conservative,
            _ => panic!("Unknown PCR Error Model"),
        };

        return pcr_error_model;
    }
}

impl PairHMMLikelihoodCalculationEngine {
    const DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR: f64 = 1.0;
    const MAX_STR_UNIT_LENGTH: usize = 8;
    const MAX_REPEAT_LENGTH: usize = 20;
    const MIN_ADJUSTED_QSCORE: usize = 10;
    const MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ: usize = 40;
    const DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH: usize = 3;
    const DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET: usize = 1;

    const INITIAL_QSCORE: f64 = 40.0;
    pub const HMM_BASE_QUALITIES_TAG: &'static str = "HMMQuals";
    //
    // The expected rate of random sequencing errors for a read originating from its true haplotype.
    //
    // For example, if this is 0.01, then we'd expect 1 error per 100 bp.
    //
    pub const DEFAULT_EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;

    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param pcrErrorModel model to correct for PCR indel artifacts
     */
    pub fn new(
        constant_gcp: u8,
        args: &clap::ArgMatches,
        log10_global_read_mismapping_rate: f64,
        pcr_error_model: PCRErrorModel,
        base_quality_score_threshold: u8,
        dynamic_read_disqualification: bool,
        read_disqualification_scale: f64,
        expected_error_rate_per_base: f64,
        symmetrically_normalize_alleles_to_reference: bool,
        disable_cap_read_qualities_to_mapq: bool,
        modify_soft_clipped_bases: bool,
    ) -> PairHMMLikelihoodCalculationEngine {
        assert!(
            base_quality_score_threshold >= QualityUtils::MIN_USABLE_Q_SCORE,
            "base_quality_score_threshold must be greater than or equal to 6"
        );

        let mut result = PairHMMLikelihoodCalculationEngine {
            constant_gcp,
            log10_global_read_mismapping_rate,
            pair_hmm: PairHMM::initialize(0, 0),
            dynamic_disqualification: dynamic_read_disqualification,
            read_disqualification_scale,
            expected_error_rate_per_base,
            disable_cap_read_qualities_to_mapq,
            symmetrically_normalize_alleles_to_reference,
            modify_soft_clipped_bases,
            pcr_error_model,
            base_quality_score_threshold,
            pcr_indel_error_model_cache: Vec::new(),
            input_score_imputator: PairHMMInputScoreImputator::new(constant_gcp),
        };

        result.initialize_pcr_error_model();

        return result;
    }

    fn initialize_pcr_error_model(&mut self) {
        self.pcr_indel_error_model_cache = vec![0; Self::MAX_REPEAT_LENGTH + 1];

        match self.pcr_error_model {
            PCRErrorModel::None => {
                // do nothing
            }
            _ => {
                let rate_factor = self.pcr_error_model as i64 as f64;
                for i in 0..=Self::MAX_REPEAT_LENGTH {
                    self.pcr_indel_error_model_cache[i] =
                        Self::get_error_model_adjusted_qual(i, rate_factor);
                }
            }
        }
    }

    fn get_error_model_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
        return max(
            Self::MIN_ADJUSTED_QSCORE,
            (Self::INITIAL_QSCORE
                - (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp()
                + 1.0) as usize,
        ) as u8;
    }

    pub fn compute_read_likelihoods<'a, 'b, A: AbstractReadThreadingGraph>(
        &mut self,
        assembly_result_set: &'b mut AssemblyResultSet<'a, A>,
        samples: Vec<String>,
        per_sample_read_list: HashMap<usize, Vec<BirdToolRead>>,
    ) -> AlleleLikelihoods<'a> {
        let haplotypes: Vec<Haplotype<SimpleInterval>> = assembly_result_set
            .haplotypes
            .iter()
            .cloned()
            .collect::<Vec<Haplotype<SimpleInterval>>>();
        self.initialize_pair_hmm(&haplotypes, &per_sample_read_list);
        // Add likelihoods for each sample's reads to our result
        let sample_count = samples.len();
        let mut result = AlleleLikelihoods::new(haplotypes, samples, per_sample_read_list);

        for i in 0..sample_count {
            self.compute_read_likelihoods_in_matrix(i, &mut result);
        }

        result.normalize_likelihoods(
            self.log10_global_read_mismapping_rate,
            self.symmetrically_normalize_alleles_to_reference,
        );

        if self.dynamic_disqualification {
            result.filter_poorly_modeled_evidence(Self::dynamic_log10_min_likelihood_model(
                self.read_disqualification_scale,
                Self::log10_min_true_likelihood(self.expected_error_rate_per_base, false),
            ));
        } else {
            result.filter_poorly_modeled_evidence(Self::log10_min_true_likelihood(
                self.expected_error_rate_per_base,
                true,
            ));
        };

        return result;
    }

    fn dynamic_log10_min_likelihood_model(
        dynamic_read_qual_constant: f64,
        log10_min_true_likelihood: Box<dyn Fn(&BirdToolRead) -> f64>,
    ) -> Box<dyn Fn(&BirdToolRead) -> f64> {
        Box::new(move |read| {
            let dynamic_threshold =
                Self::calculate_log10_dynamic_read_qual_threshold(read, dynamic_read_qual_constant);
            let log10_max_likelihoods_for_true_allele = (log10_min_true_likelihood)(read);
            if dynamic_threshold < log10_max_likelihoods_for_true_allele {
                debug!(
                    "For read {:?} replacing old threshold ({}) with new threshold: {}",
                    std::str::from_utf8(read.read.qname()),
                    log10_max_likelihoods_for_true_allele,
                    dynamic_threshold
                );
                return dynamic_threshold;
            } else {
                log10_max_likelihoods_for_true_allele
            }
        })
    }

    fn calculate_log10_dynamic_read_qual_threshold(
        read: &BirdToolRead,
        dynamic_read_qual_constant: f64,
    ) -> f64 {
        let mut sum_mean: f64 = 0.0;
        let mut sum_variance: f64 = 0.0;

        let base_qualities = match read.transient_attributes.get("HMMQuals") {
            None => read.read.qual(),
            Some(read_qual) => read_qual,
        };

        for qual_byte in base_qualities.iter() {
            let bq = *qual_byte as usize;
            // bound the base qualities for lookup between 1 and 40
            let entry_index = if bq <= 1 {
                0
            } else {
                min(Self::MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ, bq) - 1
            };

            let mean_offset = entry_index * Self::DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH
                + Self::DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET;
            let var_offset = mean_offset + 1;
            sum_mean += dynamic_read_qual_thresh_lookup_table[mean_offset];
            sum_variance += dynamic_read_qual_thresh_lookup_table[var_offset];
        }

        let threshold = sum_mean + dynamic_read_qual_constant * sum_variance.sqrt();
        return threshold * -0.1;
    }

    fn log10_min_true_likelihood(
        maximum_error_per_base: f64,
        cap_likelihoods: bool,
    ) -> Box<dyn Fn(&BirdToolRead) -> f64> {
        Box::new(move |read| {
            // TODO this might be replaced by an explicit calculation
            let qualified_read_length = match read.transient_attributes.get("HMMQuals") {
                None => read.len(),
                Some(value) => value.len(),
            };
            let max_errors_for_read = if cap_likelihoods {
                min(
                    OrderedFloat(2.0),
                    OrderedFloat(f64::ceil(
                        qualified_read_length as f64 * maximum_error_per_base,
                    )),
                )
                .into_inner()
            } else {
                f64::ceil(qualified_read_length as f64 * maximum_error_per_base)
            };

            let log10_qual_per_base = -4.0;

            return max_errors_for_read * log10_qual_per_base;
        })
    }

    fn compute_read_likelihoods_in_matrix(
        &mut self,
        sample_index: usize,
        likelihoods: &mut AlleleLikelihoods,
    ) {
        // Modify the read qualities by applying the PCR error model and capping the minimum base,
        // insertion,deletion qualities
        let processed_reads = self.modify_read_qualities(
            likelihoods
                .evidence_by_sample_index
                .get_mut(&sample_index)
                .unwrap(),
        );

        self.pair_hmm.compute_log10_likelihoods(
            sample_index,
            likelihoods,
            processed_reads,
            &self.input_score_imputator,
        );
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    fn modify_read_qualities(&self, reads: &mut Vec<BirdToolRead>) -> Vec<BirdToolRead> {
        let mut result = Vec::new();

        for read in reads.iter_mut() {
            if self.modify_soft_clipped_bases {
                let bases = read.read.seq().as_bytes();
                let mut read_quals = read.read.qual().to_vec();
                let mut read_ins_quals = ReadUtils::get_base_insertion_qualities(read);
                let mut read_del_quals = ReadUtils::get_base_deletion_qualities(read);
                self.apply_pcr_error_model(&bases, &mut read_ins_quals, &mut read_del_quals);

                Self::cap_minimum_read_qualities(
                    read,
                    &mut read_quals,
                    &mut read_ins_quals,
                    &mut read_del_quals,
                    self.base_quality_score_threshold,
                    self.disable_cap_read_qualities_to_mapq,
                );

                // Store the actual qualities
                read.set_transient_attribute(format!("HMM_BASE_QUALITIES_TAG"), read_quals.clone());
                // Create a new copy of the read and sets its base qualities to the modified versions.
                result.push(Self::create_quality_modified_read(
                    &read,
                    bases,
                    read_quals,
                    read_ins_quals,
                    read_del_quals,
                ));
            } else {
                let maybe_unclipped = ReadClipper::new(&read).hard_clip_soft_clipped_bases();
                let bases = maybe_unclipped.read.seq().as_bytes();
                let mut read_quals = maybe_unclipped.read.qual().to_vec();
                let mut read_ins_quals = ReadUtils::get_base_insertion_qualities(&maybe_unclipped);
                let mut read_del_quals = ReadUtils::get_base_deletion_qualities(&maybe_unclipped);
                self.apply_pcr_error_model(&bases, &mut read_ins_quals, &mut read_del_quals);

                Self::cap_minimum_read_qualities(
                    &maybe_unclipped,
                    &mut read_quals,
                    &mut read_ins_quals,
                    &mut read_del_quals,
                    self.base_quality_score_threshold,
                    self.disable_cap_read_qualities_to_mapq,
                );

                // Store the actual qualities
                read.set_transient_attribute(format!("HMM_BASE_QUALITIES_TAG"), read_quals.clone());
                // Create a new copy of the read and sets its base qualities to the modified versions.
                result.push(Self::create_quality_modified_read(
                    &read,
                    bases,
                    read_quals,
                    read_ins_quals,
                    read_del_quals,
                ));
            }
        }

        return result;
    }

    fn cap_minimum_read_qualities(
        read: &BirdToolRead,
        read_quals: &mut Vec<u8>,
        read_ins_quals: &mut Vec<u8>,
        read_del_quals: &mut Vec<u8>,
        base_quality_score_threshold: u8,
        disable_cap_read_qualities_to_mapq: bool,
    ) {
        for i in 0..read_quals.len() {
            if !disable_cap_read_qualities_to_mapq {
                read_quals[i] = min(read_quals[i], read.read.mapq());
            };
            read_quals[i] = Self::set_to_fixed_value_if_too_low(
                read_quals[i],
                base_quality_score_threshold,
                QualityUtils::MIN_USABLE_Q_SCORE,
            );

            read_ins_quals[i] = Self::set_to_fixed_value_if_too_low(
                read_ins_quals[i],
                QualityUtils::MIN_USABLE_Q_SCORE,
                QualityUtils::MIN_USABLE_Q_SCORE,
            );

            read_del_quals[i] = Self::set_to_fixed_value_if_too_low(
                read_del_quals[i],
                QualityUtils::MIN_USABLE_Q_SCORE,
                QualityUtils::MIN_USABLE_Q_SCORE,
            );
        }
    }

    fn set_to_fixed_value_if_too_low(current_val: u8, min_qual: u8, fixed_qual: u8) -> u8 {
        if current_val < min_qual {
            return fixed_qual;
        } else {
            current_val
        }
    }

    /**
     * Creates a new GATKRead with the source read's header, read group and mate
     * information, but with the following fields set to user-supplied values:
     *  - Read Bases
     *  - Base Qualities
     *  - Base Insertion Qualities
     *  - Base Deletion Qualities
     *
     *  Cigar string is empty (not-null)
     *
     * Use this method if you want to create a new GATKRead based on
     * another GATKRead, but with modified bases and qualities
     *
     * @param read a read to copy the header from
     * @param readBases an array containing the new bases you wish use in place of the originals
     * @param baseQualities an array containing the new base qualities you wish use in place of the originals
     * @param baseInsertionQualities an array containing the new base insertion qaulities
     * @param baseDeletionQualities an array containing the new base deletion qualities
     * @return a read with modified bases and qualities, safe for the GATK
     */
    fn create_quality_modified_read(
        read: &BirdToolRead,
        read_bases: Vec<u8>,
        base_qualities: Vec<u8>,
        base_insertion_qualities: Vec<u8>,
        base_deletion_qualities: Vec<u8>,
    ) -> BirdToolRead {
        let mut processed_read =
            ReadUtils::empty_read_with_quals_and_bases(&read, &base_qualities, &read_bases);
        ReadUtils::set_insertion_base_qualities(&mut processed_read, &base_insertion_qualities);
        ReadUtils::set_deletion_base_qualities(&mut processed_read, &base_deletion_qualities);
        return processed_read;
    }

    fn apply_pcr_error_model(
        &self,
        read_bases: &[u8],
        read_ins_quals: &mut Vec<u8>,
        read_del_quals: &mut Vec<u8>,
    ) {
        match self.pcr_error_model {
            PCRErrorModel::None => {
                // pass
            }
            _ => {
                for i in 1..read_bases.len() {
                    let repeat_length = Self::find_tandem_repeat_units(read_bases, i - 1).1;
                    read_ins_quals[i - 1] = min(
                        read_ins_quals[i - 1],
                        self.pcr_indel_error_model_cache[repeat_length],
                    );
                    read_del_quals[i - 1] = min(
                        read_del_quals[i - 1],
                        self.pcr_indel_error_model_cache[repeat_length],
                    );
                }
            }
        }
    }

    fn find_tandem_repeat_units(read_bases: &[u8], offset: usize) -> (Vec<u8>, usize) {
        let mut max_bw = 0;
        let mut best_bw_repeat_unit = vec![read_bases[offset]];

        for str in 1..=Self::MAX_STR_UNIT_LENGTH {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset + 1).checked_sub(str).is_none() {
                break;
            };

            // get backward repeat unit and # repeats
            max_bw = VariantContextUtils::find_number_of_repetitions_main(
                read_bases,
                offset + 1 - str,
                str,
                read_bases,
                0,
                offset + 1,
                false,
            );
            if max_bw > 1 {
                best_bw_repeat_unit = read_bases[(offset + 1 - str)..(offset + 1)].to_vec();
                break;
            };
        }

        let mut best_repeat_unit = best_bw_repeat_unit;
        let mut max_rl = max_bw;

        if offset < read_bases.len() - 1 {
            let mut best_fw_repeat_unit = vec![read_bases[offset + 1]];
            let mut max_fw = 0;

            for str in 1..=Self::MAX_STR_UNIT_LENGTH {
                // fix repeat unit length
                //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
                if offset + str + 1 > read_bases.len() {
                    break;
                };

                // get forward repeat unit and # repeats
                max_fw = VariantContextUtils::find_number_of_repetitions_main(
                    read_bases,
                    offset + 1,
                    str,
                    read_bases,
                    offset + 1,
                    read_bases.len() - offset,
                    true,
                );
                if max_fw > 1 {
                    best_fw_repeat_unit = read_bases[(offset + 1)..(offset + str + 1)].to_vec();
                    break;
                };
            }

            // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
            if best_fw_repeat_unit == best_repeat_unit {
                max_rl = max_bw + max_fw;
            } else {
                // tandem repeat starting forward from current offset.
                // It could be the case that best BW unit was different from FW unit, but that BW still contains FW unit.
                // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
                // but correct representation at that place might be (C)4.
                // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
                // representations to total
                let test_string = &read_bases[0..(offset + 1)];
                max_bw = VariantContextUtils::find_number_of_repetitions(
                    &best_fw_repeat_unit,
                    test_string,
                    false,
                );
                max_rl = max_fw + max_bw;
                best_repeat_unit = best_fw_repeat_unit;
            }
        }

        if max_rl > Self::MAX_REPEAT_LENGTH {
            max_rl = Self::MAX_REPEAT_LENGTH;
        };

        return (best_repeat_unit, max_rl);
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */

    fn initialize_pair_hmm(
        &mut self,
        haplotypes: &Vec<Haplotype<SimpleInterval>>,
        per_sample_read_list: &HashMap<usize, Vec<BirdToolRead>>,
    ) {
        let read_max_length = per_sample_read_list
            .values()
            .par_bridge()
            .flat_map(|e| e.par_iter())
            .map(|e| e.len())
            .max()
            .unwrap_or(0);
        let max_haplotype_length: usize = haplotypes
            .par_iter()
            .map(|hap| hap.allele.len())
            .max()
            .unwrap_or(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        self.pair_hmm = PairHMM::initialize(read_max_length, max_haplotype_length);
    }
}

pub struct PairHMMInputScoreImputator {
    constant_gcp: u8,
}

impl PairHMMInputScoreImputator {
    pub fn new(gcp: u8) -> Self {
        PairHMMInputScoreImputator { constant_gcp: gcp }
    }

    pub fn del_open_penalties(&self, read: &BirdToolRead) -> Vec<u8> {
        ReadUtils::get_base_deletion_qualities(read)
    }

    pub fn ins_open_penalties(&self, read: &BirdToolRead) -> Vec<u8> {
        ReadUtils::get_base_insertion_qualities(read)
    }

    pub fn gap_continuation_penalties(&self, read: &BirdToolRead) -> Vec<u8> {
        vec![self.constant_gcp; read.len()]
    }
}
