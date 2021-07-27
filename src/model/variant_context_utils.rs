use model::variant_context::VariantContext;
use rayon::prelude::*;

pub struct VariantContextUtils {}

impl VariantContextUtils {
    /**
     *
     * @param vc
     * @param refBasesStartingAtVCWithoutPad    Ref bases excluding the initial base of the variant context where the alt matches the ref.
     *                                          For example, if the reference sequence is GATCCACCACCAGTCGA and we have a deletion
     *                                          of one STR unit CCA, it is represented as a variant context TCCA -> T, where the 'T' is
     *                                          the padding base.  In this case, {@code refBasesStartingAtVCWithoutPad} is CCACCACCAGTCGA.
     * @return
     */
    pub fn get_num_tandem_repeat_units(
        vc: &mut VariantContext,
        ref_bases_starting_at_vc_without_pad: &[u8],
    ) -> Option<(Vec<usize>, Vec<u8>)> {
        if !vc.is_indel() {
            // only indels are tandem repeats
            return None;
        };

        let ref_allele = vc.get_reference();
        let ref_allele_bases_whole = ref_allele.get_bases();
        let ref_allele_bases = &ref_allele_bases_whole[1..ref_allele.length()];

        let mut repeat_unit = Vec::new();
        let mut lengths = Vec::new();
        for allele in vc.get_alternate_alleles() {
            let allele_bases = allele.get_bases();
            if allele_bases.len() > 1 {
                let result = Self::get_num_tandem_repeat_units_main(
                    ref_allele_bases,
                    &allele_bases[1..allele_bases.len()],
                    ref_bases_starting_at_vc_without_pad,
                );
                match result {
                    Some(result) => {
                        if result.0[0] == 0 || result.0[1] == 0 {
                            return None;
                        }

                        if lengths.is_empty() {
                            lengths.push(result.0[0])
                        };

                        lengths.push(result.0[1]);
                        repeat_unit = result.1.to_vec();
                    }
                    None => return None,
                }
            } else {
                return None;
            }
        }
        if lengths.is_empty() {
            return Some((lengths, repeat_unit));
        } else {
            return None;
        }
    }

    pub fn get_num_tandem_repeat_units_main(
        ref_bases: &[u8],
        alt_bases: &[u8],
        remain_ref_context: &[u8],
    ) -> Option<(Vec<usize>, Vec<u8>)> {
        /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
          Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)2.
        */
        let mut long_b;

        // find first repeat unit based on either ref or alt, whichever is longer
        if alt_bases.len() > ref_bases.len() {
            long_b = alt_bases;
        } else {
            long_b = ref_bases;
        };

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        let repeat_unit_length = Self::find_repeated_substring(long_b);
        if repeat_unit_length > 0 {
            let repeat_unit = long_b[0..repeat_unit_length].to_vec();

            let mut repetition_count = vec![0; 2];
            // look for repetitions forward on the ref bases (i.e. starting at beginning of ref bases)
            let repetitions_in_ref =
                Self::find_number_of_repetitions(&repeat_unit, ref_bases, true);
            let mut extended_ref = ref_bases.to_vec();
            extended_ref.par_extend(remain_ref_context);
            let mut extended_alt = alt_bases.to_vec();
            extended_alt.par_extend(remain_ref_context);
            repetition_count[0] =
                Self::find_number_of_repetitions(&repeat_unit, &extended_ref, true)
                    .checked_sub(repetitions_in_ref)
                    .unwrap_or(0);
            repetition_count[1] =
                Self::find_number_of_repetitions(&repeat_unit, &extended_alt, true)
                    .checked_sub(repetitions_in_ref)
                    .unwrap_or(0);
            return Some((repetition_count, repeat_unit));
        } else {
            return None;
        }
    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    pub fn find_repeated_substring(bases: &[u8]) -> usize {
        for rep_length in 1..=bases.len() {
            let candidate_repeat_unit = &bases[0..rep_length];
            let mut all_bases_match = true;
            for start in rep_length..bases.len() {
                // check that remaining of string is exactly equal to repeat unit
                let base_piece = &bases[start..start + candidate_repeat_unit.len()];
                if candidate_repeat_unit != base_piece {
                    all_bases_match = false;
                    break;
                }
            }

            if all_bases_match {
                return rep_length;
            }
        }

        return 0;
    }

    /**
     * Finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Non-empty substring represented by byte array
     * @param testString             String to test (represented by byte array), may be empty
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     *
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    pub fn find_number_of_repetitions(
        repeat_unit: &[u8],
        test_string: &[u8],
        leading_repeats: bool,
    ) -> usize {
        if test_string.len() == 0 {
            return 0;
        }
        return Self::find_number_of_repetitions_main(
            repeat_unit,
            0,
            repeat_unit.len(),
            test_string,
            0,
            test_string.len(),
            leading_repeats,
        );
    }

    /**
     * Finds number of repetitions a string consists of.
     * Same as {@link #findNumberOfRepetitions} but operates on subarrays of a bigger array to save on copying.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnitFull             Non-empty substring represented by byte array
     * @param offsetInRepeatUnitFull     the offset in repeatUnitFull from which to read the repeat unit
     * @param repeatUnitLength           length of the repeat unit
     * @param testStringFull             string to test (represented by byte array), may be empty
     * @param offsetInTestStringFull     the offset in offsetInRepeatUnitFull from which to read the test string
     * @param testStringLength           length of the test string
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    pub fn find_number_of_repetitions_main(
        repeat_unit_full: &[u8],
        offset_in_repeat_unit_full: usize,
        repeat_unit_length: usize,
        test_string_full: &[u8],
        offset_in_test_string_full: usize,
        test_string_length: usize,
        leading_repeats: bool,
    ) -> usize {
        assert!(
            repeat_unit_length <= repeat_unit_full.len(),
            "repeat unit lengths do not match"
        );

        if test_string_length == 0 {
            return 0;
        }
        assert!(
            test_string_length <= test_string_full.len(),
            "repeat unit lengths do not match"
        );
        assert!(
            test_string_length >= repeat_unit_length,
            "test string too short"
        );
        let length_difference = test_string_length - repeat_unit_length;
        if leading_repeats {
            let mut num_repeats = 0;
            // look forward on the test string
            for start in (0..=length_difference).step_by(repeat_unit_length) {
                if &repeat_unit_full[start + offset_in_repeat_unit_full..]
                    == &test_string_full[start + offset_in_test_string_full..]
                {
                    num_repeats += 1;
                } else {
                    return num_repeats;
                }
            }
            return num_repeats;
        } else {
            // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
            let mut num_repeats = 0;
            // look backward on the test string
            for start in (0..=length_difference).rev().step_by(repeat_unit_length) {
                if &repeat_unit_full[start + offset_in_repeat_unit_full..]
                    == &test_string_full[start + offset_in_test_string_full..]
                {
                    num_repeats += 1;
                } else {
                    return num_repeats;
                }
            }
            return num_repeats;
        }
    }
}
