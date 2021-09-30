use genotype::genotype_builder::{AttributeObject, GenotypesContext};
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use hashlink::LinkedHashMap;
use linked_hash_set::LinkedHashSet;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::VariantContext;
use model::variants::SPAN_DEL_ALLELE;
use rayon::prelude::*;
use reads::alignment_utils::AlignmentUtils;
use std::collections::HashMap;
use std::collections::{HashSet, VecDeque};
use std::ops::Range;
use utils::simple_interval::{Locatable, SimpleInterval};
use utils::vcf_constants::*;

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
        let ref_allele_bases = &ref_allele_bases_whole[1..ref_allele.len()];

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

        let length_difference = test_string_length as i64 - repeat_unit_length as i64;
        if leading_repeats {
            let mut num_repeats = 0;
            // look forward on the test string
            for start in (0..=length_difference).step_by(repeat_unit_length) {
                if Self::test_equal_slice(
                    test_string_full,
                    start as usize + offset_in_test_string_full,
                    repeat_unit_full,
                    offset_in_repeat_unit_full,
                    repeat_unit_length,
                ) {
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
                if Self::test_equal_slice(
                    test_string_full,
                    start as usize + offset_in_test_string_full,
                    repeat_unit_full,
                    offset_in_repeat_unit_full,
                    repeat_unit_length,
                ) {
                    num_repeats += 1;
                } else {
                    return num_repeats;
                }
            }
            return num_repeats;
        }
    }

    /// Tests whether two slices are equal.
    /// left and right offset represent first index to test, and length is the length of the unit
    /// to be tested
    pub fn test_equal_slice<T: Eq + PartialEq>(
        left: &[T],
        left_offset: usize,
        right: &[T],
        right_offset: usize,
        length: usize,
    ) -> bool {
        left[left_offset..(left_offset + length)] == right[right_offset..(right_offset + length)]
    }

    /**
     * Returns a {@link Allele#NO_CALL NO_CALL} allele list provided the ploidy.
     *
     * @param ploidy the required ploidy.
     *
     * @return never {@code null}, but an empty list if {@code ploidy} is equal or less than 0. The returned list
     *   might or might not be mutable.
     */
    pub fn no_call_alleles(ploidy: usize) -> Vec<ByteArrayAllele> {
        let no_call_string = ".";
        vec![ByteArrayAllele::new(no_call_string.as_bytes(), false); ploidy]
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name.
     * simpleMerge does not verify any more unique sample names EVEN if genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE. One should use
     * SampleUtils.verifyUniqueSamplesNames to check that before using simpleMerge.
     *
     * For more information on this method see: http://www.thedistractionnetwork.com/programmer-problem/
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param filteredAreUncalled       are filtered records uncalled?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    pub fn simple_merge(
        unsorted_vcs: Vec<VariantContext>,
        priority_list_of_vcs: Option<Vec<String>>,
        original_num_of_vcs: usize,
        filtered_record_merge_type: FilteredRecordMergeType,
        genotype_merge_options: GenotypeMergeType,
        filtered_are_uncalled: bool,
    ) -> Option<VariantContext> {
        if unsorted_vcs.is_empty() {
            return None;
        }

        if priority_list_of_vcs.is_some() {
            if priority_list_of_vcs.as_ref().unwrap().len() != original_num_of_vcs {
                panic!("the number of the original VariantContexts must be the same as the number of VariantContexts in the priority list")
            };
        };

        let pre_filtered_vcs = Self::sort_variant_contexts_by_priority(
            unsorted_vcs,
            priority_list_of_vcs,
            &genotype_merge_options,
        );

        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        let mut VCs: Vec<VariantContext> = pre_filtered_vcs
            .into_par_iter()
            .filter(|vc| !filtered_are_uncalled || vc.is_not_filtered())
            .collect::<Vec<VariantContext>>();

        if VCs.is_empty() {
            return None;
        }

        // establish the baseline info from the first VC
        let VCs_count = VCs.len();
        // let mut first: &mut VariantContext;
        // let name: &String;
        let ref_allele = Self::determine_reference_allele(&VCs, None).cloned();
        debug!("Proposed ref allele: {:?}", &ref_allele);

        let mut alleles = LinkedHashSet::new();
        let mut filters = HashSet::new();
        let mut attributes: LinkedHashMap<String, AttributeObject> = LinkedHashMap::new();
        let mut inconsistent_attributes = LinkedHashSet::new();
        let mut variant_sources = LinkedHashSet::new(); // contains the set of sources we found in our set of VCs that are variant
                                                        // let mut rs_IDs = HashSet::new(); // most of the time there's one id

        let mut longest_vc = VCs[0].loc.clone();
        let mut depth = 0.0;
        let mut log10_p_error = 1.0;
        let mut any_vc_had_filters_applied = false;
        let mut genotypes = GenotypesContext::create(10);

        // counting the number of filtered and variant VCs
        let mut n_filtered = 0;

        // cycle through and add info from the other VCs, making sure the loc/reference matches
        for vc in VCs.iter_mut() {
            assert!(
                longest_vc.get_start() == vc.loc.get_start(),
                "BUG: Attempting to merge VariantContexts with different start sites: first = {:?}, second = {:?}", longest_vc, vc.loc
            );
            if vc.loc.size() > longest_vc.size() {
                longest_vc = vc.loc.clone();
            };

            n_filtered += if vc.is_filtered() { 1 } else { 0 };
            if vc.is_variant() {
                variant_sources.insert(vc.source.clone());
            };
            let allele_mapping =
                Self::resolve_incompatible_alleles(ref_allele.as_ref().unwrap(), &vc);
            alleles.extend(allele_mapping.values().into_iter().cloned());

            Self::merge_genotypes(
                &mut genotypes,
                &vc,
                allele_mapping,
                genotype_merge_options == GenotypeMergeType::Uniquify,
            );
            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if log10_p_error == 1.0 {
                log10_p_error = vc.get_log10_p_error();
            };

            filters.extend(vc.filters.clone());
            any_vc_had_filters_applied = any_vc_had_filters_applied | !vc.filters.is_empty();

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if vc.attributes.contains_key(&*DEPTH_KEY) {
                depth += match vc
                    .attributes
                    .get(&*DEPTH_KEY)
                    .unwrap_or(&AttributeObject::f64(0.0))
                {
                    AttributeObject::f64(val) => *val,
                    _ => panic!(
                        "Incorrect AttributeObject type {:?}",
                        vc.attributes.get(&*DEPTH_KEY)
                    ),
                };
            };

            for (key, value) in vc.attributes.iter() {
                // only output annotations that have the same value in every input VC
                // if we don't like the key already, don't go anywhere
                if !inconsistent_attributes.contains(key) {
                    let already_found = attributes.contains_key(key);
                    if already_found {
                        let bound_value = attributes.get(key).unwrap();
                        let bound_value_is_missing = bound_value.is_empty();
                        if !bound_value_is_missing && bound_value != value {
                            // we found the value but we're inconsistent, put it in the exclude list
                            inconsistent_attributes.insert(key.clone());
                            attributes.remove(key);
                        } else if !already_found || bound_value_is_missing {
                            // no value
                            attributes.insert(key.clone(), value.clone());
                        }
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate),
        // as well as allele-dependent attributes like AC,AF, and AD
        let alleles = alleles.into_iter().collect::<Vec<ByteArrayAllele>>();
        for vc in VCs.iter_mut() {
            if vc.get_alleles_ref().len() == 1 {
                continue;
            };

            if Self::has_pl_incompatibilities(&alleles, &vc.alleles) {
                Self::strip_pls_and_ad(&mut genotypes);
                // TODO: Remove potentially stale AC and AF tags on the VC here
                break;
            }
        }

        if (&filtered_record_merge_type == &FilteredRecordMergeType::KeepIfAnyUnfiltered
            && n_filtered != VCs_count)
            || &filtered_record_merge_type == &FilteredRecordMergeType::KeepUnconditional
        {
            filters.clear();
        }

        if depth > 0.0 {
            attributes.insert(DEPTH_KEY.to_string(), AttributeObject::f64(depth));
        }

        let mut builder =
            VariantContext::build(longest_vc.tid, longest_vc.start, longest_vc.end, alleles);
        builder.genotypes = genotypes;
        builder.log10_p_error = log10_p_error;
        if any_vc_had_filters_applied {
            builder.filters = filters;
        };
        builder.attributes = attributes;

        return Some(builder);
    }

    pub fn strip_pls_and_ad(genotypes: &mut GenotypesContext) {
        genotypes.genotypes_mut().par_iter_mut().for_each(|g| {
            g.pl = Vec::new();
            g.ad = Vec::new();
        })
    }

    fn has_pl_incompatibilities<'a, I>(allele_set_1: I, allele_set_2: I) -> bool
    where
        I: IntoIterator<Item = &'a ByteArrayAllele>,
    {
        let mut it1 = allele_set_1.into_iter().peekable();
        let mut it2 = allele_set_2.into_iter().peekable();

        while it1.peek().is_some() && it2.peek().is_some() {
            let a1 = it1.next().unwrap();
            let a2 = it2.next().unwrap();

            if a1 != a2 {
                return true;
            }
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.peek().is_some() || it2.peek().is_some();
    }

    fn merge_genotypes<'b>(
        merged_genotypes: &'b mut GenotypesContext,
        one_vc: &'b VariantContext,
        mut allele_mapping: AlleleMapper<'b>,
        uniquify_samples: bool,
    ) {
        //TODO: should we add a check for cases when the genotypeMergeOption is REQUIRE_UNIQUE
        for g in one_vc.get_genotypes().genotypes() {
            let name = Self::merged_sample_name(&one_vc.source, &g.sample_name, uniquify_samples);
            if !merged_genotypes.contains_sample(&name) {
                // only add if the name is new
                let mut new_g = g.clone();

                if uniquify_samples || allele_mapping.needs_remapping() {
                    if allele_mapping.needs_remapping() {
                        new_g.alleles = allele_mapping.remap_allele_list(&g.alleles);
                    };
                    new_g.sample_name = name;
                }

                merged_genotypes.add(new_g);
            }
        }
    }

    pub fn merged_sample_name(track_name: &str, sample_name: &str, uniquify: bool) -> String {
        if uniquify {
            return format!("{}.{}", sample_name, track_name);
        } else {
            return format!("{}", sample_name);
        }
    }

    pub fn resolve_incompatible_alleles<'b>(
        ref_allele: &'b ByteArrayAllele,
        vc: &'b VariantContext,
    ) -> AlleleMapper<'b> {
        if ref_allele == vc.get_reference() {
            let mut am = AlleleMapper::default();
            am.set_vc(Some(vc));
            return am;
        } else {
            let mut map = Self::create_allele_mapping(
                ref_allele,
                vc.get_reference_and_index(),
                vc.get_alleles_with_index(),
            );
            map.insert(0, ref_allele.clone());
            let mut am = AlleleMapper::default();
            am.set_map(Some(map));
            return am;
        }
    }

    //TODO as part of a larger refactoring effort {@link #createAlleleMapping} can be merged with {@link ReferenceConfidenceVariantContextMerger#remapAlleles}.
    /**
     * Create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
     *
     * The refAllele is the longest reference allele seen at this start site.
     * So imagine it is:
     * refAllele: ACGTGA
     * myRef:     ACGT
     * myAlt:     A
     *
     * We need to remap all of the alleles in vc to include the extra GA so that
     * myRef => refAllele and myAlt => AGA
     *
     * @param refAllele          the new (extended) reference allele
     * @param inputRef           the reference allele that may need to be extended
     * @param inputAlts          the alternate alleles that may need to be extended
     * @return a non-null mapping of original alleles to new (extended) ones
     */
    pub fn create_allele_mapping<'a>(
        ref_allele: &ByteArrayAllele,
        input_ref: (usize, &ByteArrayAllele),
        input_alleles: Vec<(usize, &'a ByteArrayAllele)>,
    ) -> HashMap<usize, ByteArrayAllele> {
        assert!(
            ref_allele.len() > input_ref.1.len(),
            "BUG: Input ref is longer than ref_allele"
        );
        let extra_bases = &ref_allele.get_bases()[input_ref.1.len()..];

        let mut map = HashMap::new();
        for (idx, a) in input_alleles.into_iter() {
            if Self::is_non_symbolic_extendable_allele(a) {
                let extended = ByteArrayAllele::extend(a, extra_bases);
                map.insert(idx, extended);
            } else if a == &*SPAN_DEL_ALLELE {
                map.insert(idx, a.clone());
            };
        }

        return map;
    }

    fn is_non_symbolic_extendable_allele(allele: &ByteArrayAllele) -> bool {
        return !(allele.is_ref
            || allele.is_symbolic
            || allele == &ByteArrayAllele::new(&[ByteArrayAllele::SPAN_DEL as u8], false));
    }

    pub fn get_size(vc: &VariantContext) -> usize {
        return vc.loc.get_end() - vc.loc.get_start() + 1;
    }

    /**
     * Determines the common reference allele
     *
     * @param VCs    the list of VariantContexts
     * @param loc    if not null, ignore records that do not begin at this start location
     * @return possibly null Allele
     */
    pub fn determine_reference_allele<'b>(
        VCs: &'b Vec<VariantContext>,
        loc: Option<SimpleInterval>,
    ) -> Option<&'b ByteArrayAllele> {
        let mut ref_allele = None;
        for vc in VCs {
            if Self::context_matches_loc(vc, &loc) {
                let my_ref = vc.get_reference();
                ref_allele = Some(Self::determine_reference_allele_from_alleles(
                    ref_allele,
                    Some(my_ref),
                ));
            };
        }

        return ref_allele;
    }

    pub fn determine_reference_allele_from_alleles<'a>(
        ref_1: Option<&'a ByteArrayAllele>,
        ref_2: Option<&'a ByteArrayAllele>,
    ) -> &'a ByteArrayAllele {
        match ref_1 {
            None => match ref_2 {
                None => {
                    panic!("Both reference alleles were None, m8 don't do that.");
                }
                Some(ref_2) => return ref_2,
            },
            Some(ref_1) => match ref_2 {
                None => return ref_1,
                Some(ref_2) => {
                    if ref_1.len() < ref_2.len() {
                        return ref_2;
                    } else if ref_2.len() < ref_1.len() {
                        return ref_1;
                    } else if ref_1.len() == ref_2.len() && ref_1 != ref_2 {
                        panic!("The provided reference alleles do not appear to represent the same position, {:?} vs. {:?}", ref_1, ref_2)
                    } else {
                        return ref_1;
                    }
                }
            },
        }
    }

    pub fn context_matches_loc(vc: &VariantContext, loc: &Option<SimpleInterval>) -> bool {
        match loc {
            None => return true,
            Some(loc) => loc.get_start() == vc.loc.get_start(),
        }
    }

    pub fn sort_variant_contexts_by_priority(
        mut unsorted_vcs: Vec<VariantContext>,
        priority_list_of_vcs: Option<Vec<String>>,
        merge_option: &GenotypeMergeType,
    ) -> Vec<VariantContext> {
        if merge_option == &GenotypeMergeType::Prioritize && priority_list_of_vcs.is_none() {
            panic!("cannot merge calls by priority with a None priority list")
        };

        if priority_list_of_vcs.is_none() || merge_option == &GenotypeMergeType::Unsorted {
            return unsorted_vcs;
        } else {
            let priority_list_of_vcs = priority_list_of_vcs.unwrap();
            unsorted_vcs.par_sort_unstable_by_key(|k| {
                priority_list_of_vcs
                    .iter()
                    .position(|x| x == &k.source)
                    .unwrap()
            });
            return unsorted_vcs;
        }
    }

    /**
     * Trim the alleles in inputVC forward and reverse, as requested
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @param trimForward should we trim up the alleles from the forward direction?
     * @param trimReverse should we trim up the alleles from the reverse direction?
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    pub fn trim_alleles(
        input_vc: &VariantContext,
        trim_forward: bool,
        trim_reverse: bool,
    ) -> VariantContext {
        if input_vc.get_n_alleles() <= 1 || input_vc.get_alleles().iter().any(|a| a.len() == 1) {
            return input_vc.clone();
        }

        let sequences = input_vc
            .get_alleles()
            .iter()
            .filter(|a| !a.is_symbolic)
            .map(|a| a.bases.as_slice())
            .collect::<Vec<&[u8]>>();
        let mut ranges = input_vc
            .get_alleles()
            .iter()
            .filter(|a| !a.is_symbolic)
            .map(|a| 0..a.len() as i32)
            .collect::<Vec<Range<i32>>>();

        let shifts = AlignmentUtils::normalize_alleles(
            &sequences,
            ranges.iter_mut().collect::<Vec<&mut Range<i32>>>(),
            0,
            true,
        );

        let end_trim = shifts.1;
        let start_trim = -shifts.0;

        let empty_allele = ranges.iter().any(|r| r.len() == 0);
        let restore_one_base_at_end = empty_allele && start_trim == 0;
        let restore_one_base_at_start = empty_allele && start_trim > 0;

        // if the end trimming consumed all the bases, leave one base
        let end_bases_to_clip = if restore_one_base_at_end {
            end_trim - 1
        } else {
            end_trim
        };
        let start_bases_to_clip = if restore_one_base_at_start {
            start_trim - 1
        } else {
            start_trim
        };

        return Self::do_trim_alleles(
            input_vc,
            if trim_forward { start_bases_to_clip } else { 0 } - 1,
            if trim_reverse { end_bases_to_clip } else { 0 },
        );
    }

    /**
     * Trim up alleles in inputVC, cutting out all bases up to fwdTrimEnd inclusive and
     * the last revTrim bases from the end
     *
     * @param inputVC a non-null input VC
     * @param fwdTrimEnd bases up to this index (can be -1) will be removed from the start of all alleles
     * @param revTrim the last revTrim bases of each allele will be clipped off as well
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    fn do_trim_alleles<'b>(
        input_vc: &'b VariantContext,
        fwd_trim_end: i32,
        rev_trim: i32,
    ) -> VariantContext {
        if fwd_trim_end == -1 && rev_trim == 0 {
            // nothing to do
            return input_vc.clone();
        };

        let mut alleles = Vec::new();
        let mut original_to_trimmed_allele_map = HashMap::new(); // usize to usize

        let mut allele_index = 0;
        for (i, a) in input_vc.get_alleles().iter().enumerate() {
            if a.is_symbolic {
                alleles.push(a.clone());
                original_to_trimmed_allele_map.insert(a, allele_index);
                allele_index += 1;
            } else {
                // get bases for current allele and create a new one with trimmed bases
                let new_bases =
                    &a.bases[(fwd_trim_end + 1) as usize..(a.len() - rev_trim as usize)];
                let trimmed_allele = ByteArrayAllele::new(new_bases, a.is_ref);
                alleles.push(trimmed_allele);
                original_to_trimmed_allele_map.insert(a, allele_index);
                allele_index += 1;
            }
        }

        // let mut allele_mapper = AlleleMapper::default();
        // allele_mapper.
        let genotypes = Self::update_genotypes_with_mapped_alleles(
            input_vc.get_genotypes(),
            original_to_trimmed_allele_map,
            &alleles,
        );

        let start = input_vc.loc.get_start() + (fwd_trim_end + 1) as usize;
        let mut builder = VariantContext::build_from_vc(input_vc);
        builder.loc.start = start;
        builder.loc.end = start + alleles[0].len() - 1;
        builder.alleles = alleles;
        builder.genotypes = genotypes;
        return builder;
    }

    fn update_genotypes_with_mapped_alleles<'b>(
        original_genotypes: &'b GenotypesContext,
        allele_mapper: HashMap<&'b ByteArrayAllele, usize>,
        new_alleles: &'b Vec<ByteArrayAllele>,
    ) -> GenotypesContext {
        let mut updated_genotypes = GenotypesContext::create(original_genotypes.len());

        for genotype in original_genotypes.genotypes() {
            let updated_alleles = genotype
                .alleles
                .par_iter()
                .map(|a| {
                    let new_allele_index = allele_mapper.get(&a);
                    match new_allele_index {
                        None => a.clone(),
                        Some(new_allele_index) => new_alleles[*new_allele_index].clone(),
                    }
                })
                .collect::<Vec<ByteArrayAllele>>();

            let mut genotype_builder = genotype.clone();
            genotype_builder.alleles = updated_alleles;
            updated_genotypes.add(genotype_builder);
        }

        return updated_genotypes;
    }

    /**
     * Trim the alleles in inputVC from the reverse direction
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @return a non-null VariantContext (may be == to inputVC) with alleles trimmed up
     */
    pub fn reverse_trim_alleles(input_vc: &VariantContext) -> VariantContext {
        Self::trim_alleles(input_vc, false, true)
    }

    /**
     * For testing purposes only.  Create a site-only VariantContext at contig:start containing alleles
     *
     * @param name the name of the VC
     * @param contig the contig for the VC
     * @param start the start of the VC
     * @param alleleStrings a non-null, non-empty list of strings for the alleles.  The first will be the ref allele, and others the
     *                      alt.  Will compute the stop of the VC from the length of the reference allele
     * @return a non-null VariantContext
     */
    pub fn make_from_alleles(
        name: String,
        contig: usize,
        start: usize,
        allele_strings: Vec<&str>,
    ) -> VariantContext {
        let mut alleles = Vec::new();
        let length = allele_strings[0].len();

        let mut first = true;
        for allele_string in allele_strings {
            alleles.push(ByteArrayAllele::new(allele_string.as_bytes(), first));
            first = false;
        }

        let mut b = VariantContext::build(contig, start, start + length - 1, alleles);
        return b;
    }
}

#[derive(Debug, Clone, Ord, PartialOrd, PartialEq, Eq)]
pub enum GenotypeMergeType {
    /**
     * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
     */
    Uniquify,
    /**
     * Take genotypes in priority order (see the priority argument).
     */
    Prioritize,
    /**
     * Take the genotypes in any order.
     */
    Unsorted,
    /**
     * Require that all samples/genotypes be unique between all inputs.
     */
    RequireUnique,
}

#[derive(Debug, Clone, Eq, Ord, PartialOrd, PartialEq)]
pub enum FilteredRecordMergeType {
    /**
     * Union - leaves the record if any record is unfiltered.
     */
    KeepIfAnyUnfiltered,
    /**
     * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
     */
    KeepIfAllUnfiltered,
    /**
     * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
     */
    KeepUnconditional,
}

pub struct AlleleMapper<'b> {
    vc: Option<&'b VariantContext>,
    map: Option<HashMap<usize, ByteArrayAllele>>,
}

impl<'b> AlleleMapper<'b> {
    pub fn default() -> Self {
        Self {
            vc: None,
            map: None,
        }
    }

    pub fn set_vc(&mut self, vc: Option<&'b VariantContext>) {
        self.vc = vc;
    }

    pub fn set_map(&mut self, map: Option<HashMap<usize, ByteArrayAllele>>) {
        self.map = map;
    }

    pub fn needs_remapping(&self) -> bool {
        self.map.is_some()
    }

    pub fn values(&self) -> Vec<&ByteArrayAllele> {
        if self.map.is_some() {
            self.map
                .as_ref()
                .unwrap()
                .values()
                .collect::<Vec<&ByteArrayAllele>>()
        } else {
            self.vc.as_ref().unwrap().get_alleles_ref()
        }
    }

    pub fn remap_allele(&'b self, idx: usize, a: &'b ByteArrayAllele) -> &'b ByteArrayAllele {
        match &self.map {
            None => return a,
            Some(map) => {
                if map.contains_key(&idx) {
                    return map.get(&idx).unwrap();
                } else {
                    return a;
                }
            }
        }
    }

    pub fn remap_allele_list(&self, alleles: &Vec<ByteArrayAllele>) -> Vec<ByteArrayAllele> {
        return alleles
            .par_iter()
            .enumerate()
            .map(|(idx, a)| self.remap_allele(idx, a).clone())
            .collect::<Vec<ByteArrayAllele>>();
    }
}
