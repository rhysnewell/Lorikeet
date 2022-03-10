use annotator::variant_annotation::Annotation;
use annotator::variant_annotator_engine::VariantAnnotationEngine;
use assembly::assembly_based_caller_utils::AssemblyBasedCallerUtils;
use genotype::genotype_builder::{Genotype, GenotypesContext};
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use genotype::genotyping_engine::GenotypingEngine;
use haplotype::called_haplotypes::CalledHaplotypes;
use haplotype::event_map::EventMap;
use haplotype::haplotype::Haplotype;
use haplotype::homogenous_ploidy_model::HomogeneousPloidyModel;
use haplotype::independent_samples_genotype_model::IndependentSamplesGenotypesModel;
use hashlink::linked_hash_map::LinkedHashMap;
use hashlink::LinkedHashSet;
use model::allele_likelihoods::AlleleLikelihoods;
use model::allele_list::AlleleList;
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::variant_context::VariantContext;
use model::variant_context_utils::VariantContextUtils;
use model::variants::SPAN_DEL_ALLELE;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use reference::reference_reader::ReferenceReader;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use utils::simple_interval::{Locatable, SimpleInterval};
use utils::errors::BirdToolError;

#[derive(Debug, Clone)]
pub struct HaplotypeCallerGenotypingEngine {
    genotyping_engine: GenotypingEngine,
    genotyping_model: IndependentSamplesGenotypesModel,
    ploidy_model: HomogeneousPloidyModel,
    snp_heterozygosity: f64,
    indel_heterozygosity: f64,
    max_genotype_count_to_enumerate: usize,
    practical_allele_count_for_ploidy: HashMap<usize, usize>,
    do_physical_phasing: bool,
}

impl HaplotypeCallerGenotypingEngine {
    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param doAlleleSpecificCalcs Whether the AS_QUAL key should be calculated and added to newly genotyped variants.
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} is {@code null}.
     */
    pub fn new(
        args: &clap::ArgMatches,
        samples: Vec<String>,
        do_physical_phasing: bool,
        sample_ploidy: usize,
        // apply_bqd: bool, This is a DRAGEN-GATK param, I ain't dealing with that
    ) -> Self {
        let genotyping_engine = GenotypingEngine::make(args, samples.clone(), false, sample_ploidy);
        Self {
            genotyping_engine,
            do_physical_phasing,
            genotyping_model: IndependentSamplesGenotypesModel::default(),
            ploidy_model: HomogeneousPloidyModel::new(samples, sample_ploidy),
            max_genotype_count_to_enumerate: 1024,
            snp_heterozygosity: args
                .value_of("snp-heterozygosity")
                .unwrap()
                .parse()
                .unwrap(),
            indel_heterozygosity: args
                .value_of("indel-heterozygosity")
                .unwrap()
                .parse()
                .unwrap(),
            practical_allele_count_for_ploidy: HashMap::new(),
        }
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     *
     * @param haplotypes                             Haplotypes to assign likelihoods to
     * @param readLikelihoods                        Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     * @param givenAlleles                Alleles to genotype
     * @param emitReferenceConfidence whether we should add a &lt;NON_REF&gt; alternative allele to the result variation contexts.
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @param withBamOut whether to annotate reads in readLikelihoods for future writing to bamout
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     */
    pub fn assign_genotype_likelihoods<'b>(
        &mut self,
        mut haplotypes: LinkedHashSet<Haplotype<SimpleInterval>>,
        mut read_likelihoods: AlleleLikelihoods<Haplotype<SimpleInterval>>,
        per_sample_filtered_read_list: HashMap<usize, Vec<BirdToolRead>>,
        ref_bases: &'b [u8],
        ref_loc: &'b SimpleInterval,
        active_region_window: &'b SimpleInterval,
        // tracker: FeatureContext,
        given_alleles: Vec<VariantContext>,
        emit_reference_confidence: bool,
        max_mnp_distance: usize,
        header: &'b [String],
        ploidy: usize,
        args: &'b clap::ArgMatches,
        reference_reader: &'b ReferenceReader,
        stand_min_confidence: f64,
        // with_bam_out: bool,
    ) -> Result<CalledHaplotypes, BirdToolError> {
        let mut haplotypes = haplotypes.into_iter().collect::<Vec<Haplotype<SimpleInterval>>>();
        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        let mut start_pos_key_set = match EventMap::build_event_maps_for_haplotypes(
            &mut haplotypes,
            ref_bases,
            &ref_loc,
            max_mnp_distance,
        ) {
            Ok(result) => result,
            Err(error) => return Err(error)
        };

        // Walk along each position in the key set and create each event to be outputted
        let mut called_haplotypes = HashSet::new();
        let mut return_calls = Vec::new();
        let mut no_call_alleles = VariantContextUtils::no_call_alleles(ploidy);
        let read_qualifies_for_genotyping_predicate =
            Self::compose_read_qualifies_for_genotyping_predicate();
        debug!("haplotypes at assignment {:?}", &haplotypes);

        for loc in start_pos_key_set {
            if loc < active_region_window.get_start() || loc > active_region_window.get_end() {
                continue;
            };

            debug!("Loc {}", loc);
            let events_at_this_loc =
                AssemblyBasedCallerUtils::get_variant_contexts_from_active_haplotypes(
                    loc,
                    &haplotypes,
                    !args.is_present("disable-spanning-event-genotyping"),
                );
            debug!("events at this loc {:?}", &events_at_this_loc);
            let events_at_this_loc_with_span_dels_replaced = Self::replace_span_dels(
                events_at_this_loc,
                &ByteArrayAllele::new(
                    &ref_bases[loc - ref_loc.start..(loc - ref_loc.start + 1)],
                    true,
                ),
                loc,
            );
            debug!(
                "events at this loc replaced {:?}",
                &events_at_this_loc_with_span_dels_replaced
            );

            let mut merged_vc = AssemblyBasedCallerUtils::make_merged_variant_context(
                events_at_this_loc_with_span_dels_replaced,
            );
            debug!("merged vc {:?}", &merged_vc);

            match merged_vc {
                None => continue,
                Some(mut merged_vc) => {
                    let merged_alleles_list_size_before_possible_trimming =
                        merged_vc.get_alleles_ref().len();

                    let mut allele_mapper = AssemblyBasedCallerUtils::create_allele_mapper(
                        &merged_vc,
                        loc,
                        &haplotypes,
                        !args.is_present("disable-spanning-event-genotyping"),
                    );

                    debug!("allele mapper {:?}", &allele_mapper);

                    debug!(
                        "Genotyping event at {} with alleles {:?} and genotypes {:?}",
                        loc,
                        merged_vc.get_alleles_ref(),
                        merged_vc.get_genotypes()
                    );
                    self.remove_alt_alleles_if_too_many_genotypes(
                        ploidy,
                        &mut allele_mapper,
                        &mut merged_vc,
                    );

                    debug!("Alleles in read likelihoods {:?}", read_likelihoods.alleles);
                    let mut read_allele_likelihoods = read_likelihoods.marginalize(&allele_mapper);
                    let mut variant_calling_relevant_overlap = SimpleInterval::new(
                        merged_vc.loc.tid,
                        merged_vc.loc.start,
                        merged_vc.loc.end,
                    )
                    .expand_within_contig(
                        args.value_of("allele-informative-reads-overlap-margin")
                            .unwrap()
                            .parse()
                            .unwrap(),
                        *reference_reader
                            .target_lens
                            .get(&merged_vc.loc.tid)
                            .unwrap() as usize,
                    );

                    // We want to retain evidence that overlaps within its softclipping edges.
                    read_allele_likelihoods.retain_evidence(
                        &read_qualifies_for_genotyping_predicate,
                        &variant_calling_relevant_overlap,
                    );

                    read_allele_likelihoods
                        .set_variant_calling_subset_used(&variant_calling_relevant_overlap);

                    // TODO: sample contamination downsampling occurs here. Won't worry about this for nmow
                    //      as it would require a clone of read_likelihoods
                    debug!(
                        "======================================================================="
                    );
                    debug!(
                        "Event at: {:?} with {} reads and {} disqualified",
                        &merged_vc.loc,
                        read_allele_likelihoods.evidence_count(),
                        read_allele_likelihoods
                            .filtered_evidence_by_sample_index
                            .get(&0)
                            .unwrap()
                            .len()
                    );
                    debug!("Genotypes {:?}", &merged_vc.genotypes);
                    debug!(
                        "======================================================================="
                    );

                    if emit_reference_confidence {
                        // TODO: Deletes alleles and replaces with symbolic non ref?
                        // Not sure we care about this
                    }

                    let mut genotypes = self.calculate_gls_for_this_event(
                        &read_allele_likelihoods,
                        &merged_vc,
                        &no_call_alleles,
                        ref_bases,
                        loc - ref_loc.get_start(),
                    );
                    debug!("New genotypes {:?}", &genotypes);

                    // TODO: Some extra DRAGEN parameterization is possible here
                    let mut gpc = self.resolve_genotype_prior_calculator(
                        loc - ref_loc.get_start() + 1,
                        self.snp_heterozygosity,
                        self.indel_heterozygosity,
                    );

                    let mut variant_context_builder = VariantContext::build_from_vc(&merged_vc);
                    variant_context_builder.genotypes = genotypes;
                    debug!(
                        "Variant context allele values {:?}",
                        &variant_context_builder.alleles
                    );
                    let mut call = self.genotyping_engine.calculate_genotypes(
                        variant_context_builder,
                        self.ploidy_model.ploidy,
                        &gpc,
                        &given_alleles,
                        stand_min_confidence,
                    );

                    match call {
                        None => continue, // pass,
                        Some(mut call) => {
                            // read_allele_likelihoods = read_likelihoods.marginalize(&allele_mapper);
                            // read_allele_likelihoods.retain_evidence(
                            //     &Self::compose_read_qualifies_for_genotyping_predicate(),
                            //     &variant_calling_relevant_overlap,
                            // );


                            let mut annotated_call;
                            if call.get_alleles().len() != read_allele_likelihoods.number_of_alleles() {
                                let mut read_allele_likelihoods = AlleleLikelihoods::consume_likelihoods(
                                    AlleleList::new(&call.alleles),
                                    read_allele_likelihoods
                                );

                                let non_ref_index = if call.get_alternate_alleles_with_index().len() > 0 {
                                    Some(call.get_alternate_alleles_with_index()[0].0)
                                } else {
                                    None
                                };
                                read_allele_likelihoods.update_non_ref_allele_likelihoods(
                                    AlleleList::new(&call.alleles),
                                    None
                                );

                                // Skim the filtered map based on the location so that we do not add filtered read that are going to be removed
                                // right after a few lines of code below.
                                debug!("Called allele {:?}", &call.alleles);
                                debug!("Called genotypes {:?}", &call.genotypes);
                                let overlapping_filtered_reads = Self::overlapping_filtered_reads(
                                    &per_sample_filtered_read_list,
                                    variant_calling_relevant_overlap,
                                );
                                debug!("Overlapping filtered reads {:?}", overlapping_filtered_reads.iter().map(|(_, v)| v.len()).collect::<Vec<usize>>());
                                read_allele_likelihoods.add_evidence(overlapping_filtered_reads, 0.0);

                                annotated_call = self.make_annotated_call(
                                    merged_alleles_list_size_before_possible_trimming,
                                    &mut read_allele_likelihoods,
                                    &mut call,
                                );
                            } else {

                                // Skim the filtered map based on the location so that we do not add filtered read that are going to be removed
                                // right after a few lines of code below.
                                debug!("Called allele {:?}", &call.alleles);
                                debug!("Called genotypes {:?}", &call.genotypes);
                                let overlapping_filtered_reads = Self::overlapping_filtered_reads(
                                    &per_sample_filtered_read_list,
                                    variant_calling_relevant_overlap,
                                );
                                debug!("Overlapping filtered reads {:?}", overlapping_filtered_reads.iter().map(|(_, v)| v.len()).collect::<Vec<usize>>());
                                read_allele_likelihoods.add_evidence(overlapping_filtered_reads, 0.0);

                                annotated_call = self.make_annotated_call(
                                    merged_alleles_list_size_before_possible_trimming,
                                    &mut read_allele_likelihoods,
                                    &mut call,
                                );
                            };

                            debug!("Annotated call {:?}", &annotated_call);

                            if annotated_call
                                .get_genotypes()
                                .genotypes()
                                .iter()
                                .map(|g| g.dp - g.ad[0])
                                .sum::<i64>()
                                >= 2
                            {
                                // at least two supporting reads
                                debug!(">= 2 supporting reads");
                                return_calls.push(annotated_call);
                                call.alleles
                                    .into_iter()
                                    .enumerate()
                                    .map(|(idx, a)| allele_mapper.remove(&idx))
                                    .for_each(|a| {
                                        match a {
                                            None => {
                                                // do nothing
                                            }
                                            Some(a) => {
                                                called_haplotypes.extend(a.into_iter().cloned())
                                            }
                                        }
                                    });
                            }
                        }
                    }
                }
            }
        }

        debug!(
            "Potential return calls {:?} and called haplotypes {:?}",
            &return_calls, &called_haplotypes
        );
        let phased_calls = if self.do_physical_phasing {
            AssemblyBasedCallerUtils::phase_calls(return_calls, &called_haplotypes)
        } else {
            return_calls
        };

        return Ok(CalledHaplotypes::new(phased_calls, called_haplotypes));
    }

    fn overlapping_filtered_reads(
        per_sample_filtered_read_list: &HashMap<usize, Vec<BirdToolRead>>,
        loc: SimpleInterval,
    ) -> HashMap<usize, Vec<BirdToolRead>> {
        let mut overlapping_filtered_reads =
            HashMap::with_capacity(per_sample_filtered_read_list.len());

        for (sample_index, original_list) in per_sample_filtered_read_list {
            if original_list.is_empty() {
                continue;
            };
            let new_list = original_list
                .into_iter()
                .filter(|r| r.overlaps(&loc))
                .map(|r| r.clone())
                .collect::<Vec<BirdToolRead>>();

            if !new_list.is_empty() {
                overlapping_filtered_reads.insert(*sample_index, new_list);
            }
        }

        return overlapping_filtered_reads;
    }

    fn make_annotated_call<'b, A: Allele>(
        &self,
        // ref_seq: &[u8],
        // ref_loc: &SimpleInterval,
        // tracker: FeatureContext,
        // reference_reader: &ReferenceReader,
        // samples: &Vec<String>,
        // merged_vc: &VariantContext,
        merged_alleles_list_size_before_possible_trimming: usize,
        read_allele_likelihoods: &mut AlleleLikelihoods<A>,
        call: &'b mut VariantContext,
        // annotation_engine: VariantAnnotationEnging
    ) -> VariantContext {
        // let locus = merged_vc.loc.clone();
        // let ref_loc_interval = ref_loc.clone();

        // TODO: This function does a bunch of annotation that I'm not sure we need to worry about
        //       Can revisit if it is causing issues. So we will skipf or not
        debug!("Call length {}", call.get_alleles().len());
        let mut untrimmed_result = VariantAnnotationEngine::annotate_context(
            call,
            read_allele_likelihoods,
            Box::new(|a: &Annotation| true),
        );

        debug!(
            "untrimmed len {} vs {}",
            untrimmed_result.get_alleles().len(),
            merged_alleles_list_size_before_possible_trimming
        );

        if untrimmed_result.get_alleles_ref().len()
            == merged_alleles_list_size_before_possible_trimming
        {
            return untrimmed_result;
        } else {
            return VariantContextUtils::reverse_trim_alleles(&untrimmed_result);
        }
    }

    fn resolve_genotype_prior_calculator(
        &self,
        pos: usize,
        snp_heterozygosity: f64,
        indel_heterozygosity: f64,
    ) -> GenotypePriorCalculator {
        return GenotypePriorCalculator::assuming_hw(
            snp_heterozygosity.log10(),
            indel_heterozygosity.log10(),
            None,
        );
    }

    /**
     * For a particular event described in inputVC, form PL vector for each sample by looking into allele read map and filling likelihood matrix for each allele
     * @param readLikelihoods          Allele map describing mapping from reads to alleles and corresponding likelihoods
     * @param mergedVC               Input VC with event to genotype
     * @return                       GenotypesContext object wrapping genotype objects with PLs
     */
    fn calculate_gls_for_this_event<'b>(
        &'b mut self,
        read_likelihoods: &'b AlleleLikelihoods<Haplotype<SimpleInterval>>,
        merged_vc: &'b VariantContext,
        no_call_alleles: &'b Vec<ByteArrayAllele>,
        padded_reference: &'b [u8],
        offset_for_ref_into_event: usize,
    ) -> GenotypesContext {
        let vc_alleles = &merged_vc.alleles;
        let allele_list: AlleleList<ByteArrayAllele> = if read_likelihoods.number_of_alleles() == vc_alleles.len() {
            read_likelihoods.get_allele_list_byte_array()
        } else {
            AlleleList::new(vc_alleles)
        };
        // let allele_list = read_likelihoods.alleles.clone();
        debug!(
            "Read likelihood {:#?}",
            read_likelihoods.values_by_sample_index
        );
        let likelihoods = self.genotyping_model.calculate_likelihoods(
            &allele_list,
            read_likelihoods.get_allele_list_byte_array(),
            read_likelihoods,
            &self.ploidy_model,
            padded_reference,
            offset_for_ref_into_event,
        );

        debug!("Proposed likelihoods {:#?}", &likelihoods);

        let sample_count = self.genotyping_engine.samples.len();
        let mut result = GenotypesContext::create(sample_count);
        for (s, likelihood) in likelihoods.into_iter().enumerate() {
            let mut genotype_builder = Genotype::build_from_likelihoods(
                self.ploidy_model.ploidy,
                likelihood,
                self.genotyping_engine.samples[s].clone(),
            );
            genotype_builder.alleles = no_call_alleles.clone();
            debug!("Adding genotype {:#?}", &genotype_builder);
            result.add(genotype_builder);
        }

        return result;
    }

    /**
     * If the number of alleles is so high that enumerating all possible genotypes is impractical, as determined by
     * {@link #maxGenotypeCountToEnumerate}, remove alt alleles from the input {@code alleleMapper} that are
     * not well supported by good-scored haplotypes.
     * Otherwise do nothing.
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * After the remove operation, entries in map are guaranteed to have the same relative order as they were in the input map,
     * that is, entries will be only be removed but not not shifted relative to each other.
     *  @param ploidy        ploidy of the sample
     * @param alleleMapper  original allele to haplotype map
     */
    fn remove_alt_alleles_if_too_many_genotypes<'b>(
        &mut self,
        ploidy: usize,
        allele_mapper: &mut LinkedHashMap<usize, Vec<&'b Haplotype<SimpleInterval>>>,
        merged_vc: &mut VariantContext,
    ) {
        let original_allele_count = allele_mapper.len();
        let max_genotype_count_to_enumerate = self.max_genotype_count_to_enumerate;
        let practical_allele_count = self
            .practical_allele_count_for_ploidy
            .entry(ploidy)
            .or_insert_with(|| {
                GenotypeLikelihoodCalculators::compute_max_acceptable_allele_count(
                    ploidy,
                    max_genotype_count_to_enumerate,
                )
            });

        if original_allele_count > *practical_allele_count {
            let mut alleles_to_keep = Self::which_alleles_to_keep_based_on_hap_scores(
                allele_mapper,
                merged_vc,
                *practical_allele_count,
            );
            allele_mapper.retain(|allele, _| alleles_to_keep.contains(&allele));

            debug!(
                "At position {:?} removed alt alleles where ploidy is {} and original allele count \
                is {}, whereas after trimming the allele count becomes {}. Alleles kept are: {:?}",
                merged_vc.loc, ploidy, original_allele_count, practical_allele_count, &alleles_to_keep
            );

            Self::remove_excess_alt_alleles_from_vc(merged_vc, alleles_to_keep)
        }
    }

    /**
     * Returns an VC that is similar to {@code inputVC} in every aspect except that alleles not in {@code allelesToKeep}
     * are removed in the returned VC.
     * @throws IllegalArgumentException if 1) {@code allelesToKeep} is null or contains null elements; or
     *                                     2) {@code allelesToKeep} doesn't contain a reference allele; or
     *                                     3) {@code allelesToKeep} is not a subset of {@code inputVC.getAlleles()}
     */
    fn remove_excess_alt_alleles_from_vc(
        input_vc: &mut VariantContext,
        alleles_to_keep: Vec<usize>,
    ) {
        let input_len = input_vc.alleles.len();
        assert!(
            alleles_to_keep
                .iter()
                .any(|a| input_vc.alleles[*a].is_reference()),
            "alleles to keep doesn't contain reference allele!"
        );
        assert!(
            alleles_to_keep.iter().all(|a| a < &input_len),
            "alleles to keep is not a subset of input VC alleles"
        );

        if input_vc.alleles.len() == alleles_to_keep.len() {
            // do nothing
        } else {
            // input_vc.alleles.retain(|a| alleles_to_keep.contains(&a));
            let mut index = 0;
            input_vc.alleles.retain(|allele| {
                let retain = alleles_to_keep.contains(&index);
                index += 1;
                retain
            });
        }
    }

    /**
     * Returns a list of alleles that is a subset of the key set of input map {@code alleleMapper}.
     * The size of the returned list is min({@code desiredNumOfAlleles}, alleleMapper.size()).
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * Entries in the returned list are guaranteed to have the same relative order as they were in the input map.
     *
     * @param alleleMapper          original allele to haplotype map
     * @param desiredNumOfAlleles   desired allele count, including ref allele
     */
    fn which_alleles_to_keep_based_on_hap_scores<'b>(
        allele_mapper: &mut LinkedHashMap<usize, Vec<&'b Haplotype<SimpleInterval>>>,
        merged_vc: &mut VariantContext,
        desired_num_of_alleles: usize,
    ) -> Vec<usize> {
        if allele_mapper.len() <= desired_num_of_alleles {
            return allele_mapper.keys().map(|a| *a).collect::<Vec<usize>>();
        }

        let mut allele_max_priority_q = BinaryHeap::new();
        for allele in allele_mapper.keys() {
            let mut hap_scores = allele_mapper
                .get(allele)
                .unwrap()
                .iter()
                .map(|h| h.score)
                .collect::<Vec<OrderedFloat<f64>>>();
            hap_scores.sort();
            let highest_score = if hap_scores.len() > 0 {
                hap_scores[hap_scores.len() - 1].into()
            } else {
                f64::NEG_INFINITY
            };

            let second_highest_score = if hap_scores.len() > 1 {
                hap_scores[hap_scores.len() - 2].into()
            } else {
                f64::NEG_INFINITY
            };

            allele_max_priority_q.push(AlleleScoredByHaplotype::new(
                &merged_vc.alleles[*allele],
                highest_score,
                second_highest_score,
                *allele,
            ));
        }

        let mut alleles_to_retain = LinkedHashSet::new();
        let mut current_allele;
        while alleles_to_retain.len() < desired_num_of_alleles && allele_max_priority_q.len() > 0 {
            current_allele = allele_max_priority_q.pop().unwrap().get_allele_index();
            alleles_to_retain.insert(current_allele);
        }

        return allele_mapper
            .keys()
            .filter(|a| alleles_to_retain.contains(a))
            .map(|a| *a)
            .collect::<Vec<usize>>();
    }

    fn replace_span_dels(
        events_at_this_loc: Vec<&VariantContext>,
        ref_allele: &ByteArrayAllele,
        loc: usize,
    ) -> Vec<VariantContext> {
        return events_at_this_loc
            .into_iter()
            .map(|vc| Self::replace_with_span_del_vc(vc.clone(), ref_allele, loc))
            .collect::<Vec<VariantContext>>();
    }

    fn replace_with_span_del_vc(
        variant_context: VariantContext,
        ref_allele: &ByteArrayAllele,
        loc: usize,
    ) -> VariantContext {
        if variant_context.loc.get_start() == loc {
            return variant_context;
        } else {
            let mut builder = VariantContext::build_from_vc(&variant_context);
            builder.loc.start = loc;
            builder.loc.end = loc;
            builder.alleles = vec![ref_allele.clone(), SPAN_DEL_ALLELE.clone()];
            return builder;
        }
    }

    /**
     * Composes the appropriate test to determine if a read is to be retained for evidence/likelihood calculation for variants
     * located in a target region.
     * @param hcArgs configuration that may affect the criteria use to retain or filter-out reads.
     * @return never {@code null}.
     */
    pub fn compose_read_qualifies_for_genotyping_predicate(
    ) -> Box<dyn Fn(&BirdToolRead, &SimpleInterval) -> bool> {
        // TODO: DRAGEN has a check here usign args
        return Box::new(|read: &BirdToolRead, target: &SimpleInterval| -> bool {
            // NOTE: we must make this comparison in target -> read order because occasionally realignment/assembly produces
            // reads that consume no reference bases and this can cause them to overlap adjacent
            target.overlaps(read)
        });
    }
}

/**
 * A utility class that provides ordering information, given best and second best haplotype scores.
 * If there's a tie between the two alleles when comparing their best haplotype score, the second best haplotype score
 * is used for breaking the tie. In the case that one allele doesn't have a second best allele, i.e. it has only one
 * supportive haplotype, its second best score is set as {@link Double#NEGATIVE_INFINITY}.
 * In the extremely unlikely cases that two alleles, having the same best haplotype score, neither have a second
 * best haplotype score, or the same second best haplotype score, the order is exactly the same as determined by
 * {@link Allele#compareTo(Allele)}.
 */
#[derive(Debug, Clone)]
struct AlleleScoredByHaplotype<'a> {
    allele: &'a ByteArrayAllele,
    best_haplotype_score: f64,
    second_best_haplotype_score: f64,
    allele_index: usize,
}

impl<'a> AlleleScoredByHaplotype<'a> {
    pub fn new(
        allele: &'a ByteArrayAllele,
        best_haplotype_score: f64,
        second_best_haplotype_score: f64,
        allele_index: usize,
    ) -> AlleleScoredByHaplotype<'a> {
        Self {
            allele,
            best_haplotype_score,
            second_best_haplotype_score,
            allele_index,
        }
    }

    pub fn get_allele(self) -> &'a ByteArrayAllele {
        self.allele
    }

    pub fn get_allele_index(&self) -> usize {
        self.allele_index
    }
}

impl<'a> Ord for AlleleScoredByHaplotype<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.allele.is_ref && !other.allele.is_ref {
            return Ordering::Less;
        } else if !self.allele.is_ref && other.allele.is_ref {
            return Ordering::Greater;
        } else if self.best_haplotype_score > other.best_haplotype_score {
            return Ordering::Less;
        } else if self.best_haplotype_score < other.best_haplotype_score {
            return Ordering::Greater;
        } else if (self.second_best_haplotype_score - other.second_best_haplotype_score).abs()
            > f64::EPSILON
        {
            if self.second_best_haplotype_score > other.second_best_haplotype_score {
                return Ordering::Less;
            } else {
                return Ordering::Greater;
            }
        } else {
            return self.allele.cmp(other.allele);
        }
    }
}

impl<'a> PartialOrd for AlleleScoredByHaplotype<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> PartialEq for AlleleScoredByHaplotype<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.allele == other.allele
    }
}

impl<'a> Eq for AlleleScoredByHaplotype<'a> {}
