use assembly::assembly_region::AssemblyRegion;
use assembly::assembly_result_set::AssemblyResultSet;
use genotype::genotype_builder::AttributeObject;
use haplotype::haplotype::Haplotype;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use haplotype::location_and_alleles::LocationsAndAlleles;
use haplotype::reference_confidence_model::ReferenceConfidenceModel;
use hashlink::linked_hash_map::LinkedHashMap;
use itertools::Itertools;
use model::allele_likelihoods::{AlleleLikelihoods, ReadIndexer};
use model::byte_array_allele::{Allele, ByteArrayAllele};
use model::location_and_alleles::LocationAndAlleles;
use model::variant_context::VariantContext;
use model::variant_context_utils::{
    FilteredRecordMergeType, GenotypeMergeType, VariantContextUtils,
};
use model::variants::*;
use pair_hmm::pair_hmm_likelihood_calculation_engine::{
    PCRErrorModel, PairHMMLikelihoodCalculationEngine,
};
use rayon::prelude::*;
use read_error_corrector::nearby_kmer_error_corrector::NearbyKmerErrorCorrector;
use read_threading::abstract_read_threading_graph::AbstractReadThreadingGraph;
use read_threading::read_threading_assembler::ReadThreadingAssembler;
use read_threading::read_threading_graph::ReadThreadingGraph;
use reads::alignment_utils::AlignmentUtils;
use reads::bird_tool_reads::BirdToolRead;
use reads::cigar_utils::CigarUtils;
use reads::read_clipper::ReadClipper;
use reads::read_utils::ReadUtils;
use reference::reference_reader::ReferenceReader;
use rust_htslib::bam::ext::BamRecordExtensions;
use smith_waterman::bindings::{SWOverhangStrategy, SWParameters};
use smith_waterman::smith_waterman_aligner::{NEW_SW_PARAMETERS, STANDARD_NGS};
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use utils::math_utils::MathUtils;
use utils::quality_utils::QualityUtils;
use utils::simple_interval::{Locatable, SimpleInterval};
use utils::vcf_constants::*;

lazy_static! {
    static ref PHASE_01: PhaseGroup = PhaseGroup::new("0|1".to_string(), 1);
    static ref PHASE_10: PhaseGroup = PhaseGroup::new("1|0".to_string(), 0);
    // pub static ref HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY: Box<dyn Fn(&Haplotype<SimpleInterval>) -> i32> = Box::new(
    //         |h: &Haplotype<SimpleInterval>| {
    //             let cigar = &h.cigar;
    //             let reference_term = if h.is_ref() {
    //                 1
    //             } else {
    //                 0
    //             };
    //             let cigar_term = 1 - cigar.0.len() as i32;
    //             return reference_term + cigar_term
    //         }
    //     );
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PhaseGroup {
    description: String,
    alt_allele_index: usize,
}

impl PhaseGroup {
    pub fn new(description: String, alt_allele_index: usize) -> PhaseGroup {
        PhaseGroup {
            description,
            alt_allele_index,
        }
    }

    pub fn get_description(&self) -> String {
        self.description.clone()
    }

    pub fn get_alt_allele_index(&self) -> usize {
        self.alt_allele_index
    }
}

/// Collection of functions used during various assembly stages
///
/// @author Rhys Newell <rhys.newell@.hdr.qut.edu.au>
pub struct AssemblyBasedCallerUtils {}

impl AssemblyBasedCallerUtils {
    const REFERENCE_PADDING_FOR_ASSEMBLY: usize = 500;
    // After trimming to fit the assembly window, throw away read stubs shorter than this length
    // if we don't, the several bases left of reads that end just within the assembly window can
    // get realigned incorrectly.  See https://github.com/broadinstitute/gatk/issues/5060
    pub const MINIMUM_READ_LENGTH_AFTER_TRIMMING: usize = 10;
    pub const NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO: usize = 5;

    pub fn finalize_regions(
        region: &mut AssemblyRegion,
        error_correct_reads: bool,
        dont_use_soft_clipped_bases: bool,
        min_tail_quality: u8,
        correct_overlapping_base_qualities: bool,
        soft_clip_low_quality_ends: bool,
    ) {
        if !region.is_finalized() {
            let min_tail_quality_to_use = if error_correct_reads {
                HaplotypeCallerEngine::MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION as u8
            } else {
                min_tail_quality
            };

            let mut reads_to_use = Vec::new();
            for original_read in region.get_reads().iter() {
                // TODO unclipping soft clips may introduce bases that aren't in the extended region if the unclipped bases
                // TODO include a deletion w.r.t. the reference.  We must remove kmers that occur before the reference haplotype start
                let mut read = if dont_use_soft_clipped_bases
                    || !ReadUtils::has_well_defined_fragment_size(&original_read)
                {
                    ReadClipper::new(original_read).hard_clip_soft_clipped_bases()
                } else {
                    ReadClipper::new(original_read).revert_soft_clipped_bases()
                };

                if read.get_start() <= read.get_end() {
                    read = if read.read.is_unmapped() {
                        read
                    } else {
                        ReadClipper::new(&read).hard_clip_adaptor_sequence()
                    };

                    if !read.is_empty() && read.read.seq_len_from_cigar(false) > 0 {
                        read = ReadClipper::hard_clip_to_region(
                            read,
                            region.get_padded_span().get_start(),
                            region.get_padded_span().get_end(),
                        );

                        if read.get_start() <= read.get_end()
                            && read.len() > 0
                            && read.overlaps(&region.get_padded_span())
                        {
                            // NOTE: here we make a defensive copy of the read if it has not been modified by the above operations
                            // which might only make copies in the case that the read is actually clipped
                            reads_to_use.push(read);
                        }
                    }
                }
            }

            reads_to_use.par_sort_unstable();

            // handle overlapping read pairs from the same fragment
            // if correct_overlapping_base_qualities {
            //
            // }
            region.clear_reads();
            region.add_all(reads_to_use);
            region.set_finalized(true);
        }
    }

    pub fn haplotype_alignment_tiebreaking_priority(
    ) -> Box<dyn Fn(&Haplotype<SimpleInterval>) -> i32> {
        Box::new(|h: &Haplotype<SimpleInterval>| {
            let cigar = &h.cigar;
            let reference_term = if h.is_ref() { 1 } else { 0 };
            let cigar_term = 1 - cigar.0.len() as i32;
            return reference_term + cigar_term;
        })
    }

    pub fn reference_tiebreaking_priority<A: Allele>() -> Box<dyn Fn(&A) -> i32> {
        Box::new(|h: &A| if h.is_reference() { 1 } else { 0 })
    }

    /**
     * Returns a map with the original read indices (sample_index, evidence_index) as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    pub fn realign_reads_to_their_best_haplotype(
        original_read_likelihoods: &AlleleLikelihoods<Haplotype<SimpleInterval>>,
        ref_haplotype: &Haplotype<SimpleInterval>,
        padded_reference_loc: &SimpleInterval,
    ) -> HashMap<ReadIndexer, BirdToolRead> {
        let best_alleles = original_read_likelihoods
            .best_alleles_breaking_ties_main(Self::haplotype_alignment_tiebreaking_priority());

        return best_alleles
            .iter()
            .map(|best_allele| {
                let best_haplotype = &original_read_likelihoods
                    .alleles
                    .get_allele(best_allele.allele_index.unwrap());
                let original_read = &original_read_likelihoods
                    .evidence_by_sample_index
                    .get(&best_allele.sample_index)
                    .unwrap()[best_allele.evidence_index];
                let is_informative = best_allele.is_informative();
                let realigned_read = AlignmentUtils::create_read_aligned_to_ref(
                    original_read,
                    best_haplotype,
                    ref_haplotype,
                    padded_reference_loc.get_start(),
                    is_informative,
                );
                (
                    ReadIndexer::new(best_allele.sample_index, best_allele.evidence_index),
                    realigned_read,
                )
            })
            .collect::<HashMap<ReadIndexer, BirdToolRead>>();
    }

    // /**
    //  *  Modify base qualities when paired reads overlap to account for the possibility of PCR error.
    //  *
    //  *  Overlapping mates provded independent evidence as far as sequencing error is concerned, but their PCR errors
    //  *  are correlated.  The base qualities are thus limited by the sequencing base quality as well as half of the PCR
    //  *  quality.  We use half of the PCR quality because downstream we treat read pairs as independent, and summing two halves
    //  *  effectively gives the PCR quality of the pairs when taken together.
    //  *
    //  * @param reads the list of reads to consider
    //  * @param samplesList   list of samples
    //  * @param readsHeader   bam header of reads' source
    //  * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
    //  * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR
    //  * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR
    //  */
    // TODO: Need to be able to split apart reads by sample name for this to work, not high priority
    // pub fn clean_overlapping_read_pairs(
    //     read: &mut Vec<BirdToolRead>, set_conflicting_to_zero: bool,
    //     half_of_pcr_snv_qual: Option<usize>, half_of_pcr_indel_qual: Option<usize>
    // ) {
    //     for
    // }

    /**
     * High-level function that runs the assembler on the given region's reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     */
    pub fn assemble_reads<'b>(
        mut region: AssemblyRegion,
        given_alleles: &Vec<VariantContext>,
        args: &clap::ArgMatches,
        reference_reader: &'b mut ReferenceReader,
        assembly_engine: &mut ReadThreadingAssembler,
        correct_overlapping_base_qualities: bool,
        sample_names: &Vec<String>,
    ) -> AssemblyResultSet<ReadThreadingGraph> {
        Self::finalize_regions(
            &mut region,
            args.is_present("error-correct-reads"),
            args.is_present("dont-use-soft-clipped-bases"),
            args.value_of("min-base-quality")
                .unwrap()
                .parse::<u8>()
                .unwrap()
                - 1,
            correct_overlapping_base_qualities,
            args.is_present("soft-clip-low-quality-ends"),
        );
        debug!(
            "Assembling {:?} with {} reads:    (with overlap region = {:?})",
            region.get_span(),
            region.get_reads().len(),
            region.get_padded_span()
        );

        let full_reference_with_padding = region
            .get_assembly_region_reference(reference_reader, Self::REFERENCE_PADDING_FOR_ASSEMBLY);
        let padded_reference_loc = Self::get_padded_reference_loc(
            &region,
            Self::REFERENCE_PADDING_FOR_ASSEMBLY,
            &reference_reader,
        );
        debug!("Padded reference location {:?}", &padded_reference_loc);
        let mut ref_haplotype =
            Self::create_reference_haplotype(&region, &padded_reference_loc, reference_reader);

        let mut read_error_corrector;

        if args.is_present("error-correct-reads") {
            read_error_corrector = Some(NearbyKmerErrorCorrector::default(
                args.value_of("kmer-length-for-read-error-correction")
                    .unwrap()
                    .parse::<usize>()
                    .unwrap(),
                HaplotypeCallerEngine::MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
                args.value_of("min-observations-for-kmers-to-be-solid")
                    .unwrap()
                    .parse::<usize>()
                    .unwrap(),
                full_reference_with_padding.as_slice(),
            ))
        } else {
            read_error_corrector = None
        }

        let region_padded_start = region.get_padded_span().get_start();
        let mut assembly_result_set = assembly_engine.run_local_assembly(
            region,
            &mut ref_haplotype,
            full_reference_with_padding,
            padded_reference_loc,
            read_error_corrector,
            sample_names,
            *STANDARD_NGS,
            *NEW_SW_PARAMETERS,
        );

        if !given_alleles.is_empty() {
            Self::add_given_alleles(
                region_padded_start,
                given_alleles,
                args.value_of("max-mnp-distance").unwrap().parse().unwrap(),
                *NEW_SW_PARAMETERS,
                &ref_haplotype,
                &mut assembly_result_set,
            );
        }

        return assembly_result_set;
    }

    fn add_given_alleles<A: AbstractReadThreadingGraph>(
        assembly_region_start: usize,
        given_alleles: &Vec<VariantContext>,
        max_mnp_distance: usize,
        haplotype_to_reference_sw_parameters: SWParameters,
        ref_haplotype: &Haplotype<SimpleInterval>,
        assembly_result_set: &mut AssemblyResultSet<A>,
    ) {
        let active_region_start = ref_haplotype.alignment_start_hap_wrt_ref;
        let mut grouped_by = HashMap::new(); // vcs grouped by start
        assembly_result_set
            .get_variation_events(max_mnp_distance)
            .into_iter()
            .for_each(|vc| {
                let pos = grouped_by.entry(vc.loc.get_start()).or_insert(Vec::new());
                pos.push(vc);
            });
        let mut assembled_variants = grouped_by
            .into_iter()
            .map(|(i, vcs)| (i, Self::make_merged_variant_context(vcs).unwrap()))
            .collect::<HashMap<usize, VariantContext>>();

        for given_vc in given_alleles {
            let mut assembled_haplotypes = assembly_result_set.get_haplotype_list();
            let mut assembled_vc = assembled_variants.get(&given_vc.loc.get_start());
            let given_vc_ref_length = given_vc.get_reference().len();
            let mut unassembled_given_alleles = Vec::new();
            let mut longer_ref;
            match &assembled_vc {
                None => {
                    longer_ref = given_vc.get_reference();
                }
                Some(assembled_vc) => {
                    if given_vc_ref_length > assembled_vc.get_reference().len() {
                        longer_ref = given_vc.get_reference();
                    } else {
                        longer_ref = assembled_vc.get_reference();
                    }
                }
            };

            match assembled_vc {
                None => {
                    unassembled_given_alleles = given_vc
                        .get_alternate_alleles()
                        .into_iter()
                        .map(|vc| vc.clone())
                        .collect::<Vec<ByteArrayAllele>>();
                }
                Some(assembled_vc) => {
                    // map all alleles to the longest common reference
                    let mut assembled_allele_set = HashSet::new();
                    if longer_ref.len() == assembled_vc.get_reference().len() {
                        assembled_allele_set.extend(
                            assembled_vc
                                .get_alternate_alleles()
                                .into_iter()
                                .map(|vc| vc.clone()),
                        );
                    } else {
                        assembled_allele_set
                            .extend(VariantContextUtils::remap_alleles(assembled_vc, longer_ref));
                    };

                    let mut given_allele_set = HashSet::new();
                    if longer_ref.len() == given_vc_ref_length {
                        given_allele_set.extend(
                            given_vc
                                .get_alternate_alleles()
                                .into_iter()
                                .map(|vc| vc.clone()),
                        );
                    } else {
                        given_allele_set
                            .extend(VariantContextUtils::remap_alleles(given_vc, longer_ref));
                    };

                    unassembled_given_alleles = given_allele_set
                        .into_iter()
                        .filter(|a| !assembled_allele_set.contains(a))
                        .collect::<Vec<ByteArrayAllele>>();
                }
            }

            let unassembled_non_symbolic_alleles = unassembled_given_alleles
                .into_iter()
                .filter(|a| {
                    let bases = a.get_bases();
                    !(ByteArrayAllele::would_be_no_call_allele(bases)
                        || ByteArrayAllele::would_be_null_allele(bases)
                        || ByteArrayAllele::would_be_star_allele(bases)
                        || ByteArrayAllele::would_be_symbolic_allele(bases))
                })
                .collect::<Vec<ByteArrayAllele>>();

            // choose the highest-scoring haplotypes along with the reference for building force-calling haplotypes
            let base_haplotypes = if unassembled_non_symbolic_alleles.is_empty() {
                Vec::new()
            } else {
                assembled_haplotypes.sort_by(|h1, h2| {
                    h1.is_reference()
                        .cmp(&h2.is_reference())
                        .then_with(|| h1.score.cmp(&h2.score))
                });
                assembled_haplotypes.reverse();
                assembled_haplotypes
                    .into_iter()
                    .take(Self::NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO)
                    .collect::<Vec<Haplotype<SimpleInterval>>>()
            };

            for given_allele in unassembled_non_symbolic_alleles {
                for base_haplotype in base_haplotypes.iter() {
                    match &base_haplotype.event_map {
                        Some(event) => {
                            if event
                                .get_variant_contexts()
                                .iter()
                                .any(|vc| vc.loc.overlaps(&given_vc.loc))
                            {
                                continue;
                            };
                        }
                        _ => {
                            // pass
                        }
                    };

                    let mut inserted_haplotype = base_haplotype.insert_allele(
                        longer_ref,
                        &given_allele,
                        active_region_start + given_vc.loc.get_start() - assembly_region_start,
                    );

                    if let Some(mut inserted_haplotype) = inserted_haplotype {
                        let cigar = CigarUtils::calculate_cigar(
                            ref_haplotype.get_bases(),
                            inserted_haplotype.get_bases(),
                            SWOverhangStrategy::Indel,
                            &haplotype_to_reference_sw_parameters,
                        );
                        inserted_haplotype.set_cigar(cigar.unwrap().0);
                        inserted_haplotype.set_genome_location(
                            ref_haplotype.get_genome_location().unwrap().clone(),
                        );
                        inserted_haplotype.set_alignment_start_hap_wrt_ref(active_region_start);
                        assembly_result_set.add_haplotype(inserted_haplotype);
                    }
                }
            }
        }

        assembly_result_set.regenerate_variation_events(max_mnp_distance);
    }

    // Contract: the List<Allele> alleles of the resulting VariantContext is the ref allele followed by alt alleles in the
    // same order as in the input vcs
    pub fn make_merged_variant_context(vcs: Vec<VariantContext>) -> Option<VariantContext> {
        if vcs.is_empty() {
            return None;
        }

        let haplotype_sources = vcs
            .iter()
            .map(|vc| vc.source.clone())
            .collect::<Vec<String>>();
        let orignal_num_of_vcs = haplotype_sources.len();

        return VariantContextUtils::simple_merge(
            vcs,
            Some(haplotype_sources),
            orignal_num_of_vcs,
            FilteredRecordMergeType::KeepIfAnyUnfiltered,
            GenotypeMergeType::Prioritize,
            false,
        );
    }

    pub fn get_variant_contexts_from_given_alleles<'b>(
        loc: usize,
        active_alleles_to_genotype: &'b Vec<VariantContext>,
        include_spanning_events: bool,
    ) -> Vec<VariantContext> {
        let mut unique_locations_and_alleles = HashSet::new();
        let mut results = Vec::new();

        let mut given_allele_source_count = 0;
        for given_allele_vc in active_alleles_to_genotype.iter() {
            if given_allele_vc.loc.get_start() <= loc && given_allele_vc.loc.get_end() >= loc {
                if !(include_spanning_events || given_allele_vc.loc.get_start() == loc) {
                    continue;
                }
                let mut allele_count = 0;
                for given_alt_allele in given_allele_vc.get_alternate_alleles() {
                    let allele_set = vec![
                        given_allele_vc.get_reference().clone(),
                        given_alt_allele.clone(),
                    ];

                    let vc_source_name =
                        format!("Comp{}Allele{}", given_allele_source_count, allele_count);
                    // check if this event is already in the list of events due to a repeat in the input alleles track
                    let mut candidate_event_to_add = VariantContext::build_from_vc(given_allele_vc);
                    candidate_event_to_add.add_source(vc_source_name);
                    candidate_event_to_add.add_alleles(allele_set);

                    let location_and_alleles = LocationAndAlleles::new(
                        candidate_event_to_add.loc.get_start(),
                        candidate_event_to_add.get_alleles().clone(),
                    );

                    if !unique_locations_and_alleles.contains(&location_and_alleles) {
                        unique_locations_and_alleles.insert(location_and_alleles);
                        results.push(candidate_event_to_add);
                    }
                    allele_count += 1;
                }
            }
            given_allele_source_count += 1;
        }
        return results;
    }

    /**
     * Returns the list of events discovered in assembled haplotypes that are active at this location. The results will
     * include events that span the current location if includeSpanningEvents is set to true; otherwise it will only
     * include events that have loc as their start position.
     * @param loc The start position we are genotyping
     * @param haplotypes list of active haplotypes at the current location
     * @param includeSpanningEvents If true, will also return events that span loc
     */
    pub fn get_variant_contexts_from_active_haplotypes<'b>(
        loc: usize,
        haplotypes: &'b Vec<Haplotype<SimpleInterval>>,
        include_spanning_event: bool,
    ) -> Vec<&'b VariantContext> {
        let mut results = Vec::new();
        let mut unique_locations_and_alleles = HashSet::new();

        haplotypes
            .iter()
            .flat_map(|h| h.event_map.as_ref().unwrap().get_overlapping_events(loc))
            .filter(|v| (include_spanning_event || v.loc.get_start() == loc))
            .for_each(|v| {
                let location_and_alleles =
                    LocationsAndAlleles::new(v.loc.get_start(), v.get_alleles());
                if !unique_locations_and_alleles.contains(&location_and_alleles) {
                    unique_locations_and_alleles.insert(location_and_alleles);
                    results.push(&*v)
                }
            });

        return results;
    }

    pub fn get_padded_reference_loc(
        region: &AssemblyRegion,
        reference_padding: usize,
        reference_reader: &ReferenceReader,
    ) -> SimpleInterval {
        let pad_left = region
            .get_padded_span()
            .get_start()
            .checked_sub(reference_padding)
            .unwrap_or(0);
        let pad_right = min(
            region.get_padded_span().get_end() + reference_padding,
            reference_reader.get_contig_length(region.get_contig()) as usize - 1, // 0-based inclusive hence -1 to avoid creating size that is +1 bp larger than reference
        );

        return SimpleInterval::new(region.get_contig(), pad_left, pad_right);
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    pub fn create_reference_haplotype<L: Locatable>(
        region: &AssemblyRegion,
        padded_reference_loc: &SimpleInterval,
        reference_reader: &mut ReferenceReader,
    ) -> Haplotype<L> {
        return ReferenceConfidenceModel::create_reference_haplotype(
            region,
            region
                .get_assembly_region_reference(reference_reader, 0)
                .as_slice(),
            padded_reference_loc,
        );
    }

    pub fn create_reference_haplotype_from_bytes<L: Locatable>(
        region: &AssemblyRegion,
        padded_reference_loc: &SimpleInterval,
        reference: &[u8],
    ) -> Haplotype<L> {
        return ReferenceConfidenceModel::create_reference_haplotype(
            region,
            reference,
            padded_reference_loc,
        );
    }

    /**
     * Returns a mapping from Allele in the mergedVC, which represents all of the alleles being genotyped at loc,
     * to a list of Haplotypes that support that allele. If the mergedVC includes a spanning deletion allele, all haplotypes
     * that support spanning deletions will be assigned to that allele in the map.
     *
     * @param mergedVC The merged variant context for the locus, which includes all active alternate alleles merged to a single reference allele
     * @param loc The active locus being genotyped
     * @param haplotypes Haplotypes for the current active region
     * @param emitSpanningDels If true will map spanning events to a * allele instead of reference // TODO add a test for this behavior
     * @return
     */
    pub fn create_allele_mapper<'b>(
        merged_vc: &VariantContext,
        loc: usize,
        haplotypes: &'b Vec<Haplotype<SimpleInterval>>,
        emit_spanning_dels: bool,
    ) -> LinkedHashMap<usize, Vec<&'b Haplotype<SimpleInterval>>> {
        let mut result = LinkedHashMap::new();

        let (ref_index, ref_allele) = merged_vc.get_reference_and_index();
        result.insert(ref_index, Vec::new());

        //Note: we can't use the alleles implied by eventsAtThisLoc because they are not yet merged to a common reference
        //For example, a homopolymer AAAAA reference with a single and double deletion would yield (i) AA* A and (ii) AAA* A
        //in eventsAtThisLoc, when in mergedVC it would yield AAA* AA A
        merged_vc
            .get_alleles()
            .iter()
            .enumerate()
            .filter(|(idx, a)| !a.is_symbolic && !a.is_ref)
            .for_each(|(idx, a)| {
                result.insert(idx, Vec::new());
            });

        for h in haplotypes {
            let spanning_events = h.event_map.as_ref().unwrap().get_overlapping_events(loc);

            if spanning_events.is_empty() {
                let allele_vec = result.entry(ref_index).or_insert(Vec::new());
                allele_vec.push(h);
            } else {
                for spanning_event in spanning_events {
                    if spanning_event.loc.get_start() == loc {
                        // the event starts at the current location
                        if spanning_event.get_reference().len() == ref_allele.len() {
                            // reference allele lengths are equal; we can just use the spanning event's alt allele
                            // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                            let index = merged_vc
                                .alleles
                                .iter()
                                .position(|a| a == spanning_event.get_alternate_alleles()[0]);
                            match index {
                                None => {
                                    // pass
                                }
                                Some(index) => {
                                    if result.contains_key(&index) {
                                        // variant contexts derived from the event map have only one alt allele each, so we can just
                                        // grab the first one (we're not assuming that the sample is biallelic)
                                        let mut allele_vec = result.get_mut(&index).unwrap();
                                        allele_vec.push(h);
                                    }
                                }
                            }
                        } else if spanning_event.get_reference().len()
                            < merged_vc.get_reference().len()
                        {
                            // spanning event has shorter ref allele than merged VC; we need to pad out its alt allele
                            let spanning_event_allele_mapping_to_merge_vc =
                                VariantContextUtils::create_allele_mapping(
                                    merged_vc.get_reference(),
                                    spanning_event.get_reference_and_index(),
                                    spanning_event.get_alleles_with_index(),
                                );

                            let remapped_spanning_event_alt_allele =
                                spanning_event_allele_mapping_to_merge_vc
                                    .get(&spanning_event.get_alternate_alleles_with_index()[0].0)
                                    .unwrap();
                            // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                            let index = merged_vc
                                .alleles
                                .iter()
                                .position(|a| a == remapped_spanning_event_alt_allele);
                            match index {
                                None => {
                                    // pass
                                }
                                Some(index) => {
                                    if result.contains_key(&index) {
                                        // variant contexts derived from the event map have only one alt allele each, so we can just
                                        // grab the first one (we're not assuming that the sample is biallelic)
                                        let mut allele_vec = result.get_mut(&index).unwrap();
                                        allele_vec.push(h);
                                    }
                                }
                            }
                            // if result.contains_key(&remapped_spanning_event_alt_allele) {
                            //     let mut allele_vec = result.get_mut(&remapped_spanning_event_alt_allele).unwrap();
                            //     allele_vec.push(h);
                            // }
                        } else {
                            // the process of creating the merged VC in AssemblyBasedCallerUtils::makeMergedVariantContext should have
                            // already padded out the reference allele, therefore this spanning VC must not be in events at this site
                            // because we're in GGA mode and it's not an allele we want
                            continue;
                        }
                    } else {
                        if emit_spanning_dels {
                            // the event starts prior to the current location, so it's a spanning deletion
                            let index = merged_vc
                                .alleles
                                .iter()
                                .position(|a| a == &*SPAN_DEL_ALLELE);
                            match index {
                                None => {
                                    // I guess pass? We can't exactly insert it into merged_vc
                                    result.get_mut(&ref_index).unwrap().push(h);
                                    break;
                                }
                                Some(index) => {
                                    result.get_mut(&index).unwrap().push(h);
                                    break;
                                }
                            }
                            // if !result.contains_key(&*SPAN_DEL_ALLELE) {
                            //     result.insert(*SPAN_DEL_ALLELE.clone(), Vec::new());
                            // }
                            //
                            // result.get_mut(&*SPAN_DEL_ALLELE).unwrap().push(h);
                            // break
                        } else {
                            result.get_mut(&ref_index).unwrap().push(h);
                            break;
                        }
                    }
                }
            }
        }

        return result;
    }

    pub fn get_alleles_consistent_with_given_alleles<'b>(
        given_alleles: &'b Vec<VariantContext>,
        merged_vc: &'b VariantContext,
    ) -> HashSet<&'b ByteArrayAllele> {
        if given_alleles.is_empty() {
            return HashSet::new();
        }

        let vcs = AssemblyBasedCallerUtils::get_variant_contexts_from_given_alleles(
            merged_vc.loc.get_start(),
            given_alleles,
            false,
        );

        let given_alt_and_ref_alleles_in_original_context = vcs
            .iter()
            .flat_map(|vc| {
                let refr = vc.get_reference();
                let alt = vc.get_alternate_alleles();
                alt.into_iter()
                    .map(|allele| (allele, refr))
                    .collect::<Vec<(&ByteArrayAllele, &ByteArrayAllele)>>()
            })
            .collect::<Vec<(&ByteArrayAllele, &ByteArrayAllele)>>();

        let result = merged_vc
            .get_alternate_alleles()
            .into_iter()
            .map(|allele| (allele, merged_vc.get_reference()))
            .filter(|alt_and_ref| {
                given_alt_and_ref_alleles_in_original_context
                    .iter()
                    .any(|given_alt_and_ref| {
                        AssemblyBasedCallerUtils::alleles_are_consistent(
                            given_alt_and_ref,
                            alt_and_ref,
                        )
                    })
            })
            .map(|alt_and_ref_pair| alt_and_ref_pair.0)
            .collect::<HashSet<&ByteArrayAllele>>();

        return result;
    }

    fn alleles_are_consistent<A: Allele>(
        alt_and_ref_1: &(&A, &A),
        alt_and_ref_2: &(&A, &A),
    ) -> bool {
        let alt1 = &alt_and_ref_1.0;
        let alt2 = &alt_and_ref_2.0;

        if alt1.is_symbolic() || alt2.is_symbolic() {
            return false;
        } else {
            // let size_diff_1 = alt_1.length() - alt_and_ref_1.1.length();
            // let size_diff_2 = alt_2.length() - alt_and_ref_2.1.length();

            alt_and_ref_1 == alt_and_ref_2
        }
    }

    pub fn split_reads_by_sample(
        reads: Vec<BirdToolRead>,
        n_samples: usize,
    ) -> HashMap<usize, Vec<BirdToolRead>> {
        let mut return_map = HashMap::new();
        (0..n_samples).into_iter().for_each(|sample_index| {
            return_map.entry(sample_index).or_insert(Vec::new());
        });

        for read in reads.into_iter() {
            let read_vec = return_map.entry(read.sample_index).or_insert(Vec::new());
            read_vec.push(read);
        }

        return return_map;
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    pub fn create_likelihood_calculation_engine(
        args: &clap::ArgMatches,
        handle_soft_clips: bool,
    ) -> PairHMMLikelihoodCalculationEngine {
        //AlleleLikelihoods::normalizeLikelihoods uses std::f64::NEGATIVE_INFINITY as a flag to disable capping
        let log10_global_read_mismapping_rate = if args
            .value_of("phred-scaled-global-read-mismapping-rate")
            .unwrap()
            .parse::<f64>()
            .unwrap()
            < 0.0
        {
            std::f64::NEG_INFINITY
        } else {
            QualityUtils::qual_to_error_prob_log10(
                args.value_of("phred-scaled-global-read-mismapping-rate")
                    .unwrap()
                    .parse::<u8>()
                    .unwrap(),
            )
        };

        PairHMMLikelihoodCalculationEngine::new(
            args.value_of("pair-hmm-gap-continuation-penalty")
                .unwrap()
                .parse()
                .unwrap(),
            log10_global_read_mismapping_rate,
            PCRErrorModel::new(args),
            args.value_of("base-quality-score-threshold")
                .unwrap()
                .parse()
                .unwrap(),
            args.is_present("enable-dynamic-read-disqualification-for-genotyping"),
            args.value_of("dynamic-read-disqualification-threshold")
                .unwrap()
                .parse()
                .unwrap(),
            args.value_of("expected-mismatch-rate-for-read-disqualification")
                .unwrap()
                .parse()
                .unwrap(),
            !args.is_present("disable-symmetric-hmm-normalizing"),
            args.is_present("disable-cap-base-qualities-to-map-quality"),
            handle_soft_clips,
        )
    }

    /**
     * Tries to phase the individual alleles based on pairwise comparisons to the other alleles based on all called haplotypes
     *
     * @param calls             the list of called alleles
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return a non-null list which represents the possibly phased version of the calls
     */
    pub fn phase_calls(
        mut calls: Vec<VariantContext>,
        called_haplotypes: &HashSet<Haplotype<SimpleInterval>>,
    ) -> Vec<VariantContext> {
        // construct a mapping from alternate allele to the set of haplotypes that contain that allele
        let haplotype_map = Self::construct_haplotype_mapping(&calls, called_haplotypes);

        // construct a mapping from call to phase set ID
        let phase_set_mapping = Self::construct_phase_set_mapping(&calls, haplotype_map);
        let unique_counter_end_value = phase_set_mapping.values().map(|p| p.0).unique().count();

        // we want to establish (potential) *groups* of phased variants, so we need to be smart when looking at pairwise phasing partners
        Self::construct_phase_groups(&mut calls, phase_set_mapping, unique_counter_end_value);

        return calls;
    }

    /**
     * Assemble the phase groups together and update the original calls accordingly
     *
     * @param originalCalls    the original unphased calls
     * @param phaseSetMapping  mapping from call (variant context) to phase group ID
     * @param indexTo          last index (exclusive) of phase group IDs
     * @return a non-null list which represents the possibly phased version of the calls
     */
    fn construct_phase_groups(
        original_calls: &mut Vec<VariantContext>,
        phase_set_mapping: HashMap<usize, (usize, PhaseGroup)>,
        index_to: usize,
    ) {
        // let mut phased_calls = vec![VariantContext::empty(0, 0, 0); original_calls.len()];

        // if we managed to find any phased groups, update the VariantContexts
        for count in 0..index_to {
            // get all of the (indexes of the) calls that belong in this group (keeping them in the original order)
            let mut indexes = Vec::new();
            for index in 0..original_calls.len() {
                let call = &original_calls[index];
                if phase_set_mapping.contains_key(&index) {
                    if phase_set_mapping.get(&index).unwrap().0 == count {
                        indexes.push(index)
                    }
                }
            }

            if indexes.len() < 2 {
                panic!("Somehow we have a group of phased variants that has fewer than 2 members")
            }

            // create a unique ID based on the leftmost one
            let unique_id = Self::create_unique_id(&original_calls[indexes[0]]);

            // create the phase set identifier, which is the position of the first variant in the set
            let phase_set_id = original_calls[indexes[0]].loc.start;

            // udpate the vcs
            for index in indexes {
                Self::phase_vc(
                    &mut original_calls[index],
                    unique_id.clone(),
                    &phase_set_mapping.get(&index).unwrap().1,
                    phase_set_id,
                )
            }
        }
    }

    /**
     * Add physical phase information to the provided variant context
     *
     * @param vc   the variant context
     * @param ID   the ID to use
     * @param phaseGT the phase GT string to use
     * @return phased non-null variant context
     */
    fn phase_vc<'b>(
        vc: &'b mut VariantContext,
        id: String,
        phase_gt: &'b PhaseGroup,
        phase_set_id: usize,
    ) {
        // let mut phased_genotypes = Vec::new();

        for g in vc.genotypes.genotypes_mut() {
            // let alleles = &mut g.alleles;
            // let mut new_alleles = alleles.clone();
            let phased_alt_allele_index = phase_gt.alt_allele_index;

            if g.is_het() && !Self::is_site_specific_alt_allele(&g.alleles[phased_alt_allele_index])
            {
                g.alleles.reverse();
            };

            g.is_phased = true;
            g.attribute(
                HAPLOTYPE_CALLER_PHASING_ID_KEY.to_string(),
                AttributeObject::String(id.clone()),
            );
            g.attribute(
                HAPLOTYPE_CALLER_PHASING_GT_KEY.to_string(),
                AttributeObject::String(phase_gt.get_description().clone()),
            );
            g.attribute(
                PHASE_SET_KEY.to_string(),
                AttributeObject::UnsizedInteger(phase_set_id),
            );
        }
    }

    /**
     * Create a unique identifier given the variant context
     *
     * @param vc   the variant context
     * @return non-null String
     */
    fn create_unique_id(vc: &VariantContext) -> String {
        return format!(
            "{}_{}_{}",
            vc.loc.start,
            std::str::from_utf8(vc.get_reference().bases.as_slice()).unwrap(),
            std::str::from_utf8(vc.get_alternate_alleles()[0].bases.as_slice()).unwrap()
        );
    }

    /**
     * Construct the mapping from call (variant context) to phase set ID
     *
     * @param totalAvailableHaplotypes the total number haplotypes that support alternate alleles of called variants
     * @param originalCalls    the original unphased calls
     * @param haplotypeMap     mapping from alternate allele to the set of haplotypes that contain that allele
     * @return a map from each variant context to a pair with the phase set ID and phase group of the alt allele
     *  note this may be empty in impossible-to-phase situations
     */
    fn construct_phase_set_mapping(
        original_calls: &Vec<VariantContext>,
        haplotype_map: HashMap<usize, HashSet<&Haplotype<SimpleInterval>>>,
    ) -> HashMap<usize, (usize, PhaseGroup)> {
        // count the total number of alternate haplotypes
        let mut haplotypes_with_called_variants: HashSet<&Haplotype<SimpleInterval>> =
            HashSet::new();
        haplotype_map.values().for_each(|h| {
            haplotypes_with_called_variants.extend(h);
        });

        let total_available_haplotypes = haplotypes_with_called_variants.len();

        let mut phase_set_mapping = HashMap::new();

        let num_calls = original_calls.len();
        let mut unique_counter = 0;

        // use the haplotype mapping to connect variants that are always/never present on the same haplotypes
        for i in 0..(num_calls).checked_sub(1).unwrap_or(0) {
            let call = &original_calls[i];
            let haplotypes_with_call = haplotype_map.get(&i);
            match haplotypes_with_call {
                None => continue,
                Some(haplotypes_with_call) => {
                    if haplotypes_with_call.is_empty() {
                        continue;
                    } else {
                        // NB: callIsOnAllAltHaps does not necessarily mean a homozygous genotype, since this method does
                        // not consider the reference haplotype
                        let call_is_on_all_alt_haps =
                            haplotypes_with_call.len() == total_available_haplotypes;

                        // if the call is on all haplotypes but we only use some of them to phase it with another variant
                        // we need to keep track of which ones are still active for downstream phasing.
                        // ie if the call is on all alt haps but we phase it with a site that is only on one of the alt haps,
                        // we remove the haplotypes with ref at the comp site for the purposes of phasing with additional variants
                        // downstream. This set keeps track of what call haplotypes are available for phasing downstream for "callIsOnAllAltHaps" variants.
                        let mut call_haplotypes_available_for_phasing = HashSet::new();

                        for j in (i + 1)..num_calls {
                            let comp = &original_calls[j];
                            let haplotypes_with_comp = haplotype_map.get(&j);
                            match haplotypes_with_comp {
                                None => continue,
                                Some(haplotypes_with_comp) => {
                                    if haplotypes_with_comp.is_empty() {
                                        continue;
                                    } else {
                                        // if the variants are together on all alt haplotypes, record that fact. NB that this does not mean that the
                                        // genotype for the variant is homozygous since this method does not consider the ref haplotype
                                        let comp_is_on_all_alt_haps = haplotypes_with_comp.len()
                                            == total_available_haplotypes;

                                        if ((haplotypes_with_call.len()
                                            == haplotypes_with_comp.len())
                                            && haplotypes_with_comp
                                                .iter()
                                                .all(|h| haplotypes_with_call.contains(h)))
                                            || (call_is_on_all_alt_haps
                                                && haplotypes_with_comp.iter().all(|h| {
                                                    call_haplotypes_available_for_phasing
                                                        .contains(*h)
                                                }))
                                            || comp_is_on_all_alt_haps
                                        {
                                            // create a new group if these are the first entries
                                            if !phase_set_mapping.contains_key(&i) {
                                                // note that if the comp is already in the map then that is very bad because it means that there is
                                                // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                                                // situation, we should abort if we encounter it.
                                                if phase_set_mapping.contains_key(&j) {
                                                    phase_set_mapping.clear();
                                                    return phase_set_mapping;
                                                };

                                                // An important note: even for homozygous variants we are setting the phase as "0|1" here. We
                                                // don't know at this point whether the genotype is homozygous vs the variant being on all alt
                                                // haplotypes, for one thing. Also, we cannot possibly know for sure at this time that the genotype for this
                                                // sample will actually be homozygous downstream: there are steps in the pipeline that are liable
                                                // to change the genotypes.  Because we can't make those assumptions here, we have decided to output
                                                // the phase as if the call is heterozygous and then "fix" it downstream as needed.
                                                phase_set_mapping
                                                    .insert(i, (unique_counter, PHASE_01.clone()));
                                                phase_set_mapping
                                                    .insert(j, (unique_counter, PHASE_01.clone()));

                                                // if the call was on all alternate haps but the comp isn't, we need to narrow down the set of
                                                // alt haps we'll consider for further phasing other variants with the call
                                                call_haplotypes_available_for_phasing
                                                    .retain(|h| haplotypes_with_comp.contains(&h));

                                                unique_counter += 1;
                                            } else if !phase_set_mapping.contains_key(&j) {
                                                // otherwise it's part of an existing group so use that group's unique ID
                                                let call_phase =
                                                    phase_set_mapping.get(&i).unwrap().clone();
                                                phase_set_mapping.insert(j, call_phase);
                                            }
                                        } else if (haplotypes_with_call.len()
                                            + haplotypes_with_comp.len())
                                            == total_available_haplotypes
                                        {
                                            // if the variants are apart on *all* alternate haplotypes, record that fact

                                            let intersection = haplotypes_with_call
                                                .intersection(haplotypes_with_comp)
                                                .collect::<HashSet<&&Haplotype<SimpleInterval>>>();
                                            if intersection.is_empty() {
                                                // create a new group if these are the first entries
                                                if !phase_set_mapping.contains_key(&i) {
                                                    // note that if the comp is already in the map then that is very bad because it means that there is
                                                    // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                                                    // situation, we should abort if we encounter it.
                                                    if phase_set_mapping.contains_key(&j) {
                                                        phase_set_mapping.clear();
                                                        return phase_set_mapping;
                                                    };

                                                    phase_set_mapping.insert(
                                                        i,
                                                        (unique_counter, PHASE_01.clone()),
                                                    );
                                                    phase_set_mapping.insert(
                                                        j,
                                                        (unique_counter, PHASE_10.clone()),
                                                    );
                                                    unique_counter += 1;
                                                } else if !phase_set_mapping.contains_key(&j) {
                                                    // otherwise it's part of an existing group so use that group's unique ID
                                                    let call_phase =
                                                        phase_set_mapping.get(&i).unwrap().clone();
                                                    phase_set_mapping.insert(
                                                        j,
                                                        (
                                                            call_phase.0,
                                                            if call_phase.1 == *PHASE_01 {
                                                                PHASE_10.clone()
                                                            } else {
                                                                PHASE_01.clone()
                                                            },
                                                        ),
                                                    );
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return phase_set_mapping;
    }

    /**
     * Construct the mapping from alternate allele to the set of haplotypes that contain that allele
     *
     * @param originalCalls    the original unphased calls
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return non-null Map
     */
    fn construct_haplotype_mapping<'b>(
        original_calls: &Vec<VariantContext>,
        called_haplotypes: &'b HashSet<Haplotype<SimpleInterval>>,
    ) -> HashMap<usize, HashSet<&'b Haplotype<SimpleInterval>>> {
        let mut haplotype_map = HashMap::new();

        for (index, call) in original_calls.iter().enumerate() {
            // don't try to phase if there is not exactly 1 alternate allele
            if !Self::is_biallelic_with_one_site_specific_alternate_allele(call) {
                haplotype_map.insert(index, HashSet::new());
            } else {
                // keep track of the haplotypes that contain this particular alternate allele
                let alt = Self::get_site_specific_alt_allele(call);
                match alt {
                    None => {
                        haplotype_map.insert(index, HashSet::new());
                    }
                    Some(alt) => {
                        let haps_with_allele = called_haplotypes
                            .iter()
                            .filter(|h| match &h.event_map {
                                None => false,
                                Some(event_map) => event_map.map.values().any(|vc| {
                                    vc.loc.get_start() == call.loc.get_start()
                                        && vc.get_alternate_alleles().contains(&alt)
                                }),
                            })
                            .collect::<HashSet<&Haplotype<SimpleInterval>>>();

                        haplotype_map.insert(index, haps_with_allele);
                    }
                }
            }
        }

        return haplotype_map;
    }

    /**
     * If at least one exists, returns a concrete (not NONREF) site-specific (starting at the current POS) alternate allele
     * from within the current variant context.
     */
    fn get_site_specific_alt_allele<'b>(call: &'b VariantContext) -> Option<&'b ByteArrayAllele> {
        let allele = call
            .get_alternate_alleles()
            .into_iter()
            .find(|a| Self::is_site_specific_alt_allele(a));
        return allele;
    }

    /**
     * Is this variant bi-allelic?  This implementation is very much specific to this class so shouldn't be pulled out into a generalized place.
     *
     * @param vc the variant context
     * @return true if this variant context is bi-allelic, ignoring the NON-REF symbolic allele and '*' symbolic allele, false otherwise
     */
    fn is_biallelic_with_one_site_specific_alternate_allele(vc: &VariantContext) -> bool {
        return vc
            .get_alternate_alleles()
            .iter()
            .filter(|a| Self::is_site_specific_alt_allele(a))
            .count()
            == 1;
    }

    /**
     * A site-specific alternate allele is one that represents concrete (i.e. not NONREF) variation that begins at the
     * site (i.e. not '*', which represents a concrete alternate allele that begins upstream of the current site).
     */
    fn is_site_specific_alt_allele(a: &ByteArrayAllele) -> bool {
        return !(a.is_ref
            || a.bases.as_slice() == "<NON_REF>".as_bytes()
            || a.bases.as_slice() == "<*>".as_bytes()
            || a.bases.as_slice() == "*".as_bytes());
    }
}
