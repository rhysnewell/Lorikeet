use assembly::assembly_region::AssemblyRegion;
use graphs::chain_pruner::ChainPruner;
use graphs::multi_sample_edge::MultiSampleEdge;
use haplotype::haplotype::Haplotype;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use haplotype::reference_confidence_model::ReferenceConfidenceModel;
use model::location_and_alleles::LocationAndAlleles;
use model::variant_context::VariantContext;
use model::variants::Allele;
use rayon::prelude::*;
use read_error_corrector::nearby_kmer_error_corrector::NearbyKmerErrorCorrector;
use read_threading::multi_debruijn_vertex::MultiDeBruijnVertex;
use read_threading::read_threading_assembler::ReadThreadingAssembler;
use reads::read_clipper::ReadClipper;
use reads::read_utils::ReadUtils;
use reference::reference_reader::ReferenceReader;
use rust_htslib::bam::ext::BamRecordExtensions;
use std::cmp::{max, min};
use std::collections::HashSet;
use utils::simple_interval::{Locatable, SimpleInterval};

lazy_static! {
    static ref PHASE_01: PhaseGroup = PhaseGroup::new("0|1".to_string(), 1);
    static ref PHASE_10: PhaseGroup = PhaseGroup::new("1|0".to_string(), 0);
}

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

pub struct AssemblyBasedCallerUtils {}

impl AssemblyBasedCallerUtils {
    const REFERENCE_PADDING_FOR_ASSEMBLY: usize = 500;

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

    /**
     *  Modify base qualities when paired reads overlap to account for the possibility of PCR error.
     *
     *  Overlapping mates provded independent evidence as far as sequencing error is concerned, but their PCR errors
     *  are correlated.  The base qualities are thus limited by the sequencing base quality as well as half of the PCR
     *  quality.  We use half of the PCR quality because downstream we treat read pairs as independent, and summing two halves
     *  effectively gives the PCR quality of the pairs when taken together.
     *
     * @param reads the list of reads to consider
     * @param samplesList   list of samples
     * @param readsHeader   bam header of reads' source
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR
     */
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
    pub fn assemble_reads(
        mut region: AssemblyRegion,
        given_alleles: &Vec<VariantContext>,
        args: &clap::ArgMatches,
        reference_reader: &mut ReferenceReader,
        assembly_engine: &mut ReadThreadingAssembler,
        correct_overlapping_base_qualities: bool,
        sample_names: &Vec<String>,
    ) {
        Self::finalize_regions(
            &mut region,
            args.is_present("error-correct-reads"),
            args.is_present("dont-use-soft-clipped-bases"),
            args.value_of("base-quality-threshold")
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
        let mut ref_haplotype =
            Self::create_reference_haplotype(&region, &padded_reference_loc, reference_reader);

        let pileup_error_correction_log_odds = args
            .value_of("pileup-correction-log-odds")
            .unwrap()
            .parse::<f64>()
            .unwrap();
        let mut read_error_corrector;
        if pileup_error_correction_log_odds == std::f64::NEG_INFINITY {
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
        } else {
            read_error_corrector = None
        };

        let assembly_result_set = assembly_engine.run_local_assembly(
            region,
            &mut ref_haplotype,
            full_reference_with_padding.as_slice(),
            padded_reference_loc,
            read_error_corrector,
            sample_names,
        );
    }

    pub fn get_variant_contexts_from_given_alleles(
        loc: usize,
        active_alleles_to_genotype: Vec<VariantContext>,
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

    pub fn get_padded_reference_loc(
        region: &AssemblyRegion,
        reference_padding: usize,
        reference_reader: &ReferenceReader,
    ) -> SimpleInterval {
        let pad_left = max(
            region
                .get_padded_span()
                .get_start()
                .checked_sub(reference_padding)
                .unwrap_or(0),
            1,
        );
        let pad_right = min(
            region.get_padded_span().get_end() + reference_padding,
            reference_reader.get_contig_length(region.get_contig()) as usize,
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

    pub fn get_alleles_consistent_with_given_alleles(
        given_alleles: Vec<VariantContext>,
        merged_vc: &VariantContext,
    ) -> HashSet<Allele> {
        if given_alleles.is_empty() {
            return HashSet::new();
        }

        let given_alt_and_ref_alleles_in_original_context =
            AssemblyBasedCallerUtils::get_variant_contexts_from_given_alleles(
                merged_vc.loc.get_start(),
                given_alleles,
                false,
            )
            .into_par_iter()
            .flat_map(|vc| {
                let refr = vc.get_reference().clone();
                let alt = vc.get_alternate_alleles();
                alt.into_par_iter()
                    .map(move |allele| (allele, refr.clone()))
            })
            .collect::<Vec<(Allele, Allele)>>();

        let result = merged_vc
            .get_alternate_alleles()
            .par_iter()
            .map(|allele| (allele.clone(), merged_vc.get_reference().clone()))
            .filter(|alt_and_ref| {
                given_alt_and_ref_alleles_in_original_context
                    .par_iter()
                    .any(|given_alt_and_ref| {
                        AssemblyBasedCallerUtils::alleles_are_consistent(
                            given_alt_and_ref,
                            alt_and_ref,
                        )
                    })
            })
            .map(|alt_and_ref_pair| alt_and_ref_pair.0)
            .collect::<HashSet<Allele>>();

        return result;
    }

    fn alleles_are_consistent(
        alt_and_ref_1: &(Allele, Allele),
        alt_and_ref_2: &(Allele, Allele),
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
}
