use model::variants::{Allele, NON_REF_ALLELE, Filter, Variant};
use genotype::genotype_builder::{GenotypesContext, GenotypeAssignmentMethod, Genotype};
use genotype::genotype_prior_calculator::GenotypePriorCalculator;
use utils::math_utils::MathUtils;
use rayon::prelude::*;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use genotype::genotype_likelihood_calculators::GenotypeLikelihoodCalculators;
use utils::vcf_constants::*;
use utils::simple_interval::{Locatable, SimpleInterval};
use std::collections::{HashMap, HashSet};
use ordered_float::OrderedFloat;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{IndexedReader, Record, Reader, Read};
use rust_htslib::bcf::header::HeaderView;
use reference::reference_reader::ReferenceReader;
use std::ops::Range;
use itertools::Itertools;

#[derive(Debug, Clone)]
pub struct VariantContext {
    pub loc: SimpleInterval,
    // variant alleles
    pub alleles: Vec<Allele>,
    // per sample likelihoods
    pub genotypes: GenotypesContext,
    pub source: String,
    pub log10_p_error: f64,
    pub filters: HashSet<Filter>,
    pub attributes: HashMap<String, Vec<f64>>
}

impl VariantContext {
    pub const MAX_ALTERNATE_ALLELES: usize = 180;
    pub const SUM_GL_THRESH_NOCALL: f64 = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    pub fn empty(tid: usize, start: usize, end: usize) -> VariantContext {
        VariantContext {
            loc: SimpleInterval::new(tid, start, end),
            alleles: Vec::new(),
            genotypes: GenotypesContext::empty(),
            source: "".to_string(),
            log10_p_error: 0.0,
            filters: HashSet::new(),
            attributes: HashMap::new(),
        }
    }

    pub fn build(tid: usize, start: usize, end: usize, alleles: Vec<Allele>) -> VariantContext {
        VariantContext {
            loc: SimpleInterval::new(tid, start, end),
            alleles,
            genotypes: GenotypesContext::empty(),
            source: "".to_string(),
            log10_p_error: 0.0,
            filters: HashSet::new(),
            attributes: HashMap::new(),
        }
    }

    pub fn build_from_vc(vc: &VariantContext) -> VariantContext {
        VariantContext {
            loc: vc.loc.clone(),
            alleles: Vec::new(),
            genotypes: GenotypesContext::empty(),
            source: "".to_string(),
            log10_p_error: 0.0,
            filters: HashSet::new(),
            attributes: HashMap::new(),
        }
    }

    pub fn attributes(&mut self, attributes: HashMap<String, Vec<f64>>) {
        self.attributes.par_extend(attributes);
    }

    pub fn set_attribute<S: Into<String>>(&mut self, tag: S, value: Vec<f64>) {
        self.attributes.insert(tag.into(), value);
    }

    pub fn filter(&mut self, filter: Filter) {
        self.filters.insert(filter);
    }

    pub fn get_log10_p_error(&self) -> f64 {
        self.log10_p_error
    }

    pub fn get_phred_scaled_qual(&self) -> f64 {
        (self.log10_p_error * (-10.0)) + 0.0
    }

    pub fn log10_p_error(&mut self, error: f64) {
        self.log10_p_error = error
    }

    pub fn add_alleles(&mut self, alleles: Vec<Allele>) {
        self.alleles = alleles
    }

    pub fn add_source(&mut self, source: String) {
        self.source = source
    }

    pub fn get_n_samples(&self) -> usize {
        self.genotypes.len()
    }

    pub fn get_n_alleles(&self) -> usize {
        self.alleles.len()
    }

    pub fn add_genotypes(
        &mut self,
        genotypes: Vec<Genotype>
    ) {
        for genotype in genotypes.iter() {
            if genotype.pl.len() != self.alleles.len() {
                panic!(
                    "Number of likelihoods does not match number of alleles at position {} on tid {}",
                    self.loc.get_start(), self.loc.get_contig()
                )
            }
        }
        self.genotypes = GenotypesContext::new(genotypes);
    }

    pub fn has_non_ref_allele(&self) -> bool {
        self.alleles.iter().any(|a| !a.is_reference())
    }

    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get into trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    pub fn has_too_many_alternative_alleles(&self) -> bool {
        if self.get_n_alleles() <= GenotypeLikelihoods::MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED {
            false
        } else {
            true
        }
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested {@link GenotypeAssignmentMethod}
     *
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param genotypeLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles with respect to which the likelihoods are defined
     */
    pub fn make_genotype_call(
        ploidy: usize,
        gb: &mut Genotype,
        assignment_method: &GenotypeAssignmentMethod,
        genotype_likelihoods: Option<Vec<f64>>,
        alleles_to_use: &Vec<Allele>,
        original_gt: &Vec<Allele>,
        gpc: &GenotypePriorCalculator
    ) {
        match genotype_likelihoods {
            Some(genotype_likelihoods) => {
                match assignment_method {
                    &GenotypeAssignmentMethod::SetToNoCall => {
                        gb.no_call_alleles(ploidy);
                        gb.no_qg();
                    },
                    &GenotypeAssignmentMethod::UsePLsToAssign => {
                        if !VariantContext::is_informative(&genotype_likelihoods) {
                            gb.no_call_alleles(ploidy);
                            gb.no_qg();
                        } else {
                            let max_likelihood_index = MathUtils::max_element_index(&genotype_likelihoods, 0, genotype_likelihoods.len());
                            let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles_to_use.len());
                            let allele_counts = gl_calc.genotype_allele_counts_at(max_likelihood_index);
                            let final_alleles = allele_counts.as_allele_list(alleles_to_use);
                            if final_alleles.contains(&*NON_REF_ALLELE) {
                                gb.no_call_alleles(ploidy);
                                gb.pl = GenotypeLikelihoods::from_log10_likelihoods(vec![0.0; genotype_likelihoods.len()]);
                            } else {
                                gb.alleles = final_alleles;
                            }
                            let num_alt_alleles = alleles_to_use.len() - 1;
                            if num_alt_alleles > 0 {
                                gb.log10_p_error(
                                    GenotypeLikelihoods::get_gq_log10_from_likelihoods_on_the_fly(max_likelihood_index, &genotype_likelihoods)
                                );
                            }
                        }
                    },
                    &GenotypeAssignmentMethod::SetToNoCallNoAnnotations => {
                        gb.no_call_alleles(ploidy);
                        gb.no_annotations();
                    },
                    &GenotypeAssignmentMethod::BestMatchToOriginal => {
                        let mut best = Vec::new();
                        let refr = &alleles_to_use[0];
                        for original_allele in original_gt.iter() {
                            if alleles_to_use.contains(original_allele) || original_allele.is_no_call() {
                                best.push(original_allele.clone())
                            } else {
                                best.push(refr.clone())
                            }
                        }
                        gb.alleles = best;
                    },
                    &GenotypeAssignmentMethod::UsePosteriorProbabilities => {
                        // Calculate posteriors.
                        let mut gl_calc = GenotypeLikelihoodCalculators::get_instance(ploidy, alleles_to_use.len());
                        let log10_priors = gpc.get_log10_priors(&mut gl_calc, alleles_to_use);
                        let log10_posteriors = MathUtils::ebe_add(&log10_priors, &genotype_likelihoods);
                        let normalized_log10_posteriors = MathUtils::scale_log_space_array_for_numeric_stability(
                            &log10_posteriors
                        );
                        // Update GP and PG annotations:
                        gb.attribute(
                            (&(*GENOTYPE_POSTERIORS_KEY)).to_string(),
                            normalized_log10_posteriors.par_iter().map(|v| {
                                if *v == 0.0 { // the reason for the == 0.0 is to avoid a signed 0 output "-0.0"
                                    0.0
                                } else {
                                    v * -10.
                                }
                            }).collect::<Vec<f64>>()
                        );

                        gb.attribute(
                            (&(*GENOTYPE_PRIOR_KEY)).to_string(),
                            log10_priors.par_iter().map(|v| {
                                if *v == 0.0 {
                                    0.0
                                } else {
                                    v * -10.
                                }
                            }).collect::<Vec<f64>>()
                        );
                        // Set the GQ accordingly
                        let max_posterior_index = MathUtils::max_element_index(&log10_posteriors, 0, log10_posteriors.len());
                        if alleles_to_use.len() > 0 {
                            gb.log10_p_error(VariantContext::get_gq_log10_from_posteriors(
                                max_posterior_index, &normalized_log10_posteriors,
                            ));
                        }
                        // Finally we update the genotype alleles.
                        gb.alleles(gl_calc.genotype_allele_counts_at(max_posterior_index).as_allele_list(alleles_to_use));
                    },
                    _ => panic!("Unknown GenotypeAssignmentMethod")
                }
            },
            _ => {
                gb.no_call_alleles(ploidy);
                gb.no_qg();
            }
        }
    }

    fn get_gq_log10_from_posteriors(best_genotype_index: usize, log10_posteriors: &[f64]) -> f64 {
        match log10_posteriors.len() {
            0 | 1 => {
                1.0
            },
            2 => {
                if best_genotype_index == 0 {
                    log10_posteriors[1]
                } else {
                    log10_posteriors[0]
                }
            },
            3 => {
                std::cmp::min(OrderedFloat(0.0), OrderedFloat(MathUtils::log10_sum_log10_two_values(
                    log10_posteriors[
                            if best_genotype_index == 0 {
                                2
                            } else {
                                best_genotype_index - 1
                            }
                        ],
                    log10_posteriors[
                            if best_genotype_index == 2 {
                                0
                            } else {
                                best_genotype_index + 1
                            }
                        ]
                    )
                )).into_inner()
            },
            _ => {
                if best_genotype_index == 0 {
                    MathUtils::log10_sum_log10(log10_posteriors, 1, log10_posteriors.len())
                } else if best_genotype_index == (log10_posteriors.len() - 1) {
                    MathUtils::log10_sum_log10(log10_posteriors, 0, best_genotype_index)
                } else {
                    std::cmp::min(
                        OrderedFloat(0.0),
                        OrderedFloat(MathUtils::log10_sum_log10_two_values(
                            MathUtils::log10_sum_log10(
                                log10_posteriors, 0, best_genotype_index
                            ),
                            MathUtils::log10_sum_log10(
                                log10_posteriors, best_genotype_index + 1, log10_posteriors.len()
                            )
                        ))
                    ).into_inner()
                }
            }
        }
    }

    pub fn is_informative(gls: &[f64]) -> bool {
        gls.par_iter().sum::<f64>() < VariantContext::SUM_GL_THRESH_NOCALL
    }

    /**
     * Subset the samples in VC to reference only information with ref call alleles
     *
     * Preserves DP if present
     *
     * @param vc the variant context to subset down to
     * @param defaultPloidy defaultPloidy to use if a genotype doesn't have any alleles
     * @return a GenotypesContext
     */
    pub fn subset_to_ref_only(vc: &mut VariantContext, default_ploidy: usize) -> GenotypesContext {
        if default_ploidy < 1 {
            panic!("default_ploidy must be >= 1, got {}", default_ploidy)
        } else {
            let old_gts = vc.get_genotypes().clone();
            if old_gts.is_empty() { return old_gts }

            let mut new_gts = GenotypesContext::create(old_gts.size());

            let ref_allele = vc.get_reference();
            let diploid_ref_alleles = vec![ref_allele.clone(); 2];

            // create the new genotypes
            for g in old_gts.genotypes() {
                let g_ploidy = if g.get_ploidy() == 0 { default_ploidy } else { g.get_ploidy() };
                let ref_alleles = if g_ploidy == 2 { diploid_ref_alleles.clone() } else { vec![ref_allele.clone(); g_ploidy] };

                let genotype = Genotype::build_from_alleles(ref_alleles);
                new_gts.add(genotype);
            }

            return new_gts
        }
    }

    pub fn get_genotypes(&self) -> &GenotypesContext {
        &self.genotypes
    }

    pub fn get_genotypes_mut(&mut self) -> &mut GenotypesContext {
        &mut self.genotypes
    }

    pub fn get_alleles(&self) -> &Vec<Allele> {
        &self.alleles
    }

    pub fn get_reference(&self) -> &Allele {
        &self.alleles[0]
    }

    pub fn get_dp(&self) -> i64 {
        self.genotypes.get_dp()
    }

    pub fn get_alternate_alleles(&self) -> Vec<Allele> {
        self.alleles[1..].to_vec()
    }

    pub fn process_vcf_in_region(
        indexed_vcf: &mut IndexedReader,
        tid: u32,
        start: u64,
        end: u64
    ) -> Vec<VariantContext> {
        indexed_vcf.fetch(tid, start, Some(end));

        let variant_contexts = indexed_vcf.records().into_iter().map(|record|{
            let mut vcf_record = record.unwrap();
            let vc = Self::from_vcf_record(
                &mut vcf_record,
            );
            vc.unwrap()
        }).collect::<Vec<VariantContext>>();

        return variant_contexts
    }

    pub fn process_vcf_from_path(vcf_path: &str) -> Vec<VariantContext> {

        let vcf_reader = Reader::from_path(vcf_path);
        match vcf_reader {
            Ok(ref mut reader) => {
                let variant_contexts = reader.records().into_iter().map(|record|{
                    let mut vcf_record = record.unwrap();
                    let vc = Self::from_vcf_record(
                        &mut vcf_record,
                    );
                    vc.unwrap()
                }).collect::<Vec<VariantContext>>();

                return variant_contexts
            }
            Err(_) => {
                debug!("No VCF records found for {}", vcf_path);
                return Vec::new()
            }
        }
    }

    pub fn from_vcf_record(
        record: &mut Record,
    ) -> Option<VariantContext> {
        let variants = Self::collect_variants(record, false, false, None);
        if variants.len() > 0 {

            // Get elements from record
            let mut vc = Self::build(
                record.rid().unwrap() as usize,
                record.pos() as usize,
                record.pos() as usize,
                variants,
            );
            vc.log10_p_error(record.qual() as f64);

            // Separate filters into hashset of filter struct
            let filters = record.filters();
            let header = record.header();
            for filter in filters {
                vc.filter(Filter::from_result(std::str::from_utf8(
                    &header.id_to_name(filter)[..],
                )));
            }

            // let mut depths = record.format(b"AD").float().unwrap_or(vec![0.0; variants.len()]);
            // let depths = depths.par_iter().map(|dp| dp[0] as f64).collect::<Vec<f64>>();
            //
            // vc.set_attribute("AD", depths);
            Some(vc)
        } else {
            None
        }
    }

    /// Collect variants from a given ´bcf::Record`.
    pub fn collect_variants(
        record: &mut Record,
        omit_snvs: bool,
        omit_indels: bool,
        indel_len_range: Option<Range<u32>>,
    ) -> Vec<Allele> {
        let pos = record.pos();
        let svlens = match record.info(b"SVLEN").integer() {
            // Gets value from SVLEN tag in VCF record
            Ok(Some(svlens)) => Some(
                svlens
                    .into_iter()
                    .map(|l| {
                        if !l.is_missing() {
                            Some(l.abs() as u32)
                        } else {
                            None
                        }
                    })
                    .collect_vec(),
            ),
            _ => None,
        };
        let end = match record.info(b"END").integer() {
            Ok(Some(end)) => {
                let end = end[0] as u32 - 1;
                Some(end)
            }
            _ => None,
        };

        // check if len is within the given range
        let is_valid_len = |svlen| {
            if let Some(ref len_range) = indel_len_range {
                // TODO replace with Range::contains once stabilized
                if svlen < len_range.start || svlen >= len_range.end {
                    return false;
                }
            }
            true
        };

        let is_valid_insertion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
            alt_allele == b"<INS>"
                || (ref_allele.len() < alt_allele.len()
                && ref_allele == &alt_allele[..ref_allele.len()])
        };

        let is_valid_deletion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
            alt_allele == b"<DEL>"
                || (ref_allele.len() > alt_allele.len()
                && &ref_allele[..alt_allele.len()] == alt_allele)
        };

        let is_valid_inversion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
            alt_allele == b"<INV>" || (ref_allele.len() == alt_allele.len())
        };

        let is_valid_mnv = |ref_allele: &[u8], alt_allele: &[u8]| {
            alt_allele == b"<MNV>" || (ref_allele.len() == alt_allele.len())
        };

        let variants = {
            let alleles = record.alleles();
            let ref_allele = alleles[0];
            let mut variant_vec = vec![];
            alleles
                .iter()
                .enumerate()
                .for_each(|(i, alt_allele)| {
                    let is_reference = i == 0;
                    if alt_allele == b"<*>" {
                        // dummy non-ref allele, signifying potential homozygous reference site
                        if omit_snvs {
                            variant_vec.push(Allele::new(Variant::None, is_reference))
                        } else {
                            variant_vec.push(Allele::new(Variant::None, is_reference))
                        }
                    } else if alt_allele == b"<DEL>" {
                        if let Some(ref svlens) = svlens {
                            if let Some(svlen) = svlens[i] {
                                variant_vec.push(Allele::new(Variant::Deletion(svlen), is_reference))
                            } else {
                                // TODO fail with an error in this case
                                variant_vec.push(Allele::new(Variant::None, is_reference))
                            }
                        } else {
                            // TODO fail with an error in this case
                            variant_vec.push(Allele::new(Variant::None, is_reference))
                        }
                    } else if alt_allele[0] == b'<' {
                        // TODO Catch <DUP> structural variants here
                        // skip any other special alleles
                        variant_vec.push(Allele::new(Variant::None, is_reference))
                    } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        // SNV
                        if omit_snvs {
                            variant_vec.push(Allele::new(Variant::None, is_reference))
                        } else if alt_allele == &ref_allele {
                            variant_vec.push(Allele::new(Variant::None, is_reference))
                        } else {
                            variant_vec.push(Allele::new(Variant::SNV(alt_allele[0]), is_reference))
                        }
                    } else {
                        let indel_len =
                            (alt_allele.len() as i32 - ref_allele.len() as i32).abs() as u32;
                        // TODO fix position if variant is like this: cttt -> ct

                        if (omit_indels || !is_valid_len(indel_len))
                            && is_valid_mnv(ref_allele, alt_allele)
                        {
                            // println!("MNV 1 {:?} {:?}", ref_allele, alt_allele);
                            variant_vec.push(Allele::new(Variant::MNV(alt_allele.to_vec()), is_reference))
                        } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                            variant_vec.push(Allele::new(Variant::Deletion(
                                (ref_allele.len() - alt_allele.len()) as u32,
                            ), is_reference))
                        } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                            variant_vec.push(Allele::new(Variant::Insertion(
                                alt_allele[ref_allele.len()..].to_owned(),
                            ), is_reference))
                        } else if is_valid_mnv(ref_allele, alt_allele) {
                            // println!("MNV 2 {:?} {:?}", ref_allele, alt_allele);
                            variant_vec.push(Allele::new(Variant::MNV(alt_allele.to_vec()), is_reference))
                        }
                    }
                });
            variant_vec
        };

        variants
    }

    pub fn get_contig_vcf_tid(vcf_header: &HeaderView, contig_name: &[u8]) -> Option<u32> {
        match vcf_header.name2rid(contig_name) {
            Ok(rid) => {
                return Some(rid)
            },
            Err(_) => {
                // Remove leading reference stem
                match vcf_header.name2rid(ReferenceReader::split_contig_name(contig_name, '~' as u8)) {
                    Ok(rid) => {
                        return Some(rid)
                    },
                    Err(_) => {
                        return None
                    }
                }
            }
        }
    }
}