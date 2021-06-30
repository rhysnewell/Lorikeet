use model::variant_context::VariantContext;
use model::location_and_alleles::LocationAndAlleles;
use model::variants::Allele;
use std::collections::HashSet;
use utils::simple_interval::Locatable;
use rayon::prelude::*;
use assembly::assembly_region::AssemblyRegion;

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
            alt_allele_index
        }
    }

    pub fn get_description(&self) -> String {
        self.description.clone()
    }

    pub fn get_alt_allele_index(&self) -> usize {
        self.alt_allele_index
    }
}

pub struct AssemblyBasedCallerUtils {

}

impl AssemblyBasedCallerUtils {

    /**
     * High-level function that runs the assembler on the given region's reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     */
    pub fn assemble_reads(
        region: AssemblyRegion,
        given_alleles: Vec<VariantContext>,
        args: &clap::ArgMatches,
        reference_reader: &
    )

    pub fn get_variant_contexts_from_given_alleles(
        loc: usize,
        active_alleles_to_genotype: Vec<VariantContext>,
        include_spanning_events: bool
    ) -> Vec<VariantContext> {
        let mut unique_locations_and_alleles = HashSet::new();
        let mut results = Vec::new();

        let mut given_allele_source_count = 0;
        for given_allele_vc in active_alleles_to_genotype.iter() {
            if given_allele_vc.loc.get_start() <= loc && given_allele_vc.loc.get_end() >= loc {
                if !(include_spanning_events || given_allele_vc.loc.get_start() == loc) {
                    continue
                }
                let mut allele_count = 0;
                for given_alt_allele in given_allele_vc.get_alternate_alleles() {
                    let allele_set = vec![given_allele_vc.get_reference().clone(), given_alt_allele.clone()];

                    let vc_source_name = format!("Comp{}Allele{}", given_allele_source_count, allele_count);
                    // check if this event is already in the list of events due to a repeat in the input alleles track
                    let mut candidate_event_to_add = VariantContext::build_from_vc(given_allele_vc);
                    candidate_event_to_add.add_source(vc_source_name);
                    candidate_event_to_add.add_alleles(allele_set);

                    let location_and_alleles = LocationAndAlleles::new(candidate_event_to_add.loc.get_start(), candidate_event_to_add.get_alleles().clone());

                    if !unique_locations_and_alleles.contains(&location_and_alleles) {
                        unique_locations_and_alleles.insert(location_and_alleles);
                        results.push(candidate_event_to_add);
                    }
                    allele_count += 1;
                }
            }
            given_allele_source_count += 1;
        }
        return results
    }

    pub fn get_alleles_consistent_with_given_alleles(given_alleles: Vec<VariantContext>, merged_vc: &VariantContext) -> HashSet<Allele> {
        if given_alleles.is_empty() {
            return HashSet::new();
        }

        let given_alt_and_ref_alleles_in_original_context =
            AssemblyBasedCallerUtils::get_variant_contexts_from_given_alleles(merged_vc.loc.get_start(), given_alleles, false)
                .into_par_iter()
                .flat_map(|vc| {
                    let refr = vc.get_reference().clone();
                    let alt = vc.get_alternate_alleles();
                    alt.into_par_iter().map(move |allele| {
                        (allele, refr.clone())
                    })
                }).collect::<Vec<(Allele, Allele)>>();

        let result = merged_vc.get_alternate_alleles().par_iter()
            .map(|allele| {
                (allele.clone(), merged_vc.get_reference().clone())
            })
            .filter(|alt_and_ref| {
                given_alt_and_ref_alleles_in_original_context
                    .par_iter().any(|given_alt_and_ref| {
                    AssemblyBasedCallerUtils::alleles_are_consistent(given_alt_and_ref, alt_and_ref)
                })
            })
            .map(|alt_and_ref_pair| alt_and_ref_pair.0)
            .collect::<HashSet<Allele>>();

        return result
    }

    fn alleles_are_consistent(alt_and_ref_1: &(Allele, Allele), alt_and_ref_2: &(Allele, Allele)) -> bool {
        let alt1 = &alt_and_ref_1.0;
        let alt2 = &alt_and_ref_2.0;

        if alt1.is_symbolic() || alt2.is_symbolic() {
            return false
        } else {
            // let size_diff_1 = alt_1.length() - alt_and_ref_1.1.length();
            // let size_diff_2 = alt_2.length() - alt_and_ref_2.1.length();

            alt_and_ref_1 == alt_and_ref_2
        }
    }
}