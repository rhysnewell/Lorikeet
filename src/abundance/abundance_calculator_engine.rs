use abundance::strain_abundances_calculator::StrainAbundanceCalculator;
use annotator::variant_annotation::VariantAnnotations;
use genotype::genotype_builder::AttributeObject;
use hashlink::LinkedHashMap;
use model::variant_context::VariantContext;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Calculates the per sample strain abundance for a list of variant contexts
/// that have been annotated with their potential strain assignments
///
/// @author Rhys Newell <rhys.newell@.hdr.qut.edu.au>
pub struct AbundanceCalculatorEngine<'a> {
    variant_contexts: Vec<VariantContext>,
    reference_name: &'a str,
    output_prefix: &'a str,
    sample_names: &'a Vec<&'a str>,
}

impl<'a> AbundanceCalculatorEngine<'a> {
    // Welcome. To the house... of Abundance

    pub fn new(
        variant_contexts: Vec<VariantContext>,
        reference_name: &'a str,
        output_prefix: &'a str,
        sample_names: &'a Vec<&'a str>,
    ) -> Self {
        Self {
            variant_contexts,
            reference_name,
            output_prefix,
            sample_names,
        }
    }

    pub fn run_abundance_calculator(
        mut self,
        mut n_strains: usize,
        n_samples: usize,
    ) -> Vec<VariantContext> {
        // The initialization vector for the EM algorithm
        let mut abundance_vectors = Vec::with_capacity(n_samples);
        let reference_present = self.reference_strain_potentially_present(n_samples);
        let mut per_sample_reference_presence = vec![reference_present; n_samples];

        n_strains += if reference_present { 1 } else { 0 };

        let mut strain_ids = (0..n_strains).into_iter().collect::<Vec<usize>>();
        let mut abundance_key = HashMap::with_capacity(n_strains);
        let mut strain_id_key = HashMap::with_capacity(n_strains);

        loop {
            for (abundance_index, strain_index) in strain_ids.iter().enumerate() {
                abundance_key.insert(*strain_index, abundance_index);
                strain_id_key.insert(abundance_index, *strain_index);
            }

            debug!(
                "Populating sample_vec with genotypes... {:?}",
                &abundance_key,
            );

            abundance_vectors = (0..n_samples)
                .into_iter()
                .map(|sample_vec_index| {
                    strain_ids
                        .iter()
                        .map(|strain_index| StrainAbundanceCalculator::new(*strain_index))
                        .collect::<Vec<StrainAbundanceCalculator>>()
                })
                .collect::<Vec<Vec<StrainAbundanceCalculator>>>();

            debug!(
                "Number of genotypes {} and genotype vectors {:?}",
                n_strains, abundance_vectors
            );

            for (sample_index, sample_vector) in abundance_vectors.iter_mut().enumerate() {
                for vc in self.variant_contexts.iter_mut() {
                    let strains = match vc
                        .attributes
                        .get_mut(VariantAnnotations::Strain.to_key())
                        .unwrap()
                    {
                        AttributeObject::VecUnsize(vec) => vec,
                        _ => panic!("Strain annotation was not the correct AttributeObject"),
                    };
                    // We divide the total depth of variant here by the total amount of strains that
                    // variant occurs in. E.g. if a variant had a depth of 6
                    // and occurred in 3 genotypes, then for each genotype its initialization value would be 2
                    let variant_depth = vc.genotypes.genotypes()[sample_index].ad[1] as f64;
                    let weight = variant_depth; // strains.len() as f64;

                    let mut strains_to_remove = HashSet::new();
                    for strain in strains.iter() {
                        let abundance_index = match abundance_key.get(strain) {
                            Some(idx) => *idx,
                            None => {
                                strains_to_remove.insert(*strain);
                                continue;
                            }
                        };

                        sample_vector[abundance_index].variant_weights.push(weight);

                        // Collect the other abundance indices that this variant is associated with
                        let other_abundance_indices = strains
                            .iter()
                            .filter_map(|idx| match abundance_key.get(idx) {
                                Some(a_idx) => Some(*a_idx),
                                None => None,
                            })
                            .collect();

                        sample_vector[abundance_index]
                            .variant_genotype_ids
                            .push(other_abundance_indices);
                    }

                    // Remove the strains that weren't present
                    for to_remove in strains_to_remove.into_iter() {
                        let strain_index = strains.iter().position(|i| *i == to_remove).unwrap();
                        strains.remove(strain_index);
                    }

                    let mut reference_present = &mut per_sample_reference_presence[sample_index];
                    if *reference_present {
                        // We divide the total depth of variant here by the total amount of strains that
                        // variant occurs in. E.g. if a variant had a depth of 6
                        // and occurred in 3 genotypes, then for each genotype its initialization value would be 2
                        let reference_depth = vc.genotypes.genotypes()[sample_index].ad[0] as f64;
                        let weight = reference_depth; // (n_strains - strains.len()) as f64;

                        let reference_strain_index = n_strains - 1;
                        let abundance_index = match abundance_key.get(&reference_strain_index) {
                            Some(idx) => *idx,
                            None => {
                                *reference_present = false;
                                continue;
                            }
                        };

                        sample_vector[abundance_index].variant_weights.push(weight);

                        // Collect the other abundance indices that this variant is associated with
                        // For the refernce variant this is the strain ids that aren't currently in
                        // strains vector
                        let other_abundance_indices = (0..n_strains)
                            .into_iter()
                            .filter_map(|idx| {
                                if !strains.contains(&idx) {
                                    match abundance_key.get(&idx) {
                                        Some(a_idx) => Some(*a_idx),
                                        None => None,
                                    }
                                } else {
                                    None
                                }
                            })
                            .collect();

                        sample_vector[abundance_index]
                            .variant_genotype_ids
                            .push(other_abundance_indices);
                    }
                }
            }

            debug!("Calculating abundances...");

            abundance_vectors
                .par_iter_mut()
                .enumerate()
                .for_each(|(idx, sample_calculators)| {
                    debug!("Genotype Vector before EM {} {:?}", idx, sample_calculators);
                    StrainAbundanceCalculator::calculate_abundances(sample_calculators);
                    debug!("Genotype Vector after EM {} {:?}", idx, sample_calculators);
                });

            // Vector of counters for each genotype
            // If the counter reaches the same value as the number of samples
            // Then that genotype had an abundance weighting of 0 in
            // every sample and therefore does not exist so remove it.
            let mut strains_to_remove = vec![0; n_strains];
            for sample_vector in abundance_vectors.iter() {
                for abundance_calculator in sample_vector.iter() {
                    if abundance_calculator.abundance_weight == 0.0 {
                        strains_to_remove
                            [*abundance_key.get(&abundance_calculator.index).unwrap()] += 1;
                    }
                }
            }

            debug!("strain removal counts {:?}", &strains_to_remove);
            let mut something_removed = false;
            for (strain_index, count) in strains_to_remove.into_iter().enumerate() {
                if count == n_samples {
                    let strain_id_to_remove = strain_id_key.get(&strain_index).unwrap();
                    let position_of_strain_id = strain_ids
                        .iter()
                        .position(|s| s == strain_id_to_remove)
                        .unwrap();
                    strain_ids.remove(position_of_strain_id);
                    something_removed = true
                }
            }

            if !something_removed {
                break;
            } else {
                debug!("Genotype removed, rerunning abundance calculations");
            }
        }

        self.print_strain_coverages(abundance_vectors);

        self.variant_contexts
    }

    fn print_strain_coverages(&self, abundance_vectors: Vec<Vec<StrainAbundanceCalculator>>) {
        debug!("Printing strain coverages {}", self.reference_name);
        let file_name = format!(
            "{}/{}_strain_coverages.tsv",
            self.output_prefix, self.reference_name,
        );

        let file_path = Path::new(&file_name);

        let mut file_open = match File::create(file_path) {
            Ok(coverage_file) => coverage_file,
            Err(e) => {
                println!("Cannot create file {:?}", e);
                std::process::exit(1)
            }
        };

        // rearrange the genotype vector for better printing
        // Just free genotype struct from memory but keep the abundance weight
        let mut printing_genotype: LinkedHashMap<usize, Vec<f64>> = LinkedHashMap::new();
        for (sample_idx, abundance_vector) in abundance_vectors.into_iter().enumerate() {
            for abundance_calculator in abundance_vector.into_iter() {
                let genotype_info = printing_genotype
                    .entry(abundance_calculator.index)
                    .or_insert(vec![0.; self.sample_names.len()]);
                genotype_info[sample_idx] = abundance_calculator.abundance_weight
            }
        }

        writeln!(
            file_open,
            "##source=lorikeet-v{}",
            env!("CARGO_PKG_VERSION")
        );
        for (sample_idx, sample_name) in self.sample_names.iter().enumerate() {
            // remove tmp file name from sample id
            writeln!(
                file_open,
                "##sample=<ID={}, name={}>",
                sample_idx + 1,
                sample_name
            );
        }

        // Print header line
        write!(file_open, "{: <10}", "strainID").unwrap();
        for sample in 0..self.sample_names.len() {
            write!(file_open, "\t{: <10}", sample + 1).unwrap();
        }
        write!(file_open, "\n").unwrap();

        for (strain_id, abundances) in printing_genotype.iter() {
            write!(file_open, "strain_{: <10}", strain_id,).unwrap();

            for coverage in abundances.iter() {
                write!(file_open, "\t{: <10}", coverage).unwrap();
            }
            write!(file_open, "\n").unwrap();
        }
    }

    fn reference_strain_potentially_present(&self, n_samples: usize) -> bool {
        let mut reference_presence_counters = vec![0; n_samples];
        for vc in self.variant_contexts.iter() {
            for (sample_index, genotype) in vc.get_genotypes().genotypes().into_iter().enumerate() {
                if genotype.ad[0] > 0 {
                    reference_presence_counters[sample_index] += 1;
                }
            }
        }

        reference_presence_counters
            .iter()
            .any(|count| *count == self.variant_contexts.len())
    }
}
