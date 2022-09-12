use annotator::variant_annotation::{Annotation, VariantAnnotations};
use genotype::genotype_builder::AttributeObject;
use itertools::Itertools;
use mathru::algebra::abstr::Sign;
use model::variant_context::VariantContext;
use model::variant_context_utils::VariantContextUtils;
use ndarray::Array2;
use std::cmp::min;
use std::fs::File;
use std::io::Write;
use std::path::Path;
// use polars;

/// Holds the population and consensus ANI & Fst arrays
/// Compares ANI between samples (Non-diagonal cells) and the ANI of sample compared to reference (Diagonals)
///
/// Defined by InStrain: https://www.nature.com/articles/s41587-020-00797-0?proof=t%3B#Sec9
/// Consensus ANI:
///     1 - number_of_consensus_allele_diff / compared_position_count
/// Population ANI:
///     1 - population_alleles_diff / compared_position_count
///
/// Consensus Allele == Allele with highest depth at position
///     - If this is the same between two samples, then do not incremenet number_of_consensus_alleles_diff
///
/// Population ANI. If you have the reference and another allele
/// at a position then this is, for all intents and purposes, the same as the reference/other sample
/// for popANI.
///
/// Subpopulation ANI is an extension of PopulationANI, where a difference is measured if the two
/// populations being compared don't share all the same alleles
///
/// To summarize:
/// Consensus ANI - Share the same consensus Allele
/// Population ANI - Share at least one common allele
/// Subpopulation ANI - Share all alleles
///
/// Since lorikeet calls Indels, we compare the length of the allele rather than just the position
/// So the rather than alleles different, it is bases different.
pub struct ANICalculator {
    popANI: Array2<f64>,
    subpopANI: Array2<f64>,
    conANI: Array2<f64>,
    // fst: Array2<f64>
}

impl ANICalculator {
    pub fn new(n_samples: usize) -> Self {
        Self {
            popANI: Array2::default((n_samples, n_samples)),
            subpopANI: Array2::default((n_samples, n_samples)),
            conANI: Array2::default((n_samples, n_samples)),
            // fst: Array2::default((n_samples, n_samples)),
        }
    }

    pub fn run_calculator(
        &mut self,
        contexts: &mut [VariantContext],
        output_prefix: &str,
        sample_names: &[&str],
        reference_name: &str,
        genome_size: u64,
        passing_sites: Option<Vec<Vec<i32>>>,
        qual_by_depth_filter: f64,
        qual_threshold: f64,
        depth_per_sample_filter: i64,
    ) {
        let compared_bases =
            self.calculate_compared_bases(passing_sites, genome_size, sample_names.len());

        self.calculate_from_contexts(
            contexts,
            genome_size,
            qual_by_depth_filter,
            qual_threshold,
            depth_per_sample_filter,
            compared_bases,
        );

        Self::write_ani_tables(
            output_prefix,
            sample_names,
            reference_name,
            &self.conANI,
            "consensus_ani",
        );
        Self::write_ani_tables(
            output_prefix,
            sample_names,
            reference_name,
            &self.popANI,
            "population_ani",
        );
        Self::write_ani_tables(
            output_prefix,
            sample_names,
            reference_name,
            &self.subpopANI,
            "subpopulation_ani",
        );
    }

    fn calculate_compared_bases(
        &mut self,
        passing_sites: Option<Vec<Vec<i32>>>,
        genome_size: u64,
        n_samples: usize,
    ) -> Array2<f64> {
        let mut compared_bases = Array2::default((n_samples, n_samples));
        match passing_sites {
            Some(passing_sites) => {
                for s1_ind in 0..passing_sites.len() {
                    let s1 = &passing_sites[s1_ind];
                    for s2_ind in (s1_ind + 1)..passing_sites.len() {
                        let s2 = &passing_sites[s2_ind];

                        let mut i1 = 0; // index in each sample array
                        let mut i2 = 0;

                        let mut used1 = 0; // a buffer containing the number of used bases at the current index
                        let mut used2 = 0;
                        let mut differing_bases = 0; // number of bases found to be incomparable

                        while i1 < s1.len() && i2 < s2.len() {
                            // the current values
                            let val1 = s1[i1];
                            let val2 = s2[i2];

                            let abs1 = val1.abs() - used1;
                            let abs2 = val2.abs() - used2;

                            if val1.sign() != val2.sign() {
                                differing_bases += min(abs1, abs2);
                            }
                            used1 += min(abs1, abs2);
                            used2 += min(abs1, abs2);

                            if used1 >= val1.abs() && used2 >= val2.abs() {
                                i1 += 1;
                                i2 += 2;
                                used1 -= val1.abs();
                                used2 -= val2.abs();
                            } else if used1 >= val1.abs() {
                                i1 += 1;
                                used1 -= val1.abs();
                            } else {
                                i2 += 1;
                                used2 -= val2.abs();
                            }
                        }

                        let comparable_bases = (genome_size - differing_bases as u64) as f64;
                        compared_bases[[s1_ind, s2_ind]] = comparable_bases;
                        compared_bases[[s2_ind, s1_ind]] = comparable_bases;
                    }
                }
            }
            _ => {
                compared_bases.iter_mut().for_each(|val| {
                    *val = genome_size as f64;
                });
            }
        }
        compared_bases
    }

    /// Takes refernce to a vec of variant contexts and compares the consensus and population
    /// ANI between each sample. The input contexts need to be non split i.e. prior to being
    /// put through genotyping pipeline
    fn calculate_from_contexts(
        &mut self,
        contexts: &mut [VariantContext],
        genome_size: u64,
        qual_by_depth_filter: f64,
        qual_threshold: f64,
        depth_per_sample_filter: i64,
        compared_bases: Array2<f64>,
    ) {
        let n_samples = self.conANI.ncols();

        for context in contexts {
            let n_alleles = context.get_n_alleles();

            let mut consenus_allele_indices = Vec::with_capacity(n_samples);
            let mut present_alleles = Vec::with_capacity(n_samples);
            if !VariantContextUtils::passes_thresholds(
                context,
                qual_by_depth_filter,
                qual_threshold,
            ) {
                continue;
            }

            // println!("Context passes {} {}", context.log10_p_error, context.get_dp());
            // don't consider poor quality variant sites
            for sample_idx_1 in 0..n_samples {
                if consenus_allele_indices.len() == sample_idx_1 {
                    // get consensus of first sample
                    consenus_allele_indices.push(
                        context
                            .get_consensus_allele_index(sample_idx_1)
                            .unwrap_or_default(),
                    );
                    // which alleles are present in first sample
                    let mut which_are_present =
                        context.alleles_present_in_sample(sample_idx_1, depth_per_sample_filter);

                    present_alleles.push(which_are_present);
                }

                if !present_alleles[sample_idx_1].iter().any(|val| *val) {
                    continue; // nothing present here
                }

                for sample_idx_2 in 0..n_samples {
                    if consenus_allele_indices.len() == sample_idx_2 {
                        // get consensus of first sample, default to ref
                        consenus_allele_indices.push(
                            context
                                .get_consensus_allele_index(sample_idx_2)
                                .unwrap_or_default(),
                        );
                        // which alleles are present in first sample with at least two supporting reads
                        let mut which_are_present = context
                            .alleles_present_in_sample(sample_idx_2, depth_per_sample_filter);
                        present_alleles.push(which_are_present);
                    }

                    if !present_alleles[sample_idx_2].iter().any(|val| *val) {
                        continue; // nothing present here
                    }

                    if sample_idx_1 != sample_idx_2 {
                        let consensus_1 = &consenus_allele_indices[sample_idx_1];
                        let allele_present_1 = &present_alleles[sample_idx_1];

                        let consensus_2 = &consenus_allele_indices[sample_idx_2];
                        let allele_present_2 = &present_alleles[sample_idx_2];

                        if consensus_1 != consensus_2 {
                            if context.alleles[*consensus_1].len() > 1
                                || context.alleles[*consensus_2].len() > 1
                            {
                                let bases_different = (context.alleles[*consensus_1].len() as f64
                                    - context.alleles[*consensus_2].len() as f64)
                                    .abs();
                                self.conANI[[sample_idx_1, sample_idx_2]] += bases_different;
                            } else {
                                self.conANI[[sample_idx_1, sample_idx_2]] += 1.0;
                            }
                        }

                        let mut bases_different = 0.0;
                        let mut divisor = 0.0;
                        for (idx, (present_1, present_2)) in allele_present_1
                            .iter()
                            .zip(allele_present_2.iter())
                            .enumerate()
                        {
                            if present_1 != present_2 {
                                bases_different += context.alleles[idx].len() as f64;
                                divisor += 1.0;
                            }
                        }

                        bases_different =
                            bases_different / if divisor > 0.0 { divisor } else { 1.0 };

                        if !allele_present_1
                            .iter()
                            .zip(allele_present_2.iter())
                            .any(|(present_1, present_2)| *present_1 && *present_2)
                        {
                            // if they share ANY alleles, then popANI does not change

                            self.popANI[[sample_idx_1, sample_idx_2]] += bases_different;
                        }

                        if allele_present_1 != allele_present_2 {
                            self.subpopANI[[sample_idx_1, sample_idx_2]] += bases_different;
                        }
                    } else {
                        let consensus_1 = &consenus_allele_indices[sample_idx_1];
                        let allele_present_1 = &present_alleles[sample_idx_1];

                        if *consensus_1 != 0 {
                            if context.alleles[*consensus_1].len() > 1
                                || context.alleles[0].len() > 1
                            {
                                let bases_different = (context.alleles[*consensus_1].len() as f64
                                    - context.alleles[0].len() as f64)
                                    .abs();
                                self.conANI[[sample_idx_1, sample_idx_2]] += bases_different;
                            } else {
                                self.conANI[[sample_idx_1, sample_idx_2]] += 1.0;
                            }
                        }

                        if !allele_present_1[0] {
                            // reference not present
                            let mut bases_different = 0.0;
                            let mut divisor = 0.0;
                            for (idx, present) in allele_present_1.iter().enumerate() {
                                if *present {
                                    bases_different += context.alleles[idx].len() as f64;
                                    divisor += 1.0;
                                }
                            }
                            bases_different =
                                bases_different / if divisor > 0.0 { divisor } else { 1.0 };

                            self.popANI[[sample_idx_1, sample_idx_2]] += bases_different;
                            self.subpopANI[[sample_idx_1, sample_idx_2]] += bases_different;
                        }
                    }
                }
            }
        }

        // let length = genome_size as f64;
        self.popANI
            .iter_mut()
            .zip(compared_bases.iter())
            .for_each(|(val, length)| {
                *val = 1.0 - (*val / length);
            });

        self.subpopANI
            .iter_mut()
            .zip(compared_bases.iter())
            .for_each(|(val, length)| {
                *val = 1.0 - (*val / length);
            });

        self.conANI
            .iter_mut()
            .zip(compared_bases.iter())
            .for_each(|(val, length)| {
                *val = 1.0 - (*val / length);
            })
    }

    fn write_ani_tables(
        output_prefix: &str,
        sample_names: &[&str],
        reference_name: &str,
        table: &Array2<f64>,
        table_name: &str,
    ) {
        debug!("Printing ani calculations {}", reference_name);
        let file_name = format!("{}/{}_{}.tsv", output_prefix, reference_name, table_name);

        let file_path = Path::new(&file_name);

        let mut file_open = match File::create(file_path) {
            Ok(file) => file,
            Err(e) => {
                println!("Cannot create file {:?}", e);
                std::process::exit(1)
            }
        };

        writeln!(
            file_open,
            "##source=lorikeet-v{}",
            env!("CARGO_PKG_VERSION")
        );

        for (sample_idx, sample_name) in sample_names.iter().enumerate() {
            // remove tmp file name from sample id
            writeln!(
                file_open,
                "##sample=<ID={}, name={}>",
                sample_idx + 1,
                sample_name
            );
        }

        // Print header line
        write!(file_open, "{: <10}", "SampleID").unwrap();
        for sample in 0..sample_names.len() {
            write!(file_open, "\t{: <8}", sample + 1).unwrap();
        }
        writeln!(file_open).unwrap();

        let mut counter = 0;
        for ani_vals in table.rows() {
            write!(file_open, "{}", counter + 1).unwrap();
            counter += 1;
            for ani in ani_vals {
                write!(file_open, "\t{:.8}", ani).unwrap();
            }
            writeln!(file_open).unwrap();
        }
    }
}
