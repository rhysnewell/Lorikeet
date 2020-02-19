use ndarray::{Array2, Array1};
use rayon::prelude::*;
use ndarray_npy::write_npy;
use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex};
use std::process;

// Public enum to handle different Array dimensions as return of function
#[derive(Debug, Clone)]
pub enum VariantMatrix {
    Array1(Array1<f64>),
    Array2(Array2<f64>),
    Constraints(Array1<f64>),
}

impl VariantMatrix {

    pub fn write_npy(&self, tmp_path: &str) {
        match self {
            VariantMatrix::Array1(array1) | VariantMatrix::Constraints(array1) => {
                write_npy(&tmp_path, array1.to_owned()).expect("Unable to write to tempfile");
            },
            VariantMatrix::Array2(array2) => {
                write_npy(&tmp_path, array2.to_owned()).expect("Unable to write to tempfile");
            },
        }
    }

    fn index(&mut self, row_index: usize, col_index: usize, n: usize, distance1: f64, distance2: Option<f64>) {
        match self {
            VariantMatrix::Array1(array1) | VariantMatrix::Constraints(array1) => {
                let condensed_index = get_condensed_index(row_index, col_index, n).unwrap();
                array1[[condensed_index]] = distance1;
            },
            VariantMatrix::Array2(array2) => {
                array2[[row_index, col_index]] = distance1;
                array2[[col_index, row_index]] = distance1;
            },
        }
    }

    pub fn get_array2(&self) -> Array2<f64> {
        match self {
            VariantMatrix::Array2(array) => {
                return array.clone()
            },
            _ => {
                error!("Method cannot be used on Array1 VariantMatrix");
                process::exit(1)
            }
        }
    }
}

pub fn get_condensed_distances(variant_info_all: &[(&i32, String, (Vec<f64>, Vec<f64>), &i32)],
                               indels_map: &mut HashMap<i32, HashMap<i32, BTreeMap<String, BTreeSet<i64>>>>,
                               snps_map: &mut HashMap<i32, HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>,
                               geom_means_var: &[f64],
                               geom_means_dep: &[f64],
                               sample_count: i32) -> Arc<Mutex<VariantMatrix>> {

    let variant_distances = Arc::new(
                Mutex::new(
                    VariantMatrix::Array2(Array2::<f64>::zeros(
                        (variant_info_all.len(), variant_info_all.len())))));

//    let mut constraints = Arc::new(Mutex::new(
//       VariantMatrix::Constraints(Array1::<f64>::zeros((variant_info_all.len().pow(2) - variant_info_all.len())/2)))
//    );

    debug!("Filling matrix of size {}",
           (variant_info_all.len().pow(2) - variant_info_all.len())/2);

    // Create variable to store mean of abundance if only one sample
    let vector_mean: f64 = 0.;

    let n = variant_info_all.len();
    // produced condensed pairwise distances
    // described here: https://docs.rs/kodama/0.2.2/kodama/
    (0..variant_info_all.len()-1)
        .into_par_iter().for_each(|row_index| {
        let mut row_variant_set = &BTreeSet::new();
        let row_info = &variant_info_all[row_index];
        // lazily get the row variant read id set
        if indels_map[&row_info.3].contains_key(&row_info.0) {
            if indels_map[&row_info.3][&row_info.0].contains_key(&row_info.1) {
                row_variant_set = &indels_map[&row_info.3][&row_info.0][&row_info.1];
            } else if snps_map[&row_info.3].contains_key(&row_info.0) {
                let var_char = row_info.1.as_bytes()[0] as char;
                if snps_map[&row_info.3][&row_info.0].contains_key(&var_char) {
                    row_variant_set = &snps_map[&row_info.3][&row_info.0][&var_char];
                }
            }
        } else if snps_map[&row_info.3].contains_key(&row_info.0) {
            let var_char = row_info.1.as_bytes()[0] as char;
            if snps_map[&row_info.3][&row_info.0].contains_key(&var_char) {
                row_variant_set = &snps_map[&row_info.3][&row_info.0][&var_char];
            }
        }

        let row_start = *row_info.0 as usize;
        let row_end = row_start + row_info.1.len() - 1;

        (row_index+1..variant_info_all.len())
            .into_par_iter().for_each(|col_index| {
            if row_index == col_index {
                variant_distances.lock().unwrap().index(row_index,
                                                        col_index,
                                                        n,
                                                        0.,
                                                        None);
            } else {
                let mut col_variant_set = &BTreeSet::new();
                let col_info = &variant_info_all[col_index];
                if indels_map[&col_info.3].contains_key(&col_info.0) {
                    if indels_map[&col_info.3][&col_info.0].contains_key(&col_info.1) {
                        col_variant_set = &indels_map[&col_info.3][&col_info.0][&col_info.1];
                    } else if snps_map[&col_info.3].contains_key(&col_info.0) {
                        let var_char = col_info.1.as_bytes()[0] as char;
                        if snps_map[&col_info.3][&col_info.0].contains_key(&var_char) {
                            col_variant_set = &snps_map[&col_info.3][&col_info.0][&var_char];
                        }
                    }
                } else if snps_map[&col_info.3].contains_key(&col_info.0) {
                    let var_char = col_info.1.as_bytes()[0] as char;
                    if snps_map[&col_info.3][&col_info.0].contains_key(&var_char) {
                        col_variant_set = &snps_map[&col_info.3][&col_info.0][&var_char];
                    }
                }

//                let col_start = *col_info.0 as usize;
//                let col_end = col_start + col_info.1.len() - 1;

                // Calculate constraints value for Pmfcc
                // positive means they can't occur together i.e. two variants one position
                // negative means they must occur together i.e. share reads
                // Jaccard Similarity Modified for total read depth
                // Calculates the observed frequency of two variants together
                // |A (inter) B| / ((depth(A) + depth(B) - |A (inter) B|)
                let mut constraint: f64 = 1.;
                if row_info.0 == col_info.0 && row_info.3 == col_info.3 {
                    constraint = -1.;
//                    constraints.lock().unwrap().index(row_index,
//                                                      col_index,
//                                                      n,
//                                                      constraint,
//                                                      None);
                } else {
                    let intersection_len = row_variant_set
                        .intersection(&col_variant_set).collect::<HashSet<_>>().len() as f64;
                    constraint = 1. - (intersection_len + 1.) /
                        ((row_info.2).0.iter().sum::<f64>() +
                            (col_info.2).0.iter().sum::<f64>() - intersection_len + 1.);
//                    if constraint > 0. {
//                        constraints.lock().unwrap().index(row_index,
//                                                          col_index,
//                                                          n,
//                                                          constraint,
//                                                          None);
//                    }
                }

                let mut distance: f64 = 0.;

                // If the variants share positions, then instantly they can't be in the same
                // gentoype so max distance
//                if row_start <= col_end && col_start <= row_end {
//                    distance = 1.;
                    {

                        if (row_info.2).1.len() > 1 {

                            // Calculate the log-ratio variance across compositions
                            // Essentially analogous to correlation
                            let mut log_vec = Arc::new(
                                Mutex::new(Vec::new()));
                            (row_info.2).1.par_iter()
                                .zip((col_info.2).1.par_iter()).for_each(|(r_freq, c_freq)|{
                                let mut log_vec = log_vec.lock().unwrap();
                                log_vec.push(((r_freq)/ (c_freq)).ln() as f64);
                            });
                            let log_vec = log_vec.lock().unwrap();

                            let clr = |input: &Vec<f64>| -> Vec<f64> {
                                let output = input.iter().enumerate().map(|(i,v)| {
                                    (v / geom_means_var[i] as f64).ln()
                                }).collect();
                                return output
                            };

                            let get_mean = |input: &Vec<f64>| -> f64 {
                                let sum = input.iter().sum::<f64>();
                                sum / input.len() as f64
                            };

                            let row_vals: Vec<f64> = clr(&(row_info.2).1);

                            let col_vals: Vec<f64> = clr(&(col_info.2).1);

                            let mean_row = get_mean(&row_vals);

                            let mean_col = get_mean(&col_vals);

                            let mean = get_mean(&log_vec);


                            // calculate the variance of the log vector
                            let log_var = log_vec.iter().map(|&value|{
                                let diff = mean - value;
                                diff * diff
                            }).sum::<f64>() / log_vec.len() as f64;

                            // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4870310/ eq. 2
                            // p = 2*cov(Ai, Aj) / (Var(Ai) + Var(Aj))

                            // Distance as https://en.wikipedia.org/wiki/Cosine_similarity
                            let mut row_var = 0.;
                            let mut col_var = 0.;
                            let mut covar = 0.;

                            row_vals.iter()
                                .zip(col_vals.iter()).for_each(|(r_freq, c_freq)| {
                                row_var += (r_freq - mean_row).powf(2.);
                                col_var += (c_freq - mean_col).powf(2.);
                                covar += (r_freq - mean_row) * (c_freq - mean_col)
                            });

                            row_var = row_var / row_vals.len() as f64;
                            col_var = col_var / col_vals.len() as f64;
                            covar = covar / row_vals.len() as f64;

                            // technically correlation -1 to 1
                            distance = (2. * covar) / (row_var + col_var);

                            // 0 to 2
                            distance += 1.;
//                            distance = 1. - (-log_var.powf(1. / 2.)).exp();
                            if constraint < 0. {
                                distance = 0.
                            } else {
                                distance += (1. - constraint);
                            }
//                            if distance > 2. {
//                                distance = 2.;
//                            }
                            if distance < 0. {
                                distance = 0.;
                            } else {
                                variant_distances.lock().unwrap().index(row_index,
                                                                        col_index,
                                                                        n,
                                                                        distance,
                                                                        None);
                            }
                        } else {

                            let row_freq = (row_info.2).1[0];
                            let col_freq = (col_info.2).1[0];

                            let row_depth = (row_info.2).0[0];
                            let col_depth = (col_info.2).0[0];

                            distance = ((row_freq / geom_means_var[0] as f64 - col_freq / geom_means_var[0] as f64).powf(2.)
                                + (row_depth / geom_means_dep[0] as f64 - col_depth / geom_means_dep[0] as f64).powf(2.)).powf(1. / 2.);

                            if constraint < 0. {
                                distance = 1.
                            } else {
                                distance *= constraint
                            }

                            if distance < 0. {
                                error!("Negative value encountered {} invalid for NMF", distance);
                                process::exit(1);
                            } else {
                                variant_distances.lock().unwrap().index(row_index,
                                                                        col_index,
                                                                        n,
                                                                        distance,
                                                                        None);
                            }
                        }
                    }
                }
            });
        });
//    return variant_distances
    debug!("Distances {:?}", variant_distances);
//    variant_distances.lock().unwrap().write_npy(dist_file);
//    constraints.lock().unwrap().write_npy(cons_file);
    return variant_distances
}

// helper function to get the index of condensed matrix from it square form
fn get_condensed_index(i: usize, j: usize, n: usize) -> Option<usize>{
    if i == j {
        return None
    } else {
        return Some(n*i - i*(i+1)/2 + j - 1 - i)
//        return Some(n*(n-1)/2 - (n - row_i)*(n - row_i - 1)/2 + col_j - row_i - 1)
    }
}