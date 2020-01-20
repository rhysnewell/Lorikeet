use ndarray::{Array2, Array1, ArrayView};
use rayon::prelude::*;
use ndarray_npy::write_npy;
use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex, MutexGuard};

// Public enum to handle different Array dimensions as return of function
#[derive(Debug)]
pub enum VariantMatrix {
    Array1(Array1<f32>),
    Array2(Array2<f32>),
    Constraints(Array1<f32>),
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

    fn index(&mut self, row_index: usize, col_index: usize, n: usize, distance1: f32, distance2: Option<f32>) {
        match self {
            VariantMatrix::Array1(array1) | VariantMatrix::Constraints(array1) => {
                let condensed_index = get_condensed_index(row_index, col_index, n).unwrap();
                array1[[condensed_index]] = distance1;
            },
            VariantMatrix::Array2(array2) => {
                array2[[row_index, col_index]] = distance1;
                array2[[col_index, row_index]] = distance2.unwrap();
            },
        }
    }
}

pub fn get_condensed_distances(variant_info_all: &[(&i32, String, (f32, Vec<f32>), &i32)],
                               indels_map: &mut HashMap<i32, HashMap<i32, BTreeMap<String, BTreeSet<i64>>>>,
                               snps_map: &mut HashMap<i32, HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>,
                               geom_means: &[f32],
                               sample_count: i32,
                               dist_file: &str,
                               cons_file: &str) {

    let mut variant_distances = match sample_count {
        1 => {
            Arc::new(
                Mutex::new(
                    VariantMatrix::Array1(Array1::<f32>::zeros((variant_info_all.len().pow(2) - variant_info_all.len())/2))))
        },
        _ => {
            Arc::new(
                Mutex::new(
                    VariantMatrix::Array1(Array1::<f32>::zeros((variant_info_all.len().pow(2) - variant_info_all.len())/2))))
        }
    };

//    let mut constraints = Arc::new(Mutex::new(
//       VariantMatrix::Constraints(Array1::<f32>::zeros((variant_info_all.len().pow(2) - variant_info_all.len())/2)))
//    );

    debug!("Filling matrix of size {}",
           (variant_info_all.len().pow(2) - variant_info_all.len())/2);

    // Create variable to store mean of abundance if only one sample
    let mut vector_mean: f32 = 0.;

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
//                let mut variant_distances = variant_distances.lock().unwrap();
//                variant_distances[[row_index, col_index]] = 0.;
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
//                let mut constraint: f32 = 0.;
//                if row_info.0 == col_info.0 && row_info.3 == col_info.3 {
//                    constraint = 1.;
//                    constraints.lock().unwrap().index(row_index,
//                                                      col_index,
//                                                      n,
//                                                      constraint,
//                                                      None);
//                } else {
//                    let intersection_len = row_variant_set
//                        .intersection(&col_variant_set).collect::<HashSet<_>>().len() as f32;
//                    constraint = -((intersection_len + 1.) /
//                        ((row_info.2).0 + (col_info.2).0 - intersection_len + 1.));
//                    if constraint < 0. {
//                        constraints.lock().unwrap().index(row_index,
//                                                          col_index,
//                                                          n,
//                                                          constraint,
//                                                          None);
//                    }
//                }

                let mut distance: f32 = 0.;

                // If the variants share positions, then instantly they can't be in the same
                // gentoype so max distance
//                if row_start <= col_end && col_start <= row_end {
//                    distance = 1.;
                {

                    // Distance will be defined as the mean between jaccard dist
                    // and dist_f
                    let mut corr = 0.;
                    let mut w = 0;
                    if (row_info.2).1.len() > 1 {
                        // Calculate the log-ratio variance across compositions
                        // Essentially analogous to correlation
                        let mut log_vec = Arc::new(
                            Mutex::new(Vec::new()));
                        (row_info.2).1.par_iter()
                            .zip((col_info.2).1.par_iter()).for_each(|(r_freq, c_freq)|{
                            let mut log_vec = log_vec.lock().unwrap();
                            log_vec.push(((r_freq + 1.)/ (c_freq + 1.)).ln() as f32);
                        });
                        let log_vec = log_vec.lock().unwrap();
                        let sum = log_vec.iter().sum::<f32>();
                        let mean = sum / log_vec.len() as f32;
                        // calculate the variance of the log vector
                        let variance = log_vec.iter().map(|&value|{
                            let diff = mean - value;
                            diff * diff
                        }).sum::<f32>() / log_vec.len() as f32;


                        variant_distances.lock().unwrap().index(row_index,
                                                                col_index,
                                                                n,
                                                                distance,
                                                                None);

                    } else {
//                        if vector_mean == 0. {
//                            // loop through infos, should only happen once
//                            vector_mean = vector_info_all.par_iter().map(|info_tup|{
//                                (info_tup.2).1[0]
//                            }).sum::<f32>() / vector_info_all.len();
//                        }
                        let mut d_kl_a: f32 = 0.;
                        let mut d_kl_b: f32 = 0.;
                        let row_freq = (row_info.2).1[0];
                        let col_freq = (col_info.2).1[0];
                        if row_freq == 1. || col_freq == 1. {
                            // since the lim x->0 of xln(x) = 0
                            d_kl_a = row_freq * (row_freq / col_freq).ln();
                            d_kl_b = col_freq * (col_freq / row_freq).ln();
                        } else {
                            d_kl_a = row_freq * (row_freq / col_freq).ln()
                                + (1. - row_freq) * ((1. - row_freq) / (1. - col_freq)).ln();

                            d_kl_b = col_freq * (col_freq / row_freq).ln()
                                + (1. - col_freq) * ((1. - col_freq) / (1. - row_freq)).ln();
                        }
//                        let intersection_len = (row_variant_set
//                            .intersection(&col_variant_set).collect::<HashSet<_>>().len()) as f32;


                        if d_kl_a < 0. {
                            d_kl_a = 0.
                        }
                        if d_kl_b < 0. {
                            d_kl_b = 0.
                        }

                        variant_distances.lock().unwrap().index(row_index,
                                                                col_index,
                                                                n,
                                                                d_kl_a,
                                                                Some(d_kl_b));

                    }
                }
            }
        });
    });
//    return variant_distances
    variant_distances.lock().unwrap().write_npy(dist_file);
//    constraints.lock().unwrap().write_npy(cons_file);
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