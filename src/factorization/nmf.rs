use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use rayon::prelude::*;

use ndarray::{Array, Array2, Axis};
use crate::factorization::{seeding::Seed, nmf_std};
use factorization::seeding::SeedFunctions;
//use crate::matrix_handling;

pub enum Factorization {
    NMF {
        v: Array2<f32>,
        v1: Option<Array2<f32>>,
        h1: Option<Array2<f32>>,
        w1: Option<Array2<f32>>,
        h: Option<Array2<f32>>,
        w: Option<Array2<f32>>,
        seed: Seed,
        final_obj: f32,
        rank: usize,
        n_run: usize,
        update: String,
        objective: String,
        conn_change: usize,
        max_iter: usize,
        min_residuals: f32,
    }
}

pub trait RunFactorization {
    fn initialize(&mut self,
                  input: Array2<f32>,
                  r: usize,
                  nrun: usize,
                  updatemethod: &str,
                  objectivemethod: &str,
                  connchange: usize,
                  miter: usize,
                  minresiduals: f32);

    fn factorize(&mut self);
}

impl RunFactorization for Factorization {
    fn initialize(&mut self,
                  input: Array2<f32>,
                  r: usize,
                  nrun: usize,
                  updatemethod: &str,
                  objectivemethod: &str,
                  connchange: usize,
                  miter: usize,
                  minresiduals: f32) {
        match self {
            Factorization::NMF {
                ref mut v,
                ref mut v1,
                ref mut h1,
                ref mut w1,
                ref mut h,
                ref mut w,
                ref mut seed,
                ref mut final_obj,
                ref mut rank,
                ref mut n_run,
                ref mut update,
                ref mut objective,
                ref mut conn_change,
                ref mut max_iter,
                ref mut min_residuals,
            } => {
                *v = input;
                *v1 = None;
                *h1 = None;
                *w1 = None;
                *h = None;
                *w = None;
                *seed = Seed::new_nnsvd(rank, &input);
                *final_obj = 0.;
                *rank = r;
                *n_run = nrun;
                *update = updatemethod.to_string();
                *objective = objectivemethod.to_string();
                *conn_change = connchange;
                *max_iter = miter;
                *min_residuals = minresiduals;
            }
        }
    }

    fn factorize(&mut self) {
        match self {
            Factorization::NMF {
                ref mut v,
                ref mut v1,
                ref mut h1,
                ref mut w1,
                ref mut h,
                ref mut w,
                ref mut seed,
                ref mut final_obj,
                ref mut rank,
                ref mut n_run,
                ref mut update,
                ref mut objective,
                ref mut conn_change,
                ref mut max_iter,
                ref mut min_residuals,
            } => {
                //        Compute matrix factorization.
                //
                //        Return fitted factorization model.
                let is_satisfied = |p_obj: f32, c_obj: f32, run: usize| -> bool {
                    if self.max_iter <= run {
                        false
                    } else if run > 0 && (p_obj - c_obj) < self.min_residuals {
                        false
                    } else if run > 0 && c_obj > p_obj {
                        false
                    } else {
                        true
                    }
                };

                let update_wh = || {
                    match update {
                        String::from("euclidean") => {
                            // Update basis and mixture matrix based on
                            // Euclidean distance multiplicative update rules.
                            let h_unwrap = h.unwrap();
                            let w_unwrap = w.unwrap();
                            *h = Some(h_unwrap * (w_unwrap.t().dot(v)
                                / w_unwrap.t().dot(w_unwrap.dot(h_unwrap))));
                            *w = Some(w_unwrap * (v.dot(h_unwrap.t()
                                / w_unwrap.dot(h_unwrap.dot(h_unwrap.t())))));
                        },
                        String::from("divergence") => {
                            // Update basis and mixture matrix based on
                            // Divergence distance multiplicative update rules.
                            let h_unwrap = h.unwrap();
                            let w_unwrap = w.unwrap();
                            let h1 = Array::from_elem(
                                (1, v.shape()[1]), w_unwrap.sum_axis(Axis(0)));
                            *h = Some(h_unwrap * (
                                w_unwrap.t().dot(
                                    v / (w_unwrap.dot(
                                        h_unwrap))))) / h1;

                            let w1 = Array::from_elem(
                                (v.shape()[0], 1), h_unwrap.sum_axis(Axis(1)));
                            *w = Some(w_unwrap * ((v / (w_unwrap.dot(h_unwrap))).dot(h_unwrap.t()) / w1));

                        },
                        _ => {},
                    }
                };


                let mut best_obj = Arc::new(Mutex::new(0.))

                (0..self.n_run).into_par_iter().for_each(|run|{
                    seed.initialize(self.v);
                    let (wsvd, hsvd) = seed.get_wh();
                    *w = Some(*wsvd);
                    *h = Some(*hsvd);

                    let mut p_obj = std::f32::MAX;
                    let mut c_obj = std::f32::MAX;
                    if run == 0 {
                        let mut best_obj = best_obj.lock().unwrap();
                        best_obj = c_obj;
                    }
                    while is_satisfied(p_obj, c_obj, run) {
                        p_obj = c_obj;
                        // nimfa adjusts small values here tow avoid underflows, but not sure
                        // if necessary

                    }
                });
            }
        }
    }
}
