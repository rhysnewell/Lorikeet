use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use rayon::prelude::*;

use ndarray::{Array, Array2, Axis};
use crate::factorization::{seeding::Seed, nmf_std};
use factorization::seeding::SeedFunctions;
//use crate::matrix_handling;

pub enum Update {
    Euclidean,
    Divergence,
    None
}

impl Update {
    pub fn from_str(to_match: &str) -> Self {
        match to_match {
            "euclidean" | "Euclidean" => {
                return Update::Euclidean
            },
            "divergence" | "Divergence" => {
                return Update::Divergence
            },
            _ => {
                return Update::None
            }
        };
    }
}

pub enum Objective {
    Fro,
    Div,
    Conn,
    None,
}

impl Objective {
    pub fn from_str(to_match: &str) -> Self {
        match to_match {
            "fro" | "frobenius" | "f" => {
                return Objective::Fro
            },
            "divergence" | "Divergence" | "div" | "d" => {
                return Objective::Div
            },
            "conn" | "connectivity" => {
                return Objective::Conn
            },
            _ => {
                return Objective::None
            }
        };
    }
}

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
        update: Update,
        objective: Objective,
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
                *seed = Seed::new_nndsvd(*rank, &input);
                *v = input;
                *v1 = None;
                *h1 = None;
                *w1 = None;
                *h = None;
                *w = None;
                *final_obj = 0.;
                *rank = r;
                *n_run = nrun;
                *update = Update::from_str(updatemethod);
                *objective = Objective::from_str(objectivemethod);
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
                    if *max_iter <= run {
                        false
                    } else if run > 0 && (p_obj - c_obj) < *min_residuals {
                        false
                    } else if run > 0 && c_obj > p_obj {
                        false
                    } else {
                        true
                    }
                };

                let update_wh = || {
                    match update {
                        Update::Euclidean => {
                            // Update basis and mixture matrix based on
                            // Euclidean distance multiplicative update rules.
                            let h_unwrap = h.as_ref().unwrap();
                            let w_unwrap = w.as_ref().unwrap();

                            let lower_dot: Array2<f32> = w_unwrap.t()
                                .dot(&w_unwrap
                                    .dot(h_unwrap));
                            let upper_dot: Array2<f32> = w_unwrap.t().dot(v);

                            let inner_dot: Array2<f32> = w_unwrap
                                .dot(&h_unwrap
                                    .dot(&h_unwrap.t()));

                            *w = Some(w_unwrap * &(v.dot(&h_unwrap.t()) / inner_dot));
                            *h = Some(h_unwrap * &(upper_dot / lower_dot));
                        },
                        Update::Divergence => {
                            // Update basis and mixture matrix based on
                            // Divergence distance multiplicative update rules.
                            let h_unwrap = h.as_ref().unwrap();
                            let w_unwrap = w.as_ref().unwrap();
                            let h1: Array2<f32> = Array::from_elem(
                                (1, v.shape()[1]), w_unwrap.sum_axis(Axis(0))[0]);

                            let inner_dot: Array2<f32> = w_unwrap.dot(h_unwrap);
                            let inner_elop: Array2<f32> = v.clone() / inner_dot;
                            let h_inner: Array2<f32> = w_unwrap.t()
                                .dot(&(inner_elop));

                            let w1: Array2<f32> = Array::from_elem(
                                (v.shape()[0], 1), h_unwrap.sum_axis(Axis(1))[0]);

                            let mut inner_dot: Array2<f32> = w_unwrap.dot(h_unwrap);
                            let inner_elop: Array2<f32> = v.clone() / inner_dot;
                            let w_inner: Array2<f32> = inner_elop.dot(&h_unwrap.t());

                            *h = Some(h_unwrap * &(h_inner / h1));
                            *w = Some(w_unwrap * &(w_inner / w1));

                        },
                        _ => {},
                    }
                };


                let mut best_obj = Arc::new(Mutex::new(0.));

                (0..*n_run).into_iter().for_each(|run|{
                    seed.initialize(v);
                    let (wsvd, hsvd) = seed.get_wh();
                    *w = Some(wsvd.to_owned());
                    *h = Some(hsvd.to_owned());

                    let mut p_obj = std::f32::MAX;
                    let mut c_obj = std::f32::MAX;
                    if run == 0 {
                        let mut best_obj = best_obj.lock().unwrap();
                        *best_obj = c_obj;
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
