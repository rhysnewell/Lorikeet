use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use rayon::prelude::*;

use ndarray::{Array, Array2, Axis, Zip};
use crate::factorization::{seeding::Seed, nmf_std};
use factorization::seeding::SeedFunctions;
use std::process;

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
        consecutive_conn: f32,
        max_iter: usize,
        min_residuals: f32,
        cons: Array2<f32>,
        old_cons: Array2<f32>,
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

    fn is_satisfied(&mut self, p_obj: &f32,
                    c_obj: &f32,
                    run: &usize) -> bool;

    fn update_wh(&mut self);

    fn objective_update(&mut self) -> f32;
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
                ref mut consecutive_conn,
                ref mut max_iter,
                ref mut min_residuals,
                ref mut cons,
                ref mut old_cons,
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
                *consecutive_conn = 0.;
                *max_iter = miter;
                *min_residuals = minresiduals;
                *cons = Array::zeros((v.shape()[0], v.shape()[1]));
                *old_cons = Array::zeros((v.shape()[0], v.shape()[1]));
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
                ref mut consecutive_conn,
                ref mut max_iter,
                ref mut min_residuals,
                ref mut cons,
                ref mut old_cons,
            } => {
                //        Compute matrix factorization.
                //
                //        Return fitted factorization model.

                // Defined all variables first so they can be used inside closures
                let mut best_obj = Arc::new(Mutex::new(0.));
                let mut p_obj = std::f32::MAX;
                let mut c_obj = std::f32::MAX;

                for run in (0..*n_run).into_iter() {
                    seed.initialize(v);
                    let (wsvd, hsvd) = seed.get_wh();
                    *w = Some(wsvd.to_owned());
                    *h = Some(hsvd.to_owned());

                    p_obj = std::f32::MAX;
                    c_obj = std::f32::MAX;
                    if run == 0 {
                        let mut best_obj = best_obj.lock().unwrap();
                        *best_obj = c_obj;
                    }
                    while self.is_satisfied(&p_obj,
                                       &c_obj,
                                       &run) {
                        p_obj = c_obj;
                        self.update_wh();
                        c_obj = self.objective_update();
                        // nimfa adjusts small values here tow avoid underflows, but not sure
                        // if necessary
                    }

                };
            }
        }
    }

    fn is_satisfied(&mut self, p_obj: &f32,
                    c_obj: &f32,
                    run: &usize) -> bool {
        match self {
            Factorization::NMF {
                ref mut consecutive_conn,
                objective,
                min_residuals,
                conn_change,
                max_iter,
                ..
            } => {
                let conn = match objective {
                    Objective::Conn => true,
                    _ => false
                };
                if *max_iter <= *run {
                    false
                } else if conn {
                    if c_obj >= &1. {
                        *consecutive_conn *= 0.;
                    } else {
                        *consecutive_conn += 1.;
                    }
                    if *consecutive_conn >= *conn_change as f32 {
                        false
                    } else {
                        true
                    }
                } else if run > &0 && (p_obj - c_obj) < *min_residuals {
                    false
                } else if run > &0 && c_obj > p_obj {
                    false
                } else {
                    true
                }
            }
        }
    }

    fn update_wh(&mut self){
        match self {
            Factorization::NMF{
                update,
                ref mut h,
                ref mut w,
                ref mut v,
                ..
            } => {
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
            }
        }

    }

    fn objective_update(&mut self) -> f32 {
        match self {
            Factorization::NMF{
                objective,
                w,
                h,
                ref mut v,
                ref mut cons,
                ref mut old_cons,
                ..
            } => {
                match objective {
                    Objective::Fro => {
                        // Compute squared Frobenius norm of a target matrix and its NMF estimate.
                        let w_unwrap = w.as_ref().unwrap();
                        let h_unwrap = h.as_ref().unwrap();

                        let r = v.clone() - w_unwrap.dot(h_unwrap);
                        return (r.clone() * r).sum()
                    },
                    Objective::Div => {
                        // Compute divergence of target matrix from its NMF estimate.
                        let w_unwrap = w.as_ref().unwrap();
                        let h_unwrap = h.as_ref().unwrap();

                        let va = w_unwrap.dot(h_unwrap);
                        let inner_elop = (v.clone() / va.clone()).mapv(|x|{x.ln()});

                        return ((v.clone() * inner_elop) - v.clone() + va).sum()
                    },
                    Objective::Conn => {
                        // Compute connectivity matrix and compare it to connectivity matrix
                        // from previous iteration.
                        // Return logical value denoting whether connectivity matrix has changed
                        // from previous iteration.
                        let w_unwrap = w.as_ref().unwrap();
                        let h_unwrap = h.as_ref().unwrap();
                        let mut idx = Array::zeros(
                            (h_unwrap.shape()[1]));
                        h_unwrap.outer_iter().enumerate()
                            .for_each(|(col_idx, row)|{
                                let argmax = row.iter()
                                    .fold(None,|m,&x|
                                        m.map_or(Some(x), |mv|
                                            Some(if x > mv {x} else {mv})));
                                idx[col_idx] = argmax.unwrap();

                            });
                        let mat1 = Array::from_elem(
                            (v.shape()[1], 1), idx.clone());
                        let mat2 = Array::from_elem(
                            (1, v.shape()[1]), idx.t());

                        let mut new_cons = Array::zeros(
                            (v.shape()[0], v.shape()[1]));

                        Zip::from(&mut new_cons)
                            .and(&mat1)
                            .and(&mat2)
                            .apply(|a, b, c| {
                                if b == c {
                                    *a = 1.
                                } else {
                                    *a = 0.
                                }
                            });

                        *old_cons = cons.clone();
                        *cons = new_cons;

                        let mut connectivity_change = 0.;
                        Zip::from(cons)
                            .and(old_cons)
                            .apply(|a, b| {
                                if a != b {
                                    connectivity_change += 1.;
                                }
                            });

                        return connectivity_change
                    },
                    Objective::None => {
                        process::exit(1);
                    },
                }
            }
        }
    }
}
