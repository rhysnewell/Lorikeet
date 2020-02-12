use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use rayon::prelude::*;
use ordered_float::NotNan;

use ndarray::{Array, Array2, Axis, Zip};
use crate::factorization::{seeding::Seed, nmf_std};
use factorization::seeding::SeedFunctions;
use std::process;

//use crate::matrix_handling;
#[derive(Debug, Clone, Copy)]
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

#[derive(Debug, Clone, Copy)]
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

    fn is_satisfied(p_obj: &f32,
                    c_obj: &f32,
                    run: &usize,
                    min_residuals: &f32,
                    max_iter: &usize) -> bool;

    fn update_wh(v: &Array2<f32>,
                 w: Option<Array2<f32>>,
                 h: Option<Array2<f32>>,
                 update: &Update) -> (Option<Array2<f32>>, Option<Array2<f32>>);

    fn objective_update(v: &Array2<f32>,
                        w: &Option<Array2<f32>>,
                        h: &Option<Array2<f32>>,
                        cons: &Array2<f32>,
                        old_cons: &Array2<f32>,
                        objective: &Objective) -> (f32, Option<Array2<f32>>);

    fn basis(&self) -> &Array2<f32>;

    fn target(&self) -> &Array2<f32>;

    fn coef(&self) -> &Array2<f32>;

    fn fitted(&self) -> Array2<f32>;

    fn distance(&self, metric: Objective) -> f32;

    fn residuals(&self) -> Array2<f32>;

    fn predict(&self, what: &str) -> Array2<NotNan<f32>>;
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
                v,
                v1,
                h1,
                w1,
                ref mut h,
                ref mut w,
                seed,
                final_obj,
                rank,
                n_run,
                update,
                objective,
                conn_change,
                max_iter,
                min_residuals,
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
                let mut w_ret = None;
                let mut h_ret = None;
                let mut cons_ret = Array::zeros(
                    (v.shape()[0], v.shape()[1]));
                let mut old_cons_ret = Array::zeros(
                    (v.shape()[0], v.shape()[1]));
                let mut consecutive_conn = 0.;

                for run in (0..*n_run).into_iter() {
                    let (wsvd, hsvd) = seed.initialize(&v);
                    w_ret = Some(wsvd);
                    h_ret = Some(hsvd);

                    p_obj = std::f32::MAX;
                    c_obj = std::f32::MAX;
                    if run == 0 {
                        let mut best_obj = best_obj.lock().unwrap();
                        *best_obj = c_obj;
                    }
                    while Factorization::is_satisfied(&p_obj,
                                            &c_obj,
                                            &run,
                                            &min_residuals,
                                            &max_iter) {

                        // to satisfy borrow checker, connectivity has to be checked here
                        match objective {
                            Objective::Conn => {
                                if c_obj >= 1. {
                                    consecutive_conn *= 0.;
                                } else {
                                    consecutive_conn += 1.;
                                }
                                if consecutive_conn >= *conn_change as f32 {
                                    break
                                }
                            },
                            _ => {}
                        };

                        p_obj = c_obj;
                        let (mut w_update, h_update)
                            = Factorization::update_wh(&v, w_ret, h_ret, &update);
                        w_ret = w_update;
                        h_ret = h_update;
                        let (c_obj, new_cons) =
                            Factorization::objective_update(&v,
                                                          &w_ret,
                                                          &h_ret,
                                                          &cons_ret,
                                                          &old_cons_ret,
                                                            &objective);
                        match new_cons {
                            Some(array) => {
                                old_cons_ret = cons_ret.clone();
                                cons_ret = array;
                            },
                            None => {},
                        }
                        // nimfa adjusts small values here tow avoid underflows, but not sure
                        // if necessary
                    }

                };
                *w = w_ret;
                *h = h_ret;
                *cons = cons_ret;
                *old_cons = old_cons_ret;
            }
        }
    }

    fn is_satisfied(p_obj: &f32,
                    c_obj: &f32,
                    run: &usize,
                    min_residuals: &f32,
                    max_iter: &usize) -> bool {
            if *max_iter <= *run {
                false
            } else if run > &0 && (p_obj - c_obj) < *min_residuals {
                false
            } else if run > &0 && c_obj > p_obj {
                false
            } else {
                true
            }
        }

    fn update_wh(v: &Array2<f32>, w: Option<Array2<f32>>, h: Option<Array2<f32>>, update: &Update) -> (Option<Array2<f32>>, Option<Array2<f32>>){
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

                let w_return = Some(w_unwrap * &(v.dot(&h_unwrap.t()) / inner_dot));
                let h_return = Some(h_unwrap * &(upper_dot / lower_dot));

                return (w_return, h_return)
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

                let h_return = Some(h_unwrap * &(h_inner / h1));
                let w_return = Some(w_unwrap * &(w_inner / w1));
                return (w_return, h_return)
            },
            _ => {
                process::exit(1)
            },
        }
    }

    fn objective_update(v: &Array2<f32>,
                        w: &Option<Array2<f32>>,
                        h: &Option<Array2<f32>>,
                        cons: &Array2<f32>,
                        old_cons: &Array2<f32>,
                        objective: &Objective) -> (f32, Option<Array2<f32>>) {
        match objective {
            Objective::Fro => {
                // Compute squared Frobenius norm of a target matrix and its NMF estimate.
                let w_unwrap = w.as_ref().unwrap();
                let h_unwrap = h.as_ref().unwrap();

                let r = v.clone() - w_unwrap.dot(h_unwrap);
                return ((r.clone() * r).sum(), None)
            },
            Objective::Div => {
                // Compute divergence of target matrix from its NMF estimate.
                let w_unwrap = w.as_ref().unwrap();
                let h_unwrap = h.as_ref().unwrap();

                let va = w_unwrap.dot(h_unwrap);
                let inner_elop = (v.clone() / va.clone()).mapv(|x|{x.ln()});

                return (((v.clone() * inner_elop) - v.clone() + va).sum(), None)
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
                        let notnan_row: Vec<NotNan<f32>> = row
                            .iter()
                            .cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();

                        let max = notnan_row.iter().max().unwrap();
                        let argmax = notnan_row.iter().position(|element| element == max).unwrap();
                        idx[col_idx] = argmax;
                    });
                let mat1 = Array::from_elem(
                    (v.shape()[1], 1), idx.clone());
                let mat2 = Array::from_elem(
                    (1, v.shape()[1]), idx.t());

                let mut new_cons: Array2<f32> = Array::zeros(
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


                let mut connectivity_change = 0.;
                Zip::from(&mut new_cons)
                    .and(cons)
                    .apply(|a, b| {
                        if a != b {
                            connectivity_change += 1.;
                        }
                    });

                return (connectivity_change, Some(new_cons))
            },
            Objective::None => {
                process::exit(1);
            },
        }
    }

    fn basis(&self) -> &Array2<f32> {
        match self {
            Factorization::NMF {
                w,
                ..
            } => {
                return w.as_ref().unwrap()
            }
        }
    }

    fn target(&self) -> &Array2<f32> {
        match self {
            Factorization::NMF {
                v,
                ..
            } => {
                return v
            }
        }
    }

    fn coef(&self) -> &Array2<f32> {
        match self {
            Factorization::NMF {
                h,
                ..
            } => {
                return h.as_ref().unwrap()
            }
        }
    }

    fn fitted(&self) -> Array2<f32> {
        match self {
            Factorization::NMF {
                w,
                h,
                ..
            } => {
                let w_unwrap = w.as_ref().unwrap();
                let h_unwrap = h.as_ref().unwrap();

                return w_unwrap.dot(h_unwrap)
            }
        }
    }

    fn distance(&self, metric: Objective) -> f32 {
        match self {
            Factorization::NMF {
                v,
                w,
                h,
                ..
            } => {
                match metric {
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
                        let inner_elop = (v.clone() / va.clone()).mapv(|x| { x.ln() });

                        return ((v.clone() * inner_elop) - v.clone() + va).sum()
                    },
                    _ => {
                        process::exit(1)
                    }
                }
            }
        }
    }

    fn residuals(&self) -> Array2<f32> {
        match self {
            Factorization::NMF {
                v,
                w,
                h,
                ..
            } => {
                let w_unwrap = w.as_ref().unwrap();
                let h_unwrap = h.as_ref().unwrap();

                return v.clone() - w_unwrap.dot(h_unwrap)
            }
        }
    }

    fn predict(&self, what: &str) -> Array2<NotNan<f32>> {
        // Compute the dominant basis components. The dominant basis component is
        // computed as the row index for which the entry is the maximum within the column.
        match self {
            Factorization::NMF {
                v,
                w,
                h,
                ..
            } => {
                let mut x = Array::zeros((2, 2));
                match what {
                    "samples" => {
                        x = h.as_ref().unwrap().clone();
                    },
                    "features" => {
                        x = w.as_ref().unwrap().clone();
                    },
                    _ => {
                        error!("Invalid option for prediction parameter {}", what);
                        process::exit(1)
                    }
                };


                let mut idx = Array::zeros(
                    (x.shape()[1], 2));

                x.outer_iter().enumerate()
                    .for_each(|(col_idx, row)|{
                        let notnan_row: Vec<NotNan<f32>> = row
                            .iter()
                            .cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();

                        let max = notnan_row.iter().max().unwrap();
                        let argmax = notnan_row.iter().position(|element| element == max).unwrap();
                        idx[[col_idx, 0]] = *max;
                        idx[[col_idx, 1]] = NotNan::from(argmax as f32);
                    });

                let sums = x.sum_axis(Axis(0));

                let mut prob = Array::zeros((x.shape()[1]));

                Zip::from(&mut prob)
                    .and(&sums)
                    .and(&idx.slice(s![.., 0]))
                    .apply(|prob, sum, e|{
                        *prob = *e / NotNan::from(sum + 1f32.powf(-5.))
                    });

                stack![Axis(1), idx, prob.insert_axis(Axis(1))]

            }
        }
    }
}
