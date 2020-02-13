use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use rayon::prelude::*;
use ordered_float::NotNan;

use ndarray::{Array, Array2, Axis, Zip};
use ndarray::prelude::*;
use crate::factorization::{seeding::Seed, nmf_std};
use factorization::seeding::SeedFunctions;
use std::process;
use std::f32;

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

impl Factorization {
    pub fn new_nmf(input: Array2<f32>,
               r: usize,
               nrun: usize,
               updatemethod: &str,
               objectivemethod: &str,
               connchange: usize,
               miter: usize,
               minresiduals: f32) -> Factorization {
        Factorization::NMF {
            seed: Seed::new_nndsvd(r, &input),
            v: input,
            v1: None,
            h1: None,
            w1: None,
            h: None,
            w: None,
            final_obj: 0.,
            rank: r,
            n_run: nrun,
            update: Update::from_str(updatemethod),
            objective: Objective::from_str(objectivemethod),
            conn_change: connchange,
            max_iter: miter,
            min_residuals: minresiduals,
            cons: Array::zeros((1, 1)),
            old_cons: Array::zeros((1, 1)),
        }
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
                    max_iter: &usize,
                    objective: &Objective) -> bool;

    fn update_wh(v: &Array2<f32>,
                 w: Option<Array2<f32>>,
                 h: Option<Array2<f32>>,
                 update: &Update) -> (Array2<f32>, Array2<f32>);

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

    fn rss (&self) -> f32;

    fn evar(&self) -> f32;
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
//                    debug!("Initialized H: {:?}", h_ret);

                    p_obj = std::f32::MAX;
                    c_obj = std::f32::MAX;

                    let mut iteration = 0;
                    while Factorization::is_satisfied(&p_obj,
                                            &c_obj,
                                            &iteration,
                                            &min_residuals,
                                            &max_iter,
                                            &objective) {

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
                        let (mut w_update, mut h_update)
                            = Factorization::update_wh(&v, w_ret, h_ret, &update);


                        // Adjust small values to prevent numerical underflow
                        w_update.par_mapv_inplace(|x| {
                            if x > f32::EPSILON {
                                x
                            } else {
                                f32::EPSILON
                            }
                        });

                        h_update.par_mapv_inplace(|x| {
                            if x > f32::EPSILON {
                                x
                            } else {
                                f32::EPSILON
                            }
                        });

                        w_ret = Some(w_update);
                        h_ret = Some(h_update);
                        let (new_c_obj, new_cons) =
                            Factorization::objective_update(&v,
                                                          &w_ret,
                                                          &h_ret,
                                                          &cons_ret,
                                                          &old_cons_ret,
                                                            &objective);
                        c_obj = new_c_obj;

                        match new_cons {
                            Some(array) => {
                                old_cons_ret = cons_ret.clone();
                                cons_ret = array;
                            },
                            None => {},
                        }
                        debug!("Consecutive Conn: {} c_obj {} p_obj {}", consecutive_conn,
                                 c_obj, p_obj);
                        iteration += 1;

                    }
                    if c_obj < *best_obj.lock().unwrap() || run == 0 {
                        let mut best_obj = best_obj.lock().unwrap();
                        *best_obj = c_obj;
                    }
                };
                debug!("H: {:?}", h_ret);
                debug!("W: {:?}", w_ret);
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
                    max_iter: &usize,
                    objective: &Objective) -> bool {
            if *max_iter <= *run {
                false
            } else if match objective {
                Objective::Conn => true,
                _ => false,
            }{
                // connectivity test outside this function
                true
            } else if run > &0 && (p_obj - c_obj) < *min_residuals {
                false
            } else if run > &0 && c_obj > p_obj {
                false
            } else {
                true
            }
        }

    fn update_wh(v: &Array2<f32>,
                 w: Option<Array2<f32>>,
                 h: Option<Array2<f32>>,
                 update: &Update) -> (Array2<f32>, Array2<f32>){
        match update {
            Update::Euclidean => {
                // Update basis and mixture matrix based on
                // Euclidean distance multiplicative update rules.
                let h_unwrap = h.as_ref().unwrap();
                let w_unwrap = w.as_ref().unwrap();

                // Build individual parts of functions
                // Function 1
                let w_dot_h = w_unwrap.dot(h_unwrap);
                let w_t = w_unwrap.t();

                let w_t_dot_w_dot_h: Array2<f32> = w_t.dot(&w_dot_h);
                let w_t_dot_v: Array2<f32> = w_t.dot(v);

                let upper_div_lower_a = w_t_dot_v / w_t_dot_w_dot_h ;
                //

                // Function 2
                let h_t = h_unwrap.t();
                let h_dot_h_t = h_unwrap.dot(&h_t);
                let w_dot_h_dot_h_t = w_unwrap.dot(&h_dot_h_t);
                let v_dot_h_t = v.dot(&h_t);

                let upper_div_lower_b = v_dot_h_t / w_dot_h_dot_h_t;

                let h_return = h_unwrap * &(upper_div_lower_a);
                let w_return = w_unwrap * &(upper_div_lower_b);

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

                let h_return = h_unwrap * &(h_inner / h1);
                let w_return = w_unwrap * &(w_inner / w1);
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
                let mut idx = Arc::new(Mutex::new(Array1::<usize>::zeros(
                    (h_unwrap.shape()[1]))));

                h_unwrap.axis_iter(Axis(1)).into_par_iter().enumerate()
                    .for_each(|(col_idx, row)|{
                        let notnan_row: Vec<NotNan<f32>> = row.into_par_iter().cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();

                        let max = notnan_row.par_iter().max().unwrap();
                        let argmax = notnan_row.par_iter().position(|element| element == max).unwrap();
                        let mut idx = idx.lock().unwrap();
                        idx[col_idx] = argmax;
                    });
                let mut idx = idx.lock().unwrap();
                let mat1 = idx.broadcast((v.shape()[1], v.shape()[1])).unwrap();
                let mat2 = mat1.t();

                let mut new_cons: Array2<f32> = Array::zeros(
                    (v.shape()[0], v.shape()[1]));
//                let mut mat1 = mat1.lock().unwrap();
//                let mut mat2 = mat2.lock().unwrap();

                Zip::from(&mut new_cons)
                    .and(mat1)
                    .and(mat2)
                    .par_apply(|a, b, c| {
                        if b == c {
                            *a = 1.
                        } else {
                            *a = 0.
                        }
                    });


                let connectivity_change = Arc::new(Mutex::new(0.));
                Zip::from(&mut new_cons)
                    .and(cons)
                    .par_apply(|a, b| {
                        if a != b {
                            let mut connectivity_change = connectivity_change.lock().unwrap();
                            *connectivity_change += 1.;
                        }
                    });
                let connectivity_change = connectivity_change.lock().unwrap().clone();
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

                x.axis_iter(Axis(1)).enumerate()
                    .for_each(|(col_idx, row)|{
                        let notnan_row: Vec<NotNan<f32>> = row
                            .into_par_iter()
                            .cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();

                        let max = notnan_row.par_iter().max().unwrap();
                        let argmax = notnan_row.par_iter().position(|element| element == max).unwrap();
                        idx[[col_idx, 0]] = *max;
                        idx[[col_idx, 1]] = NotNan::from(argmax as f32);
                    });

                let sums = x.sum_axis(Axis(0));

                let mut prob = Array::zeros((x.shape()[1]));

                Zip::from(&mut prob)
                    .and(&sums)
                    .and(&idx.slice(s![.., 0]))
                    .par_apply(|prob, sum, e|{
                        *prob = *e / NotNan::from(sum + 1e-5)
                    });

                stack![Axis(1), idx, prob.insert_axis(Axis(1))]

            }
        }
    }

    fn rss(&self) -> f32 {
        match self {
            Factorization::NMF {
                ..
            } => {
                let x = self.residuals();
                let rss = (&x * &x).sum();
                return rss
            }
        }
    }

    fn evar(&self) -> f32 {
        match self {
            Factorization::NMF {
                v,
                ..
            } => {
                let x = self.rss();
                let var = (v * v).sum();
                debug!("RSS {} VAR {}", x, var);
                return 1. - x / var
            }
        }
    }
}
