use std::error::Error;
use std::sync::{Arc, Mutex};
use std::prelude::*;
use ordered_float::NotNan;
use approx;

use ndarray::{Array, Array2, Axis, Zip, stack};
use ndarray::prelude::*;
use rayon::prelude::*;
use crate::factorization::{seeding::Seed};
use factorization::seeding::SeedFunctions;
use std::process;
use std::f64;


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
        v: Array2<f64>,
        v1: Option<Array2<f64>>,
        h1: Option<Array2<f64>>,
        w1: Option<Array2<f64>>,
        h: Array2<f64>,
        w: Array2<f64>,
        seed: Seed,
        final_obj: f64,
        rank: usize,
        n_run: usize,
        update: Update,
        objective: Objective,
        conn_change: usize,
        max_iter: usize,
        min_residuals: f64,
        cons: Array2<f64>,
        old_cons: Array2<f64>,
    }
}

impl Factorization {
    pub fn new_nmf(input: Array2<f64>,
               r: usize,
               nrun: usize,
               updatemethod: &str,
               objectivemethod: &str,
               seedmethod: Seed,
               connchange: usize,
               miter: usize,
               minresiduals: f64) -> Factorization {
        Factorization::NMF {
            seed: seedmethod,
            v: input,
            v1: None,
            h1: None,
            w1: None,
            h: Array2::zeros((1, 1)),
            w: Array2::zeros((1, 1)),
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
                  input: Array2<f64>,
                  r: usize,
                  nrun: usize,
                  updatemethod: &str,
                  objectivemethod: &str,
                  seedmethod: Seed,
                  connchange: usize,
                  miter: usize,
                  minresiduals: f64);

    fn estimate_rank(&mut self);

    fn factorize(v: &Array2<f64>,
                seed: Seed,
                final_obj: f64,
                rank: usize,
                update: Update,
                objective: Objective,
                conn_change: usize,
                max_iter: usize,
                min_residuals: f64) -> (Array2<f64>, Array2<f64>);

    fn is_satisfied(p_obj: &f64,
                    c_obj: &f64,
                    run: &usize,
                    min_residuals: &f64,
                    max_iter: &usize,
                    objective: &Objective) -> bool;

    fn update_wh(v: &Array2<f64>,
                 w: Array2<f64>,
                 h: Array2<f64>,
                 update: &Update) -> (Array2<f64>, Array2<f64>);

    fn objective_update(v: &Array2<f64>,
                        w: &Array2<f64>,
                        h: &Array2<f64>,
                        cons: &Array2<f64>,
                        old_cons: &Array2<f64>,
                        objective: &Objective) -> (f64, Option<Array2<f64>>);

    fn adjustment(input: &mut Array2<f64>) -> Array2<f64>;

    fn basis(&self) -> &Array2<f64>;

    fn target(&self) -> &Array2<f64>;

    fn coef(&self) -> &Array2<f64>;

    fn fitted(&self) -> Array2<f64>;

    fn distance(&self, metric: Objective) -> f64;

    fn residuals(&self) -> Array2<f64>;

    fn predict(&self, what: &str) -> Array2<NotNan<f64>>;

    fn probabilities(&self, what: &str) -> Array2<NotNan<f64>>;

    fn rss (&self) -> f64;

    fn evar(&self) -> f64;
}

impl RunFactorization for Factorization {
    fn initialize(&mut self,
                  input: Array2<f64>,
                  r: usize,
                  nrun: usize,
                  updatemethod: &str,
                  objectivemethod: &str,
                  seedmethod: Seed,
                  connchange: usize,
                  miter: usize,
                  minresiduals: f64) {
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
                *seed = seedmethod;
                *v = input;
                *v1 = None;
                *h1 = None;
                *w1 = None;
                *h = Array2::zeros((1, 1));
                *w = Array2::zeros((1, 1));
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

    fn estimate_rank(&mut self) {
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
                let best_obj = Arc::new(Mutex::new(vec![0.; *rank]));
                let mut best_rank = 0;


                (2..*rank+1).into_iter().for_each(|r| {
                    let (w_ret, h_ret) = Factorization::factorize(v, *seed, *final_obj,
                                                                 r, *update, Objective::Fro, *conn_change,
                                                                 50, *min_residuals);
                    let (c_obj, new_cons) =
                        Factorization::objective_update(&v,
                                                        &w_ret,
                                                        &h_ret,
                                                        &Array2::zeros((1, 1)),
                                                        &Array2::zeros((1, 1)),
                                                        &Objective::Fro);
                    let mut best_obj = best_obj.lock().unwrap();
                    best_obj[r-1] = c_obj;

                });
                let best_obj = best_obj.lock().unwrap();
                let mut prev = std::f64::MAX;

                for (r, curr) in best_obj.iter().enumerate() {
                    debug!("PREV {} CURR {}", prev, curr);
                    if (&prev - curr) > 1e-3 {
                        prev = *curr;
                        best_rank = r;
                    } else {
                        break
                    }
                }
                info!("Best NMF rank: {}", best_rank);

                let (w_ret, h_ret) = Factorization::factorize(v, *seed, *final_obj,
                                                              best_rank, *update, *objective, *conn_change,
                                                              *max_iter, *min_residuals);

                *w = w_ret.clone();
                *h = h_ret.clone();

            }
        }
    }

    fn factorize(v: &Array2<f64>,
                 seed: Seed,
                 final_obj: f64,
                 rank: usize,
                 update: Update,
                 objective: Objective,
                 conn_change: usize,
                 max_iter: usize,
                 min_residuals: f64) -> (Array2<f64>, Array2<f64>){
        let mut p_obj = std::f64::MAX;
        let mut c_obj = std::f64::MAX;
        let mut cons_ret = Array::zeros(
            (v.shape()[0], v.shape()[1]));
        let mut old_cons_ret = Array::zeros(
            (v.shape()[0], v.shape()[1]));
        let mut consecutive_conn = 0.;

        let (mut w_ret, mut h_ret) = seed.initialize(&v, &rank);

//                    debug!("Initialized H: {:?}", h_ret);

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
                    if consecutive_conn >= conn_change as f64 {
                        break
                    }
                },
                _ => {}
            };

            p_obj = c_obj;
            let (mut w_update, mut h_update)
                = Factorization::update_wh(&v,  w_ret, h_ret, &update);

            // This code does not work and it seems to function fine without it
//                        // Adjust small values to prevent numerical underflow
            w_ret = Factorization::adjustment(&mut w_update);
            h_ret = Factorization::adjustment(&mut h_update);

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
//            debug!("Consecutive Conn: {} c_obj {} p_obj {}", consecutive_conn,
//                   c_obj, p_obj);
            iteration += 1;
        }
        info!("NMF using rank {} objective function value {}", rank, c_obj);
        return (w_ret, h_ret)
    }

    fn adjustment(input: &mut Array2<f64>) -> Array2<f64> {
        input.par_mapv_inplace(|x| {
            if x > f64::EPSILON {
                x
            } else if x.is_nan() {
                f64::EPSILON
            } else {
                f64::EPSILON
            }
        });

        return input.to_owned()
    }

    fn is_satisfied(p_obj: &f64,
                    c_obj: &f64,
                    run: &usize,
                    min_residuals: &f64,
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

    fn update_wh(v: &Array2<f64>,
                 w: Array2<f64>,
                 h: Array2<f64>,
                 update: &Update) -> (Array2<f64>, Array2<f64>){
        match update {
            Update::Euclidean => {
                // Update basis and mixture matrix based on
                // Euclidean distance multiplicative update rules.
                // Build individual parts of functions
                // H is updated first, and then used to update W
                // Function 1
                let w_t = w.t();
                let mut lower_1 = w_t.dot(v);
                let lower_2 = w_t.dot(&w.dot(&h));
                let mut lower_1 = w_t.dot(v) / w_t.dot(&w.dot(&h));

                let h = h * &(lower_1);

                // Function 2
                let mut h_dot_h_t = h.dot(&h.t());
                h_dot_h_t = Factorization::adjustment(&mut h_dot_h_t);
                let mut w_dot_h_dot_h_t = w.dot(&h_dot_h_t);
                w_dot_h_dot_h_t = Factorization::adjustment(&mut w_dot_h_dot_h_t);
                let mut v_dot_h_t = v.dot(&h.t());
                v_dot_h_t = Factorization::adjustment(&mut v_dot_h_t);
                let upper_div_lower_b = v_dot_h_t / w_dot_h_dot_h_t;
                let w = w * upper_div_lower_b;

                return (w, h)
            },
            Update::Divergence => {
                // Update basis and mixture matrix based on
                // Divergence distance multiplicative update rules.
                let h1: Array2<f64> = Array::from_elem(
                    (1, v.shape()[1]), w.sum_axis(Axis(0))[0]);

                let inner_dot: Array2<f64> = w.dot(&h);
                let inner_elop: Array2<f64> = v.clone() / inner_dot;
                let h_inner: Array2<f64> = w.t()
                    .dot(&(inner_elop));

                let w1: Array2<f64> = Array::from_elem(
                    (v.shape()[0], 1), h.sum_axis(Axis(1))[0]);

                let inner_dot: Array2<f64> = w.dot(&h);
                let inner_elop: Array2<f64> = v.clone() / inner_dot;
                let w_inner: Array2<f64> = inner_elop.dot(&h.t());

                let h = h * &(h_inner / h1);
                let w = w * &(w_inner / w1);

                return (w, h)
            },
            _ => {
                process::exit(1)
            },
        }
    }

    fn objective_update(v: &Array2<f64>,
                        w: &Array2<f64>,
                        h: &Array2<f64>,
                        cons: &Array2<f64>,
                        old_cons: &Array2<f64>,
                        objective: &Objective) -> (f64, Option<Array2<f64>>) {
        match objective {
            Objective::Fro => {
                // Compute squared Frobenius norm of a target matrix and its NMF estimate.
                let r = v.clone() - w.dot(h);
                return ((r.clone() * r).sum(), None)
            },
            Objective::Div => {
                // Compute divergence of target matrix from its NMF estimate.
                let va = w.dot(h);
                let mut inner_elop = v.clone() / va.clone();
                inner_elop.par_mapv_inplace(|x|{x.ln()});

                return (((v.clone() * inner_elop) - v.clone() + va).sum(), None)
            },
            Objective::Conn => {
                // Compute connectivity matrix and compare it to connectivity matrix
                // from previous iteration.
                // Return logical value denoting whether connectivity matrix has changed
                // from previous iteration.
                let idx = Arc::new(Mutex::new(Array1::<usize>::zeros(
                    h.shape()[1])));

                h.axis_iter(Axis(1)).into_par_iter().enumerate()
                    .for_each(|(col_idx, col)|{
                        let notnan_row: Vec<NotNan<f64>> = col.into_par_iter().cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();
                        let max = notnan_row.par_iter().max().expect("No maximum found");
                        let argmax = notnan_row.par_iter().position(|element| element == max).unwrap();
                        let mut idx = idx.lock().unwrap();
                        idx[col_idx] = argmax;
                    });
                let idx = idx.lock().unwrap();
//                debug!("IDX: {:?}", idx);

                let mat1 = idx.broadcast((v.shape()[1], v.shape()[1])).unwrap();
                let mat2 = mat1.t();


                let mut new_cons: Array2<f64> = Array::zeros(
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
//                debug!("Eqaulity: {}", &new_cons==cons);
                let connectivity_change = connectivity_change.lock().unwrap().clone();
                return (connectivity_change, Some(new_cons))
            },
            Objective::None => {
                process::exit(1);
            },
        }
    }

    fn basis(&self) -> &Array2<f64> {
        match self {
            Factorization::NMF {
                w,
                ..
            } => {
                return w
            }
        }
    }

    fn target(&self) -> &Array2<f64> {
        match self {
            Factorization::NMF {
                v,
                ..
            } => {
                return v
            }
        }
    }

    fn coef(&self) -> &Array2<f64> {
        match self {
            Factorization::NMF {
                h,
                ..
            } => {
                return h
            }
        }
    }

    fn fitted(&self) -> Array2<f64> {
        match self {
            Factorization::NMF {
                w,
                h,
                ..
            } => {
                return w.dot(h)
            }
        }
    }

    fn distance(&self, metric: Objective) -> f64 {
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
                        let r = v.clone() - w.dot(h);
                        return (r.clone() * r).sum()
                    },
                    Objective::Div => {
                        // Compute divergence of target matrix from its NMF estimate.
                        let va = w.dot(h);
                        let mut inner_elop = v.clone() / va.clone();
                        inner_elop.par_mapv_inplace(|x| { x.ln() });

                        return ((v.clone() * inner_elop) - v.clone() + va).sum()
                    },
                    _ => {
                        process::exit(1)
                    }
                }
            }
        }
    }

    fn residuals(&self) -> Array2<f64> {
        match self {
            Factorization::NMF {
                v,
                w,
                h,
                ..
            } => {

                return v.clone() - w.dot(h)
            }
        }
    }

    fn predict(&self, what: &str) -> Array2<NotNan<f64>> {
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
                        x = h.to_owned()
                    },
                    "features" => {
                        x = w.to_owned()
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
                        let notnan_row: Vec<NotNan<f64>> = row
                            .into_par_iter()
                            .cloned()
                            .map(NotNan::new)
                            .filter_map(Result::ok)
                            .collect();

                        let max = notnan_row.par_iter().max().unwrap();
                        let argmax = notnan_row.par_iter().position(|element| element == max).unwrap();
                        idx[[col_idx, 0]] = *max;
                        idx[[col_idx, 1]] = NotNan::from(argmax as f64);
                    });

                let sums = x.sum_axis(Axis(0));

                let mut prob = Array::zeros(x.shape()[1]);

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

    fn probabilities(&self, what: &str) -> Array2<NotNan<f64>> {
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
                        x = h.to_owned()
                    },
                    "features" => {
                        x = w.to_owned()
                    },
                    _ => {
                        error!("Invalid option for prediction parameter {}", what);
                        process::exit(1)
                    }
                };

                let mut probabilities = Array::zeros(
                    (x.shape()[1], x.shape()[0]));

                let sums = x.sum_axis(Axis(0));

                return probabilities
            }
        }
    }

    fn rss(&self) -> f64 {
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

    fn evar(&self) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array2};
    use ndarray_linalg::{norm::*};
    use std::f64;

    #[test]
    fn test_update_euclidean() {
        let v = Array2::from_shape_vec((5, 5),
                                         vec![0.0, 0.5, 1.0, 1.5, 2.0,
                                                 0.5, 0.0, 0.5, 0.5, 1.5,
                                                 1.0, 0.5, 0.0, 0.5, 1.0,
                                                 1.5, 0.5, 0.5, 0.0, 0.5,
                                                 2.0, 1.5, 1.0, 0.5, 0.0]).unwrap();
        assert_eq!(v.norm(), 4.847679857416329);

        let placeholder = Array2::zeros((1, 1));


        let mut w = Array2::from_shape_vec((5, 2),
                                       vec![0.75, 0.85,
                                               0.65, 0.95,
                                               0.55, 0.05,
                                               0.45, 0.15,
                                               0.35, 0.25]).unwrap();
        let mut h = Array2::from_shape_vec((2, 5),
                                       vec![0.15, 0.25, 0.35, 0.45, 0.55,
                                            0.05, 0.95, 0.85, 0.75, 0.65]).unwrap();

        let w_t = w.t().to_owned();

        assert_eq!(w_t.dot(&v), Array2::from_shape_vec((2, 5),
                                                         vec![2.25, 1.4 , 1.65, 1.9000000000000001 , 3.2500000000000004,
                                                              1.25, 0.9 , 1.65, 1.9 , 3.25]).unwrap());
//        assert_eq!((h.clone() * w_t.dot(&v) / w_t.dot(&w.dot(&h))), Array2::from_shape_vec((2, 5),
//                                                                                    vec![1.07569721, 0.19787986, 0.32330301, 0.47401247, 0.98146877,
//                                                                                         0.20746888, 0.43045941, 0.71601787, 0.73786408, 1.109652  ]).unwrap());
        let (p_obj, _unused) = Factorization::objective_update(&v,
                                                                    &w,
                                                                    &h,
                                                               &placeholder,
                                                               &placeholder,
                                                               &Objective::Fro);


        let (mut w_ret, mut h_ret) = Factorization::update_wh(&v,
                                                                      w,
                                                                      h,
                                                                      &Update::Euclidean);

        assert_eq!(w_ret, Array2::from_shape_vec((5, 2), vec![0.6544092289292585, 0.9987459184909812,
                                                                        0.4476568628600942, 0.6446373448603425,
                                                                        0.8967482904791207, 0.07829888125304496,
                                                                        0.7511493120741213, 0.1702699613712722,
                                                                        0.7709963753667608, 0.4058245858289912]).unwrap());

        assert_eq!(h_ret, Array2::from_shape_vec((2, 5), vec![1.0756972111553786, 0.19787985865724383, 0.3233030090972708, 0.4740124740124741, 0.9814687714481813,
                                                                        0.20746887966804983, 0.43045940843297675, 0.7160178685386088, 0.7378640776699029, 1.1096520026263952]).unwrap());
        let (c_obj, _unused) = Factorization::objective_update(&v,
                                                                 &w_ret,
                                                                    &h_ret,
                                                               &placeholder,
                                                            &placeholder,
                                                            &Objective::Fro);
        let (mut w_ret, mut h_ret) = Factorization::update_wh(&v,
                                                              w_ret,
                                                              h_ret,
                                                              &Update::Euclidean);

        assert_eq!(w_ret, Array2::from_shape_vec((5, 2), vec![0.48035360626494766, 1.1090065732365098,
                                                                         0.399928424846599, 0.6137075300859243,
                                                                         0.8457693737152376, 0.07891922698649473,
                                                                         0.851590657857285, 0.1311655618922494,
                                                                         1.0057552021993001, 0.36825246357154434]).unwrap());

        assert_eq!(h_ret, Array2::from_shape_vec((2, 5), vec![1.3195125019368716, 0.40114707497220353, 0.3485167493278629, 0.4201726819390238, 0.7681395421952684,
                                                                        0.16044803005161773, 0.5403153136083633, 0.7985953508543852, 0.8099673349643689, 1.0792614312672162]).unwrap());
        assert_eq!(p_obj, 13.150625);
        assert_eq!(c_obj, 6.8164628097490425);

//        let mut x = 0.;
//        assert_eq!(f64::EPSILON, x)
    }
}