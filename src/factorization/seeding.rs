use ndarray::{Array2, Array1, Axis, ArrayView, Ix1};
use ndarray_linalg::{SVD, convert::*, diagonal::*, Norm};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

pub enum Seed {
    Nndsvd {
        rank: usize,
        w: Array2<f32>,
        h: Array2<f32>,
    }
}

impl Seed {
    pub fn new_nndsvd(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::Nndsvd {
            rank,
            w: Array2::zeros((v.shape()[0], rank)),
            h: Array2::zeros((rank, v.shape()[1])),
        }
    }
}

pub trait SeedFunctions {
    fn initialize(&mut self, v: &Array2<f32>);

    fn get_wh(&mut self) -> (&Array2<f32>, &Array2<f32>);
}

impl SeedFunctions for Seed {
    fn initialize(&mut self, v: &Array2<f32>) {
        match self {
            Seed::Nndsvd {
                ref mut rank,
                ref mut w,
                ref mut h,
            } => {
                let (u, s, e)
                    = v.svd(true, true).unwrap();
                let e = e.unwrap();
                let e = e.t();
                let u = u.unwrap();

                // choose the first singular triplet to be nonnegative
                let s = s.into_diag();
                w.slice_mut(s![.., 0]).assign(
                    &(s[0].powf(1. / 2.) * u.slice(s![.., 0]).mapv(|x| x.abs())));
                h.slice_mut(s![0, ..]).assign(
                    &(s[0].powf(1. / 2.) * e.slice(s![.., 0]).t().mapv(|x| x.abs())));

                // generate mutex guards around w and h
                let w_guard = Arc::new(Mutex::new(w.clone()));
                let h_guard = Arc::new(Mutex::new(h.clone()));

                // second svd for the other factors
                (1..*rank).into_par_iter().for_each(|i|{
                    let uu = u.slice(s![.., i]);
                    let vv = e.slice(s![.., i]);
                    let uup = pos(&uu);
                    let uun = neg(&uu);
                    let vvp = pos(&vv);
                    let vvn = neg(&vv);
                    let n_uup = uup.norm();
                    let n_uun = uun.norm();
                    let n_vvp = vvp.norm();
                    let n_vvn = vvn.norm();
                    let termp = n_uup * n_vvp;
                    let termn = n_uun * n_vvn;
                    if termp >= termn {
                        let mut w_guard = w_guard.lock().unwrap();
                        let mut h_guard = h_guard.lock().unwrap();
                        w_guard.slice_mut(s![.., i]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (uup.mapv(|x| x * n_uup))));
                        h_guard.slice_mut(s![i, ..]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (vvp.t().mapv(|x| x * n_vvp))));;
                    } else {
                        let mut w_guard = w_guard.lock().unwrap();
                        let mut h_guard = h_guard.lock().unwrap();
                        w_guard.slice_mut(s![.., i]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (uun.mapv(|x| x * n_uun))));
                        h_guard.slice_mut(s![i, ..]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (vvn.t().mapv(|x| x * n_vvn))));;
                    }
                });
                let w_guard = w_guard.lock().unwrap();
                let h_guard = h_guard.lock().unwrap();
                *w = w_guard.mapv(|x|{
                    if x < 1.0_f32.powf(-11.) {
                        0.
                    } else {
                        x
                    }
                });

                *h = h_guard.mapv(|x|{
                    if x < 1.0_f32.powf(-11.) {
                        0.
                    } else {
                        x
                    }
                });

            }
        }
    }

    fn get_wh(&mut self) -> (&Array2<f32>, &Array2<f32>) {
        match self {
            Seed::Nndsvd {
                w,
                h,
                ..
            } => {
                (w, h)
            }
        }
    }
}

fn pos(matrix: &ArrayView<f32, Ix1>) -> Array1<f32> {
    matrix.mapv(|x| {
        if x > 0. {
            1.
        } else {
            0.
        }
    }) * matrix
}

fn neg(matrix: &ArrayView<f32, Ix1>) -> Array1<f32> {
    matrix.mapv(|x| {
        if x < 0. {
            1.
        } else {
            0.
        }
    }) * matrix
}
