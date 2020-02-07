use ndarray;
use ndarray_linalg::{SVD, convert::*, diagonal::*};
use rayon::prelude::*;

pub enum Seed {
    Nndsvd {
        rank: usize,
        w: ndarray::Array2<f32>,
        h: ndarray::Array2<f32>,
    }
}

impl Seed {
    pub fn new_nndsvd(rank: usize, v: &ndarray::Array2<f32>) -> Seed {
        Seed::Nndsvd {
            rank,
            w: ndarray::Array2::zeros((v.shape()[0], rank)),
            h: ndarray::Array2::zeros((rank, v.shape()[1])),
        }
    }
}

pub trait SeedFunctions {
    fn initialize(&mut self, v: &ndarray::Array2<f32>);

    fn pos(&self, matrix: &ndarray::Array2<f32>) -> ndarray::Array2<f32>;

    fn neg(&self, matrix: &ndarray::Array2<f32>) -> ndarray::Array2<f32>;
}

impl SeedFunctions for Seed {
    fn initialize(&mut self, v: &ndarray::Array2<f32>) {
        match self {
            Seed::Nndsvd {
                ref mut rank,
                ref mut w,
                ref mut h,
            } => {
                let (mut u, mut s, mut e)
                    = v.svd(calc_u: True, calc_vt: True).unwrap();
                e = e.t();
                s = s.diag();
                w.slice_mut(s![.., 0]).assign(s[0].powf(1. / 2.) * u.slice(s![.., 0].mapv_inplace(|x| x.abs())));
                h.slice_mut(s![0, ..]).assign(s[0].powf(1. / 2.) * e.slice(s![.., 0]).t().mapv_inplace(|x| x.aabs()));
            }
        }
    }

    fn pos(&self, matrix: &ndarray::Array2<f32>) -> ndarray::Array2<f32> {
        match self {
            Seed::Nndsvd {
                ref mut rank,
                ref mut w,
                ref mut h,
            } => {
                matrix.mapv(|x| x > 0.) * matrix
            }
        }
    }

    fn neg(&self, matrix: &ndarray::Array2<f32>) -> ndarray::Array2<f32> {
        match self {
            Seed::Nndsvd {
                ref mut rank,
                ref mut w,
                ref mut h,
            } => {
                matrix.mapv(|x| x < 0.) * matrix
            }
        }
    }
}