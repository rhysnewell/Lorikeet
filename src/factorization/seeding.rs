use ndarray;
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
            w: ndarray::Array2::zeros((v.shape[0], rank)),
            h: ndarray::Array2::zeros((rank, v.shape[1])),
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