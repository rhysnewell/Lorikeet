use std::error::Error;
use rayon::prelude::*;

use ndarray;
use crate::nmf_std;
//use crate::matrix_handling;

pub struct Nmf {
    v: ndarray::Array2<f32>,
    v1: Some(ndarray::Array2<f32>),
    h1: Some(ndarray::Array2<f32>),
    w: Some(ndarray::Array2<f32>),
    seed: String,
    rank: usize,
    n_run: usize,
    update: String,
    objective: String,
    conn_change: usize,
    max_iter: usize,
    min_residuals: usize,
}

impl Nmf {
    pub fn initialize(input: ndarray::Array2<f32>,
                      seed: &str,
                      rank: usize,
                      n_run: usize,
                      update: &str,
                      objective: &str,
                      conn_change: usize,
                      max_iter: usize,
                      min_residuals: usize) -> Self {
        Nmf {
            v: input,
            v1: None,
            h1: None,
            w: None,
            seed: seed.to_string(),
            rank,
            n_run,
            update: update.to_string(),
            objective: objective.to_string(),
            conn_change,
            max_iter,
            min_residuals,
        }
    }

    pub fn factorize(self) -> Result<Self, Box<dyn Error>> {

//        Compute matrix factorization.
//
//        Return fitted factorization model.
        (0..self.n_run).into_par_iter().for_each(|run|{
            
        });

        return Some(self)

    }
}