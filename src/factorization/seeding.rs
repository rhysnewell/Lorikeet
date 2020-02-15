use ndarray::{Array2, Array1, Axis, ArrayView, Ix1, prelude::*};
use ndarray_linalg::{SVD, convert::*, diagonal::*, Norm};
use std::sync::{Arc, Mutex};
use std::process;
use rayon::prelude::*;
use rand::{thread_rng, Rng};


#[derive(Debug, Clone, Copy)]
pub enum Seed {
    Nndsvd {
        rank: usize,
    },
    RandomVcol {
        rank: usize,
    },
    RandomC {
        rank: usize,
    },
    Random {
        rank: usize,
    },
    Fixed {
        rank: usize,
    },
    None,
}

impl Seed {
    pub fn new_nndsvd(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::Nndsvd {
            rank,
        }
    }

    pub fn new_random_vcol(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::RandomVcol {
            rank,
        }
    }

    pub fn new_random_c(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::RandomC {
            rank,
        }
    }

    pub fn new_random(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::Random {
            rank,
        }
    }

    pub fn new_fixed(rank: usize, v: &Array2<f32>) -> Seed {
        Seed::Fixed {
            rank,
        }
    }
}

pub trait SeedFunctions {
    fn initialize(&self, v: &Array2<f32>) -> (Array2<f32>, Array2<f32>);
}

impl SeedFunctions for Seed {
    fn initialize(&self, v: &Array2<f32>) -> (Array2<f32>, Array2<f32>) {
        match self {
            Seed::Nndsvd {
                rank,
            } => {
                let (u, s, e)
                    = v.svd(true, true).unwrap();
                info!("SVD calculation finished");
                let e = e.unwrap();
                let e = e.t();
                let u = u.unwrap();

                let mut w = Array2::zeros((v.shape()[0], *rank));
                let mut h = Array2::zeros((*rank, v.shape()[1]));

                // choose the first singular triplet to be nonnegative
//                let s = s.into_diag();
                debug!("S: {:?}", s);
                let mut u_slice = u.slice(s![.., 0]).to_owned();
                u_slice.par_mapv_inplace(|x| x.abs());
                let mut e_slice = e.slice(s![.., 0]).to_owned();
                e_slice = e_slice.t().to_owned();
                e_slice.par_mapv_inplace(|x| x.abs());

                w.slice_mut(s![.., 0]).assign(
                    &(s[0].powf(1. / 2.) * u_slice));
                h.slice_mut(s![0, ..]).assign(
                    &(s[0].powf(1. / 2.) * e_slice));

                // generate mutex guards around w and h
                let w_guard = Arc::new(Mutex::new(w.clone()));
                let h_guard = Arc::new(Mutex::new(h.clone()));


                // Update other factors based on associated svd factor
                (1..*rank).into_par_iter().for_each(|i|{
                    debug!("Inside Loop");

                    let uu = u.slice(s![.., i]).to_owned();
                    let vv = e.slice(s![.., i]).to_owned();
                    let mut uup = pos(&uu);
                    let mut uun = neg(&uu);
                    let vvp = pos(&vv);
                    let vvn = neg(&vv);
                    let n_uup = uup.norm();
                    let n_uun = uun.norm();
                    let n_vvp = vvp.norm();
                    let n_vvn = vvn.norm();
                    let termp = n_uup * n_vvp;
                    let termn = n_uun * n_vvn;

                    if termp >= termn {
                        debug!("First statement");
                        let mut w_guard = w_guard.lock().unwrap();
                        let mut h_guard = h_guard.lock().unwrap();

                        uup.par_mapv_inplace(|x| x * n_uup);
                        let mut vvp_t = vvp.t().to_owned();
                        vvp_t.par_mapv_inplace(|x| x * n_vvp);

                        w_guard.slice_mut(s![.., i]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (uup)));
                        h_guard.slice_mut(s![i, ..]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (vvp_t)));;
                    } else {
                        debug!("Second statement");
                        let mut w_guard = w_guard.lock().unwrap();
                        let mut h_guard = h_guard.lock().unwrap();

                        uun.par_mapv_inplace(|x| x * n_uun);
                        let mut vvn_t = vvn.t().to_owned();
                        vvn_t.par_mapv_inplace(|x| x * n_vvn);

                        w_guard.slice_mut(s![.., i]).assign(
                            &((s[i] * termn).powf(1. / 2.) / (uun)));
                        h_guard.slice_mut(s![i, ..]).assign(
                            &((s[i] * termn).powf(1. / 2.) / (vvn_t)));;
                    }
                });
                debug!("Outside Loop");
                let mut w_guard = w_guard.lock().unwrap();
                let mut h_guard = h_guard.lock().unwrap();

//                w_guard.par_mapv_inplace(|x|{
//                    if x < 1e-11 {
//                        0.
//                    } else {
//                        x
//                    }
//                });
//
//                h_guard.par_mapv_inplace(|x|{
//                    if x < 1e-11 {
//                        0.
//                    } else {
//                        x
//                    }
//                });

                let w = w_guard.clone();
                let h = h_guard.clone();

                debug!("H: {:?}", h);
                debug!("W: {:?}", w);
                return (w, h)

            },
            Seed::RandomVcol {
                rank
            } => {
                let p_c = (1. / 5. * v.shape()[1] as f32) as i32;
                let p_r = (1. / 5. * v.shape()[0] as f32) as i32;
                let mut h_guard = Arc::new(Mutex::new(Array2::zeros((*rank, v.shape()[1]))));
                let mut w_guard = Arc::new(Mutex::new(Array2::zeros((v.shape()[0], *rank))));

                let mut rng = thread_rng();

                (0..*rank).into_iter().for_each(|i| {
                    // retrieve random columns from input matrix
                    let mut random_cols = Array2::zeros((v.shape()[0], p_c as usize));
                    (0..p_c).into_iter().for_each(|idx| {
                        let col_id = rng.gen_range(0, v.shape()[1]);
                        let v_slice = v.slice(s![.., col_id]).to_owned();

                        random_cols.slice_mut(s![.., idx]).assign(&v_slice)
                    });
                    let mut w_guard = w_guard.lock().unwrap();

                    w_guard.slice_mut(s![.., i]).assign(
                        &random_cols.mean_axis(Axis(1)).unwrap());

                    // retrieve random rows from input matrix
                    let mut random_rows = Array2::zeros((p_r as usize, v.shape()[1]));
                    (0..p_r).into_iter().for_each(|idx| {
                        let row_id = rng.gen_range(0, v.shape()[0]);
                        let v_slice = v.slice(s![row_id, ..]).to_owned();
                        random_rows.slice_mut(s![idx, ..]).assign(&v_slice)
                    });
                    let mut h_guard = h_guard.lock().unwrap();

                    h_guard.slice_mut(s![i, ..]).assign(
                        &random_rows.mean_axis(Axis(0)).unwrap());;
                });
                let mut w_guard = w_guard.lock().unwrap();
                let mut h_guard = h_guard.lock().unwrap();

                let w = w_guard.clone();
                let h = h_guard.clone();

                return (w, h)
            },
            _ => process::exit(1)
        }
    }
}

fn pos(matrix: &Array1<f32>) -> Array1<f32> {
    let mut pos_mat = matrix.to_owned();
    pos_mat.par_mapv_inplace(|x| {
        if x >= 0. {
            1.
        } else {
            0.
        }
    });
    pos_mat * matrix
}

fn neg(matrix: &Array1<f32>) -> Array1<f32> {
    let mut neg_mat = matrix.to_owned();
    neg_mat.par_mapv_inplace(|x| {
        if x < 0. {
            1.
        } else {
            0.
        }
    });
    neg_mat * -matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1};
    use ndarray_linalg::svd;

    #[test]
    fn test_pos_and_neg() {
        let array = Array1::from_shape_vec((7),
                                           vec![23., -23., 0., 1., -1., -2., 2.]).unwrap();
        let pos_ret = pos(&array);
        let neg_ret = neg(&array);

        assert_eq!(pos_ret,
                   Array1::from_shape_vec((7),
                                          vec![23., 0., 0., 1., 0., 0., 2.]).unwrap());

        assert_eq!(neg_ret,
                   Array1::from_shape_vec((7),
                                          vec![0., 23., 0., 0., 1., 2., 0.]).unwrap())
    }

    #[test]
    fn test_svd() {
        let v = Array2::from_shape_vec((5, 5),
                                       vec![0.5, 0.5, 0.5, 0.5, 0.5,
                                            0.5, 0.5, 0.5, 0.5, 0.5,
                                            0.5, 0.5, 0.5, 0.5, 0.5,
                                            0.5, 0.5, 0.5, 0.5, 0.5,
                                            0.5, 0.5, 0.5, 0.5, 0.5]).unwrap();

        let (u, s, e)
            = v.svd(true, true).unwrap();
        println!("U: {:?}", u);
        println!("S: {:?}", s);
        println!("E: {:?}", e);

//        assert_eq!(1, 2);

    }
}