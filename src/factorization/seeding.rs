use ndarray::{Array2, Array1, Axis, prelude::*};
use ndarray_linalg::{SVD, Norm};
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
    RandomSym {
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
    pub fn new_nndsvd(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::Nndsvd {
            rank,
        }
    }

    pub fn new_random_vcol(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::RandomVcol {
            rank,
        }
    }

    pub fn new_random_sym(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::RandomSym {
            rank,
        }
    }

    pub fn new_random_c(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::RandomC {
            rank,
        }
    }

    pub fn new_random(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::Random {
            rank,
        }
    }

    pub fn new_fixed(rank: usize, v: &Array2<f64>) -> Seed {
        Seed::Fixed {
            rank,
        }
    }
}

pub trait SeedFunctions {
    fn initialize(&self, v: &Array2<f64>) -> (Array2<f64>, Array2<f64>);
}

impl SeedFunctions for Seed {
    fn initialize(&self, v: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        match self {
            Seed::Nndsvd {
                rank,
            } => {
                let (u, s, e)
                    = v.svd(true, true).unwrap();
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
                        debug!("Map in place 1");
                        uup.par_mapv_inplace(|x| x * n_uup);
                        debug!("vvp_t");
                        let mut vvp_t = vvp.t().to_owned();
                        debug!("Map in place 2");
                        vvp_t.par_mapv_inplace(|x| x * n_vvp);

                        debug!("First slice");
                        w_guard.slice_mut(s![.., i]).assign(
                            &((s[i] * termp).powf(1. / 2.) / (uup)));

                        debug!("Second slice");
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

                w_guard.par_mapv_inplace(|x|{
                    if x < 1e-11 {
                        0.
                    } else {
                        x
                    }
                });

                h_guard.par_mapv_inplace(|x|{
                    if x < 1e-11 {
                        0.
                    } else {
                        x
                    }
                });

                let w = w_guard.clone();
                let h = h_guard.clone();

                debug!("H: {:?}", h);
                debug!("W: {:?}", w);
                return (w, h)

            },
            Seed::RandomVcol {
                rank
            } => {
                let p_c = (1. / 2. * v.shape()[1] as f64) as usize;
                let p_r = (1. / 2. * v.shape()[0] as f64) as usize;
                let h_guard = Arc::new(Mutex::new(Array2::zeros((*rank, v.shape()[1]))));
                let w_guard = Arc::new(Mutex::new(Array2::zeros((v.shape()[0], *rank))));

                let mut rng = thread_rng();
                let mut rc_indices = Array2::zeros((*rank, p_c));
                let mut rr_indices = Array2::zeros((p_r, *rank));
                for r in (0..*rank).into_iter() {
                    // generate random column indices
                    let mut random_indices_c = Array1::zeros(p_c);
                    for c in (0..p_c).into_iter() {
                        random_indices_c[c] = rng.gen_range(0, v.shape()[1]);
                    }
                    rc_indices.slice_mut(s![r, ..]).assign(&random_indices_c);

                    // Generate random row indices
                    let mut random_indices_r = Array1::zeros(p_r);
                    for r in (0..p_r).into_iter() {
                        random_indices_r[r] = rng.gen_range(0, v.shape()[0]);
                    }
                    rr_indices.slice_mut(s![.., r]).assign(&random_indices_r);
                }
                debug!("Random Col Ids {}", rc_indices);
                debug!("Random Row Ids {}", rr_indices);


                (0..*rank).into_par_iter().for_each(|i| {
                    // retrieve random columns from input matrix
                    let random_cols = Arc::new(
                        Mutex::new(
                            Array2::zeros((v.shape()[0], p_c))));

                    (0..p_c).into_iter().for_each(|idx| {
                        debug!("inner column loop");
                        let col_id = rc_indices[[i, idx]];
                        let v_slice = v.slice(s![.., col_id]).to_owned();

                        let mut random_cols = random_cols.lock().unwrap();
                        debug!("Slicing on idx {}", idx);
                        random_cols.slice_mut(s![.., idx]).assign(&v_slice)
                    });
                    let mut w_guard = w_guard.lock().unwrap();
                    let random_cols = random_cols.lock().unwrap();
                    debug!("Slicing in column means");
                    w_guard.slice_mut(s![.., i]).assign(
                        &random_cols.mean_axis(Axis(1)).unwrap());

                    // retrieve random rows from input matrix
                    let random_rows = Arc::new(
                        Mutex::new(
                            Array2::zeros((p_r, v.shape()[1]))));

                    (0..p_r).into_iter().for_each(|idx| {
                        debug!("inner row loop");
                        let row_id = rr_indices[[idx, i]];
                        let v_slice = v.slice(s![row_id, ..]).to_owned();

                        let mut random_rows = random_rows.lock().unwrap();
                        random_rows.slice_mut(s![idx, ..]).assign(&v_slice);
                    });
                    let mut h_guard = h_guard.lock().unwrap();
                    let random_rows = random_rows.lock().unwrap();
                    debug!("Slicing in row means");
                    h_guard.slice_mut(s![i, ..]).assign(
                        &random_rows.mean_axis(Axis(0)).unwrap());;
                });
                let w_guard = w_guard.lock().unwrap();
                let h_guard = h_guard.lock().unwrap();

                let w = w_guard.clone();
                let h = h_guard.clone();

                return (w, h)
            },
            Seed::RandomSym {
                rank
            } => {
                let p_c = (1. / 2. * v.shape()[1] as f64) as usize;
                let p_r = (1. / 2. * v.shape()[0] as f64) as usize;
                let mut h_guard = Arc::new(Mutex::new(Array2::zeros((*rank, v.shape()[1]))));
                let mut w_guard = Arc::new(Mutex::new(Array2::zeros((v.shape()[0], *rank))));

                // TODO: Since we are dealing with a symmetrical pairwise matrix,
                //       would it make sense to initialize with the same rows/columns for both
                //       W and H? i.e. rr_indices becomes the transpose matrix of rc_indices.

                let mut rng = thread_rng();
                let mut rc_indices = Array2::zeros((*rank, p_c));
//                let mut rr_indices = Array2::zeros((p_r, rank));
                for r in (0..*rank).into_iter() {
                    let mut random_indices = Array1::zeros(p_c);
                    for c in (0..p_c).into_iter() {
                        random_indices[c] = rng.gen_range(0, v.shape()[1]);
                    }
                    rc_indices.slice_mut(s![r, ..]).assign(&random_indices);
                }
                let rr_indices = rc_indices.t();
                debug!("Random Col Ids {}", rc_indices);


                (0..*rank).into_par_iter().for_each(|i| {
                    // retrieve random columns from input matrix
                    let random_cols = Arc::new(
                        Mutex::new(
                            Array2::zeros((v.shape()[0], p_c))));

                    (0..p_c).into_iter().for_each(|idx| {
                        debug!("inner column loop");
                        let col_id = rc_indices[[i, idx]];
                        let v_slice = v.slice(s![.., col_id]).to_owned();

                        let mut random_cols = random_cols.lock().unwrap();
                        debug!("Slicing on idx {}", idx);
                        random_cols.slice_mut(s![.., idx]).assign(&v_slice)
                    });
                    let mut w_guard = w_guard.lock().unwrap();
                    let random_cols = random_cols.lock().unwrap();
                    debug!("Slicing in column means");
                    w_guard.slice_mut(s![.., i]).assign(
                        &random_cols.mean_axis(Axis(1)).unwrap());

                    // retrieve random rows from input matrix
                    let random_rows = Arc::new(
                        Mutex::new(
                            Array2::zeros((p_r, v.shape()[1]))));

                    (0..p_r).into_iter().for_each(|idx| {
                        debug!("inner row loop");
                        let row_id = rr_indices[[idx, i]];
                        let v_slice = v.slice(s![row_id, ..]).to_owned();

                        let mut random_rows = random_rows.lock().unwrap();
                        random_rows.slice_mut(s![idx, ..]).assign(&v_slice);
                    });
                    let mut h_guard = h_guard.lock().unwrap();
                    let random_rows = random_rows.lock().unwrap();
                    debug!("Slicing in row means");
                    h_guard.slice_mut(s![i, ..]).assign(
                        &random_rows.mean_axis(Axis(0)).unwrap());;
                });
                let w_guard = w_guard.lock().unwrap();
                let h_guard = h_guard.lock().unwrap();

                let w = w_guard.clone();
                let h = h_guard.clone();

                return (w, h)
            },
            _ => process::exit(1)
        }
    }
}

fn pos(matrix: &Array1<f64>) -> Array1<f64> {
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

fn neg(matrix: &Array1<f64>) -> Array1<f64> {
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