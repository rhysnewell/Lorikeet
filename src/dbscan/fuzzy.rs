use std::collections::HashSet;
use std::hash::Hash;
use std::f64;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use ndarray::{Array2};

fn take_arbitrary<T: Hash + Eq + Copy>(set: &mut HashSet<T>) -> Option<T> {
    let key_copy = if let Some(key_ref) = set.iter().next() {
        Some(*key_ref)
    } else {
        None
    };
    if let Some(key) = key_copy {
        set.take(&key)
    } else {
        None
    }
}

/// A trait to compute distances between points.
pub trait MetricSpace: Sized + Send + Sync {
    /// Returns the distance between `self` and `other`.
    fn distance(&self, other: &Self) -> f64;
}

#[derive(Debug, Clone)]
pub struct Point {
    pub pos: i32,
    pub var: String,
    pub geom_var: Vec<f64>,
    pub geom_dep: Vec<f64>,
    pub deps: Vec<f64>,
    pub vars: Vec<f64>,
    pub tid: i32,
}

pub fn dist_mat(input: &Vec<f64>) -> Array2<f64> {
    let mut output = Arc::new(Mutex::new(Array2::<f64>::zeros((input.len(), input.len()))));
    (0..input.len()-1).into_par_iter().for_each(|row_index|{
        (row_index+1..input.len()).into_par_iter().for_each(|col_index|{
            let euc_dist: f64 = (input[row_index] - input[col_index]).powf(2.).powf(1. / 2.);
            let mut output = output.lock().unwrap();
            output[[row_index, col_index]] = euc_dist;
            output[[col_index, row_index]] = euc_dist;
        });
    });
    let output = output.lock().unwrap().clone();
    return output

}

impl MetricSpace for Point {
    fn distance(&self, other: &Self) -> f64 {

        if self.vars.len() > 1 {
            if self.pos == other.pos && self.tid == other.tid {
                return 10.
            } else {
                let mut distance = 0.;
                let clr = |input: &Vec<f64>, geom_mean_var: &Vec<f64>| -> Vec<f64> {
                    let output = input.par_iter().enumerate().map(|(i, v)| {
                        ((v + 1.) / geom_mean_var[i] as f64).ln()
                    }).collect();
                    return output
                };


                let get_mean = |input: &Vec<f64>| -> f64 {
                    let sum = input.par_iter().sum::<f64>();
                    sum / input.len() as f64
                };

                let row_vals: Vec<f64> = clr(&self.vars, &self.geom_var);

                let col_vals: Vec<f64> = clr(&other.vars, &other.geom_var);

//                let mean_row = get_mean(&row_vals);
//
//                let mean_col = get_mean(&col_vals);
//
//                // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4870310/ eq. 2
//                // p = 2*cov(Ai, Aj) / (Var(Ai) + Var(Aj))
//
//                // Distance as https://en.wikipedia.org/wiki/Cosine_similarity
//                let mut row_var = 0.;
//                let mut col_var = 0.;
//                let mut covar = 0.;
//
//                row_vals.iter()
//                    .zip(col_vals.iter()).for_each(|(r_freq, c_freq)| {
//                    row_var += (r_freq - mean_row).powf(2.);
//                    col_var += (c_freq - mean_col).powf(2.);
//                    covar += (r_freq - mean_row) * (c_freq - mean_col)
//                });
//
//                row_var = row_var / row_vals.len() as f64;
//                col_var = col_var / col_vals.len() as f64;
//                covar = covar / row_vals.len() as f64;
//
//                if row_var == 0. && col_var == 0. {
//                    distance = 0.;
//                } else {
//                    distance = (2. * covar) / (row_var + col_var);
//                }

                let mut sum_of_diff = 0.;
                for (r, c) in row_vals.iter().zip(col_vals.iter()) {
                    sum_of_diff += (r - c).powf(2.)
                }
                distance = sum_of_diff.powf(1. / 2.);
//                debug!("Distance {}", distance);
//                 swap signs
//                distance *= -1.;
//                 move between 0 and 2
//                distance += 1.;
//                debug!("Distance {} Self {:?} Other {:?}", distance, &self, other);

//                let a_jk = dist_mat(&clr(&self.vars, &self.geom_var));
//                let b_jk = dist_mat(&clr(&other.vars, &other.geom_var));
//                let mut a_mean = 0.;
//                let mut b_mean = 0.;
//                a_jk.iter().for_each(|a| {a_mean += a});
//                b_jk.iter().for_each(|b| {b_mean += b});
//                a_mean = a_mean / a_jk.len() as f64;
//                b_mean = b_mean / b_jk.len() as f64;
//                let d_cov = Arc::new(Mutex::new(0.));
//                let d_var_a = Arc::new(Mutex::new(0.));
//                let d_var_b = Arc::new(Mutex::new(0.));
//                (0..a_jk.shape()[0]).into_par_iter().for_each(|j|{
//                    let a_j_mean: f64 = a_jk.slice(s![j, ..]).iter().sum::<f64>() / a_jk.slice(s![j, ..]).len() as f64;
//                    let b_j_mean: f64 = b_jk.slice(s![j, ..]).iter().sum::<f64>() / b_jk.slice(s![j, ..]).len() as f64;
//
//                    (0..a_jk.shape()[1]).into_par_iter().for_each(|k|{
//                        let a_k_mean: f64 = a_jk.slice(s![.., k]).iter().sum::<f64>() / a_jk.slice(s![.., k]).len() as f64;
//                        let b_k_mean: f64 = b_jk.slice(s![.., k]).iter().sum::<f64>() / b_jk.slice(s![.., k]).len() as f64;
//
//                        let A_jk: f64 = a_jk[[j, k]] - a_j_mean - a_k_mean + a_mean;
//                        let B_jk: f64 = b_jk[[j, k]] - b_j_mean - b_k_mean + b_mean;
//
//                        let mut d_cov = d_cov.lock().unwrap();
//                        *d_cov += A_jk * B_jk;
//
//                        let mut d_var_a = d_var_a.lock().unwrap();
//                        *d_var_a += A_jk.powf(2.);
//
//                        let mut d_var_b = d_var_b.lock().unwrap();
//                        *d_var_b += B_jk.powf(2.);
//
//                    })
//                });
//                let mut d_cov = d_cov.lock().unwrap().clone();
//                let mut d_var_a = d_var_a.lock().unwrap().clone();
//                let mut d_var_b = d_var_b.lock().unwrap().clone();
//                let n = self.vars.len() as f64;
//                d_cov = d_cov / n.powf(2.);
//                d_var_a = d_var_a / n.powf(2.);
//                d_var_b = d_var_b / n.powf(2.);
//
//                debug!("d_cov {} d_var_a {} d_var_b {}", d_cov, d_var_a, d_var_b);
//
//                let mut distance = d_cov / (d_var_a * d_var_b).powf(1. / 2.);
//                debug!("Distance {}", distance);



                return distance
            }
        } else {
            if self.pos == other.pos && self.tid == other.tid {
                return 20.
            } else {
                let row_freq = self.vars[0];
                let col_freq = other.vars[0];

                let row_depth = self.deps[0];
                let col_depth = other.deps[0];

                let distance = (((row_freq / self.geom_var[0] as f64).ln() - (col_freq / other.geom_var[0] as f64).ln()).powf(2.)
                    + ((row_depth / self.geom_dep[0] as f64).ln() - (col_depth / other.geom_dep[0] as f64).ln()).powf(2.)).powf(1. / 2.);

                return distance
            }
        }
    }
}

/// A high-level classification, as defined by the FuzzyDBSCAN algorithm.
#[derive(Debug, PartialEq, Serialize, Clone, Copy, Hash, Eq)]
pub enum Category {
    Core,
    Border,
    Noise,
}

/// An element of a [cluster](Cluster).
#[derive(Debug, Serialize, Clone)]
pub struct Assignment {
    /// The point index.
    pub index: usize,
    /// A (soft) label between `0.0` and `1.0`.
    pub label: f64,
    /// A high-level category.
    pub category: Category,
}

/// A group of [assigned](Assignment) points.
pub type Cluster = Vec<Assignment>;

/// An instance of the FuzzyDBSCAN algorithm.
///
/// Note that when setting `eps_min = eps_max` and `pts_min = pts_max` the algorithm will reduce to classic DBSCAN.
pub struct FuzzyDBSCAN {
    /// The minimum fuzzy local neighborhood radius.
    pub eps_min: f64,
    /// The maximum fuzzy local neighborhood radius.
    pub eps_max: f64,
    /// The minimum fuzzy neighborhood density (number of points).
    pub pts_min: f64,
    /// The maximum fuzzy neighborhood density (number of points).
    pub pts_max: f64,
}

impl FuzzyDBSCAN {
    /// Creates a new instance of the algorithm.
    pub fn new() -> Self {
        FuzzyDBSCAN {
            eps_min: f64::NAN,
            eps_max: f64::NAN,
            pts_min: f64::NAN,
            pts_max: f64::NAN,
        }
    }
}

impl FuzzyDBSCAN {
    /// Clusters a list of `points`.
    pub fn cluster<P: MetricSpace>(&self, points: &[P]) -> Vec<Cluster> {
        self.fuzzy_dbscan(points)
    }
}

impl FuzzyDBSCAN {
    fn fuzzy_dbscan<P: MetricSpace>(&self, points: &[P]) -> Vec<Cluster> {
        let mut clusters = Vec::new();
        let mut noise_cluster = Vec::new();
        let mut visited = vec![false; points.len()];
        for point_index in 0..points.len() {
            if visited[point_index] {
                continue;
            }
            visited[point_index] = true;
            let neighbor_indices = self.region_query(points, point_index);
            let point_label = self.mu_min_p(self.density(point_index, &neighbor_indices, points));
            if point_label == 0.0 {
                noise_cluster.push(Assignment {
                    index: point_index,
                    category: Category::Noise,
                    label: 1.0,
                });
            } else {
                clusters.push(self.expand_cluster_fuzzy(
                    point_label,
                    point_index,
                    neighbor_indices,
                    points,
                    &mut visited,
                ));
            }
        }
        if !noise_cluster.is_empty() {
            // we don't want noise, so we just ignore it for now
//            clusters.push(noise_cluster);
        }
        clusters
    }

    fn expand_cluster_fuzzy<P: MetricSpace>(
        &self,
        point_label: f64,
        point_index: usize,
        mut neighbor_indices: HashSet<usize>,
        points: &[P],
        visited: &mut [bool],
    ) -> Vec<Assignment> {
        let mut cluster = vec![Assignment {
            index: point_index,
            category: Category::Core,
            label: point_label,
        }];
        let mut border_points = Vec::new();
        let mut neighbor_visited = vec![false; points.len()];
        while let Some(neighbor_index) = take_arbitrary(&mut neighbor_indices) {
            neighbor_visited[neighbor_index] = true;
            visited[neighbor_index] = true;
            let neighbor_neighbor_indices = self.region_query(points, neighbor_index);
            let neighbor_label =
                self.mu_min_p(self.density(neighbor_index, &neighbor_neighbor_indices, points));
            if neighbor_label > 0.0 {
                for neighbor_neighbor_index in neighbor_neighbor_indices {
                    if !neighbor_visited[neighbor_neighbor_index] {
                        neighbor_indices.insert(neighbor_neighbor_index);
                    }
                }
                cluster.push(Assignment {
                    index: neighbor_index,
                    category: Category::Core,
                    label: neighbor_label,
                });
            } else {
                border_points.push(Assignment {
                    index: neighbor_index,
                    category: Category::Border,
                    label: f64::MAX,
                });
            }
        }
        for border_point in &mut border_points {
            for cluster_point in &cluster {
                let mu_distance =
                    self.mu_distance(&points[border_point.index], &points[cluster_point.index]);
                if mu_distance > 0.0 {
                    border_point.label =
                        cluster_point.label.min(mu_distance).min(border_point.label);
                }
            }
        }
        cluster.append(&mut border_points);
        cluster
    }

    fn region_query<P: MetricSpace>(&self, points: &[P], point_index: usize) -> HashSet<usize> {
        points
            .into_par_iter()
            .enumerate()
            .filter(|(neighbor_index, neighbor_point)| {
                *neighbor_index != point_index
                    && neighbor_point.distance(&points[point_index]) <= self.eps_max
            }).map(|(neighbor_index, _)| neighbor_index)
            .collect() //TODO: would be neat to prevent this allocation.
    }

    fn density<P: MetricSpace>(
        &self,
        point_index: usize,
        neighbor_indices: &HashSet<usize>,
        points: &[P],
    ) -> f64 {
        let mut density: f64 = neighbor_indices.par_iter().fold(|| 0.0, |sum, &neighbor_index| {
            sum + self.mu_distance(&points[point_index], &points[neighbor_index])
        }).sum();
        density + 1.0
    }

    fn mu_min_p(&self, n: f64) -> f64 {
        if n >= self.pts_max {
            1.0
        } else if n < self.pts_min {
            0.0
        } else {
            (n - self.pts_min) / (self.pts_max - self.pts_min)
        }
    }

    fn mu_distance<P: MetricSpace>(&self, a: &P, b: &P) -> f64 {
        let distance = a.distance(b);
        if distance <= self.eps_min {
            1.0
        } else if distance > self.eps_max {
            0.0
        } else {
            (self.eps_max - distance) / (self.eps_max - self.eps_min)
        }
    }
}