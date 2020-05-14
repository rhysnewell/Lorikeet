use std::collections::HashSet;
use std::hash::Hash;
use std::f64;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use model::variants::*;

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

/// An enum that decides whether a cluster should be updated or duplicated
#[derive(Debug, Clone)]
pub enum Update {
    Push,
    Clone,
}

pub fn update_clusters(original_cluster: &mut Cluster, variants: &Vec<Var>, var_to_check: &Var) -> Update {
    let skip = false;

    let to_do = Arc::new(
        Mutex::new(
            Update::Push));

    original_cluster.par_iter().for_each(|assignment|{
        if !skip {
            let var_in_cluster = &variants[assignment.index];
            if var_in_cluster.tid == var_to_check.tid && var_in_cluster.pos == var_to_check.pos {

            }
        }
    });
    let to_do = to_do.lock().unwrap().clone();
    return to_do
}

/// A trait to compute distances between points.
pub trait MetricSpace: Sized + Send + Sync {
    /// Returns the distance between `self` and `other`.
    fn distance(&self, other: &Self, geom_var: &Vec<f64>,
                geom_dep: &Vec<f64>, geom_frq: &Vec<f64>) -> f64;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Var {
    pub pos: i64,
    pub var: Variant,
    pub deps: Vec<i32>,
    pub vars: Vec<i32>,
    pub rel_abunds: Vec<f64>,
    pub tid: i32,
}

#[derive(Debug, Clone)]
pub struct Point {
    pub values: Vec<f64>,
}

impl MetricSpace for Var<> {
    fn distance(&self, other: &Self, geom_var: &Vec<f64>,
                geom_dep: &Vec<f64>, geom_frq: &Vec<f64>) -> f64 {

        if self.vars.len() > 1 {
            if self.pos == other.pos && self.tid == other.tid {
                return 10.
            } else {
                let clr = |input: &Vec<f64>, geom_means: &Vec<f64>| -> Vec<f64> {
                    let output = input.par_iter().enumerate().map(|(i, v)| {
                        ((v + 1.) / geom_means[i] as f64).ln()
                    }).collect();
                    return output
                };

                let get_mean = |input: &Vec<f64>| -> f64 {
                    let sum = input.par_iter().sum::<f64>();
                    sum / input.len() as f64
                };



                let row_vals: Vec<f64> = clr(&self.rel_abunds, geom_frq);

                let col_vals: Vec<f64> = clr(&other.rel_abunds, geom_frq);

                let mean_row = get_mean(&row_vals);

                let mean_col = get_mean(&col_vals);

                // lovell et al. Phi and Phi distance
                // https://github.com/tpq/propr/blob/17263691e9e84b778d3977fffd33b467e4205ef9/src/lr2propr.cpp
                // propr log ratios to vlr and lr2rho

                let mut row_var = 0.;
                let mut col_var = 0.;
                let mut covar = 0.;

                row_vals.iter()
                    .zip(col_vals.iter()).for_each(|(r_freq, c_freq)| {
                    row_var += (r_freq - mean_row).powf(2.);
                    col_var += (c_freq - mean_col).powf(2.);
                    covar += (r_freq - mean_row) * (c_freq - mean_col)
                });

                row_var = row_var / (row_vals.len() as f64 - 1.);
                col_var = col_var / (col_vals.len() as f64 - 1.);
                covar = covar / (row_vals.len() as f64 - 1.);

                let vlr = -2. * covar + row_var + col_var;
                let rho = 1. - vlr / (row_var + col_var);

                let phi = 1. + row_var / col_var -
                    2. * (row_var / col_var).sqrt()
                        * covar / (col_var * row_var).sqrt();

                let phi_dist = ((row_var / col_var).ln()).abs() + 2.0_f64.ln()
                    - (covar / (col_var * row_var).sqrt() + 1.).ln();

//                debug!("Phi {}, Phi-D {}, Rho {}, Rho-M {}", phi, phi_dist, rho, 1.-rho);

                return 1. - rho
            }
        } else {
            return if self.pos == other.pos && self.tid == other.tid {
                20.
            } else {
                let row_freq = self.vars[0] as f64;
                let col_freq = other.vars[0] as f64;

                let row_depth = self.deps[0] as f64;
                let col_depth = other.deps[0] as f64;

                let distance = (((row_freq / geom_var[0] as f64).ln() - (col_freq / geom_var[0] as f64).ln()).powf(2.)
                    + ((row_depth / geom_dep[0] as f64).ln() - (col_depth / geom_dep[0] as f64).ln()).powf(2.)).powf(1. / 2.);
                distance
            }
        }
    }
}

impl MetricSpace for Point {
    fn distance(&self, other: &Self, geom_var: &Vec<f64>,
                geom_dep: &Vec<f64>, geom_frq: &Vec<f64>) -> f64 {
        let mut sum_squares = Arc::new(Mutex::new(0.));
        self.values.par_iter().zip(other.values.par_iter()).for_each(|(x, y)|{
            let mut sum_squares = sum_squares.lock().unwrap();
            *sum_squares += (x - y).powf(2.)
        });
        let sum_squares: f64 = *sum_squares.lock().unwrap();
        let dist: f64 = sum_squares.sqrt();
        return dist
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
#[derive(Debug, Serialize, Clone, PartialEq)]
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
    /// The minimum threshold required for a label to become a Core point.
    pub phi: f64,
    /// The geometric mean of the depth of the variants across samples (as a vector).
    pub geom_var: Vec<f64>,
    /// The geometric mean of the total depth at each base across samples (as a vector).
    pub geom_dep: Vec<f64>,
    /// The geometric mean of the relative abundances of the variants across samples (as a vector).
    pub geom_frq: Vec<f64>,
}

impl FuzzyDBSCAN {
    /// Creates a new instance of the algorithm.
    pub fn new() -> Self {
        FuzzyDBSCAN {
            eps_min: f64::NAN,
            eps_max: f64::NAN,
            pts_min: f64::NAN,
            pts_max: f64::NAN,
            phi: f64::NAN,
            geom_var: Vec::new(),
            geom_dep: Vec::new(),
            geom_frq: Vec::new(),
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
            info!("{} Variants Clustered as noise during Fuzzy DBSCAN", noise_cluster.len());
            clusters.push(noise_cluster);
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
            if neighbor_label >= self.phi {
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
        border_points.par_iter_mut().for_each(|border_point|{
            for cluster_point in &cluster {
                let mu_distance =
                    self.mu_distance(&points[border_point.index], &points[cluster_point.index]);
                if mu_distance > 0.0 {
                    border_point.label =
                        cluster_point.label.min(mu_distance).min(border_point.label);
                }
            }
        });
        cluster.append(&mut border_points);
        cluster
    }

    fn region_query<P: MetricSpace>(&self, points: &[P], point_index: usize) -> HashSet<usize> {
        points
            .into_par_iter()
            .enumerate()
            .filter(|(neighbor_index, neighbor_point)| {
                *neighbor_index != point_index
                    && neighbor_point.distance(&points[point_index],
                                               &self.geom_var,
                                               &self.geom_dep,
                                               &self.geom_frq) <= self.eps_max
            }).map(|(neighbor_index, _)| neighbor_index)
            .collect() //TODO: would be neat to prevent this allocation.
    }

    fn density<P: MetricSpace>(
        &self,
        point_index: usize,
        neighbor_indices: &HashSet<usize>,
        points: &[P],
    ) -> f64 {
        let density: f64 = neighbor_indices.par_iter().fold(|| 0.0, |sum, &neighbor_index| {
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
        let distance = a.distance(b, &self.geom_var,
                                  &self.geom_dep, &self.geom_frq);
        if distance <= self.eps_min {
            1.0
        } else if distance > self.eps_max {
            0.0
        } else {
            (self.eps_max - distance) / (self.eps_max - self.eps_min)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

//    #[test]
//    fn test_clustering() {
//        let mut points =  vec![Var {
//
//        }];
//    }


}