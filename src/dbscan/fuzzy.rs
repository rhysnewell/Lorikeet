use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::Itertools;
use model::variants::*;
// use rand::Rng;
use rayon::prelude::*;
use std::collections::HashSet;
use std::f64;
use std::hash::Hash;
use std::sync::Arc;

#[allow(unused)]
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
    fn distance(
        &self,
        other: &Self,
        geom_var: &Vec<f64>,
        geom_dep: &Vec<f64>,
        geom_frq: &Vec<f64>,
    ) -> f64;

    fn clash(&self, other: &Self) -> bool;

    fn var(&self) -> &Variant;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Var {
    pub pos: i64,
    pub var: Variant,
    pub deps: Vec<i32>,
    pub vars: Vec<i32>,
    //    pub rel_abunds: Vec<f64>,
    pub tid: i32,
    pub reads: HashSet<Vec<u8>>,
}

#[allow(unused)]
impl MetricSpace for Var {
    fn distance(
        &self,
        other: &Self,
        geom_var: &Vec<f64>,
        geom_dep: &Vec<f64>,
        geom_frq: &Vec<f64>,
    ) -> f64 {
        if self.vars.len() > 1 {
            {
                let clr = |input: &Vec<f64>,
                           geom_means: &Vec<f64>,
                           indices: &HashSet<usize>|
                 -> Vec<f64> {
                    let mut output_clr = Vec::with_capacity(indices.len());

                    for (i, v) in input.iter().enumerate() {
                        if indices.contains(&i) {
                            let new_v = ((v + 1.) / geom_means[i]).ln();
                            output_clr.push(new_v);
                        }
                    }
                    return output_clr;
                };

                let get_indices = |row_vals: &Vec<i32>, col_vals: &Vec<i32>| -> HashSet<usize> {
                    let mut to_return = HashSet::new();
                    let mut found_one_zero = false;
                    for (index, (val_one, val_two)) in row_vals.iter().zip(col_vals).enumerate() {
                        if val_one > &0 || val_two > &0 {
                            to_return.insert(index);
                        } else if val_one == &0 && val_two == &0 && !found_one_zero {
                            to_return.insert(index);
                            found_one_zero = true;
                        }
                        to_return.insert(index);
                    }

                    return to_return;
                };
                let indices = get_indices(&self.vars, &other.vars);
                // println!("indices {:?} row_vals {:?} col vals {:?}", &indices, &self.vars, &other.vars,);
                let row_vals: Vec<f64> = clr(
                    &self.vars.iter().map(|v| *v as f64).collect(),
                    geom_var,
                    &indices,
                );

                let col_vals: Vec<f64> = clr(
                    &other.vars.iter().map(|v| *v as f64).collect(),
                    geom_var,
                    &indices,
                );
                //                debug!("row values {:?} col values {:?}", &row_vals, &col_vals);

                let mean_row = get_mean(&row_vals, &indices);

                let mean_col = get_mean(&col_vals, &indices);

                // lovell et al. Phi and Phi distance: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075
                // Rho: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4870310/
                // propr log ratios to vlr and lr2rho: https://github.com/tpq/propr/blob/master/src/lr2propr.cpp

                let mut row_var = 0.;
                let mut col_var = 0.;
                let mut covar = 0.;

                row_vals.iter().zip(col_vals.iter()).enumerate().for_each(
                    |(index, (r_freq, c_freq))| {
                        if indices.contains(&index) {
                            row_var += (r_freq - mean_row).powf(2.);
                            col_var += (c_freq - mean_col).powf(2.);
                            covar += (r_freq - mean_row) * (c_freq - mean_col)
                        }
                    },
                );

                row_var = row_var / (row_vals.len() as f64 - 1.);
                col_var = col_var / (col_vals.len() as f64 - 1.);
                covar = covar / (row_vals.len() as f64 - 1.);

                let vlr = -2. * covar + row_var + col_var;
                let rho = 1. - vlr / (row_var + col_var);

                let phi = 1. + row_var / col_var
                    - 2. * (row_var / col_var).sqrt() * covar / (col_var * row_var).sqrt();

                let phi_dist = ((row_var / col_var).ln()).abs() + 2.0_f64.ln()
                    - (covar / (col_var * row_var).sqrt() + 1.).ln();

                // reciprocal concordance correlation coeffecient between 0 and 2
                let concordance = -2. * covar / (row_var + col_var) + 1.;

                // Read ids of first variant
                let set1 = &self.reads;

                // Read ids of second variant
                let set2 = &other.reads;

                // Add the jaccard's similarity to the hashmap for the two clusters
                let intersection: HashSet<_> = set1.intersection(&set2).collect();
                let jaccard = intersection.len() as f64 / (set1.len() + set2.len() + 1) as f64;
                // if self.pos == 3439
                //     || self.pos == 948
                //     || self.pos == 8107
                //     || self.pos == 7660
                //     || other.pos == 3440
                //     || other.pos == 949
                //     || other.pos == 8108
                //     || other.pos == 7661 {
                // if self.var == Variant::None && other.var == Variant::None {
                //     println!("{}: {:?} {}: {:?} Phi {}, Phi-D {}, concordance {}, Rho {}, Rho-M {}, row_var {} col_var {} covar {} indices {:?} deps: {:?} {:?}",
                //              self.pos, self.var, other.pos, other.var, phi, phi_dist, concordance, rho, 1. - rho, row_var, col_var, covar, &indices, self.vars, other.vars);
                //     }
                // }

                return (phi * (1. - jaccard));
            }
        } else {
            return if self.pos == other.pos && self.tid == other.tid {
                // arbitrarily high distance
                20.
            } else {
                let row_freq = self.vars[0] as f64;
                let col_freq = other.vars[0] as f64;

                let row_depth = self.deps[0] as f64;
                let col_depth = other.deps[0] as f64;

                // Calculate Aitchinson Distance
                let distance = (((row_freq / geom_var[0] as f64).ln()
                    - (col_freq / geom_var[0] as f64).ln())
                .powf(2.)
                    + ((row_depth / geom_dep[0] as f64).ln()
                        - (col_depth / geom_dep[0] as f64).ln())
                    .powf(2.))
                .powf(1. / 2.);
                debug!("Distance {}", distance);
                distance
            };
        }
    }

    fn clash(&self, other: &Self) -> bool {
        if self.tid == other.tid && self.pos == other.pos && self.var != other.var {
            match &self.var {
                Variant::MNV(mnv) => match &other.var {
                    Variant::SNV(snv) => {
                        let pos_in_mnv = (self.pos as i64 - other.pos as i64).abs() as usize;
                        if &mnv[pos_in_mnv] == snv {
                            false
                        } else {
                            true
                        }
                    }
                    _ => true,
                },
                Variant::SNV(snv) => match &other.var {
                    Variant::MNV(mnv) => {
                        let pos_in_mnv = (self.pos as i64 - other.pos as i64).abs() as usize;
                        if &mnv[pos_in_mnv] == snv {
                            false
                        } else {
                            true
                        }
                    }
                    _ => true,
                },
                Variant::None => match &other.var {
                    Variant::None => false,
                    _ => true,
                },
                _ => true,
            }
        } else {
            false
        }
    }

    fn var(&self) -> &Variant {
        &self.var
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

impl Assignment {
    pub fn new() -> Assignment {
        Assignment {
            index: 0,
            label: 0.,
            category: Category::Noise,
        }
    }
}

/// A group of [assigned](Assignment) points.
pub type Cluster = Vec<Assignment>;

trait Dedup<T: PartialEq + Clone> {
    fn clear_duplicates(&mut self);
}

impl<T: PartialEq + Clone> Dedup<T> for Vec<Cluster> {
    fn clear_duplicates(&mut self) {
        let mut already_seen = vec![];
        self.retain(|item| match already_seen.contains(item) {
            true => false,
            _ => {
                already_seen.push(item.clone());
                true
            }
        })
    }
}

/// An instance of the FuzzyDBSCAN algorithm.
///
/// Note that when setting `eps_min = eps_max` and `pts_min = pts_max` the algorithm will reduce to classic DBSCAN.
pub struct FuzzyDBSCAN {
    /// The minimum fuzzy local neighbourhood radius.
    pub eps_min: f64,
    /// The maximum fuzzy local neighbourhood radius.
    pub eps_max: f64,
    /// The minimum fuzzy neighbourhood density (number of points).
    pub pts_min: f64,
    /// The maximum fuzzy neighbourhood density (number of points).
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

// impl FuzzyDBSCAN {
//     /// Creates a new instance of the algorithm.
//     pub fn new() -> Self {
//         FuzzyDBSCAN {
//             eps_min: f64::NAN,
//             eps_max: f64::NAN,
//             pts_min: f64::NAN,
//             pts_max: f64::NAN,
//             phi: f64::NAN,
//             geom_var: Vec::new(),
//             geom_dep: Vec::new(),
//             geom_frq: Vec::new(),
//         }
//     }
// }

impl FuzzyDBSCAN {
    /// Clusters a list of `points`.
    pub fn cluster<P: MetricSpace>(
        &mut self,
        points: &[P],
        _initial_clusters: Vec<Cluster>,
        ref_name: &str,
        multi: &Arc<MultiProgress>,
        ref_idx: usize,
    ) -> Vec<Cluster> {
        self.fuzzy_dbscan(points, _initial_clusters, ref_name, multi, ref_idx)
    }
}

#[allow(unused)]
impl FuzzyDBSCAN {
    fn fuzzy_dbscan<P: MetricSpace>(
        &mut self,
        points: &[P],
        _initial_clusters: Vec<Cluster>,
        ref_name: &str,
        multi: &Arc<MultiProgress>,
        ref_idx: usize,
    ) -> Vec<Cluster> {
        let mut acceptable_distribution = false;
        let mut clusters = Vec::new();

        let sty = ProgressStyle::default_bar().template(
            "[{elapsed_precise}] {bar:40.green/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]",
        );

        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(points.len() as u64));
        pb1.set_style(sty.clone());

        let multi_inner = Arc::clone(&multi);
        let max_iter = 10;
        let mut niter = 0;

        pb1.set_message(&format!("{}: Running fuzzy DBSCAN...", ref_name));
        'outer: while !acceptable_distribution && niter < max_iter {
            let mut visited = vec![false; points.len()];
            let mut noise_cluster = Vec::new();

            {
                for point_index in 0..points.len() {
                    if visited[point_index] {
                        continue;
                    }
                    visited[point_index] = true;
                    pb1.inc(1);

                    let mut neighbour_indices = self.region_query(points, point_index);
                    let point_label =
                        self.mu_min_p(self.density(point_index, &neighbour_indices, points));
                    if point_label == 0.0 {
                        if points[point_index].var() == &Variant::None {
                            noise_cluster.push(Assignment {
                                index: point_index,
                                category: Category::Noise,
                                label: 1.0,
                            });
                        }
                    } else {
                        'expand: loop {
                            pb1.set_message(&format!("{}: Expanding cluster...", ref_name,));
                            let expansion = self.expand_cluster_fuzzy(
                                point_label,
                                point_index,
                                &mut neighbour_indices,
                                points,
                                &mut visited,
                                &multi_inner,
                                &pb1,
                                ref_idx,
                            );

                            match expansion {
                                Some(expanded) => {
                                    let vars: HashSet<&Variant> = expanded
                                        .par_iter()
                                        .map(|p| points[p.index].var())
                                        .collect();
                                    if vars.len() == 1 && vars.contains(&Variant::None) {
                                        noise_cluster.par_extend(expanded);
                                    } else {
                                        clusters.push(expanded);
                                    }
                                    pb1.set_message(&format!("{}: Cluster pushed.", ref_name,));
                                    pb1.inc(1);
                                    break 'expand;
                                }
                                None => {
                                    // Cluster parameters were too slack, tighten them up
                                    // let mut rng = rand::thread_rng();
                                    // let scaler = rng.gen_range(0.5, 0.6);
                                    let scaler = 0.5;
                                    self.eps_min = self.eps_min * scaler;
                                    self.eps_max = self.eps_max * scaler;

                                    pb1.reset();
                                    niter += 1;
                                    if niter < max_iter {
                                        clusters = Vec::new();
                                    }
                                    pb1.set_message(&format!(
                                        "{}: Clash in points in cluster. Trying next core point...",
                                        ref_name,
                                    ));
                                    continue 'outer;
                                }
                            }
                        }
                    }
                }
            }

            // Sort the clusters by smallest to largest
            clusters.par_sort_by(|a, b| a.len().cmp(&b.len()));

            debug!("Clusters {:?}", clusters.len());
            pb1.set_message(&format!(
                "{}: All points visited... {} initial clusters",
                ref_name,
                clusters.len()
            ));

            // Deduplicate clusters by sorted indices
            clusters.dedup_by(|cluster1, cluster2| {
                // Remove clusters that are completely contained within another cluster
                let dedup_key_1 = cluster1
                    .par_iter()
                    .map(|var| var.index)
                    .collect::<HashSet<usize>>();

                let dedup_key_2 = cluster2
                    .par_iter()
                    .map(|var| var.index)
                    .collect::<HashSet<usize>>();

                let intersection: HashSet<_> = dedup_key_1.intersection(&dedup_key_2).collect();

                let small_similarity = intersection.len() / dedup_key_1.len();

                debug!(
                    "dedup_key_1 {:?} dedup_key_2 {:?} Similarity {}",
                    dedup_key_1, dedup_key_2, small_similarity,
                );
                small_similarity == 1
            });

            // We use RNG here so we don't end up with a situation where it recursively goes back
            // and forth between two cluster states (no clusters, 1 cluster), which is theoretically possible
            if clusters.len() == 0 {
                // Cluster params to strict, loosen them up
                // let mut rng = rand::thread_rng();
                // let scaler = rng.gen_range(2.5, 5.0);
                let scaler = 1.3;
                self.eps_min = self.eps_min * scaler;
                self.eps_max = self.eps_max * scaler;
                pb1.set_message(&format!(
                    "{}: No viable clusters. Adjusting parameters...",
                    ref_name,
                ));
                pb1.reset();
                niter += 1;
                if niter < max_iter {
                    clusters = Vec::new();
                }
            // clusters =
            //     fuzzy_scanner.fuzzy_dbscan(points, initial_clusters, ref_name, multi, ref_idx);
            } else if clusters.len() == 1 && clusters[0].len() == points.len() {
                // Cluster parameters were to slack, tighten them up
                // let mut rng = rand::thread_rng();
                // let scaler = rng.gen_range(0.5, 0.7);
                let scaler = 0.8;
                self.eps_min = self.eps_min * scaler;
                self.eps_max = self.eps_max * scaler;
                pb1.set_message(&format!(
                    "{}: Clustering too lenient. Adjusting parameters...",
                    ref_name,
                ));
                pb1.reset();
                niter += 1;
                if niter < max_iter {
                    clusters = Vec::new();
                }
            } else if clusters.len() >= (points.len() / 2) && points.len() > 10 {
                // Each point likely formed it's own cluster, have to be careful when there are
                // only a few variants though
                // let mut rng = rand::thread_rng();
                // let scaler_eps = rng.gen_range(1.2, 1.5);
                // let scaler = rng.gen_range(1.5, 2.0);
                let scaler = 1.25;
                self.eps_min = self.eps_min * scaler;
                self.eps_max = self.eps_max * scaler;
                // self.pts_min = self.pts_min;
                // self.pts_max = self.pts_max;
                pb1.set_message(&format!(
                    "{}: Too many clusters. Adjusting parameters...",
                    ref_name,
                ));
                pb1.reset();
                niter += 1;
                if niter < max_iter {
                    clusters = Vec::new();
                }
            } else {
                if !noise_cluster.is_empty() {
                    debug!(
                        "{} Variants Clustered as noise during Fuzzy DBSCAN",
                        noise_cluster.len()
                    );
                    clusters.push(noise_cluster);
                }
                acceptable_distribution = true;
            } // What about when there is far too many clusters? Maybe this is fixed by the read phasing
        }
        pb1.finish_and_clear();

        clusters
    }

    fn expand_cluster_fuzzy<P: MetricSpace>(
        &self,
        point_label: f64,
        point_index: usize,
        mut neighbour_indices: &mut HashSet<usize>,
        points: &[P],
        visited: &mut [bool],
        multi: &MultiProgress,
        progress: &ProgressBar,
        ref_idx: usize,
    ) -> Option<Vec<Assignment>> {
        let mut cluster = vec![Assignment {
            index: point_index,
            category: Category::Core,
            label: point_label,
        }];
        let mut border_points = Vec::new();
        let mut neighbour_visited = vec![false; points.len()];
        let pb4 = multi.insert(ref_idx + 4, ProgressBar::new_spinner());
        pb4.set_style(
            ProgressStyle::default_spinner().template("[{elapsed_precise}] {spinner:.red} {msg}"),
        );

        pb4.set_message(&format!("Expanding cluster...",));

        while let Some(neighbour_index) = take_arbitrary(&mut neighbour_indices) {
            neighbour_visited[neighbour_index] = true;
            let mut visited_point = &mut visited[neighbour_index];
            if visited_point == &mut false {
                progress.inc(1);
                *visited_point = true;
            }
            let neighbour_neighbour_indices = self.region_query(points, neighbour_index);
            let neighbour_label =
                self.mu_min_p(self.density(neighbour_index, &neighbour_neighbour_indices, points));
            if neighbour_label >= self.phi {
                for neighbour_neighbour_index in neighbour_neighbour_indices {
                    if !neighbour_visited[neighbour_neighbour_index] {
                        neighbour_indices.insert(neighbour_neighbour_index);
                    }
                }
                // if self.check_for_clash(points, &cluster, neighbour_index) {
                // This suggests the cluster is too lenient. Bringing in opposing variants
                // return none and then try again with update parameters.
                // visited[neighbour_index] = false;
                // pb4.finish_and_clear();
                // return None;
                // continue 'expand;
                // } else {
                // progress.inc(1);
                cluster.push(Assignment {
                    index: neighbour_index,
                    category: Category::Core,
                    label: neighbour_label,
                });
            // }
            } else {
                // if self.check_for_clash(points, &cluster, neighbour_index) {
                // This suggests the cluster is too lenient. Bringing in opposing variants
                // return none and then try again with update parameters.
                // visited[neighbour_index] = false;
                // pb4.finish_and_clear();
                // return None;
                // continue 'expand;
                // } else {
                // progress.inc(1);
                border_points.push(Assignment {
                    index: neighbour_index,
                    category: Category::Border,
                    label: f64::MAX,
                });
                // }
            }
            pb4.inc(1);
        }
        pb4.finish_and_clear();
        border_points.par_iter_mut().for_each(|border_point| {
            let border = &points[border_point.index];
            for cluster_point in &cluster {
                let core = &points[cluster_point.index];
                let mu_distance = self.mu_distance(border, core);
                if mu_distance > 0.0 {
                    border_point.label =
                        cluster_point.label.min(mu_distance).min(border_point.label);
                }
            }
        });
        cluster.append(&mut border_points);
        // if cluster.len() > points.len() / 2 {
        //     None
        // } else {
        Some(cluster)
        // }
    }

    fn check_for_clash<P: MetricSpace>(
        &self,
        points: &[P],
        cluster: &Vec<Assignment>,
        point_index: usize,
    ) -> bool {
        let point_to_check = &points[point_index];
        // let mut clash = false;
        cluster.par_iter().any(|assignment| {
            let point_in_cluster = &points[assignment.index];
            point_to_check.clash(point_in_cluster)
        })
    }

    fn region_query<P: MetricSpace>(&self, points: &[P], point_index: usize) -> HashSet<usize> {
        points
            .into_par_iter()
            .enumerate()
            .filter(|(neighbour_index, neighbour_point)| {
                *neighbour_index != point_index
                    && neighbour_point.distance(
                        &points[point_index],
                        &self.geom_var,
                        &self.geom_dep,
                        &self.geom_frq,
                    ) <= self.eps_max
            })
            .map(|(neighbour_index, _)| neighbour_index)
            .collect() //TODO: would be neat to prevent this allocation.
    }

    fn density<P: MetricSpace>(
        &self,
        point_index: usize,
        neighbour_indices: &HashSet<usize>,
        points: &[P],
    ) -> f64 {
        let density: f64 = neighbour_indices
            .par_iter()
            .fold(
                || 0.0,
                |sum, &neighbour_index| {
                    sum + self.mu_distance(&points[point_index], &points[neighbour_index])
                },
            )
            .sum();
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
        let distance = a.distance(b, &self.geom_var, &self.geom_dep, &self.geom_frq);
        if distance <= self.eps_min {
            1.0
        } else if distance > self.eps_max {
            0.0
        } else {
            (self.eps_max - distance) / (self.eps_max - self.eps_min)
        }
    }
}

fn get_mean(input: &Vec<f64>, indices: &HashSet<usize>) -> f64 {
    let sum = input.par_iter().sum::<f64>();
    sum / indices.len() as f64
}

#[allow(unused)]
fn propd(row_vals: Vec<f64>, col_vals: Vec<f64>, indices: &HashSet<usize>) {
    let mean_row = get_mean(&row_vals, indices);

    let mean_col = get_mean(&col_vals, indices);

    // lovell et al. Phi and Phi distance: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075
    // Rho: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4870310/
    // propr log ratios to vlr and lr2rho: https://github.com/tpq/propr/blob/master/src/lr2propr.cpp

    let mut row_var = 0.;
    let mut col_var = 0.;
    let mut covar = 0.;

    // Differential proportionality requires the calculation of VLR for different groupings of
    // our measured variables: https://mran.microsoft.com/snapshot/2018-03-29/web/packages/propr/vignettes/e_differential.html

    for combos in (0..row_vals.len()).into_iter().combinations(2) {
        let ind_1 = combos[0];
        let ind_2 = combos[1];
    }
    row_vals
        .iter()
        .zip(col_vals.iter())
        .for_each(|(r_freq, c_freq)| {
            row_var += (r_freq - mean_row).powf(2.);
            col_var += (c_freq - mean_col).powf(2.);
            covar += (r_freq - mean_row) * (c_freq - mean_col)
        });

    row_var = row_var / (row_vals.len() as f64 - 1.);
    col_var = col_var / (col_vals.len() as f64 - 1.);
    covar = covar / (row_vals.len() as f64 - 1.);

    let vlr = -2. * covar + row_var + col_var;
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
