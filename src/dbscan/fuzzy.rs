use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::Itertools;
use model::variants::*;
use rayon::prelude::*;
use std::collections::HashSet;
use std::f64;
use std::hash::Hash;

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
            if self.pos == other.pos && self.tid == other.tid {
                return 40.;
            } else {
                let clr = |input: &Vec<f64>, geom_means: &Vec<f64>| -> Vec<f64> {
                    let output = input
                        .par_iter()
                        .enumerate()
                        .map(|(i, v)| ((v + 1.) / geom_means[i]).ln())
                        .collect();
                    return output;
                };

                //                debug!("row values {:?} geom_var {:?} geom_frq {:?}", &self.vars, &geom_var, &geom_frq);
                let row_vals: Vec<f64> =
                    clr(&self.vars.iter().map(|v| *v as f64).collect(), geom_var);

                let col_vals: Vec<f64> =
                    clr(&other.vars.iter().map(|v| *v as f64).collect(), geom_var);
                //                debug!("row values {:?} col values {:?}", &row_vals, &col_vals);

                let mean_row = get_mean(&row_vals);

                let mean_col = get_mean(&col_vals);

                // lovell et al. Phi and Phi distance: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075
                // Rho: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4870310/
                // propr log ratios to vlr and lr2rho: https://github.com/tpq/propr/blob/master/src/lr2propr.cpp

                let mut row_var = 0.;
                let mut col_var = 0.;
                let mut covar = 0.;

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

                debug!("Phi {}, Phi-D {}, concordance {}, Rho {}, Rho-M {}, row_var {} col_var {} covar {}",
                       phi, phi_dist, concordance, rho, 1.-rho, row_var, col_var, covar);

                return (phi_dist * (1. - jaccard));
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
        &self,
        points: &[P],
        initial_clusters: Vec<Cluster>,
        ref_name: &str,
    ) -> Vec<Cluster> {
        self.fuzzy_dbscan(points, initial_clusters, ref_name)
    }
}

impl FuzzyDBSCAN {
    fn fuzzy_dbscan<P: MetricSpace>(
        &self,
        points: &[P],
        initial_clusters: Vec<Cluster>,
        ref_name: &str,
    ) -> Vec<Cluster> {
        let mut clusters = Vec::new();
        let mut noise_cluster = Vec::new();
        let mut visited = vec![false; points.len()];
        let multi = MultiProgress::new();
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-");

        let pb1 = multi.insert(2, ProgressBar::new(points.len() as u64));
        pb1.set_style(sty.clone());

        let _ = std::thread::spawn(move || {
            multi.join_and_clear().unwrap();
        });
        pb1.set_message("Running fuzzy DBSCAN...");

        for initial in initial_clusters.iter() {
            // Work out which point in our initial clusters has the highest density to other
            // points
            // Set up multi progress bars

            let mut point_label_max = 0.;
            let mut point_index = 0;
            let mut point_neighbours = HashSet::new();

            for core_point in initial.iter() {
                if visited[core_point.index] {
                    continue;
                }
                visited[core_point.index] = true;
                let neighbour_indices = self.region_query(points, core_point.index);
                let point_label =
                    self.mu_min_p(self.density(core_point.index, &neighbour_indices, points));
                // check if new label is better than max
                if point_label > point_label_max {
                    point_label_max = point_label;
                    point_index = core_point.index;
                    point_neighbours = neighbour_indices;
                }
                pb1.inc(1);
            }
            // Extend neighbour indices with known connections
            let mut initial = initial
                .par_iter()
                .map(|link| link.index)
                .collect::<HashSet<usize>>();
            initial.remove(&point_index);
            point_neighbours.par_extend(initial.par_iter());

            // Cluster based on that core point
            clusters.push(self.expand_cluster_fuzzy(
                point_label_max,
                point_index,
                point_neighbours,
                points,
                &mut visited,
                &pb1,
            ));
        }

        // Cluster any unvisited points
        if initial_clusters.len() == 0
            || visited
                .par_iter()
                .any(|has_this_point_been_seen| has_this_point_been_seen == &false)
        {
            for point_index in 0..points.len() {
                if visited[point_index] {
                    pb1.inc(1);
                    continue;
                }
                visited[point_index] = true;
                let neighbour_indices = self.region_query(points, point_index);
                let point_label =
                    self.mu_min_p(self.density(point_index, &neighbour_indices, points));
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
                        neighbour_indices,
                        points,
                        &mut visited,
                        &pb1,
                    ));
                }
                pb1.inc(1);
            }
        }
        if !noise_cluster.is_empty() {
            debug!(
                "{} Variants Clustered as noise during Fuzzy DBSCAN",
                noise_cluster.len()
            );
            // clusters.push(noise_cluster);
        }

        // Sort the clusters by smallest to largest
        clusters.par_sort_by(|a, b| a.len().cmp(&b.len()));

        debug!("Clusters {:?}", clusters.len());
        pb1.set_message(&format!(
            "All points visited... {} initial clusters",
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

        pb1.finish_with_message(&format!(
            "Genome {}: {} possible strains... ",
            ref_name,
            clusters.len()
        ));

        clusters
    }

    fn expand_cluster_fuzzy<P: MetricSpace>(
        &self,
        point_label: f64,
        point_index: usize,
        mut neighbour_indices: HashSet<usize>,
        points: &[P],
        visited: &mut [bool],
        progress: &ProgressBar,
    ) -> Vec<Assignment> {
        let mut cluster = vec![Assignment {
            index: point_index,
            category: Category::Core,
            label: point_label,
        }];
        let mut border_points = Vec::new();
        let mut neighbour_visited = vec![false; points.len()];
        let multi = MultiProgress::new();
        let pb4 = multi.insert(2, ProgressBar::new_spinner());
        pb4.set_style(
            ProgressStyle::default_spinner()
                .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
                .template("[{elapsed_precise}] {prefix:.bold.dim} {spinner} {wide_msg}"),
        );

        let _ = std::thread::spawn(move || {
            multi.join_and_clear().unwrap();
        });

        pb4.set_message(&format!("Expanding cluster...",));

        while let Some(neighbour_index) = take_arbitrary(&mut neighbour_indices) {
            neighbour_visited[neighbour_index] = true;
            visited[neighbour_index] = true;
            progress.inc(1);
            let neighbour_neighbour_indices = self.region_query(points, neighbour_index);
            let neighbour_label =
                self.mu_min_p(self.density(neighbour_index, &neighbour_neighbour_indices, points));
            if neighbour_label >= self.phi {
                for neighbour_neighbour_index in neighbour_neighbour_indices {
                    if !neighbour_visited[neighbour_neighbour_index] {
                        neighbour_indices.insert(neighbour_neighbour_index);
                    }
                }
                cluster.push(Assignment {
                    index: neighbour_index,
                    category: Category::Core,
                    label: neighbour_label,
                });
            } else {
                border_points.push(Assignment {
                    index: neighbour_index,
                    category: Category::Border,
                    label: f64::MAX,
                });
            }
            pb4.inc(1);
        }
        pb4.finish();
        border_points.par_iter_mut().for_each(|border_point| {
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

fn get_mean(input: &Vec<f64>) -> f64 {
    let sum = input.par_iter().sum::<f64>();
    sum / input.len() as f64
}

#[allow(unused)]
fn propd(row_vals: Vec<f64>, col_vals: Vec<f64>) {
    let mean_row = get_mean(&row_vals);

    let mean_col = get_mean(&col_vals);

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
