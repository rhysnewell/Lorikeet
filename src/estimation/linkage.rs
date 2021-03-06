use bird_tool_utils::command;
use dbscan::fuzzy;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::Itertools;
use model::variants::*;
use ndarray::prelude::*;
use ndarray_npy::{read_npy, write_npy};
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};

/// Connects fuzzy DBSCAN clusters based on shared read information
#[allow(unused)]
pub fn linkage_clustering_of_clusters(
    mut clusters: Vec<Vec<usize>>,
    mut cluster_separation: Array2<f64>,
    variant_info: &Vec<fuzzy::Var>,
    anchor_size: usize,
    min_cluster_distance: f64,
    minimum_reads_in_link: usize,
    multi: &Arc<MultiProgress>,
    ref_idx: usize,
    output_prefix: String,
    pts_max: f64,
    threads: usize,
) -> Vec<Vec<usize>> {
    clusters.par_sort_by(|a, b| a.len().cmp(&b.len()));
    if clusters.len() > 2 {
        // Set up multi progress bars
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}");
        let combinations = (0..clusters.len())
            .into_iter()
            .combinations(2)
            .collect::<Vec<Vec<usize>>>();
        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(combinations.len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Phasing variants within clusters...");

        let mut distances: Arc<Mutex<Array2<f64>>> = Arc::new(Mutex::new(Array::from_elem(
            (clusters.len(), clusters.len()),
            0.,
        )));

        // The average total allele depth of the variants within a cluster
        let mut depths: Arc<Mutex<HashMap<usize, f64>>> = Arc::new(Mutex::new(HashMap::new()));

        combinations.into_par_iter().for_each(|cluster_indices| {
            let cluster1_id = &cluster_indices[0];
            let cluster2_id = &cluster_indices[1];

            // boolean stating whether there is enough information to combine these clusters
            let mut extended = false;

            // Loop through each variant in cluster and expand based on shared read content
            let cluster1 = &clusters[*cluster1_id];
            let cluster2 = &clusters[*cluster2_id];

            // Reads in cluster 1 and 2 which do not clash with any variants
            let mut read_set1: HashSet<Vec<u8>> = HashSet::new();
            let mut read_set2: HashSet<Vec<u8>> = HashSet::new();

            // The reads that were found to clash in these clusters
            // let mut clash: HashSet<Vec<u8>> = HashSet::new();

            // The unique reads found to connect these clusters
            // let mut intersection_set = BTreeSet::new();
            let mut depth_1 = 0;
            let mut depth_2 = 0;

            cluster1.iter().for_each(|assigned_variant1| {
                // Get variants by index
                let var1 = &variant_info[*assigned_variant1];
                let set1 = &var1.reads;
                read_set1.par_extend(set1.clone());
                // {
                depth_1 += var1.vars.par_iter().sum::<i32>();
                // }
            });

            cluster2.iter().for_each(|assigned_variant2| {
                let var2 = &variant_info[*assigned_variant2];
                {
                    // Read ids of first variant

                    // Read ids of second variant
                    let set2 = &var2.reads;

                    read_set2.par_extend(set2.clone());

                    depth_2 += var2.vars.par_iter().sum::<i32>();

                    // pb1.inc(1);
                }
            });

            let intersection: HashSet<_> = read_set1.intersection(&read_set2).collect();
            // let mut union: HashSet<_> = read_set1.union(&read_set2).cloned().collect();
            // let union: HashSet<_> = union.union(&clash).collect();
            let cov_1 = depth_1 as f64 / cluster1.len() as f64;
            let cov_2 = depth_2 as f64 / cluster2.len() as f64;
            let mut fraction = 1.;
            // Take the fraction that results in a value < 1
            if cov_1 < cov_2 {
                fraction = cov_1 / cov_2;
            } else {
                fraction = cov_2 / cov_1;
            }

            // let jaccard = intersection.len() as f64 / union.len() as f64;
            let jaccard =
                intersection.len() as f64 / std::cmp::min(read_set1.len(), read_set2.len()) as f64;
            let jaccard_d = 1. - jaccard;
            let mut distances = distances.lock().unwrap();

            // use coverage information. If a cluster has higher coverage than another cluster
            // Then it is more likely to split across many clusters. The likelihood of two clusters
            // being put together is based on the ratio of their coverage
            distances[[*cluster1_id, *cluster2_id]] = fraction * jaccard_d;
            distances[[*cluster2_id, *cluster1_id]] = fraction * jaccard_d;
            // }

            let mut depths = depths.lock().unwrap();
            if !depths.contains_key(&cluster1_id) {
                depths.entry(*cluster1_id).or_insert(cov_1);
            }

            if !depths.contains_key(&cluster2_id) {
                depths.entry(*cluster2_id).or_insert(cov_2);
            }

            pb1.inc(1);
        });
        pb1.finish_and_clear();
        let mut distances = distances.lock().unwrap();
        let mut depths = depths.lock().unwrap();
        let mut total_cov = depths.values().sum::<f64>() / depths.len() as f64;

        let mut new_clusters: Vec<Vec<usize>> = Vec::new();
        let mut solo_clusters: Vec<Vec<usize>> = Vec::new();

        // all noise apparently
        // We will just check to see if certain clusters contain at least twice the cov
        // of another cluster
        // set of indices for clusters to verify that they have been checked previously
        let combinations = (0..clusters.len())
            .into_iter()
            .combinations(2)
            .collect::<Vec<Vec<usize>>>();
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.green/blue} {pos:>7}/{len:7} {msg}");
        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(combinations.len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Combining clusters...");

        // new_clusters = Vec::new();
        // Sets of clusters that have been joined
        let mut combined_sets = HashSet::new();
        let mut checked_single = HashSet::new();
        combinations.iter().for_each(|indices| {
            let cluster1_id = indices[0];
            let cluster2_id = indices[1];
            if true {
                // !checked_couple.contains(&vec![cluster1_id, cluster2_id]) {
                let cov_1 = depths.get(&cluster1_id).unwrap();
                let cov_2 = depths.get(&cluster2_id).unwrap();
                let distance = distances[[cluster1_id, cluster2_id]];
                let separation = cluster_separation[[cluster1_id, cluster2_id]];
                let mut combined_set = Vec::new();
                if distance <= min_cluster_distance || separation <= 0.5 {
                    let cluster1 = &clusters[cluster1_id];
                    let cluster2 = &clusters[cluster2_id];

                    let clash = check_for_clash(&variant_info, &cluster1, &cluster2);

                    if !clash {
                        // Begin generating the combined cluster
                        let mut combined = clusters[cluster1_id].clone();
                        combined.par_extend(clusters[cluster2_id].clone());
                        combined_set.push(cluster1_id);
                        combined_set.push(cluster2_id);
                        // find clusters that both have distances <= 0.5 to cluster1 and cluster2
                        for other_index in (0..clusters.len()).into_iter() {
                            let mut cluster_dists = Vec::new();
                            // Make sure we aren't checking a cluster twice
                            if combined_set.contains(&other_index) {
                                continue;
                            }

                            let mut combine = true;

                            let other_cluster = &clusters[other_index];
                            for cluster_id in combined_set.iter() {
                                let distance = distances[[other_index, *cluster_id]];
                                if distance <= min_cluster_distance {
                                    let clash = check_for_clash(
                                        &variant_info,
                                        other_cluster,
                                        &clusters[*cluster_id],
                                    );
                                    if !clash {
                                        cluster_dists.push(distance);
                                    } else {
                                        combine = false;
                                        break;
                                    }
                                }
                            }
                            if combine {
                                combined_set.push(other_index);
                                combined.par_extend(other_cluster.clone());
                            }
                        }
                        combined_set.par_sort();
                        if !combined_sets.contains(&combined_set) {
                            new_clusters.push(combined);
                            combined_sets.insert(combined_set);
                        }
                    }
                }
                // Check if a cluster should also be pushed by itself
                if !checked_single.contains(&cluster1_id) {
                    if clusters[cluster1_id].len() as f64 >= variant_info.len() as f64 * pts_max {
                        solo_clusters.push(clusters[cluster1_id].clone())
                    }
                    checked_single.insert(cluster1_id);
                }
                if !checked_single.contains(&cluster2_id) {
                    if clusters[cluster2_id].len() as f64 >= variant_info.len() as f64 * pts_max {
                        solo_clusters.push(clusters[cluster2_id].clone())
                    }
                    checked_single.insert(cluster2_id);
                }
            };
            pb1.inc(1);
        });

        pb1.finish_and_clear();

        // }

        // Sort and deduplicate clusters based on non-ref variant content
        new_clusters.par_extend(solo_clusters);
        new_clusters.par_sort_by_key(|x| x.len());
        new_clusters.par_iter_mut().for_each(|x| x.par_sort());

        // deduplicated hashset of variants stored in BTreeSet
        let mut dedup_clusters = BTreeSet::new();
        // All reference components of each cluster. Ordered.
        let mut reference = Vec::new();

        for cluster in new_clusters.iter() {
            let mut dedup = BTreeSet::new();
            let mut all_ref = BTreeSet::new();
            for ass in cluster.into_iter() {
                if variant_info[*ass].var != Variant::None {
                    dedup.insert(*ass);
                } else {
                    all_ref.insert(*ass);
                }
            }
            if !dedup_clusters.contains(&dedup) {
                dedup_clusters.insert(dedup);
                reference.push(all_ref);
            }
        }

        if dedup_clusters.len() == new_clusters.len() {
            return new_clusters;
        } else {
            // Reform the clusters, now deduplicated
            let mut new_clusters = Vec::new();
            for (variants, references) in dedup_clusters.iter().zip(reference.iter()) {
                let mut combined = Vec::new();
                combined.par_extend(variants);
                combined.par_extend(references);
                new_clusters.push(combined);
            }

            return new_clusters;
        }
    } else {
        // create placeholder jaccard hashmap when there is only one cluster to prevent nothing
        // from being returned
        // let mut jaccard_map = HashMap::new();
        // let mut input = HashMap::new();
        // input.insert(1, 1.0);
        // jaccard_map.insert(0, input);

        return clusters;
    }
}

fn check_for_clash(
    variant_info: &Vec<fuzzy::Var>,
    cluster1: &Vec<usize>,
    cluster2: &Vec<usize>,
) -> bool {
    let mut clash = false;
    for ass1 in cluster1.iter() {
        if clash {
            break;
        }
        let var1 = &variant_info[*ass1];

        for ass2 in cluster2.iter() {
            let var2 = &variant_info[*ass2];
            if var1.tid == var2.tid && var1.pos == var2.pos {
                clash = true;
                break;
            }
        }
    }
    return clash;
}

/// Connects variants into initial clusters based on shared read sets
#[allow(unused)]
pub fn linkage_clustering_of_variants(
    clusters: Vec<Vec<fuzzy::Assignment>>,
    variant_info: &Vec<fuzzy::Var>,
    anchor_size: usize,
    anchor_similarity: f64,
    minimum_reads_in_link: usize,
    multi: &Arc<MultiProgress>,
    ref_idx: usize,
) -> Vec<fuzzy::Cluster> {
    debug!("Phasing {} variants...", variant_info.len());
    if variant_info.len() > 1 {
        // Initiate the hashmap linking each variant to the variants it shares reads with
        let links = Mutex::new(HashMap::new());

        // Set up multi progress bars
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}");

        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(clusters.len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Phasing variants within clusters...");

        // let _ = std::thread::spawn(move || {
        //     multi.join_and_clear().unwrap();
        // });

        clusters.into_par_iter().for_each(|cluster| {
            let indices: BTreeSet<usize> = cluster
                .par_iter()
                .map(|assignment| assignment.index)
                .collect::<BTreeSet<usize>>();
            let extended = Arc::new(Mutex::new(false));

            // Loop through each variant in cluster and expand based on shared read content

            cluster.into_par_iter().for_each(|assigned_variant| {
                // Get variants by index
                let var1 = &variant_info[assigned_variant.index];

                variant_info
                    .par_iter()
                    .enumerate()
                    .for_each(|(index, var2)| {
                        if indices.contains(&index) {
                            if assigned_variant.index == index {
                                // do nothing as cluster already contains variant
                            } else {
                                // check for position conflict
                                // if (var1.tid == var2.tid && var1.pos == var2.pos)
                                // {
                                //     if var1.variant == Variant::MNV {
                                //
                                //     }
                                // }
                            }
                        } else {
                            // Skipping over reference variants
                            if !(var1.tid == var2.tid && var1.pos == var2.pos)
                                && !(var1.var == Variant::None || var2.var == Variant::None)
                            {
                                // Read ids of first variant
                                let set1 = &var1.reads;

                                // Read ids of second variant
                                let set2 = &var2.reads;

                                // Add the jaccard's similarity to the hashmap for the two clusters
                                let intersection: HashSet<_> = set1.intersection(&set2).collect();

                                let union: HashSet<_> = set1.union(&set2).collect();

                                if intersection.len() >= minimum_reads_in_link {
                                    // get relative frequencies of each Haplotype
                                    let pool_size = union.len() as f64;
                                    let x_11 = intersection.len() as f64 / pool_size;
                                    let p1 = set1.len() as f64 / pool_size;
                                    let q1 = set2.len() as f64 / pool_size;

                                    // Calculate Linkage D
                                    let dis = x_11 - p1 * q1;

                                    let mut links = links.lock().unwrap();
                                    // Intialize links for each indices including itself
                                    let links_out = links
                                        .entry(assigned_variant.index)
                                        .or_insert(indices.clone());
                                    links_out.insert(index);

                                    let mut extended = extended.lock().unwrap();
                                    *extended = true;
                                }
                            }
                            // pb1.inc(1);
                        }
                    });
            });
            if !*extended.lock().unwrap() {
                let mut links = links.lock().unwrap();
                // Intialize links for each indices including itself
                let links_out = links
                    .entry(*indices.iter().next().unwrap())
                    .or_insert(indices.clone());
            }
            pb1.inc(1);
        });
        pb1.finish_and_clear();

        let links = links.lock().unwrap();

        // Create condensed links sender and receiver, avoiding use of mutex
        let (condensed_links_s, condensed_links_r) = channel();
        // Set up multi progress bars
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}");

        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(links.len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Extending clusters based on read links...");

        debug!("pre filtering condensed links {:?}", &links);

        // extend the links for each anchor point by the union of all the indices
        links
            .par_iter()
            .for_each_with(condensed_links_s, |s, (_, current_links)| {
                let mut anchors = current_links.clone();

                current_links.iter().for_each(|index| {
                    match links.get(&index) {
                        Some(other_links) => {
                            // Get variants by index and check that there is no conflict in position
                            // with any established links in current links
                            let var_check_1 = &variant_info[*index];
                            let (check_s, check_r) = channel();
                            other_links.par_iter().for_each_with(
                                check_s,
                                |s_bool, index_to_check| {
                                    let var_check_2 = &variant_info[*index_to_check];
                                    s_bool
                                        .send(
                                            var_check_1.tid == var_check_2.tid
                                                && var_check_1.pos == var_check_2.pos,
                                        )
                                        .unwrap()
                                },
                            );
                            let mut checked_links: HashSet<_> = check_r.iter().collect();

                            if !checked_links.contains(&true) {
                                anchors.par_extend(other_links.par_iter())
                            }
                        }
                        _ => {}
                    }
                });

                // Filter out final links that aren't of size n
                // if anchors.len() > anchor_size {
                s.send(anchors).unwrap();
                // }
                pb1.inc(1);
            });
        pb1.set_message("Read networks established...");

        // Collect receiver into vec and sort by length
        let mut condensed_links: Vec<_> = condensed_links_r.iter().collect();
        condensed_links.par_sort_by_key(|key| key.len());

        // Filter highly similar sets
        condensed_links.dedup_by(|a, b| {
            let intersection: HashSet<_> = a.intersection(&b).collect();
            let union: HashSet<_> = a.union(&b).collect();
            let jaccard = intersection.len() as f64 / union.len() as f64;
            jaccard > anchor_similarity
        });

        pb1.set_message(&format!("{} distinct clusters...", condensed_links.len()));

        debug!("post filtering condensed links {:?}", &condensed_links);

        let (initial_clusters_s, initial_clusters_r) = channel();
        condensed_links
            .par_iter()
            .for_each_with(initial_clusters_s, |s, link_set| {
                s.send(
                    link_set
                        .par_iter()
                        .map(|link| fuzzy::Assignment {
                            index: *link,
                            label: 1.,
                            category: fuzzy::Category::Core,
                        })
                        .collect::<fuzzy::Cluster>(),
                )
                .unwrap()
            });
        let mut initial_clusters: Vec<fuzzy::Cluster> = initial_clusters_r.iter().collect();
        pb1.finish_and_clear();

        return initial_clusters;
    } else {
        Vec::new()
    }
}

/// Get all of the associated read ids for a given cluster
pub fn get_read_set(
    variants: &fuzzy::Cluster,
    variant_info: &Vec<fuzzy::Var>,
    variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
) -> HashSet<Vec<u8>> {
    let read_set = Mutex::new(HashSet::new());

    variants.par_iter().for_each(|assignment| {
        let variant = &variant_info[assignment.index];
        let base = &variant_map[&variant.tid][&variant.pos][&variant.var].reads;
        base.par_iter().for_each(|id| {
            let mut read_set = read_set.lock().unwrap();
            read_set.insert(id.clone());
        });
    });
    let read_set = read_set.lock().unwrap().clone();
    read_set
}

/// Extract the read ids associated with a particular variant
pub fn get_variant_set(
    variant: &fuzzy::Var,
    variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
) -> HashSet<Vec<u8>> {
    let variant_set = variant_map[&variant.tid][&variant.pos][&variant.var]
        .clone()
        .reads;

    return variant_set;
}

/// helper function to get the index of condensed matrix from it square form
#[allow(unused)]
fn get_condensed_index(i: usize, j: usize, n: usize) -> Option<usize> {
    if i == j {
        return None;
    } else {
        return Some(n * i - i * (i + 1) / 2 + j - 1 - i);
    }
}
