use dbscan::fuzzy;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::Itertools;
use model::variants::*;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};

/// Connects fuzzy DBSCAN clusters based on shared read information
#[allow(unused)]
pub fn linkage_clustering_of_clusters(
    mut clusters: Vec<Vec<fuzzy::Assignment>>,
    variant_info: &Vec<fuzzy::Var>,
    anchor_size: usize,
    anchor_similarity: f64,
    minimum_reads_in_link: usize,
    multi: &Arc<MultiProgress>,
    ref_idx: usize,
) -> Vec<fuzzy::Cluster> {
    clusters.par_sort_by(|a, b| a.len().cmp(&b.len()));
    if clusters.len() > 1 {
        // Initiate the hashmap linking each cluster that shares reads with another cluster
        // based on variant locations
        let mut links = HashMap::new();
        let mut original_clusters = HashMap::new();

        // Cluster IDs that seem to have all variants fully contained within another cluster
        // Based on coverage and read content
        let mut fully_contained = HashSet::new();

        // Set up multi progress bars
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}");

        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new((0..clusters.len())
            .into_iter()
            .permutations(2).collect::<Vec<Vec<usize>>>().len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Phasing variants within clusters...");

        // let _ = std::thread::spawn(move || {
        //     multi.join_and_clear().unwrap();
        // });

        (0..clusters.len())
            .into_iter()
            .permutations(2)
            .for_each(|cluster_indices| {
                let cluster1_id = &cluster_indices[0];
                let cluster2_id = &cluster_indices[1];

                // boolean stating whether there is enough information to combine these clusters
                let mut extended = false;

                // Loop through each variant in cluster and expand based on shared read content
                let cluster1 = &clusters[*cluster1_id];
                let cluster2 = &clusters[*cluster2_id];

                let mut indices: BTreeSet<usize> = cluster1
                    .par_iter()
                    .map(|assignment| assignment.index)
                    .collect::<BTreeSet<usize>>();
                // The variant indices that were found to clash in these clusters
                let mut clash = BTreeSet::new();

                // The unique reads found to connect these clusters
                let mut intersection_set = BTreeSet::new();
                let mut depth_1 = 0;
                let mut depth_2 = 0;

                cluster1.iter().for_each(|assigned_variant1| {
                    // Get variants by index
                    let var1 = &variant_info[assigned_variant1.index];
                    {
                        depth_1 +=
                            var1.vars.par_iter().sum::<i32>();
                    }
                    cluster2.iter().for_each(|assigned_variant2| {
                        let var2 = &variant_info[assigned_variant2.index];
                        {
                            if !(var1.tid == var2.tid && var1.pos == var2.pos)
                                && !(var2.var == Variant::None || var1.var == Variant::None)
                            {
                                // Read ids of first variant
                                let set1 = &var1.reads;

                                // Read ids of second variant
                                let set2 = &var2.reads;

                                // Add the jaccard's similarity to the hashmap for the two clusters
                                let intersection: BTreeSet<_> = set1.intersection(&set2).collect();

                                if intersection.len() > 0 {
                                    intersection_set =
                                        intersection_set.union(&intersection).cloned().collect();
                                    if intersection_set.len() >= minimum_reads_in_link {
                                        extended = true;
                                    }
                                }

                                depth_2 +=
                                    var2.vars.par_iter().sum::<i32>();
                                // let union: HashSet<_> = set1.union(&set2).collect();
                                //
                                // if intersection.len() >= minimum_reads_in_link {
                                //     // get relative frequencies of each Haplotype
                                //     let pool_size = union.len() as f64;
                                //     let x_11 = intersection.len() as f64 / pool_size;
                                //     let p1 = set1.len() as f64 / pool_size;
                                //     let q1 = set2.len() as f64 / pool_size;
                                //
                                //     // Calculate Linkage D
                                //     let dis = x_11 - p1 * q1;
                                // }
                            } else {
                                // Send to the clash pile
                                clash.insert(assigned_variant2.index);
                            }
                            // pb1.inc(1);
                        }
                    });
                });
                if extended {
                    let mut extend = BTreeSet::new();
                    for variant in cluster2.iter() {
                        if !clash.contains(&variant.index) {
                            extend.insert(variant.index);
                        };
                    }
                    // Intialize links for each indices including itself
                    if !fully_contained.contains(cluster1_id) {
                        let mut links_out = links.entry(*cluster1_id).or_insert(indices.clone());
                        *links_out = links_out.union(&extend).cloned().collect();
                        original_clusters.entry(*cluster1_id).or_insert(indices.clone());
                    } else if links.contains_key(cluster1_id) {
                        links.remove(cluster1_id).unwrap();
                    }

                    // Calculate mean coverage for variants
                    let coverage_1 = (depth_1 as f64) / cluster1.len() as f64;
                    let coverage_2 = (depth_2 as f64) / cluster2.len() as f64;

                    // If cluster1 has higher coverage, then it likely can appear again later
                    // Either by itself or with another.
                    if coverage_1 >= coverage_2*1.25 && cluster2.len() < cluster1.len() {
                        if original_clusters.contains_key(cluster2_id) {
                            // Remove the cluster since it seems like it can't appear by itself
                            original_clusters.remove(cluster2_id).unwrap();
                            fully_contained.insert(*cluster2_id);
                        }

                        if links.contains_key(cluster2_id) {
                            links.remove(cluster2_id).unwrap();
                        }
                        fully_contained.insert(*cluster2_id);
                    } else if original_clusters.contains_key(cluster1_id) {
                        // Remove the cluster since it seems like it can't appear by itself
                        original_clusters.remove(cluster1_id).unwrap();
                        // fully_contained.insert(*cluster1_id);
                    }
                } else {
                    // Intialize links for each indices including itself
                    if !fully_contained.contains(cluster1_id) {
                        let mut links_out = links.entry(*cluster1_id).or_insert(indices.clone());
                        // original_clusters.entry(*cluster1_id).or_insert(indices.clone());
                    } else if links.contains_key(cluster1_id) {
                        links.remove(cluster1_id).unwrap();
                    }
                }
                pb1.inc(1);
            });
        pb1.finish_and_clear();

        // Create condensed links sender and receiver, avoiding use of mutex
        let (condensed_links_s, condensed_links_r) = channel();
        // Set up multi progress bars
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}");

        let pb1 = multi.insert(ref_idx + 3, ProgressBar::new(links.len() as u64));
        pb1.set_style(sty.clone());

        pb1.set_message("Filtering extended clusters...");

        debug!("pre filtering condensed links {:?}", &links);

        // extend the links for each anchor point by the union of all the indices
        links
            .par_iter()
            .for_each_with(condensed_links_s, |s, (cluster_id, current_links)| {
                let mut anchors = current_links.clone();
                if original_clusters.contains_key(&cluster_id) {
                    let original = original_clusters.get(&cluster_id).unwrap();
                    s.send(original.clone()).unwrap();
                }

                // Filter out final links that aren't of size n
                // if anchors.len() > anchor_size {
                s.send(anchors).unwrap();
                // }
                pb1.inc(1);
            });
        pb1.set_message("Cluster networks established...");

        // Collect receiver into vec and sort by length
        let mut condensed_links: BTreeSet<_> = condensed_links_r.iter().collect();

        let mut condensed_links: Vec<_> = condensed_links.iter().collect();
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
        // create placeholder jaccard hashmap when there is only one cluster to prevent nothing
        // from being returned
        // let mut jaccard_map = HashMap::new();
        // let mut input = HashMap::new();
        // input.insert(1, 1.0);
        // jaccard_map.insert(0, input);

        return Vec::new();
    }
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
