use dbscan::fuzzy;
use dbscan::fuzzy::Cluster;
use itertools::Itertools;
use model::variants::*;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::mpsc::channel;
use std::sync::Mutex;

/// Connects fuzzy DBSCAN clusters based on shared read information
#[allow(unused)]
pub fn linkage_clustering_of_clusters(
    clusters: &mut Vec<Vec<fuzzy::Assignment>>,
    variant_info: &Vec<fuzzy::Var>,
    variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
) {
    if clusters.len() > 1 {
    } else {
        // create placeholder jaccard hashmap when there is only one cluster to prevent nothing
        // from being returned
        let mut jaccard_map = HashMap::new();
        let mut input = HashMap::new();
        input.insert(1, 1.0);
        jaccard_map.insert(0, input);
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
) -> Vec<fuzzy::Cluster> {
    debug!("Phasing {} variants...", variant_info.len());
    if variant_info.len() > 1 {
        // Initiate the hashmap linking each variant to the variants it shares reads with
        let links = Mutex::new(HashMap::new());

        clusters.into_par_iter().for_each(|cluster| {
            let indices = cluster
                .par_iter()
                .map(|assignment| assignment.index)
                .collect::<HashSet<usize>>();

            // Loop through each variant in cluster and expand based on shared read content

            cluster.into_par_iter().for_each(|assigned_variant| {
                // Get variants by index
                let var1 = &variant_info[assigned_variant.index];

                variant_info
                    .par_iter()
                    .enumerate()
                    .for_each(|(index, var2)| {
                        if indices.contains(&index) {
                            // do nothing as cluster already contains variant
                        } else {
                            // Skipping over reference variants
                            if !(var1.tid == var2.tid && var1.pos == var2.pos)
                                && !(var1.var == Variant::None || var2.var == Variant::None)
                            {
                                // Read ids of first variant
                                let set1 = &var1.reads;

                                // Read ids of second variant
                                let set2 = &var2.reads;
                                // debug!("set 1 {:?}", &set1);
                                // debug!("set 2 {:?}", &set2);

                                // Add the jaccard's similarity to the hashmap for the two clusters
                                let intersection: HashSet<_> = set1.intersection(&set2).collect();

                                let union: HashSet<_> = set1.union(&set2).collect();

                                // Scaled Jaccard Similarity Based on Minimum Set size
                                //            let jaccard = intersection.len() as f64 /
                                //                std::cmp::min(set1.len() + 1, set2.len() + 1) as f64;
                                //             Normal Jaccard's Similarity
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
                                    let links_out = links.entry(assigned_variant.index).or_insert(
                                        [assigned_variant.index]
                                            .iter()
                                            .cloned()
                                            .collect::<BTreeSet<usize>>(),
                                    );
                                    links_out.insert(index);
                                    let links_out = links.entry(index).or_insert(
                                        [index].iter().cloned().collect::<BTreeSet<usize>>(),
                                    );
                                    links_out.insert(assigned_variant.index);
                                }
                            }
                            // pb1.inc(1);
                        }
                    });
            });
        });
        // pb1.finish_with_message("Initial variant networks formed...");

        let links = links.lock().unwrap();

        // Create condensed links sender and receiver, avoiding use of mutex
        let (condensed_links_s, condensed_links_r) = channel();

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
                if anchors.len() > anchor_size {
                    s.send(anchors).unwrap();
                }
                // pb1.inc(1);
            });
        // pb1.finish_with_message("Read networks established...");

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
fn get_condensed_index(i: usize, j: usize, n: usize) -> Option<usize> {
    if i == j {
        return None;
    } else {
        return Some(n * i - i * (i + 1) / 2 + j - 1 - i);
    }
}
