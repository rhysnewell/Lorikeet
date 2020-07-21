use kodama::{Method, linkage};
use std::collections::{HashMap, HashSet, BTreeSet};
use std::sync::mpsc::channel;
use model::variants::*;
use dbscan::fuzzy;
use rayon::prelude::*;
use itertools::Itertools;
use std::sync::{Arc, Mutex};


/// Connects fuzzy DBSCAN clusters based on shared read information
#[allow(unused)]
pub fn linkage_clustering_of_clusters(
                      clusters: &Vec<Vec<fuzzy::Assignment>>,
                      variant_info: &Vec<fuzzy::Var>,
                      variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>)
                      -> (Vec<Vec<fuzzy::Assignment>>, HashMap<usize, HashMap<usize, f64>>, Vec<f64>) {

    if clusters.len() > 1 {
        let clusters_changed = Arc::new(Mutex::new(Vec::new()));
        let clusters_shared_reads = Arc::new(Mutex::new(HashMap::new()));
        let jaccard_distances = Arc::new(Mutex::new(
            vec![0.; (clusters.len().pow(2) - clusters.len()) / 2])
        );
        let changed = Arc::new(Mutex::new(0));
        let total_clusters = clusters.len();
        let clusters = Arc::new(Mutex::new(clusters.clone()));

        // Loop through each permutation of 2 clusters and observe shared variants in reads
        (0..total_clusters).into_iter()
            .combinations(2)
            .collect::<Vec<Vec<usize>>>().into_par_iter().for_each(|indices| {

            // Get clusters by index
            let clust1;
            let clust2;
            let clust1_size;
            let clust2_size;
            {
                let clusters = clusters.lock().unwrap();
                let clust1_h = clusters[indices[0]].clone();
                let clust2_h = clusters[indices[1]].clone();
                clust1_size = clust1_h.len();
                clust2_size = clust2_h.len();
                clust1 = Arc::new(Mutex::new(clust1_h));
                clust2 = Arc::new(Mutex::new(clust2_h));
            }
            // Total read set for each cluster
            let clust1_set: Arc<Mutex<HashSet<Vec<u8>>>> = Arc::new(Mutex::new(HashSet::new()));
            let clust2_set: Arc<Mutex<HashSet<Vec<u8>>>> = Arc::new(Mutex::new(HashSet::new()));

            // Cluster assignment index to avoid reappending to sets
            let clust1_index = Arc::new(Mutex::new(HashSet::new()));
            let clust2_index = Arc::new(Mutex::new(HashSet::new()));

            // Loop through first cluster
            (0..clust1_size).into_par_iter().for_each(|index1| {
                // Loop through second cluster
                let assignment1;
                {
                    let clust1 = clust1.lock().unwrap();
                    assignment1 = clust1[index1].clone();
                }
                (0..clust2_size).into_par_iter().for_each(|index2| {
                    let assignment2;
                    {
                        let clust2 = clust2.lock().unwrap();
                        assignment2 = clust2[index2].clone();
                    }
                    if assignment1.index != assignment2.index {
                        let var1 = &variant_info[assignment1.index];
                        let var2 = &variant_info[assignment2.index];

                        // Read ids of first variant
                        let set1 = &var1.reads;

                        // Read ids of second variant
                        let set2 = &var2.reads;

//                        debug!("Read IDs {:?} {:?}", set1, set2);
                        let read_intersection: HashSet<_> = set1
                            .intersection(&set2).collect();

                        // Extend each cluster read set
                        if read_intersection.len() == 0 {
                            let mut clust1_index =
                                clust1_index.lock().unwrap();
                            if !clust1_index.contains(&assignment1.index) {
                                let mut clust1_set =
                                    clust1_set.lock().unwrap();
                                clust1_set.par_extend(set1.clone());

                                clust1_index.insert(assignment1.index);
                            };
                            let mut clust2_index =
                                clust2_index.lock().unwrap();
                            if !clust2_index.contains(&assignment2.index) {
                                let mut clust2_set =
                                    clust2_set.lock().unwrap();
                                clust2_set.par_extend(set2.clone());

                                clust2_index.insert(assignment2.index);
                            };
                        } else if read_intersection.len() > 0 {
                            // Append both variant ids as we are going to extend each cluster
                            // Since the variants are connected by at least one read
//                            debug!("Read IDs {:?} {:?}", set1, set2);
                            let mut clust1_index =
                                clust1_index.lock().unwrap();
                            if !clust1_index.contains(&assignment1.index) {
                                let mut clust1_set =
                                    clust1_set.lock().unwrap();
                                clust1_set.par_extend(set1.clone());
                                clust1_index.insert(assignment1.index);
                            };

                            if !clust1_index.contains(&assignment2.index) {
                                let mut clust1_set =
                                    clust1_set.lock().unwrap();
                                clust1_set.par_extend(set2.clone());
                                clust1_index.insert(assignment2.index);
                                let mut clusters = clusters.lock().unwrap();
                                clusters[indices[0]].push(assignment2.clone());
                            };

                            let mut clust2_index =
                                clust2_index.lock().unwrap();
                            if !clust2_index.contains(&assignment2.index) {
                                let mut clust2_set =
                                    clust2_set.lock().unwrap();
                                clust2_set.par_extend(set2.clone());
                                clust2_index.insert(assignment2.index);
                            };
                            if !clust2_index.contains(&assignment1.index) {
                                let mut clust2_set =
                                    clust2_set.lock().unwrap();
                                clust2_set.par_extend(set1.clone());
                                clust2_index.insert(assignment1.index);
                                let mut clusters = clusters.lock().unwrap();
                                clusters[indices[1]].push(assignment1.clone());
                            };
                        }
                    }
                });
            });

            // Add the jaccard's similarity to the hashmap for the two clusters
            let clust1_set = clust1_set.lock().unwrap();
            let clust2_set = clust2_set.lock().unwrap();
            let intersection: HashSet<_> = clust1_set
                .intersection(&clust2_set).collect();

            let mut clusters_shared_reads = clusters_shared_reads
                .lock().unwrap();
            let cluster_map = clusters_shared_reads.entry(indices[0])
                .or_insert(HashMap::new());

            // Scaled Jaccard Similarity Based on Minimum Set size
            let jaccard = intersection.len() as f64 /
                std::cmp::min(clust1_set.len() + 1, clust2_set.len() + 1) as f64;

            debug!("Intersection Size {} {:?} {} {}", intersection.len(), indices, jaccard, 1. - jaccard);

            // Normal Jaccard's Similarity
//                let jaccard = intersection.len() as f64 /
//                    (clust1_set.len() + clust2_set.len() + 1) as f64;
            let mut jaccard_distances = jaccard_distances.lock().unwrap();
            // Check to see if we have two or more samples
            if jaccard_distances.len() > 1 {
                // Place each measure at appropriate index
                let clusters = clusters.lock().unwrap();
                let j_index = get_condensed_index(
                    indices[0], indices[1], clusters.len()).unwrap_or_else(|| 0);
                let dist = 1. - jaccard;
                jaccard_distances[j_index]
                    = dist;

            } else {
                // If only 1 sample, then just put in first index
                jaccard_distances[0] = 1. - jaccard;
            }
            // Add jaccard similarity to hashmap
            cluster_map.entry(indices[1]).or_insert(jaccard);
        });

        let clusters_shared_reads = clusters_shared_reads.lock().unwrap().clone();
        let mut jaccard_distances = jaccard_distances.lock().unwrap().to_vec();
        let clusters = clusters.lock().unwrap();
        // Perform HAC using kodama
        let dend = linkage(&mut jaccard_distances,
                           clusters.len(), Method::Single);

        debug!("Dendogram {:?} {:?} {:?}", &dend, jaccard_distances, clusters_shared_reads);
        // Step through each step in the dendrogram and combine clusters
        // that are under a certain dissimilarity
        dend.steps().into_par_iter().for_each(|step| {
            // Check to see that these are leaf clusters
            if step.cluster1 <= clusters.len() - 1 && step.cluster2 <= clusters.len() - 1 {
                // combine clusters
                if step.dissimilarity < 0.2 {
                    let mut new_cluster = Vec::new();
                    new_cluster.par_extend(clusters[step.cluster1].par_iter().cloned());
                    new_cluster.par_extend(clusters[step.cluster2].par_iter().cloned());

                    let mut clusters_changed
                        = clusters_changed.lock().unwrap();
//                        if clusters[step.cluster1].len() >= clusters[step.cluster2].len() {
//                            clusters_changed.push(clusters[step.cluster1].clone());
//                        } else {
//                            clusters_changed.push(clusters[step.cluster2].clone());
//                        }
                    let mut changed = changed.lock().unwrap();
                    *changed += 1;
                    clusters_changed.push(new_cluster);
                } else { // cluster is by itself
                    let mut clusters_changed
                        = clusters_changed.lock().unwrap();
                    clusters_changed.push(clusters[step.cluster1].clone());
                    clusters_changed.push(clusters[step.cluster2].clone());
                }
                // Check individually for leaf clusters
            } else if step.cluster1 <= clusters.len() - 1 {
                let mut clusters_changed
                    = clusters_changed.lock().unwrap();
                clusters_changed.push(clusters[step.cluster1].clone());
            } else if step.cluster2 <= clusters.len() - 1 {
                let mut clusters_changed
                    = clusters_changed.lock().unwrap();
                clusters_changed.push(clusters[step.cluster2].clone());
            }
        });

        let clusters_changed
            = clusters_changed.lock().unwrap().clone();
        let changed = changed.lock().unwrap();
        // If the number of clusters changed, then we rerun linkage clustering
        // First use of recursion properly, nice.
        if *changed > 0 {
            let (clusters_changed, clusters_shared_reads, jaccard_distances)
                = linkage_clustering_of_clusters(&clusters_changed,
                                           variant_info,
                                           variant_map);
            return (clusters_changed, clusters_shared_reads, jaccard_distances.to_vec())
        } else {
            return (clusters.clone(), clusters_shared_reads, jaccard_distances.to_vec())
        }
    } else {
        // create placeholder jaccard hashmap when there is only one cluster to prevent nothing
        // from being returned
        let mut jaccard_map = HashMap::new();
        let mut input = HashMap::new();
        input.insert(1, 1.0);
        jaccard_map.insert(0, input);

        return (clusters.clone(), jaccard_map, vec![1.0])
    }
}

/// Connects variants into initial clusters based on shared read sets
#[allow(unused)]
pub fn linkage_clustering_of_variants(variant_info: &Vec<fuzzy::Var>,
                                      anchor_size: usize,
                                      anchor_similarity: f64,
                                      minimum_reads_in_link: usize)
    -> Vec<fuzzy::Cluster>
{
    debug!("Phasing {} variants...", variant_info.len());
    if variant_info.len() > 1 {
        // Initiate the hashmap linking each variant to the variants it shares reads with
        let links = Mutex::new(HashMap::new());
        // Loop through each permutation of 2 clusters and observe shared variants in reads
        (0..variant_info.len()).into_iter()
            .permutations(2)
            .collect::<Vec<Vec<usize>>>().into_par_iter().for_each(|indices| {


            // Get variants by index
            let var1 = &variant_info[indices[0]];
            let var2 = &variant_info[indices[1]];

            if !(var1.tid == var2.tid && var1.pos == var2.pos) {

                // Read ids of first variant
                let set1 = &var1.reads;

                // Read ids of second variant
                let set2 = &var2.reads;

                // Add the jaccard's similarity to the hashmap for the two clusters
                let intersection: HashSet<_> = set1
                    .intersection(&set2).collect();

                let union: HashSet<_> = set1
                    .union(&set2).collect();

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
                    let links_out = links.entry(indices[0])
                        .or_insert([indices[0]].iter().cloned().collect::<BTreeSet<usize>>());
                    links_out.insert(indices[1]);
                    let links_out = links.entry(indices[1])
                        .or_insert([indices[1]].iter().cloned().collect::<BTreeSet<usize>>());
                    links_out.insert(indices[0]);
                }
            }
        });
        // TODO: Make sure SNPs at same location don't get put together

        let links = links.lock().unwrap();

        // Create condensed links sender and receiver, avoiding use of mutex
        let (condensed_links_s, condensed_links_r) = channel();

        // extend the links for each anchor point by the union of all the indices
        links.par_iter().for_each_with(condensed_links_s, |s, (_, current_links)|{
            let mut anchors = current_links.clone();

            current_links.iter().for_each(|index| {
                match links.get(&index) {
                    Some(other_links) => {
//                        let mut other_links: Vec<_> = other_links.par_iter().cloned().collect();
                        // Get variants by index and check that there is no conflict in position
                        // with any established links in current links
                        let var_check_1 = &variant_info[*index];
                        let (check_s, check_r) = channel();
                        other_links.par_iter().for_each_with(check_s, |s_bool, index_to_check| {
                            let var_check_2 = &variant_info[*index_to_check];
                            s_bool.send(var_check_1.tid == var_check_2.tid &&
                                var_check_1.pos == var_check_2.pos).unwrap()
                        });
                        let mut checked_links: HashSet<_> = check_r.iter().collect();

                        if !checked_links.contains(&true) {
                            anchors.par_extend(other_links.par_iter())
                        }
                    },
                    _ => {},
                }
            });

            // Filter out final links that aren't of size n
            if anchors.len() > anchor_size {
                s.send(anchors).unwrap();
            }
        });

        // Collect receiver into vec and sort by length
        let mut condensed_links: Vec<_> = condensed_links_r.iter().collect();
        condensed_links.par_sort_by_key(|key| key.len());
        debug!("pre filtering condensed links {:?}", &condensed_links);

        // Filter highly similar sets
        condensed_links.dedup_by(|a, b| {
            let intersection: HashSet<_> = a.intersection(&b).collect();
            let union: HashSet<_> = a.union(&b).collect();
            let jaccard = intersection.len() as f64 / union.len() as f64;
            jaccard > anchor_similarity
        });

        debug!("post filtering condensed links {:?}", &condensed_links);
        let (initial_clusters_s, initial_clusters_r) = channel();
        condensed_links.par_iter().for_each_with(initial_clusters_s, |s, link_set| {
            s.send(link_set.par_iter().map(|link| {
                fuzzy::Assignment {
                    index: *link,
                    label: 1.,
                    category: fuzzy::Category::Core,
                }
            }).collect::<fuzzy::Cluster>()).unwrap()
        });
        let initial_clusters: Vec<_> = initial_clusters_r.iter().collect();

        return initial_clusters
    } else {
        Vec::new()
    }
}

/// Get all of the associated read ids for a given cluster
pub fn get_read_set(variants: &fuzzy::Cluster,
                variant_info: &Vec<fuzzy::Var>,
                variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>) -> HashSet<Vec<u8>> {

    let read_set = Mutex::new(HashSet::new());

    variants.par_iter().for_each(|assignment|{
        let variant = &variant_info[assignment.index];
        let base = &variant_map[&variant.tid][&variant.pos][&variant.var].reads;
        base.par_iter().for_each(|id|{
            let mut read_set = read_set.lock().unwrap();
            read_set.insert(id.clone());
        });
    });
    let read_set = read_set.lock().unwrap().clone();
    read_set
}

/// Extract the read ids associated with a particular variant
pub fn get_variant_set(variant: &fuzzy::Var,
                   variant_map: &HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>) -> HashSet<Vec<u8>> {


    let variant_set = variant_map[&variant.tid][&variant.pos][&variant.var].clone().reads;

    return variant_set
}

/// helper function to get the index of condensed matrix from it square form
fn get_condensed_index(i: usize, j: usize, n: usize) -> Option<usize>{
    if i == j {
        return None
    } else {
        return Some(n*i - i*(i+1)/2 + j - 1 - i)
    }
}