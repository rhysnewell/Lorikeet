use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use kodama::Dendrogram;

#[derive(Debug,Clone)]
pub struct Haplotype {
    pub root_cluster_id: usize,
    pub variant_indices: HashSet<usize>,
    pub variants: HashMap<i32, BTreeMap<String, (i32, usize)>>,
    pub node_size: usize,
    pub haplotype_index: usize,
}

impl Haplotype {
    pub fn new() -> Haplotype {
        Haplotype {
            root_cluster_id: 0,
            variant_indices: HashSet::new(),
            variants: HashMap::new(),
            node_size: 0,
            haplotype_index: 0,
        }
    }

    pub fn start(node_size: usize,
                 node_id: usize,
                 index: usize) -> Haplotype {
        Haplotype {
            root_cluster_id: node_id,
            variant_indices: HashSet::new(),
            variants: HashMap::new(),
            node_size: node_size,
            haplotype_index: index,
        }
    }

    pub fn add_variants(&mut self,
                        dendrogram: &Dendrogram<f64>,
                        dendro_clusters: &HashMap<usize, BTreeMap<i32, (String, i32)>>){
        let n = dendrogram.len() + 1;
        let step_1 = &dendrogram[self.root_cluster_id];
        let mut to_search = HashSet::new();
        to_search.insert(step_1.cluster1);
        to_search.insert(step_1.cluster2);
        let mut current_step = step_1;
        while to_search.len() > 0 {
            // create new mutable hashset to insert to whilst we drain current hashset
            let mut new_search = HashSet::new();
            for cluster_label in to_search.drain() {
                // convert cluster label into index
                let step_index = cluster_label - n;
                let cur_step = &dendrogram[step_index];

                // Check to see if each cluster id corresponds to a variant index
                // If it does, then extract the variant
                // else, add cluster id to search list
                if cur_step.cluster1 < n {
                    self.variant_indices.insert(cur_step.cluster1);
                    let variant_pos =
                        dendro_clusters.get(&step_index).expect("Step ID not found");
                    for (pos, variant) in variant_pos {
                        let captured_var = self.variants
                            .entry(*pos).or_insert(BTreeMap::new());
                        captured_var.entry(variant.0.clone())
                            .or_insert((variant.1, self.haplotype_index));
                    }
                } else {
                    new_search.insert(cur_step.cluster1);
                }

                if cur_step.cluster2 < n {
                    self.variant_indices.insert(cur_step.cluster2);
                    let variant_pos =
                        dendro_clusters.get(&step_index).expect("Step ID not found");
                    for (pos, variant) in variant_pos {
                        let captured_var = self.variants
                            .entry(*pos).or_insert(BTreeMap::new());
                        captured_var.entry(variant.0.clone())
                            .or_insert((variant.1, self.haplotype_index));

                    }
                } else {
                    new_search.insert(cur_step.cluster2);
                }
            }
            to_search = new_search;
        }

    }
}


#[derive(Debug, Clone)]
pub struct Genotype {
    pub read_ids: HashSet<i32>,
    pub start_var_pos: usize,
    pub ordered_variants: HashMap<i32, String>,
    pub frequencies: Vec<f64>,
}

impl Genotype {
    pub fn start(position: usize) -> Genotype {
        Genotype {
            read_ids: HashSet::new(),
            start_var_pos: position,
            ordered_variants: HashMap::new(),
            frequencies: Vec::new(),
        }
    }

    pub fn new(&mut self, read_id: i32, position_map: &BTreeMap<i32, String>,
           variant_abundances: &HashMap<i32, BTreeMap<String, f64>>) {
        self.read_ids = HashSet::new();
        self.frequencies = Vec::new();
        self.ordered_variants = HashMap::new();
        self.read_ids.insert(read_id);
        for (pos, variant) in position_map.iter() {
            let variant_map = &variant_abundances.get(pos);
            match variant_map {
                Some(variant_map) => {
                    debug!("Fetching position: {:?} {:?} {:?}", pos, variant_map, variant);
//                self.frequencies.push(variant_map[variant].clone());
                    self.ordered_variants
                        .insert(pos.clone(), variant.to_string());
                },
                None => continue,
            }
        }
    }

    pub fn check_variants(&self, position_map: &BTreeMap<i32, String>, intersection: Vec<i32>) -> bool {
        // check variants against stored variants for a genotype along the shared positions
        let mut new_var= false;
        for check_pos in intersection.iter() {
            if self.ordered_variants.contains_key(&check_pos) {
                let current_var = match self
                    .ordered_variants
                    .get(&check_pos) {
                    Some(var) => var,
                    None => {
                        println!("Position not recorded in variant map");
                        std::process::exit(1)
                    }
                };
                let check_var = position_map.get(check_pos)
                    .expect("No variant at that position while checking genotypes");

                if current_var != check_var {
                    //Then this is a new genotype
                    new_var = true;
                    break
                } else {
                    new_var = false;
                }
            } else {
                // current variant does not contain this position, must be new genotype
                new_var = true;
                break
            }
        }
        return new_var
    }
}
