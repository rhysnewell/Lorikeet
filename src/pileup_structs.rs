use std::collections::HashMap;
use std::collections::BTreeMap;

pub enum PileupStats {
    PileupContigStats {
        tetfrequency: Vec<HashMap<char, usize>>,
        depth: Vec<usize>,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: usize
    }
}

impl PileupStats {
    pub fn new_contig_stats() -> PileupStats {
        PileupStats::PileupContigStats {
            tetfrequency: vec!(),
            depth: vec!(),
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0,
        }
    }
}

pub trait PileupFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  tet_freq: Vec<HashMap<char, usize>>,
                  depth: Vec<usize>,
                  tid: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize);

    fn calc_variants(&mut self,
                     depth_thresh: usize,
                     variant_fraction: f64);
}

impl PileupFunctions for PileupStats {
    fn setup(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut tetfrequency,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
            } => {
                *tetfrequency = vec!();
                *depth = vec!();
                *tid = 0;
                *total_indels = 0;
                *target_name = vec!();
                *target_len = 0;
            }
        }
    }

    fn add_contig(&mut self, tet_freq: Vec<HashMap<char, usize>>,
                  read_depth: Vec<usize>,
                  target_id: i32,
                  total_indels_in_contig: usize,
                  contig_name: Vec<u8>,
                  contig_len: usize) {
        match self {
            PileupStats::PileupContigStats {
                ref mut tetfrequency,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
            } => {
                *tetfrequency = tet_freq;
                *depth = read_depth;
                *tid = target_id;
                *total_indels = total_indels_in_contig;
                *target_name = contig_name;
                *target_len = contig_len;
            }
        }
    }

    fn calc_variants(&mut self, depth_thresh: usize, variant_fraction: f64){
        match self {
            PileupStats::PileupContigStats {
                ref mut tetfrequency,
                ref mut depth,
                tid,
                total_indels,
                target_name,
                target_len
            } => {
                let mut variants = BTreeMap::new();
                let mut variant_count = 0;
                let mut cursor = 0;
                let mut depth_sum = 0;
                for zipped in tetfrequency.iter().zip(depth.iter()){
                    let (tetfreq, d) = zipped;
                    let mut rel_abundance = HashMap::new();
                    if d >= &depth_thresh {
                        if tetfreq.len() > 1 {
                            variant_count += 1;
                            for (base, count) in tetfreq.iter() {
                                if *count as f64 / *d as f64 >= variant_fraction {
                                    rel_abundance.insert(base, count / d);
                                }
                            }
                        }
                    }

                    if rel_abundance.len() > 1 {
                        variants.insert(cursor, rel_abundance);
                    }

                    cursor += 1;
                    depth_sum += d;
                }

                let variant_per_base = variant_count as f32/target_len.clone() as f32;
                let contig_coverage = depth_sum as f32/target_len.clone() as f32;
                println!("{}\t{}", depth_sum, target_len);
                println!("{}\t{}\t{:?}\t{}",
                         std::str::from_utf8(&target_name[..]).unwrap(),
                         variant_per_base,
                         contig_coverage,
                         total_indels);
            }
        }
    }
}