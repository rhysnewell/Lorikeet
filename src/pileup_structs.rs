use std::collections::HashMap;
use std::collections::BTreeMap;

pub enum PileupStats {
    PileupContigStats {
        tetfrequency: Vec<HashMap<char, usize>>,
        depth: Vec<usize>,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: usize,
        contig_end_exclusion: u32,
        min_fraction_covered_bases: f32,
        min: f32,
        max: f32,
    }
}

impl PileupStats {
    pub fn new_contig_stats(min: f32, max: f32, min_fraction_covered_bases: f32,
                            contig_end_exclusion: u32) -> PileupStats {
        PileupStats::PileupContigStats {
            tetfrequency: vec!(),
            depth: vec!(),
            tid: 0,
            total_indels: 0,
            target_name: vec!(),
            target_len: 0,
            contig_end_exclusion: contig_end_exclusion,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max,
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

    fn calc_coverage(&mut self) -> f32;
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
                ..
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

    fn calc_coverage(&mut self) -> f32 {
        let total_bases = *observed_contig_length + unobserved_contig_length;
        debug!("Calculating coverage with num_covered_bases {}, observed_length {}, unobserved_length {} and counts {:?}",
               num_covered_bases, observed_contig_length, unobserved_contig_length, counts);
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
                    0.0
                } else {
                    let min_index: usize = (*min * total_bases as f32).floor() as usize;
                    let max_index: usize = (*max * total_bases as f32).ceil() as usize;
                    if *num_covered_bases == 0 {return 0.0;}
                    counts[0] += unobserved_contig_length;

                    let mut num_accounted_for: usize = 0;
                    let mut total: usize = 0;
                    let mut started = false;
                    let mut i = 0;
                    for num_covered in counts.iter() {
                        num_accounted_for += *num_covered as usize;
                        debug!("start: i {}, num_accounted_for {}, total {}, min {}, max {}", i, num_accounted_for, total, min_index, max_index);
                        if num_accounted_for >= min_index {
                            debug!("inside");
                            if started {
                                if num_accounted_for > max_index {
                                    debug!("num_accounted_for {}, *num_covered {}",
                                           num_accounted_for, *num_covered);
                                    let num_excess = num_accounted_for - *num_covered as usize;
                                    let num_wanted = match max_index >= num_excess {
                                        true => max_index - num_excess + 1,
                                        false => 0
                                    };
                                    debug!("num wanted1: {}", num_wanted);
                                    total += num_wanted * i;
                                    break;
                                } else {
                                    total += *num_covered as usize * i;
                                }
                            } else {
                                if num_accounted_for > max_index {
                                    // all coverages are the same in the trimmed set
                                    total = (max_index-min_index+1) * i;
                                    started = true
                                } else if num_accounted_for < min_index {
                                    debug!("too few on first")
                                } else {
                                    let num_wanted = num_accounted_for - min_index + 1;
                                    debug!("num wanted2: {}", num_wanted);
                                    total = num_wanted * i;
                                    started = true;
                                }
                            }
                        }
                        debug!("end i {}, num_accounted_for {}, total {}", i, num_accounted_for, total);

                        i += 1;
                    }
                    total as f32 / (max_index-min_index) as f32
                }
            }
        };
        return answer
    }
}