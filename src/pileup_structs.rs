use std::collections::HashMap;

pub enum PileupStats {
    PileupContigStats {
        tetfrequency: Vec<HashMap<char, usize>>,
        depth: Vec<usize>,
        tid: usize,
        target_name: Vec<u8>,
    }
}

impl PileupStats {
    pub fn new_contig_stats() -> PileupStats {
        PileupStats::PileupContigStats {
            tetfrequency: vec!(),
            depth: vec!(),
            tid: 0,
            target_name: vec!(),
        }
    }
}

pub trait PileupFunctions {
    fn setup(&mut self);

    fn add_contig(&mut self,
                  tet_freq: Vec<HashMap<char, usize>>,
                  depth: Vec<usize>,
                  tid: usize,
                  target_name: Vec<u8>);

    fn calc_variants(&mut self);
}

impl PileupFunctions for PileupStats {
    fn setup(&mut self) {
        match self {
            PileupStats::PileupContigStats {
                ref mut tetfrequency,
                ref mut depth,
                ref mut tid,
                ref mut target_name,
            } => {
                *tetfrequency = vec!();
                *depth = vec!();
                *tid = 0;
                *target_name = vec!();
            }
        }
    }

    fn add_contig(&mut self, tet_freq: Vec<HashMap<char, usize>>,
                  read_depth: Vec<usize>,
                  target_id: usize,
                  name: Vec<u8>) {
        match self {
            PileupStats::PileupContigStats {
                ref mut tetfrequency,
                ref mut depth,
                ref mut tid,
                ref mut target_name,
            } => {
                *tetfrequency = tet_freq;
                *depth = read_depth;
                *tid = target_id;
                *target_name = name;
            }
        }
    }
    fn calc_variants(&mut self){
        match self {
            PileupStats::PileupContigStats {ref mut tetfrequency,
                ref mut depth,
                ref mut tid,
                ref mut target_name,
            } => {

            }
        }
    }
}