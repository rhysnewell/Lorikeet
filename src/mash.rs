//use bio::io::fasta::*;
//use bio::alignment::sparse::*;
use minhashes::*;
use std::collections::HashMap;
use needletail::fastx::{fastx_bytes};
use distance::*;

#[allow(non_snake_case)]
#[derive(Debug)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    pub mashDistance: f64,
    pub commonHashes: u64,
    pub totalHashes: u64,
//    pub query: String,
//    pub reference: String,
}

#[allow(non_snake_case)]
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Sketch {
    pub name: String,
    pub seqLength: Option<u64>,
    pub numValidKmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    pub hashes: Vec<KmerCount>,
}

#[allow(non_snake_case)]
#[derive(Debug)]
pub struct MultiSketch {
    pub kmer: u8,
    pub alphabet: String,
    pub preserveCase: bool,
    pub canonical: bool,
    pub sketchSize: u32,
    pub hashType: String,
    pub hashBits: u16,
    pub hashSeed: u64,
    pub sketches: Vec<Sketch>,
}

impl Sketch {
    pub fn new(
        name: &str,
        length: u64,
        n_kmers: u64,
        kmercounts: Vec<KmerCount>,
        filters: &HashMap<String, String>,
    ) -> Self {
        Sketch {
            name: String::from(name),
            seqLength: Some(length),
            numValidKmers: Some(n_kmers),
            comment: Some(String::from("")),
            filters: Some(filters.clone()),
            hashes: kmercounts,
        }
    }

    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

//    pub fn apply_filtering(&mut self, filters: &FilterParams) -> bool {
//        let (filtered_hashes, filter_stats) = filter_sketch(&self.hashes, &filters);
//        self.hashes = filtered_hashes;
//        self.filters = Some(filter_stats);
//        true
//    }
}

pub fn mash_sequence(
    seq: &[u8],
    final_size: usize,
    kmer_length: u8,
    seed: u64, ) -> Vec<KmerCount> { // need to return Vec<KmerCount> -- Found in minhash
//    let mut sketches = Vec::with_capacity(filenames.len()); // Move outside function
//    let mut needle_seq = canonical(seq);
    let mut seq_len = 0u64;
    let mut minhash =  MinHashKmers::new(final_size, seed);
//    debug!("seq: {:?}", seq);
    match fastx_bytes(
        seq,
        |needle_seq| {
            seq_len += needle_seq.seq.len() as u64;
            debug!("seq len: {:?}", seq_len);

            for (_, kmer, is_rev_complement) in needle_seq.normalize(false).kmers(kmer_length, true) {
                debug!("kmer: {:?}", kmer);
                let rc_count = if is_rev_complement { 1u8 } else { 0u8 };
                minhash.push(kmer, rc_count);
                    }
            },
    ){
        Ok(fastx) => fastx,
        Err(e) => {
            println!("Cannot fastx bytes {:?}", e);
            std::process::exit(1)}
    };
    let _n_kmers = minhash.total_kmers() as u64;
    debug!("{:?}", minhash);
    let hashes = minhash.into_vec();

    // directory should be clipped from filename
//    let basename = path
//        .file_name()
//        .ok_or_else(|| format_err!("Couldn't get filename from path"))?;
//    let sketch = Sketch::new(
//        basename
//            .to_str()
//            .ok_or_else(|| format_err!("Couldn't make filename into string"))?,
//        seq_len,
//        n_kmers,
//        filtered_hashes,
//        &filter_stats,
//    );
//    sketches.push(sketch);
    return hashes
}

#[allow(dead_code)]
fn calc_sketch_distances(
    query_sketches: &[&Sketch],
    ref_sketches: &[Sketch],
    mash_mode: bool,
    max_distance: f64,
) -> Vec<SketchDistance> {
    let mut distances = Vec::new();
    for ref_sketch in ref_sketches.iter() {
        let rsketch = &ref_sketch.hashes;
        for query_sketch in query_sketches.iter() {
            if query_sketch == &ref_sketch {
                continue;
            }
            let qsketch = &query_sketch.hashes;
            let distance = distance(
                &qsketch,
                &rsketch,
                mash_mode,
            )
                .unwrap();
            if distance.mashDistance <= max_distance {
                distances.push(distance);
            }
        }
    }
    distances
}

