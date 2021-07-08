// /// This code mostly belongs to https://github.com/10XGenomics/longranger/blob/master/lib/pvc/src/asm_caller.rs
// /// The code is being used as it was the example code provided on the rust-debruijn docs page:
// /// https://docs.rs/debruijn/0.3.2/debruijn/
// ///
//
// use std::path::PathBuf;
// use std::fs::File;
// use std::collections::HashSet;
// use std::ops::Index;
// use std::fmt;
//
// use debruijn::{Dir, Exts};
// use debruijn::{self, Kmer};
// use debruijn::kmer;
// use debruijn::paths::{DebruijnGraph, Node};
// use debruijn::compression::{PathCompression, SimpleCompress, Simplify};
// use debruijn::filter::{self, KmerSummarizer};
// use debruijn::dna_string::DnaString;
// use debruijn::clean_graph::CleanGraph;
//
// use bio::io::fasta;
// use rust_htslib::bam::{self, Read, IndexedReader};
//
// use call;
//
// #[derive(Copy, Clone, Hash, PartialOrd, Ord, PartialEq, Eq)]
// pub struct HapCode {
//     data: u32,
// }
//
// impl HapCode {
//     fn new(hap: u8, read_num: usize) -> HapCode {
//         HapCode { data: (hap as u32) << 16 | read_num as u32 }
//     }
//
//     fn hap(&self) -> u8 {
//         (self.data >> 16) as u8
//     }
//
//     fn read(&self) -> u32 {
//         (self.data & 0xFFFF)
//     }
// }
//
//
// type Kmer1 = kmer::Kmer40;
//
// #[derive(Default)]
// struct PhasedReads {
//     // Vector containing HashSet<u32> read ids
//     // index 0 - Reads assigned to H1
//     // index 1 - Reads assigned to H2
//     // Index 2 - Unassigned reads
//     read_counts: [HashSet<u32>; 3],
// }
//
// impl fmt::Debug for PhasedReads {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f,
//                "PhasedReads {{ h1: {}, h2: {}, un: {} }}",
//                self.read_counts[0].len(),
//                self.read_counts[1].len(),
//                self.read_counts[2].len())
//     }
// }
//
// impl Clone for PhasedReads {
//     fn clone(&self) -> Self {
//         PhasedReads {
//             read_counts: [self.read_counts[0].clone(),
//                 self.read_counts[1].clone(),
//                 self.read_counts[2].clone()],
//         }
//     }
// }
//
//
// impl PhasedReads {
//     pub fn new() -> PhasedReads {
//         // Initialize PhasedReads struct
//         PhasedReads { read_counts: [HashSet::new(), HashSet::new(), HashSet::new()] }
//     }
//
//     pub fn add(&mut self, code: HapCode) {
//         // Insert read to associate HapCode
//         self.read_counts
//             .get_mut(code.hap() as usize)
//             .unwrap()
//             .insert(code.read());
//     }
//
//     pub fn reduce(mut self, other: &PhasedReads) -> PhasedReads {
//         // Combine two lists of phased reads
//         for i in 0..3 {
//             self.read_counts
//                 .get_mut(i)
//                 .unwrap()
//                 .extend(&other.read_counts[i])
//         }
//
//         self
//     }
//
//     pub fn from_hap_codes<I: Iterator<Item = HapCode>>(iter: I) -> Self {
//         // When given an iterator containing hapcodes
//         // Create a new set of phased reads from these hapcode
//         let mut pc = PhasedReads::new();
//         for code in iter {
//             pc.add(code);
//         }
//
//         pc
//     }
//
//     pub fn coverage(&self) -> usize {
//         let mut n = 0;
//         for i in 0..3 {
//             n += self.read_counts.get(i).unwrap().len();
//         }
//
//         n
//     }
// }
//
// impl Index<usize> for PhasedReads {
//     type Output = HashSet<u32>;
//     fn index(&self, index: usize) -> &HashSet<u32> {
//         self.read_counts.get(index).unwrap()
//     }
// }
//
// pub struct HapCountSummarize {
//     min_kmer_obs: usize,
// }
//
// impl KmerSummarizer<HapCode, PhasedReads> for HapCountSummarize {
//     ///The input items is an iterator over kmer observations.
//     /// Input observation is a tuple of (kmer, extensions, data).
//     /// The summarize function inspects the data and returns a tuple indicating:
//     ///     - whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation)
//     ///     - the accumulated Exts of the kmer
//     ///     - A {{PhasedReads}} data object.
//     ///
//     fn summarize<K, F: Iterator<Item = (K, Exts, HapCode)>>(&self,
//                                                             items: F)
//                                                             -> (bool, Exts, PhasedReads) {
//
//         // Store single-base extensions for a DNA Debruijn graph.
//         // 8 bits, 4 higher order ones represent extensions to the right,
//         // 4 lower order ones represent extensions to the left. For each direction the bits
//         // (from lower order to higher order) represent whether there exists an extension with
//         // each of the letters A, C, G, T. So overall the bits are: right left T G C A T G C A
//         let mut all_exts = Exts::empty();
//
//         let mut nreads = 0;
//         let mut phased_reads = PhasedReads::new();
//
//         for (_, exts, d) in items {
//             nreads += 1;
//             phased_reads.add(d);
//             all_exts = all_exts.add(exts);
//         }
//
//         (nreads >= self.min_kmer_obs, all_exts, phased_reads)
//     }
// }
//
//
// /// return Vector of sequences (seq, qual, bc, read)
// fn get_sequences(bam: &mut IndexedReader,
//                  locus: &Locus,
//                  k: usize)
//                  -> Vec<(DnaString, Exts, HapCode)> {
//     // TODO: Change this work with Base rather than Locus
//     let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
//     let mut num_reads = 0;
//
//     let mut seqs = Vec::new();
//
//     // Add phased reads to haplotype buckets
//     bam.fetch(tid, locus.start, locus.end)
//         .expect("Error fetching BAM file.");
//     for _rec in bam.records() {
//         let rec = _rec.ok().unwrap();
//
//         // Skip poor mappings & secondary alignments
//         if rec.mapq() < 10 || rec.is_secondary() {
//             continue;
//         }
//
//         let hap = call::get_haplotype(&rec).unwrap_or(3) - 1;
//         let bc = HapCode::new(hap, num_reads);
//
//         let raw_seq = call::get_seq_bounded(&rec, locus);
//
//         if raw_seq.len() < k {
//             continue;
//         }
//
//         let fwd_seq = if rec.is_reverse() {
//             event::rc_seq(&raw_seq)
//         } else {
//             raw_seq
//         };
//
//         let seq: Vec<u8> = fwd_seq.into_iter().map(debruijn::base_to_bits).collect();
//         let dna_string = DnaString::from_bytes(&seq);
//         seqs.push((dna_string, Exts::empty(), bc));
//
//         num_reads += 1;
//         if num_reads > 10_000 {
//             break;
//         }
//     }
//
//     seqs
// }
//
// #[inline(never)]
// pub fn find_hap_seqs(seqs: &Vec<(DnaString, Exts, HapCode)>,
//                      beam_search: bool,
//                      args: &Args)
//                      -> ((Vec<u8>, usize), (Vec<u8>, usize)) {
//
//     let summarizer = HapCountSummarize { min_kmer_obs: 2 };
//     let (mut valid_kmers, _) = filter::filter_kmers::<Kmer1, _, _, _, _>(seqs, summarizer, false);
//     filter::fix_exts_local(false, &mut valid_kmers);
//
//     debug!("Got valid kmers: {}", valid_kmers.len());
//
//     if valid_kmers.len() < 3 {
//         return ((Vec::new(), 0), (Vec::new(), 0));
//     }
//
//     trace!("Building graph");
//     let cmp = SimpleCompress::new(|a: PhasedReads, b: &PhasedReads| a.reduce(b));
//     let path_comp: PathCompression<Kmer1, _, _> = PathCompression::new(false, cmp);
//     let _graph = path_comp.build_nodes(&valid_kmers).finish();
//
//
//     // Now try to clean the tips.
//     let cleaner = CleanGraph::new(|node: &Node<Kmer1, PhasedReads>| {
//         node.len() < Kmer1::k() + Kmer1::k() / 3 &&
//             node.data().coverage() < 5
//     });
//     let nodes_to_censor = cleaner.find_bad_nodes(&_graph);
//
//     let cmp = SimpleCompress::new(|a: PhasedReads, b: &PhasedReads| a.reduce(b));
//     let graph = Simplify::simplify(_graph, Some(nodes_to_censor), false, cmp);
//
//     if args.flag_debug {
//         graph.print_with_data();
//     }
//
//     let (path1, path2) = if beam_search {
//
//         // Beam search of graph, with reward heuristics
//         let path1 = graph.max_path_beam(20, |e| e[0].len() as f32 - 2.0, |e| {
//             e[0].len() >= 4 ||
//                 (e[0].len() as f32 / ((e[0].len() + e[1].len() + e[2].len()) as f32) > 0.3)
//         });
//
//         let path2 = graph.max_path_beam(20, |e| e[1].len() as f32 - 2.0, |e| {
//             e[1].len() >= 4 ||
//                 (e[1].len() as f32 / ((e[0].len() + e[1].len() + e[2].len()) as f32) > 0.3)
//         });
//
//         (path1, path2)
//     } else {
//
//         // Purely greedy with stopping conditions when there's no clear winner.
//         let path1 = graph.max_path(|e| e[0].len() as f32 - 2.0, |e| {
//             e[0].len() >= 3 &&
//                 (e[0].len() as f32 / ((e[0].len() + e[1].len() + e[2].len()) as f32) > 0.3)
//         });
//
//         let path2 = graph.max_path(|e| e[1].len() as f32 - 2.0, |e| {
//             e[1].len() >= 3 &&
//                 (e[1].len() as f32 / ((e[0].len() + e[1].len() + e[2].len()) as f32) > 0.3)
//         });
//
//         (path1, path2)
//     };
//
//     let seq1 = graph.sequence_of_path(path1.iter());
//     let seq2 = graph.sequence_of_path(path2.iter());
//
//
//     if args.flag_debug {
//         graph.to_dot(PathBuf::from("db.dot"),
//                      &|d| format!("[{}, {}, {}]", d[0].len(), d[1].len(), d[2].len()));
//         graph.to_gfa_with_tags(PathBuf::from("db.gfa"),
//                                |node| format!("RC:i:{}", node.data().coverage()));
//     }
//
//     // Count reads supporting each path -- will constitute a crude QV
//     let score1 = path_hap_summary(&path1, &graph).read_counts[0].len();
//     let score2 = path_hap_summary(&path2, &graph).read_counts[1].len();
//
//     ((seq1.to_ascii_vec(), score1), (seq2.to_ascii_vec(), score2))
// }
//
// fn path_hap_summary(path: &Vec<(usize, Dir)>,
//                     graph: &DebruijnGraph<Kmer1, PhasedReads>)
//                     -> PhasedReads {
//     let mut phased_read_counts: PhasedReads = PhasedReads::new();
//
//     for &(node_id, _) in path.iter() {
//         let node = graph.get_node(node_id);
//         phased_read_counts = phased_read_counts.reduce(node.data());
//     }
//
//     phased_read_counts
// }