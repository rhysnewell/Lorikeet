use std;
use std::collections::{HashMap, HashSet};
use std::collections::BTreeMap;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use pileup_structs::*;
use mosdepth_genome_coverage_estimators::*;
use bam_generator::*;
use coverage_takers::*;
use FlagFilter;
use ReadsMapped;
use std::str;
//use rm::linalg::Matrix;
//use rm::linalg::Vector;
use std::fs::File;
use bio::io::fasta::*;
use bio::alignment::sparse::*;
use std::io::prelude::*;



pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter,
    depth_threshold: usize,
    var_fraction: f64,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32) {

    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let mut pileup_matrix = PileupMatrix::new_contig_stats();
        {
            let header = bam_generated.header().clone();
            let target_names = header.target_names();
            let bam_pileups = bam_generated.pileups();
            let mut ref_seq: Vec<u8> = Vec::new();
            let mut nuc_freq: Vec<HashMap<char, HashSet<i32>>> = Vec::new();
            let mut indels = Vec::new();
            let mut tet_freq = BTreeMap::new();
            let mut depth = Vec::new();
            let mut last_tid: i32 = -2; // no such tid in a real BAM file
            let mut total_indels_in_current_contig = 0;
            let mut read_cnt_id = 0;
            let mut read_to_id = HashMap::new();
            let mut previous_read_positions = HashMap::new();
//            let mut indel_start_sites: HashMap<i32, Vec<Indel>> = HashMap::new();
            let mut base;

            let mut process_previous_contigs = |last_tid: i32,
                                                depth: Vec<usize>,
                                                nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
                                                tet_freq,
                                                indels,
                                                total_indels_in_current_contig,
                                                ref_sequence: Vec<u8>| {
                if last_tid != -2 {

                    let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                    let contig_name = target_names[last_tid as usize].to_vec();

                    let mut pileup_struct = PileupStats::new_contig_stats(min,
                                                                          max,
                                                                          min_fraction_covered_bases,
                                                                          contig_end_exclusion);
                    debug!("INDELS: {:?}", indels);
                    debug!("nuc frequency: {:?}", nuc_freq);




                    pileup_struct.add_contig(nuc_freq,
                                               tet_freq,
                                               depth,
                                               indels,
                                               last_tid,
                                               total_indels_in_current_contig,
                                               contig_name,
                                               contig_len);

                    let coverage = pileup_struct.calc_coverage();

                    pileup_struct.calc_variants(depth_threshold,
                                                var_fraction);
                    pileup_struct.print_variants(ref_sequence);
//                    pileup_struct.generate_variant_contig(ref_sequence,
//                                                          depth_threshold.clone(),
//                                                          var_fraction.clone());

                    pileup_matrix.add_contig(pileup_struct,
                                             target_names.len() as usize);

                }
            };

            for p in bam_pileups {

                let pileup = p.unwrap();
                let tid = pileup.tid() as i32;
                // if reference has changed, print the last record
                if tid != last_tid {


                    process_previous_contigs(
                        last_tid,
                        depth,
                        nuc_freq,
                        tet_freq,
                        indels,
                        total_indels_in_current_contig,
                        ref_seq);
                    nuc_freq = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    indels = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
//                    indel_start_sites = HashMap::new();
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    total_indels_in_current_contig = 0;
                    reference.fetch_all(std::str::from_utf8(target_names[tid as usize]).unwrap());
                    ref_seq = Vec::new();
                    reference.read(&mut ref_seq);
                    let kmers = hash_kmers(&ref_seq, 4);
                    tet_freq = BTreeMap::new();
                    for (tet, loc) in kmers.iter(){
                        tet_freq.entry(tet.to_vec()).or_insert(loc.len());
                    }
                } else {
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                }

                for alignment in pileup.alignments() {
                    // Check if new read to id

                    if !read_to_id.contains_key(&alignment.record().qname().to_vec()) {
                        read_to_id.entry(alignment.record().qname().to_vec())
                            .or_insert(read_cnt_id);
                        read_cnt_id += 1;
                    }

                    let qpos = match alignment.qpos() {
                        Some(position) => {
                            previous_read_positions.insert(read_to_id[&alignment.record().qname().to_vec()], position);
                            position
                        },
                        None => {
                            let position = previous_read_positions.entry(
                                read_to_id[&alignment.record().qname().to_vec()]).or_insert(0);
                            *position
                        },
                    };
                    if !alignment.is_del() && !alignment.is_refskip() {
                        if alignment.record().seq().len() == 0 {
                            debug!("Zero length read: {:?}", alignment.record().seq());
                            break
                        }
                        base = alignment.record().seq()[qpos] as char;
                        if base != ref_seq[pileup.pos() as usize] as char {
                            let id = nuc_freq[pileup.pos() as usize].entry(base)
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                        }
                    }
                    // mark indel start
                    match alignment.indel() {

                        bam::pileup::Indel::Ins(len) => {

                            let insert = match str::from_utf8(&alignment.record().seq().as_bytes()[
                                qpos + 1..qpos + 1 + len as usize]) {
                                Ok(ins) => {ins.to_string()},
                                Err(e) => {"".to_string()},
                            };
//                            debug!("Ins len: {} cigar: {:?} id: {}", len, alignment.record().cigar(), read_to_id[&alignment.record().qname().to_vec()]);
//                            debug!("Indel: {} at {}", insert, pileup.pos());
                            let id = nuc_freq[pileup.pos() as usize].entry('I')
                                                                    .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            let id = indels[pileup.pos() as usize].entry(insert)
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            total_indels_in_current_contig += 1;
                        },

                        bam::pileup::Indel::Del(len) => {
//                            debug!("Del len: {} cigar: {:?} id: {}",
//                                   len, alignment.record().cigar(), read_to_id[&alignment.record().qname().to_vec()]);
//                            debug!("Indel: {} at {}",
//                                   std::iter::repeat("N").take(len as usize).collect::<String>(),
//                                   pileup.pos());
                            let id = nuc_freq[pileup.pos() as usize].entry('D')
                                                                        .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            let id = indels[pileup.pos() as usize].entry(
                                std::iter::repeat("N").take(len as usize).collect::<String>())
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            total_indels_in_current_contig += 1;
                        },

                        bam::pileup::Indel::None => ()
                    }
                }
            }

            process_previous_contigs(
                last_tid,
                depth,
                nuc_freq,
                tet_freq,
                indels,
                total_indels_in_current_contig,
                ref_seq);
        }
//        pileup_matrix.print_matrix();
        bam_generated.finish();
    }
}


fn print_previous_zero_coverage_contigs<T: CoverageTaker>(
    last_tid: i32,
    current_tid: i32,
    coverage_estimators: &Vec<CoverageEstimator>,
    target_names: &Vec<&[u8]>,
    coverage_taker: &mut T,
    header: &bam::HeaderView) {
    let mut my_tid = last_tid + 1;
    while my_tid < current_tid {
        debug!("printing zero coverage for tid {}", my_tid);
        coverage_taker.start_entry(
            my_tid as usize,
            std::str::from_utf8(target_names[my_tid as usize]).unwrap());
        for ref coverage_estimator in coverage_estimators.iter() {
            coverage_estimator.print_zero_coverage(
                coverage_taker, header.target_len(my_tid as u32).unwrap());
        }
        coverage_taker.finish_entry();
        my_tid += 1;
    };
}
