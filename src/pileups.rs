use std;
use std::collections::HashMap;
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
use rm::linalg::Matrix;
use rm::linalg::Vector;
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
            let mut nuc_freq = Vec::new();
            let mut depth = Vec::new();
            let mut last_tid: i32 = -2; // no such tid in a real BAM file
            let mut total_indels_in_current_contig = 0;
            let mut base;
//            let mut pileup_struct = PileupStats::new_contig_stats(min,
//                                                     max,
//                                                     min_fraction_covered_bases,
//                                                     contig_end_exclusion);

            let mut process_previous_contigs = |last_tid: i32,
                                                depth,
                                                nuc_freq,
                                                total_indels_in_current_contig| {
                if last_tid != -2 {
//                    println!("tid {} and {} indels", last_tid,
//                           total_indels_in_current_contig);

//                    println!("{} depth {:?}", std::str::from_utf8(target_names[last_tid as usize]).unwrap(), depth);
                    let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                    let contig_name = target_names[last_tid as usize].to_vec();
                    reference.fetch_all(std::str::from_utf8(target_names[last_tid as usize]).unwrap());
                    let mut seq: Vec<u8> = Vec::new();
                    reference.read(&mut seq);
                    let kmers = hash_kmers(&seq, 4);
                    let mut tet_freq = BTreeMap::new();
                    for (tet, loc) in kmers.iter(){
                        tet_freq.entry(tet.to_vec()).or_insert(loc.len());
                    }
//                    println!("{:?}", tet_freq);


                    let mut pileup_struct = PileupStats::new_contig_stats(min,
                                                                          max,
                                                                          min_fraction_covered_bases,
                                                                          contig_end_exclusion);
                    pileup_struct.add_contig(nuc_freq,
                                       tet_freq,
                                       depth,
                                       last_tid,
                                       total_indels_in_current_contig,
                                       contig_name,
                                       contig_len);

                    let coverage = pileup_struct.calc_coverage();

                    pileup_struct.calc_variants(depth_threshold,
                                                var_fraction);
                    pileup_matrix.add_contig(pileup_struct,
                                             target_names.len() as usize)
                }
            };

            for p in bam_pileups {
                // if reference has changed, print the last record

                let pileup = p.unwrap();
                let tid = pileup.tid() as i32;
                // if reference has changed, print the last record
                if tid != last_tid {
                    process_previous_contigs(
                        last_tid,
                        depth,
                        nuc_freq,
                        total_indels_in_current_contig);
                    nuc_freq = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    total_indels_in_current_contig = 0;
                } else {
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                }

                for alignment in pileup.alignments() {

                    if !alignment.is_del() && !alignment.is_refskip() {

                        base = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                        let count = nuc_freq[pileup.pos() as usize].entry(base)
                            .or_insert(vec!());
                        count.push(alignment.record().qname().to_vec());

                    } else {
                        let count = nuc_freq[pileup.pos() as usize]
                            .entry('N' as char).or_insert(vec!());
                        count.push(alignment.record().qname().to_vec());
                    }
                    // mark indel start
                    match alignment.indel() {
                        bam::pileup::Indel::Ins(len) | bam::pileup::Indel::Del(len) => {
                            total_indels_in_current_contig += 1;
                        },
                        bam::pileup::Indel::None => ()
                    }
                }
            }
//        print!("{:?}", nuc_freq);
            process_previous_contigs(
                last_tid,
                depth,
                nuc_freq,
                total_indels_in_current_contig);
        }
        pileup_matrix.print_matrix();
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
