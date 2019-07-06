use std;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use pileup_structs::*;

use mosdepth_genome_coverage_estimators::*;
use pileup_structs::*;
use bam_generator::*;
use coverage_takers::*;
use FlagFilter;
use ReadsMapped;
use std::str;


pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter) {

    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        {
            let header = bam_generated.header().clone();
            let target_names = header.target_names();
            let bam_pileups = bam_generated.pileups();
            let mut tet_freq = Vec::new();
            let mut depth = Vec::new();
            let mut last_tid: i32 = -2; // no such tid in a real BAM file
            let mut total_indels_in_current_contig = 0;
            let mut base;

            let mut process_previous_contigs = |last_tid,
                                                tid,
                                                depth,
                                                tet_freq,
                                                total_indels_in_current_contig| {
                if last_tid != -2 {
                    debug!("tid {} and {} indels", last_tid,
                           total_indels_in_current_contig);
                    println!("{} depth {:?}", std::str::from_utf8(target_names[last_tid as usize]).unwrap(), depth);
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
                        tid,
                        depth,
                        tet_freq,
                        total_indels_in_current_contig);
                    tet_freq = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    total_indels_in_current_contig = 0;
                }

                for alignment in pileup.alignments() {
                    if !alignment.is_del() && !alignment.is_refskip() {
                        base = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                        let count = tet_freq[pileup.pos() as usize].entry(base)
                            .or_insert(0);
                        *count += 1;
                        depth[pileup.pos() as usize] += 1
                    }
                    // mark indel start
                    match alignment.indel() {
                        bam::pileup::Indel::Ins(len) | bam::pileup::Indel::Del(len) => total_indels_in_current_contig += 1,
                        bam::pileup::Indel::None => ()
                    }
                }
            }
//        print!("{:?}", tet_freq);
            process_previous_contigs(
                last_tid,
                target_names.len() as i32,
                depth,
                tet_freq,
                total_indels_in_current_contig);
        }
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
