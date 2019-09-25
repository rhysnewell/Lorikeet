use std;
use std::collections::{HashMap, HashSet};
use std::collections::BTreeMap;
use rust_htslib::bam;

use pileup_structs::*;
use pileup_matrix::*;
use bam_generator::*;
use FlagFilter;
use std::str;

use std::fs;
use std::fs::File;
use std::io::prelude::*;

use bio::alignment::sparse::*;



pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    _print_zero_coverage_contigs: bool,
    _flag_filters: FlagFilter,
    depth_threshold: usize,
    var_fraction: f64,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32,
    variant_file_name: String,
    print_consensus: bool) {

//    let mut pileup_matrix = PileupMatrix::new_contig_stats();

    let consensus_variant_fasta = match File::create(variant_file_name.clone()) {
        Ok(fasta) => fasta,
        Err(e) => {
            println!("Cannot create file {:?}", e);
            std::process::exit(1)
        },
    };

    // Check to see if we are writing a consensus genome
    if !print_consensus {
        fs::remove_file(variant_file_name);
    }

    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        {

            let header = bam_generated.header().clone(); // bam header
            let target_names = header.target_names(); // contig names
            let bam_pileups = bam_generated.pileups(); // pileups for each genomic pos
            let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

            // for each genomic position, only has hashmap when variants are present. Includes read ids
            let mut nuc_freq: Vec<HashMap<char, HashSet<i32>>> = Vec::new();
            let mut indels = Vec::new();

            let mut read_starts = HashMap::new(); // read start map
            let mut depth = Vec::new(); // genomic depth
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


                    // adds contig info to pileup struct
                    pileup_struct.add_contig(nuc_freq,
                                               depth,
                                               indels,
                                               last_tid,
                                               total_indels_in_current_contig,
                                               contig_name,
                                               contig_len);

                    // calculates coverage across contig
                    pileup_struct.calc_coverage();

                    // filters variants across contig
                    pileup_struct.calc_variants(depth_threshold,
                                                var_fraction);
//                    pileup_struct.generate_variant_matrix();

                    // calculates minimum number of genotypes possible for each variant location
                    pileup_struct.generate_genotypes();

                    // prints results of variants calling
                    pileup_struct.print_variants(ref_sequence.clone(), depth_threshold);

                    if print_consensus {
                        // Write consensus contig to fasta
                        // i.e. the most abundant variant at each position from this set of reads
                        let contig_n = ">".to_owned() +
                            &str::from_utf8(target_names[last_tid as usize]).unwrap().to_string() +
                            "\n";

                        let mut consensus_clone = consensus_variant_fasta.try_clone().unwrap();
                        consensus_clone.write_all(contig_n.as_bytes()).unwrap();
                        pileup_struct.generate_variant_contig(ref_sequence.clone(),
                                                              consensus_clone);
                    }

//                    pileup_matrix.add_contig(pileup_struct,
//                                             target_names.len() as usize);

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
//                    contig_name = str::from_utf8(target_names[tid as usize]).unwrap().to_string();
                    total_indels_in_current_contig = 0;
                    match reference.fetch_all(std::str::from_utf8(target_names[tid as usize]).unwrap()) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)},
                    };
                    ref_seq = Vec::new();
                    match reference.read(&mut ref_seq) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)},
                    };

                } else {
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                }

                for alignment in pileup.alignments() {
                    // Check if new read to id

                    if !read_to_id.contains_key(&alignment.record().qname().to_vec()) {
                        read_to_id.entry(alignment.record().qname().to_vec())
                            .or_insert(read_cnt_id);
                        read_starts.entry(read_cnt_id).or_insert(pileup.pos());
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
                                Err(_e) => {"".to_string()},
                            };

                            let id = indels[pileup.pos() as usize].entry(insert)
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            total_indels_in_current_contig += 1;
                        },

                        bam::pileup::Indel::Del(len) => {

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
                indels,
                total_indels_in_current_contig,
                ref_seq);
        }
        bam_generated.finish();
    }

}

// Method for printing out contig statistics, including genotype per variant in metabat format
pub fn pileup_contigs<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    _print_zero_coverage_contigs: bool,
    _flag_filters: FlagFilter,
    depth_threshold: usize,
    var_fraction: f64,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32) {

    let mut pileup_matrix = PileupMatrix::new_matrix();

    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        let stoit_name = bam_generated.name().to_string();
        pileup_matrix.add_sample(stoit_name);
        {
            let header = bam_generated.header().clone(); // bam header
            let target_names = header.target_names(); // contig names
            let bam_pileups = bam_generated.pileups(); // pileups for each genomic pos
            let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

            // for each genomic position, only has hashmap when variants are present. Includes read ids
            let mut nuc_freq: Vec<HashMap<char, HashSet<i32>>> = Vec::new();
            let mut indels = Vec::new();

            let mut read_starts = HashMap::new(); // read start map
            let mut tet_freq = BTreeMap::new(); // kmer frequencies
            let mut depth = Vec::new(); // genomic depth
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

                    // adds contig info to pileup struct
                    pileup_struct.add_contig(nuc_freq,
                                             depth,
                                             indels,
                                             last_tid.clone(),
                                             total_indels_in_current_contig,
                                             contig_name,
                                             contig_len);

                    // calculates coverage across contig
                    pileup_struct.calc_coverage();

                    // filters variants across contig
                    pileup_struct.calc_variants(depth_threshold,
                                                var_fraction);

                    // calculates minimum number of genotypes possible for each variant location
                    pileup_struct.generate_genotypes();

                    pileup_matrix.add_kmers(last_tid,
                                            target_names.len() as usize,
                                            tet_freq);

                    pileup_matrix.add_contig(pileup_struct);

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
//                    contig_name = str::from_utf8(target_names[tid as usize]).unwrap().to_string();
                    total_indels_in_current_contig = 0;
                    match reference.fetch_all(std::str::from_utf8(target_names[tid as usize]).unwrap()) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)},
                    };
                    ref_seq = Vec::new();
                    match reference.read(&mut ref_seq) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)},
                    };
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
                        read_starts.entry(read_cnt_id).or_insert(pileup.pos());
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
                                Err(_e) => {"".to_string()},
                            };

                            let id = indels[pileup.pos() as usize].entry(insert)
                                                                  .or_insert(HashSet::new());
                            id.insert(read_to_id[&alignment.record().qname().to_vec()]);
                            total_indels_in_current_contig += 1;
                        },

                        bam::pileup::Indel::Del(len) => {

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
        bam_generated.finish();
    }
    pileup_matrix.print_stats();
    pileup_matrix.print_kmers();
}