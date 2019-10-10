use std;
use std::collections::{HashMap, HashSet};
use rust_htslib::bam;

use pileup_structs::*;
use pileup_matrix::*;
use bam_generator::*;
use FlagFilter;

use std::str;
use std::fs;
use std::fs::File;
use std::io::prelude::*;



pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    _print_zero_coverage_contigs: bool,
    _flag_filters: FlagFilter,
    depth_threshold: usize,
    mapq_threshold: u8,
    var_fraction: f64,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32,
    variant_file_name: String,
    print_consensus: bool,
    n_threads: usize) {

    let mut sample_idx = 0;
    // Print file header
    println!("tid\tpos\tvariant\treference\tabundance\tdepth\tgenotypes\tsample_id");
    // Loop through bam generators in parallel
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(n_threads);
        let stoit_name = bam_generated.name().to_string();

        let stoit_file_name = stoit_name + &variant_file_name;
        // Pre-emptively create variant fasta
        let consensus_variant_fasta = match File::create(&stoit_file_name) {
            Ok(fasta) => fasta,
            Err(e) => {
                println!("Cannot create file {:?}", e);
                std::process::exit(1)
            },
        };

        // Check to see if we are writing a consensus genome
        if !print_consensus {
            match fs::remove_file(&stoit_file_name) {
                Ok(removed) => removed,
                Err(_err) => {
                    println!("Incorrect file read/write permission");
                    std::process::exit(1)
                }
            };
        }


        {

            let header = bam_generated.header().clone(); // bam header
            let target_names = header.target_names(); // contig names
            let bam_pileups = bam_generated.pileups(); // pileups for each genomic pos
            let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

            // for each genomic position, only has hashmap when variants are present. Includes read ids
            let mut nuc_freq: Vec<HashMap<char, HashSet<i32>>> = Vec::new();
            let mut indels = Vec::new();

//            let mut read_starts = HashMap::new(); // read start map
            let mut depth = Vec::new(); // genomic depth
            let mut last_tid: i32 = -2; // no such tid in a real BAM file
            let mut total_indels_in_current_contig = 0;
            let mut read_cnt_id = 0;
            let mut read_to_id = HashMap::new();
            let mut previous_read_positions = HashMap::new();
//            let mut indel_start_sites: HashMap<i32, Vec<Indel>> = HashMap::new();
            let mut base;


            for p in bam_pileups {

                let pileup = p.unwrap();
                let tid = pileup.tid() as i32;
                // if reference has changed, print the last record
                if tid != last_tid {

                    if last_tid != -2 {
                        let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                        let contig_name = target_names[last_tid as usize].to_vec();

                        process_previous_contigs_var(
                            last_tid,
                            depth,
                            nuc_freq,
                            indels,
                            min, max,
                            total_indels_in_current_contig,
                            min_fraction_covered_bases,
                            contig_end_exclusion,
                            depth_threshold,
                            var_fraction,
                            contig_len,
                            contig_name,
                            ref_seq,
                            &consensus_variant_fasta,
                            print_consensus,
                            sample_idx);
                    }

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
                    if alignment.record().seq().len() == 0 {
                        continue
                    } else if alignment.record().mapq() < mapq_threshold {
                        continue
                    }
                    // Check if new read to id
                    if !read_to_id.contains_key(&alignment.record().qname().to_vec()) {
                        read_to_id.entry(alignment.record().qname().to_vec())
                            .or_insert(read_cnt_id);
//                        read_starts.entry(read_cnt_id).or_insert(pileup.pos());
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

            if last_tid != -2 {
                let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                let contig_name = target_names[last_tid as usize].to_vec();

                process_previous_contigs_var(
                    last_tid,
                    depth,
                    nuc_freq,
                    indels,
                    min, max,
                    total_indels_in_current_contig,
                    min_fraction_covered_bases,
                    contig_end_exclusion,
                    depth_threshold,
                    var_fraction,
                    contig_len,
                    contig_name,
                    ref_seq,
                    &consensus_variant_fasta,
                    print_consensus,
                    sample_idx);
            }
        }
        bam_generated.finish();
        sample_idx += 1;
    };
}

// Method for printing out contig statistics, including genotype per variant in metabat format
pub fn pileup_contigs<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    _print_zero_coverage_contigs: bool,
    _flag_filters: FlagFilter,
    depth_threshold: usize,
    mapq_threshold: u8,
    var_fraction: f64,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32,
    output_prefix: &str,
    n_threads: usize) {

    let mut pileup_matrix = PileupMatrix::new_matrix();
    let sample_count = bam_readers.len();
    let mut sample_idx = 0;

    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(n_threads);
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

//            let mut read_starts = HashMap::new(); // read start map
            let mut depth = Vec::new(); // genomic depth
            let mut last_tid: i32 = -2; // no such tid in a real BAM file
            let mut total_indels_in_current_contig = 0;
            let mut read_cnt_id = 0;
            let mut read_to_id = HashMap::new();
            let mut previous_read_positions = HashMap::new();
            let mut base;
            let mut read_name;

            for p in bam_pileups {

                let pileup = p.unwrap();
                let tid = pileup.tid() as i32;
                // if reference has changed, print the last record
                if tid != last_tid {
                    if last_tid != -2 {
                        let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                        let contig_name = target_names[last_tid as usize].to_vec();

                        process_previous_contigs(
                            last_tid,
                            depth,
                            nuc_freq,
                            indels,
                            min, max,
                            total_indels_in_current_contig,
                            min_fraction_covered_bases,
                            contig_end_exclusion,
                            depth_threshold,
                            var_fraction,
                            contig_len,
                            contig_name,
                            &mut pileup_matrix,
                            sample_count,
                            sample_idx);
                    };
                    nuc_freq = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    indels = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    total_indels_in_current_contig = 0;
                    match reference.fetch_all(std::str::from_utf8(target_names[tid as usize]).unwrap()) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)
                        },
                    };
                    ref_seq = Vec::new();
                    match reference.read(&mut ref_seq) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!("Cannot read sequence from reference {:?}", e);
                            std::process::exit(1)
                        },
                    };
                }  else {
                    depth[pileup.pos() as usize] = pileup.depth() as usize;
                }

                for alignment in pileup.alignments() {
                    if alignment.record().seq().len() == 0 {
                        continue
                    } else if alignment.record().mapq() < mapq_threshold {
                        continue
                    }
                    read_name = alignment.record().qname().to_vec();

                    // Check if new read to id
                    if !read_to_id.contains_key(&read_name) {
                        read_to_id.entry(read_name.clone())
                                  .or_insert(read_cnt_id);
//                        read_starts.entry(read_cnt_id).or_insert(pileup.pos());
                        read_cnt_id += 1;
                    }

                    let qpos = match alignment.qpos() {
                        Some(position) => {
                            previous_read_positions.insert(read_to_id[&read_name], position);
                            position
                        },
                        None => {
                            let position = previous_read_positions.entry(
                                read_to_id[&read_name]).or_insert(0);
                            *position
                        },
                    };
                    if !alignment.is_del() && !alignment.is_refskip() {

                        base = alignment.record().seq()[qpos] as char;
                        if base != ref_seq[pileup.pos() as usize] as char {
                            let id = nuc_freq[pileup.pos() as usize].entry(base)
                                                                    .or_insert(HashSet::new());
                            id.insert(read_to_id[&read_name]);
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
            if last_tid != -2 {
                let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
                let contig_name = target_names[last_tid as usize].to_vec();

                process_previous_contigs(
                    last_tid,
                    depth,
                    nuc_freq,
                    indels,
                    min, max,
                    total_indels_in_current_contig,
                    min_fraction_covered_bases,
                    contig_end_exclusion,
                    depth_threshold,
                    var_fraction,
                    contig_len,
                    contig_name,
                    &mut pileup_matrix,
                    sample_count,
                    sample_idx);
            };
        }
        bam_generated.finish();
        sample_idx += 1;
    }
    pileup_matrix.print_stats(output_prefix);
}

fn process_previous_contigs_var(
    last_tid: i32,
    depth: Vec<usize>,
    nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
    indels: Vec<HashMap<String, HashSet<i32>>>,
    min: f32, max: f32,
    total_indels_in_current_contig: usize,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32,
    depth_threshold: usize,
    var_fraction: f64,
    contig_len: usize,
    contig_name: Vec<u8>,
    ref_sequence: Vec<u8>,
    consensus_variant_fasta: &File,
    print_consensus: bool,
    sample_idx: i32) {

    if last_tid != -2 {

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
                                 contig_name.clone(),
                                 contig_len);

        // calculates coverage across contig
        pileup_struct.calc_coverage();

        // filters variants across contig
        pileup_struct.calc_variants(depth_threshold,
                                    var_fraction);

        // calculates minimum number of genotypes possible for each variant location
        pileup_struct.generate_genotypes();

        // prints results of variants calling
        pileup_struct.print_variants(ref_sequence.clone(), sample_idx);

        if print_consensus {
            // Write consensus contig to fasta
            // i.e. the most abundant variant at each position from this set of reads
            let contig_n = ">".to_owned() +
                &str::from_utf8(&contig_name).unwrap().to_string() +
                "\n";

            let mut consensus_clone = consensus_variant_fasta.try_clone().unwrap();
            consensus_clone.write_all(contig_n.as_bytes()).unwrap();
            pileup_struct.generate_variant_contig(ref_sequence.clone(),
                                                  consensus_clone);
        }
    }
}



fn process_previous_contigs(
        last_tid: i32,
        depth: Vec<usize>,
        nuc_freq: Vec<HashMap<char, HashSet<i32>>>,
        indels: Vec<HashMap<String, HashSet<i32>>>,
        min: f32, max: f32,
        total_indels_in_current_contig: usize,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u32,
        depth_threshold: usize,
        var_fraction: f64,
        contig_len: usize,
        contig_name: Vec<u8>,
        pileup_matrix: &mut PileupMatrix,
        sample_count: usize,
        sample_idx: usize) {

        if last_tid != -2 {

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

            pileup_matrix.add_contig(pileup_struct, sample_count, sample_idx);
    }
}
