use bio::io::gff;
use std;
use std::collections::{HashMap, HashSet};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use pileup_structs::*;
use pileup_matrix::*;
use bam_generator::*;
use FlagFilter;

use std::str;
use std::fs;
use std::fs::File;
use std::io::prelude::*;


pub fn predict_evolution<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    gff_reader: bio::io::gff::Reader<File>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter,
    mapq_threshold: u8,
    min_var_depth: usize,
    min: f32, max: f32,
    min_fraction_covered_bases: f32,
    contig_end_exclusion: u32,
    variant_file_name: String,
    print_consensus: bool,
    n_threads: usize,
    method: &str) {

    let mut sample_idx = 0;
    let include_soft_clipping = false;
    // Print file header
    println!("tid\tpos\tvariant\treference\tabundance\tdepth\tgenotypes\tsample_id");
    // Loop through bam generators in parallel
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(n_threads);
        let stoit_name = bam_generated.name().to_string();

        let stoit_file_name = stoit_name.clone() + &variant_file_name;
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

        let header = bam_generated.header().clone(); // bam header
        let target_names = header.target_names(); // contig names
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u32 = 0;

        let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        let mut nuc_freq: Vec<HashMap<char, HashSet<i32>>> = Vec::new();
        let mut indels = Vec::new();

        let mut depth = Vec::new(); // genomic depth
        let mut last_tid: i32 = -2; // no such tid in a real BAM file
        let mut total_indels_in_current_contig = 0;
        let mut read_cnt_id = 0;
        let mut read_to_id = HashMap::new();
        let mut base;

        // for record in records
        while bam_generated.read(&mut record)
            .expect("Error while reading BAM record") == true {
            debug!("Starting with a new read.. {:?}", record);
            if (!flag_filters.include_supplementary && record.is_supplementary()) ||
                (!flag_filters.include_secondary && record.is_secondary()) ||
                (!flag_filters.include_improper_pairs && !record.is_proper_pair()){
                debug!("Skipping read based on flag filtering");
                continue;
            }

            // if reference has changed, print the last record
            let tid = record.tid();
            if !record.is_unmapped() { // if mapped
                if record.seq().len() == 0 {
                    continue
                } else if record.mapq() < mapq_threshold {
                    continue
                }
                // Check if new read to id
                if !read_to_id.contains_key(&record.qname().to_vec()) {
                    read_to_id.entry(record.qname().to_vec())
                        .or_insert(read_cnt_id);
                    read_cnt_id += 1;
                }

                // if reference has changed, print the last record
                if tid != last_tid {
                    if last_tid != -2 {
                        let contig_len = header.target_len(last_tid as u32)
                            .expect("Corrupt BAM file?") as usize;
                        let contig_name = target_names[last_tid as usize].to_vec();
                        let total_mismatches = total_edit_distance_in_current_contig -
                            total_indels_in_current_contig;

                        process_previous_contigs_var(
                            last_tid,
                            depth,
                            nuc_freq,
                            indels,
                            min, max,
                            total_indels_in_current_contig as usize,
                            min_fraction_covered_bases,
                            contig_end_exclusion,
                            min_var_depth,
                            contig_len,
                            contig_name,
                            ref_seq,
                            &consensus_variant_fasta,
                            print_consensus,
                            sample_idx,
                            method,
                            total_mismatches);
                    }
                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    num_mapped_reads_total += num_mapped_reads_in_current_contig;
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    nuc_freq = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    depth = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    indels = vec![HashMap::new(); header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];

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
                }

                num_mapped_reads_in_current_contig += 1;

                // for each chunk of the cigar string
                debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
                let mut cursor: usize = record.pos() as usize;
                let mut read_cursor: usize = 0;
                for cig in record.cigar().iter() {
                    debug!("Found cigar {:} from {}", cig, cursor);
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or = increment start and decrement end index
                            debug!("Adding M, X, or = at {} and {}", cursor, cursor + cig.len() as usize);
                            ups_and_downs[cursor] += 1;
                            for qpos in read_cursor..(read_cursor+cig.len() as usize) {
                                base = record.seq()[qpos] as char;
                                if base != ref_seq[cursor as usize] as char {
                                    let id = nuc_freq[cursor as usize].entry(base)
                                        .or_insert(HashSet::new());
                                    id.insert(read_to_id[&record.qname().to_vec()]);
                                }
                                depth[cursor] += 1;
                                cursor += 1;
                            }
                            if cursor < ups_and_downs.len() { // True unless the read hits the contig end.
                                ups_and_downs[cursor] -= 1;
                            }
                            read_cursor += cig.len() as usize;
                        },
                        Cigar::Del(_) => {
                            let id = indels[cursor as usize].entry(
                                std::iter::repeat("N").take(cig.len() as usize).collect::<String>())
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&record.qname().to_vec()]);

                            cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len();
                        },
                        Cigar::RefSkip(_) => {
                            // if D or N, move the cursor
                            cursor += cig.len() as usize;
                        },
                        Cigar::Ins(_) => {
                            let insert = match str::from_utf8(&record.seq().as_bytes()[
                                read_cursor..read_cursor + cig.len() as usize]) {
                                Ok(ins) => {ins.to_string()},
                                Err(_e) => {"".to_string()},
                            };

                            let id = indels[cursor as usize].entry(insert)
                                .or_insert(HashSet::new());
                            id.insert(read_to_id[&record.qname().to_vec()]);
                            read_cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len();
                        },
                        Cigar::SoftClip(_) => {
                            // soft clipped portions of reads can be included as insertions
                            // not sure if this correct protocol or not
                            if include_soft_clipping {
                                let insert = match str::from_utf8(&record.seq().as_bytes()[
                                    read_cursor..read_cursor + cig.len() as usize]) {
                                    Ok(ins) => {ins.to_string()},
                                    Err(_e) => {"".to_string()},
                                };

                                let id = indels[cursor as usize].entry(insert)
                                    .or_insert(HashSet::new());
                                id.insert(read_to_id[&record.qname().to_vec()]);
                                total_indels_in_current_contig += cig.len();
                            }
                            read_cursor += cig.len() as usize;
                        },
                        Cigar::HardClip(_) | Cigar::Pad(_) => {}
                    }
                }
                // Determine the number of mismatching bases in this read by
                // looking at the NM tag.
                total_edit_distance_in_current_contig += match
                    record.aux("NM".as_bytes()) {
                    Some(aux) => {
                        aux.integer() as u32
                    },
                    None => {
                        panic!("Mapping record encountered that does not have an 'NM' \
                                auxiliary tag in the SAM/BAM format. This is required \
                                to work out some coverage statistics");
                    }
                };

                debug!("At end of loop")
            }
        } if last_tid != -2 {
            let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
            let contig_name = target_names[last_tid as usize].to_vec();
            let total_mismatches = total_edit_distance_in_current_contig -
                total_indels_in_current_contig;

            process_previous_contigs_var(
                last_tid,
                depth,
                nuc_freq,
                indels,
                min, max,
                total_indels_in_current_contig as usize,
                min_fraction_covered_bases,
                contig_end_exclusion,
                min_var_depth,
                contig_len,
                contig_name,
                ref_seq,
                &consensus_variant_fasta,
                print_consensus,
                sample_idx,
                method,
                total_mismatches);

            num_mapped_reads_total += num_mapped_reads_in_current_contig;
        }


        info!("In sample '{}', found {} reads mapped out of {} total ({:.*}%)",
              stoit_name, num_mapped_reads_total,
              bam_generated.num_detected_primary_alignments(), 2,
              (num_mapped_reads_total * 100) as f64 /
                  bam_generated.num_detected_primary_alignments() as f64);

        if bam_generated.num_detected_primary_alignments() == 0 {
            warn!("No primary alignments were observed for sample {} \
               - perhaps something went wrong in the mapping?",
                  stoit_name);
        }
        bam_generated.finish();
        sample_idx += 1;
    };
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
    min_var_depth: usize,
    contig_len: usize,
    contig_name: Vec<u8>,
    ref_sequence: Vec<u8>,
    consensus_variant_fasta: &File,
    print_consensus: bool,
    sample_idx: i32,
    method: &str,
    total_mismatches: u32) {

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
        pileup_struct.calc_coverage(total_mismatches, method);

        // filters variants across contig
        pileup_struct.calc_variants(
            min_var_depth);

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