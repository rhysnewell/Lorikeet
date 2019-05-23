use std;
use std::result::Result;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use std::str;
use std::collections::BTreeSet;

//use mosdepth_genome_coverage_estimators::*;
use genomes_and_contigs::GenomesAndContigs;
use bam_generator::*;
use ReadsMapped;

pub fn mosdepth_genome_coverage_with_contig_names<R: NamedBamReader,
                                                  G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    contigs_and_genomes: &GenomesAndContigs,
    proper_pairs_only: bool)
    -> Vec<ReadsMapped> {

    let mut reads_mapped_vector = vec!();
    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        // Collect reference numbers for each genome's contigs
        let mut reference_number_to_genome_index: Vec<Option<usize>> = vec![];
        let mut num_refs_in_genomes: u32 = 0;
        let mut num_refs_not_in_genomes: u32 = 0;
        // Collect reference numbers for each genome
        let mut genome_index_to_references: Vec<Vec<u32>> =
            vec![vec!(); contigs_and_genomes.genomes.len()];
        // Reads mapped are only counted when the genome has non-zero coverage.
        let mut reads_mapped_in_each_genome: Vec<u64> = vec!(
            0; contigs_and_genomes.genomes.len());
        for (tid, name) in target_names.iter().enumerate() {
            let genome_index = contigs_and_genomes.genome_index_of_contig(
                &String::from(std::str::from_utf8(name)
                              .expect("UTF8 encoding error in BAM header file")));

            match genome_index {
                Some(i) => {
                    reference_number_to_genome_index.push(Some(i));
                    num_refs_in_genomes += 1;
                    genome_index_to_references[i].push(tid as u32);
                },
                None => {
                    reference_number_to_genome_index.push(None);
                    num_refs_not_in_genomes += 1;
                }
            }
        }
        info!("Of {} reference IDs, {} were assigned to a genome and {} were not",
              num_refs_in_genomes + num_refs_not_in_genomes,
              num_refs_in_genomes, num_refs_not_in_genomes);
        debug!("Reference number to genomes: {:?}", reference_number_to_genome_index);
        if num_refs_in_genomes == 0 {
            eprintln!("Error: There are no found reference sequences that are a part of a genome");
            std::process::exit(2);
        }


        // Iterate through bam records
        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut seen_ref_ids = BTreeSet::new();
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u32 = 0;
        let mut total_indels_in_current_contig: u32 = 0;
        while bam_generated.read(&mut record).is_ok() {
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }
            if proper_pairs_only && !record.is_proper_pair() {
                continue;
            }
            let original_tid = record.tid();
            if !record.is_unmapped() { // if mapped
                let tid = original_tid as u32;
                if tid != last_tid || doing_first {
                    debug!("Came across a new tid {}", tid);
                    if doing_first == true {
                        doing_first = false;
                    } else {

                    }

                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    last_tid = tid;
                    seen_ref_ids.insert(tid);
                }

                // TODO: move below into a function for code-reuse purposes.
                // Add coverage info for the current record
                // for each chunk of the cigar string
                match reference_number_to_genome_index[tid as usize] {
                    None => {},
                    Some(genome_index) => {
                        reads_mapped_in_each_genome[genome_index] += 1;
                        num_mapped_reads_in_current_contig += 1;
                        debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
                        let mut cursor: usize = record.pos() as usize;
                        for cig in record.cigar().iter() {
                            debug!("Found cigar {:} from {}", cig, cursor);
                            match cig {
                                Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                                    // if M, X, or =, increment start and decrement end index
                                    debug!("Adding M, X, or =, at {} and {}", cursor, cursor + cig.len() as usize);
                                    ups_and_downs[cursor] += 1;
                                    let final_pos = cursor + cig.len() as usize;
                                    if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                                        ups_and_downs[final_pos] -= 1;
                                    }
                                    cursor += cig.len() as usize;
                                },
                                Cigar::Del(_) => {
                                    cursor += cig.len() as usize;
                                    total_indels_in_current_contig += cig.len() as u32;
                                },
                                Cigar::RefSkip(_) => {
                                    cursor += cig.len() as usize;
                                },
                                Cigar::Ins(_) => {
                                    total_indels_in_current_contig += cig.len() as u32;
                                },
                                Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
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
                    }
                }
            }
        }

        let mut num_mapped_reads_total: u64 = 0;
        if doing_first && bam_generated.num_detected_primary_alignments() == 0 {
            warn!("No primary alignments were observed for sample {} \
                   - perhaps something went wrong in the mapping?",
                  stoit_name);
        } else {
            // Record the last contig
            match reference_number_to_genome_index[last_tid as usize] {
                Some(genome_index) => {
                },
                None => {}
            }

            // Print the coverages of each genome
            // Calculate the unobserved length of each genome
            let mut unobserved_lengths: Vec<u32> = vec!();
            for _ in 0..contigs_and_genomes.genomes.len() {
                unobserved_lengths.push(0)
            }
            for (ref_id, genome_id_option) in reference_number_to_genome_index.iter().enumerate() {
                let ref_id_u32: u32 = ref_id as u32;
                debug!("Seen {:?}", seen_ref_ids);
                match genome_id_option {
                    Some(genome_id) => {
                        if !seen_ref_ids.contains(&ref_id_u32) {
                            debug!("Getting target #{} from header names", ref_id_u32);
                            unobserved_lengths[*genome_id] += header.target_len(ref_id_u32).unwrap()
                        }
                    },
                    None => {}
                }
            }
            // print the genomes out
            for (i, genome) in contigs_and_genomes.genomes.iter().enumerate() {
                // Determine if any coverages are non-zero
            }
        }

        let reads_mapped = ReadsMapped {
            num_mapped_reads: num_mapped_reads_total,
            num_reads: bam_generated.num_detected_primary_alignments()
        };
        info!("In sample '{}', found {} reads mapped out of {} total ({:.*}%)",
              stoit_name, reads_mapped.num_mapped_reads,
              reads_mapped.num_reads, 2,
              (reads_mapped.num_mapped_reads * 100) as f64 / reads_mapped.num_reads as f64);
        reads_mapped_vector.push(reads_mapped);

        bam_generated.finish();
    }
    return reads_mapped_vector;
}



struct UnobservedLengthAndFirstTid {
    unobserved_contig_length: u32,
    first_tid: usize
}

fn print_last_genomes(
    num_mapped_reads_in_current_contig: u64,
    last_genome: &[u8],
    unobserved_contig_length_and_first_tid: &mut UnobservedLengthAndFirstTid,
    ups_and_downs: &Vec<i32>,
    total_edit_distance_in_current_contig: u32,
    total_indels_in_current_contig: u32,
    current_genome: &[u8],
    single_genome: bool,
    target_names: &std::vec::Vec<&[u8]>,
    split_char: u8,
    header: &rust_htslib::bam::HeaderView,
    tid_to_print_zeros_to: u32) {

//    debug!("ups_and_downs {:?}", &ups_and_downs);

    // Determine coverage of previous genome

}



pub fn mosdepth_genome_coverage<R: NamedBamReader,
                                G: NamedBamReaderGenerator<R>> (
    bam_readers: Vec<G>,
    split_char: u8,
    proper_pairs_only: bool,
    single_genome: bool)
    -> Vec<ReadsMapped> {
    let mut reads_mapped_vector = vec!();
    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        let fill_genome_length_forwards = |current_tid, target_genome| {
            // Iterating reads skips over contigs with no mapped reads, but the
            // length of these contigs is required to calculate the average
            // across all contigs. This closure returns the number of bases in
            // contigs with tid > current_tid that are part of the current
            // genome.
            let mut extra: u32 = 0;
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                if single_genome ||
                    extract_genome(my_tid, &target_names, split_char) == target_genome {

                    extra += header.target_len(my_tid)
                        .expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
                } else {
                    break;
                }
            }
            return extra
        };

        let fill_genome_length_backwards_to_last = |current_tid, last_tid, target_genome| {
            if current_tid == 0 {return 0};
            let mut extra: u32 = 0;
            let mut my_tid = last_tid + 1;
            while my_tid < current_tid {
                if single_genome ||
                    extract_genome(my_tid, &target_names, split_char) == target_genome {

                    extra += header.target_len(my_tid)
                        .expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
                } else {
                    break;
                }
            }
            return extra
        };

        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: &[u8] = "error genome".as_bytes();
        let mut unobserved_contig_length_and_first_tid = UnobservedLengthAndFirstTid {
            unobserved_contig_length: 0,
            first_tid: 0
        };
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut num_mapped_reads_in_current_genome: u64 = 0;
        let mut total_edit_distance_in_current_contig: u32 = 0;
        let mut total_indels_in_current_contig: u32 = 0;
        while bam_generated.read(&mut record).is_ok() {
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }
            if proper_pairs_only && !record.is_proper_pair() {
                continue;
            }
            let original_tid = record.tid();
            if !record.is_unmapped() {
                // if reference has changed, finish a genome or not
                let tid = original_tid as u32;
                let current_genome: &[u8] = match single_genome {
                    true => "".as_bytes(),
                    false => extract_genome(tid as u32, &target_names, split_char)
                };
                if tid != last_tid || doing_first {
                    debug!("Processing a change in tid, from {} to {} (first is {}). Current \
                           unobserved_and_first: unobserved {}, first {}",
                          last_tid,
                          tid,
                          doing_first,
                          unobserved_contig_length_and_first_tid.unobserved_contig_length,
                          unobserved_contig_length_and_first_tid.first_tid);
                    if doing_first == true {
                        unobserved_contig_length_and_first_tid = fill_genome_length_backwards(
                            tid,
                            current_genome,
                            single_genome,
                            &target_names,
                            split_char,
                            &header);
                        last_genome = current_genome;
                        debug!("doing first..");
                        doing_first = false;


                    } else if current_genome == last_genome {
                        debug!("Found {} reads mapped to tid {}",
                               num_mapped_reads_in_current_contig, last_tid);
                        // Collect the length of reference sequences from this
                        // genome that had no hits that were just skipped over.
                        debug!("Filling unobserved from {} to {}", last_tid, tid);
                        unobserved_contig_length_and_first_tid.unobserved_contig_length +=
                            fill_genome_length_backwards_to_last(
                                tid, last_tid as u32, current_genome);
                    } else {
                        debug!("Found {} reads mapped to tid {}",
                               num_mapped_reads_in_current_contig, last_tid);
                        // Collect the length of refs from the end of the last genome that had no hits
                        debug!("Filling unobserved from {} to {} for {}",
                               last_tid, tid, &str::from_utf8(last_genome).unwrap());
                        unobserved_contig_length_and_first_tid.unobserved_contig_length +=
                            fill_genome_length_backwards_to_last(
                                tid, last_tid as u32, last_genome);

//                        let positive_coverage = print_last_genomes(
//                            num_mapped_reads_in_current_contig,
//                            last_genome,
//                            &mut unobserved_contig_length_and_first_tid,
//                            &ups_and_downs,
//                            total_edit_distance_in_current_contig,
//                            total_indels_in_current_contig,
//                            current_genome,
//                            single_genome,
//                            &target_names,
//                            split_char,
//                            &header,
//                            tid,
//                        );

                        num_mapped_reads_in_current_genome = 0;
                        last_genome = current_genome;

                        unobserved_contig_length_and_first_tid = fill_genome_length_backwards(
                            tid,
                            current_genome,
                            single_genome,
                            &target_names,
                            split_char,
                            &header);
                        debug!(
                            "Setting unobserved contig length to be {}",
                            unobserved_contig_length_and_first_tid.unobserved_contig_length);
                    }

                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    last_tid = tid;
                }

                // Add coverage info for the current record
                // for each chunk of the cigar string
                debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
                num_mapped_reads_in_current_contig += 1;
                num_mapped_reads_in_current_genome += 1;
                let mut cursor: usize = record.pos() as usize;
                for cig in record.cigar().iter() {
                    debug!("Found cigar {:} from {}", cig, cursor);
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or =, increment start and decrement end index
                            debug!("Adding M, X, or =, at {} and {}", cursor, cursor + cig.len() as usize);
                            ups_and_downs[cursor] += 1;
                            let final_pos = cursor + cig.len() as usize;
                            if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                                ups_and_downs[final_pos] -= 1;
                            }
                            cursor += cig.len() as usize;
                        },
                        Cigar::Del(_) => {
                            cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len() as u32;
                        },
                        Cigar::RefSkip(_) => {
                            cursor += cig.len() as usize;
                        },
                        Cigar::Ins(_) => {
                            total_indels_in_current_contig += cig.len() as u32;
                        },
                        Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
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
            }
        }

        if doing_first && bam_generated.num_detected_primary_alignments() == 0 {
            warn!("No primary alignments were observed for sample {} \
                   - perhaps something went wrong in the mapping?",
                  stoit_name);
        } else {
            // Print the last genome
            // Give the single genome a dummy name
            if single_genome {
                last_genome = "genome1".as_bytes()
            }

            debug!("Found {} reads mapped to tid {}",
                   num_mapped_reads_in_current_contig, last_tid);
            // Collect the length of refs from the end of the last genome that had no hits
            debug!("Filling unobserved from {} to end for {}",
                   last_tid, &str::from_utf8(last_genome).unwrap());
            unobserved_contig_length_and_first_tid.unobserved_contig_length +=
                fill_genome_length_forwards(last_tid, last_genome);

//            let positive_coverage = print_last_genomes(
//                num_mapped_reads_in_current_contig,
//                last_genome,
//                &mut unobserved_contig_length_and_first_tid,
//                &ups_and_downs,
//                total_edit_distance_in_current_contig,
//                total_indels_in_current_contig,
//                b"",
//                single_genome,
//                &target_names,
//                split_char,
//                &header,
//                header.target_count() - 1,
//            );
        }

        let reads_mapped = ReadsMapped {
            num_mapped_reads: num_mapped_reads_total,
            num_reads: bam_generated.num_detected_primary_alignments()
        };
        info!("In sample '{}', found {} reads mapped out of {} total ({:.*}%)",
              stoit_name, reads_mapped.num_mapped_reads,
              reads_mapped.num_reads, 2,
              (reads_mapped.num_mapped_reads * 100) as f64 / reads_mapped.num_reads as f64);
        reads_mapped_vector.push(reads_mapped);

        bam_generated.finish();
    }
    return reads_mapped_vector;
}



fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    debug!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 str::from_utf8(target_name).unwrap()));
    return &target_name[(0..offset)];
}

fn fill_genome_length_backwards(
    current_tid: u32,
    target_genome: &[u8],
    single_genome: bool,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    header: &bam::HeaderView)
    -> UnobservedLengthAndFirstTid {
    debug!("At start of fill_genome_length_backwards, found current_tid {}, target_genome {}",
          current_tid, str::from_utf8(target_genome).unwrap());
    if current_tid == 0 {
        return UnobservedLengthAndFirstTid {
            unobserved_contig_length: 0,
            first_tid: current_tid as usize
        }
    }

    let mut extra: u32 = 0;
    let mut my_tid = current_tid - 1;
    while single_genome ||
        extract_genome(my_tid, &target_names, split_char) == target_genome {
            extra += header.target_len(my_tid)
                .expect("Malformed bam header or programming error encountered");
            if my_tid == 0 {
                return UnobservedLengthAndFirstTid {
                    unobserved_contig_length: extra,
                    first_tid: 0 as usize
                }
            } else {
                my_tid -= 1;
            }
        }
    debug!("Returning UnobservedLengthAndFirstTid length {}, first_tid {}",
          extra, my_tid+1);
    return UnobservedLengthAndFirstTid {
        unobserved_contig_length: extra,
        first_tid: (my_tid+1) as usize
    }
}

// Print zero coverage for genomes that have no reads mapped. Genomes are
// detected from the header, counting backwards from the current tid until the
// last seen genome is encountered, or we reach the beginning of the tid array.
fn print_previous_zero_coverage_genomes2<'a>(
    last_genome: &[u8],
    current_genome: &[u8],
    current_tid: u32,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    header: &bam::HeaderView){

    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    let mut genomes_to_print: Vec<&[u8]> = vec![];
    let mut genome_first_tids: Vec<usize> = vec![];
    let mut genomes_unobserved_length: Vec<u32> = vec![];
    let mut unobserved_length = 0;
    // Need to record the first TID from each genome, but we are iterating down.
    // Gah.
    let mut last_first_id = None;
    loop {
        let genome = extract_genome(tid, &target_names, split_char);
        debug!("in print_previous_zero_coverage_genomes2: tid {}, genome {}", tid, str::from_utf8(genome).unwrap());
        if genome == last_genome { break; }
        else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            // Push the last
            match last_first_id {
                Some(id) => {
                    if genome != last_genome {
                        genome_first_tids.push(id as usize);
                        genomes_to_print.push(my_current_genome);
                        genomes_unobserved_length.push(unobserved_length);
                    }
                },
                None => {}
            }
            my_current_genome = genome;
            last_first_id = Some(tid);
            unobserved_length = header.target_len(tid).unwrap();
        } else if genome != current_genome {
            last_first_id = Some(tid);
            unobserved_length += header.target_len(tid).unwrap();
        }
        if tid == 0 { break; }
        tid = tid - 1;
    };
    match last_first_id {
        Some(id) => {
            genome_first_tids.push(id as usize);
            genomes_to_print.push(my_current_genome);
            genomes_unobserved_length.push(unobserved_length);},
        None => {}
    }
    debug!("genomes_to_print {:?}, genome_first_tids {:?}: unobserved: {:?}",
           genomes_to_print,
           genome_first_tids,
           genomes_unobserved_length);
}

/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
    where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}
