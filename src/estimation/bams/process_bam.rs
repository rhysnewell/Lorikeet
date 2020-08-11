use rust_htslib::bam::{self, record::Cigar};
use std;
use std::collections::{HashMap, HashSet};

use coverm::bam_generator::*;
use estimation::codon_structs::*;
use estimation::contig_variants::*;
use estimation::variant_matrix::*;
use model::variants::*;
use rayon::prelude::*;
use utils::*;

use coverm::genomes_and_contigs::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use std::str;
use tempfile::NamedTempFile;

/// Process all reads in a BAM file
#[allow(unused)]
pub fn process_bam<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bam_generator: G,
    sample_count: usize,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    variant_matrix: &mut VariantMatrix,
    split_threads: usize,
    m: &clap::ArgMatches,
    output_prefix: &str,
    coverage_fold: f32,
    codon_table: &CodonTable,
    min_var_depth: usize,
    contig_end_exclusion: u64,
    min: f32,
    max: f32,
    ani: f32,
    mode: &str,
    include_soft_clipping: bool,
    include_indels: bool,
    flag_filters: &FlagFilter,
    mapq_threshold: u8,
    method: &str,
    sample_groups: &HashMap<&str, HashSet<String>>,
    genomes_and_contigs: &GenomesAndContigs,
    reference_map: &HashMap<usize, String>,
    concatenated_genomes: &Option<NamedTempFile>,
) {
    let mut bam_generated = bam_generator.start();

    // Adjust the sample index if the bam is from long reads
    let mut longread = false;

    //    let mut bam_properties =
    //        AlignmentProperties::default(InsertSize::default());

    let stoit_name = bam_generated.name().to_string().replace("/", ".");
    // adjust sample index for longread bams

    bam_generated.set_threads(split_threads);
    debug!("Managed to set threads.");

    let header = bam_generated.header().clone(); // bam header
    let target_names = header.target_names(); // contig names

    for target_name in target_names.iter() {
        debug!("Target name {:?}", &str::from_utf8(target_name));
    }

    let mut record: bam::record::Record = bam::Record::new();
    let mut ups_and_downs: Vec<i32> = Vec::new();
    // Current contig name
    let mut contig_name = Vec::new();
    let mut contig_name_str = String::new();
    let mut num_mapped_reads_total: u64 = 0;
    let mut num_mapped_reads_in_current_contig: u64 = 0;
    let mut total_edit_distance_in_current_contig: u64 = 0;

    let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig
    let mut ref_idx = 0;
    let mut last_tid: i32 = -2; // no such tid in a real BAM file
    let mut total_indels_in_current_contig = 0;

    let sample_idx = match variant_matrix {
        VariantMatrix::VariantContigMatrix { sample_names, .. } => {
            debug!(
                "sample names {:?} and stoit_name {:?}",
                &sample_names, &stoit_name
            );
            sample_names
                .iter()
                .position(
                    |p|
                        p.contains(&stoit_name.replace("lorikeet-genome", "")))
                .unwrap()
        }
    };

    debug!("sample groups {:?}", sample_groups);
    if sample_groups.contains_key("long") {
        for longread_name in sample_groups["long"].iter() {
            if stoit_name.contains(longread_name) {
                longread = true;
                debug!("Longread {} {}", stoit_name, sample_idx);
                break;
            }
        }
    }

    let mut reference = match concatenated_genomes {
        Some(reference_path) => {
            match bio::io::fasta::IndexedReader::from_file(&reference_path.path()) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(&reference_path.path().to_str().unwrap()),
            }
        }
        None => panic!("Concatenated reference file does not exist"),
    };

    // for record in records
    let mut skipped_reads = 0;

    while bam_generated
        .read(&mut record)
        .expect("Error while reading BAM record")
        == true
    {
        if (!flag_filters.include_supplementary && record.is_supplementary() && !longread)
            || (!flag_filters.include_secondary && record.is_secondary() && !longread)
            || (!flag_filters.include_improper_pairs && !record.is_proper_pair() && !longread)
        {
            skipped_reads += 1;
            continue;
        } else if !flag_filters.include_secondary && record.is_secondary() && longread {
            skipped_reads += 1;
            continue;
        }
        // if reference has changed, print the last record
        let tid = record.tid();
        if !record.is_unmapped() {
            // if mapped
            if record.seq().len() == 0 {
                skipped_reads += 1;
                continue;
            } else if record.mapq() < mapq_threshold {
                skipped_reads += 1;
                continue;
            }

            // if reference has changed, print the last record
            if tid != last_tid {
                if tid < last_tid {
                    error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                    panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                }
                if last_tid != -2 {
                    let contig_len = header
                        .target_len(last_tid as u32)
                        .expect("Corrupt BAM file?") as usize;


                    let total_mismatches =
                        total_edit_distance_in_current_contig - total_indels_in_current_contig;

                    match variant_matrix.variants_of_contig(ref_idx, last_tid) {
                        Some(map) => {
                            debug!("Variant Matrix {:?}", map);
                        },
                        None => {
                            debug!("Ref idx {} and tid {}", ref_idx, last_tid);

                        }
                    }

                    // Retrieve the reference based on the reference index from reference_map
                    // let reference_path = reference_map.get(&ref_idx).expect("Unable to retrieve reference path");

                    match reference.fetch_all(std::str::from_utf8(&contig_name[..]).unwrap()) {
                        Ok(reference) => reference,
                        Err(e) => match reference.fetch_all(&format!(
                            "{}~{}",
                            &genomes_and_contigs.genomes[ref_idx],
                            std::str::from_utf8(&contig_name[..]).unwrap()
                        )) {
                            Ok(reference) => reference,
                            Err(e) => {
                                println!(
                                    "Cannot read sequence from reference {} {:?}",
                                    format!(
                                        "{}~{}",
                                        &genomes_and_contigs.genomes[ref_idx],
                                        std::str::from_utf8(&contig_name[..]).unwrap()
                                    ),
                                    e,
                                );
                                std::process::exit(1);
                            }
                        },
                    };
                    ref_seq = Vec::new();
                    match reference.read(&mut ref_seq) {
                        Ok(reference) => reference,
                        Err(e) => {
                            println!(
                                "Cannot read sequence from reference {} {:?}",
                                std::str::from_utf8(&contig_name[..]).unwrap(),
                                e,
                            );
                            std::process::exit(1)
                        }
                    };

                    process_previous_contigs_var(
                        mode,
                        ref_idx,
                        last_tid,
                        ups_and_downs,
                        coverage_estimators,
                        min,
                        max,
                        total_indels_in_current_contig as usize,
                        contig_end_exclusion,
                        contig_len,
                        contig_name,
                        genomes_and_contigs,
                        variant_matrix,
                        ref_seq,
                        sample_idx,
                        method,
                        total_mismatches,
                        coverage_fold,
                        num_mapped_reads_in_current_contig,
                        sample_count,
                        output_prefix,
                        &stoit_name,
                    );
                }
                contig_name = target_names[tid as usize].to_vec();
                // contig_name_str = split_contig_name(&contig_name);
                ref_idx =
                    retrieve_reference_index_from_contig(&contig_name, genomes_and_contigs);

                ups_and_downs =
                    vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];

                contig_name = target_names[tid as usize].to_vec();
                // contig_name_str = split_contig_name(&contig_name);
                // debug!("Contig name {:?}", &contig_name_str);

                // ref_idx = retrieve_reference_index_from_contig(
                //     &contig_name,
                //     genomes_and_contigs,
                // );
                last_tid = tid;
                num_mapped_reads_total += num_mapped_reads_in_current_contig;
                num_mapped_reads_in_current_contig = 0;
                total_edit_distance_in_current_contig = 0;
                total_indels_in_current_contig = 0;
                //                    nuc_freq = HashMap::new();
                //                    indels = HashMap::new();
            }

            if !record.is_supplementary() {
                num_mapped_reads_in_current_contig += 1;
            }


            // for each chunk of the cigar string
            let mut cursor: usize = record.pos() as usize;
            let quals = record.qual();
            let mut read_cursor: usize = 0;
            let read_len = record.seq().len();
            for cig in record.cigar().iter() {
                match cig {
                    Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                        // if M, X, or = increment start and decrement end index
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        // For checking against mnv
                        let mut potential_mnv = false;
                        let mut mnv_pos = 0;
                        let mut mnv = vec![];
                        let mut mnv_cursor = 0;
                        for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                            // See if read is match MNV
                            if potential_mnv && (mnv_pos < mnv.len()) {
                                let read_char = record.seq()[qpos];
                                debug!("MNV searching {} {:?} {}", &mnv_pos, &mnv, &read_char);
                                if mnv[mnv_pos] == read_char {
                                    mnv_pos += 1;
                                    debug!("pos {} length {}", &mnv_pos, &mnv.len());
                                    if mnv_pos == mnv.len() {
                                        match variant_matrix.variants(ref_idx, tid, mnv_cursor) {
                                            Some(current_variants) => {
                                                current_variants
                                                    .iter_mut()
                                                    .for_each(|(variant, base)|
                                                        match variant {
                                                            Variant::MNV(alt) => {
                                                                debug!("alt {:?} found {:?}", &alt, &mnv);
                                                                if *alt == mnv {
                                                                    base.assign_read(
                                                                        record.qname().to_vec(),
                                                                    );
                                                                    base.truedepth[sample_idx] += 1;
                                                                }
                                                            }
                                                            _ => {
                                                                debug!("Looping through non-MNV variants");
                                                            }
                                                        });
                                                mnv = vec![];
                                                mnv_pos = 0;
                                                potential_mnv = false;
                                            },
                                            None => {
                                                debug!("No associated MNV found for {:?}", &mnv);
                                                mnv = vec![];
                                                mnv_pos = 0;
                                                potential_mnv = false
                                            }
                                        };
                                        mnv = vec![];
                                        mnv_pos = 0;
                                        potential_mnv = false;
                                    }
                                } else {
                                    debug!("Read did not contain correct MNV");
                                    mnv = vec![];
                                    mnv_pos = 0;
                                    potential_mnv = false
                                }
                            }
                            match variant_matrix.variants(ref_idx, tid, cursor as i64) {
                                Some(current_variants) => {
                                    let read_char = record.seq()[qpos];
                                    current_variants.iter_mut().for_each(|(variant, base)| {
                                        match variant {
                                            Variant::SNV(alt) => {
                                                if *alt == read_char {
                                                    base.assign_read(record.qname().to_vec());
                                                    base.truedepth[sample_idx] += 1;
                                                }
                                            },
                                            // We need to check every position of the MNV
                                            Variant::MNV(alt) => {
                                                if !potential_mnv {
                                                    if alt[mnv_pos] == read_char {
                                                        mnv = alt.clone();
                                                        debug!("Potential MNV  pos {} var {:?} read {} ref {:?}", &mnv_pos, &mnv, &read_char, &base.refr);

                                                        mnv_pos += 1;
                                                        potential_mnv = true;
                                                        mnv_cursor = cursor as i64;

                                                        // Then it is automatically assigned
                                                        if mnv_pos == mnv.len() {
                                                            base.assign_read(record.qname().to_vec());
                                                            base.truedepth[sample_idx] += 1;
                                                            mnv = vec!();
                                                            mnv_pos = 0;
                                                            potential_mnv = false;
                                                        }
                                                    }
                                                }
                                            },
                                            Variant::None => {
                                                if base.refr[0] == read_char {
                                                    base.assign_read(record.qname().to_vec());
                                                    base.truedepth[sample_idx] += 1;
                                                    // info!(
                                                    //     "Reference ref {} tid {} pos {} coverages {:?}",
                                                    //     &ref_idx,
                                                    //     &tid,
                                                    //     &cursor,
                                                    //     &base.truedepth,
                                                    // )
                                                } else {
                                                    mnv = vec!();
                                                    mnv_pos = 0;
                                                    potential_mnv = false
                                                }
                                            },
                                            _ => {}
                                        }
                                    });
                                }
                                _ => {}
                            }
                            cursor += 1;
                        }

                        // debug!("CIGAR ended, resetting MNV");
                        mnv = vec![];
                        mnv_pos = 0;
                        potential_mnv = false;
                        
                        if final_pos < ups_and_downs.len() {
                            // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        read_cursor += cig.len() as usize;
                    }
                    Cigar::Del(del) => {
                        match variant_matrix.variants(ref_idx, tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Deletion(alt) => {
                                            if alt == del {
                                                base.assign_read(record.qname().to_vec());
                                                base.truedepth[sample_idx] += 1;
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }

                        cursor += cig.len() as usize;
                    }
                    Cigar::RefSkip(_) => {
                        // if D or N, move the cursor
                        cursor += cig.len() as usize;
                    }
                    Cigar::Ins(ins) => {
                        let insertion =
                            &record.seq().as_bytes()[read_cursor..read_cursor + cig.len() as usize];
                        match variant_matrix.variants(ref_idx, tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Insertion(alt) => {
                                            if String::from_utf8(alt.to_vec())
                                                .expect("Unable to convert to string")
                                                .contains(
                                                    &String::from_utf8(insertion.to_vec())
                                                        .expect("Unable to convert to string"),
                                                )
                                            {
                                                base.assign_read(record.qname().to_vec());
                                                base.truedepth[sample_idx] += 1;
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }
                        read_cursor += cig.len() as usize;
                        total_indels_in_current_contig += cig.len() as u64;
                    }
                    Cigar::SoftClip(_) => {
                        // soft clipped portions of reads can be included as structural variants
                        // not sure if this correct protocol or not
                        let insertion =
                            &record.seq().as_bytes()[read_cursor..read_cursor + cig.len() as usize];
                        match variant_matrix.variants(ref_idx, tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Insertion(alt) => {
                                            if String::from_utf8(alt.to_vec())
                                                .expect("Unable to convert to string")
                                                .contains(
                                                    &String::from_utf8(insertion.to_vec())
                                                        .expect("Unable to convert to string"),
                                                )
                                            {
                                                base.assign_read(record.qname().to_vec());
                                                base.truedepth[sample_idx] += 1;
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }
                        read_cursor += cig.len() as usize;
                    }
                    Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
            // Determine the number of mismatching bases in this read by
            // looking at the NM tag.
            total_edit_distance_in_current_contig += match record.aux("NM".as_bytes()) {
                Some(aux) => aux.integer() as u64,
                None => {
                    panic!(
                        "Mapping record encountered that does not have an 'NM' \
                            auxiliary tag in the SAM/BAM format. This is required \
                            to work out some coverage statistics"
                    );
                }
            };
        }
    }
    if last_tid != -2 {
        let contig_len = header
            .target_len(last_tid as u32)
            .expect("Corrupt BAM file?") as usize;
        contig_name = target_names[last_tid as usize].to_vec();

        ref_idx = retrieve_reference_index_from_contig(&contig_name, genomes_and_contigs);

        let total_mismatches =
            total_edit_distance_in_current_contig - total_indels_in_current_contig;

        match reference.fetch_all(std::str::from_utf8(&contig_name[..]).unwrap()) {
            Ok(reference) => reference,
            Err(_) => match reference.fetch_all(&format!(
                "{}~{}",
                &genomes_and_contigs.genomes[ref_idx],
                std::str::from_utf8(&contig_name[..]).unwrap()
            )) {
                Ok(reference) => reference,
                Err(e) => {
                    println!(
                        "Cannot read sequence from reference {} {:?}",
                        format!(
                            "{}~{}",
                            &genomes_and_contigs.genomes[ref_idx],
                            std::str::from_utf8(&contig_name[..]).unwrap()
                        ),
                        e,
                    );
                    std::process::exit(1);
                }
            },
        };
        ref_seq = Vec::new();
        match reference.read(&mut ref_seq) {
            Ok(reference) => reference,
            Err(e) => {
                println!(
                    "Cannot read sequence from reference {} {:?}",
                    std::str::from_utf8(&contig_name[..]).unwrap(),
                    e,
                );
                std::process::exit(1)
            }
        };

        process_previous_contigs_var(
            mode,
            ref_idx,
            last_tid,
            ups_and_downs,
            coverage_estimators,
            min,
            max,
            total_indels_in_current_contig as usize,
            contig_end_exclusion,
            contig_len,
            contig_name,
            &genomes_and_contigs,
            variant_matrix,
            ref_seq,
            sample_idx,
            method,
            total_mismatches,
            coverage_fold,
            num_mapped_reads_in_current_contig,
            sample_count,
            output_prefix,
            &stoit_name,
        );

        num_mapped_reads_total += num_mapped_reads_in_current_contig;
    }

    info!(
        "In sample '{}', found {} reads mapped out of {} total ({:.*}%)", // and filtered {}",
        stoit_name,
        num_mapped_reads_total,
        bam_generated.num_detected_primary_alignments(),
        2,
        (num_mapped_reads_total * 100) as f64
            / bam_generated.num_detected_primary_alignments() as f64,
        // skipped_reads
    );

    if bam_generated.num_detected_primary_alignments() == 0 {
        warn!(
            "No primary alignments were observed for sample {} \
           - perhaps something went wrong in the mapping?",
            stoit_name
        );
    }
    bam_generated.finish();
}

#[allow(unused)]
pub fn process_previous_contigs_var(
    mode: &str,
    ref_idx: usize,
    last_tid: i32,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    min: f32,
    max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u64,
    contig_len: usize,
    contig_name: Vec<u8>,
    genomes_and_contigs: &GenomesAndContigs,
    variant_matrix: &mut VariantMatrix,
    ref_sequence: Vec<u8>,
    sample_idx: usize,
    method: &str,
    total_mismatches: u64,
    coverage_fold: f32,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
    output_prefix: &str,
    stoit_name: &str,
) {
    if last_tid != -2 {
        coverage_estimators
            .par_iter_mut()
            .for_each(|estimator| estimator.setup());

        coverage_estimators.par_iter_mut().for_each(|estimator| {
            estimator.add_contig(
                &ups_and_downs,
                num_mapped_reads_in_current_contig,
                total_mismatches,
            )
        });

        let coverages: Vec<f64> = coverage_estimators
            .iter_mut()
            .map(|estimator| estimator.calculate_coverage(&vec![0]) as f64)
            .collect();

        let mut variant_struct =
            VariantStats::new_contig_stats(min as f64, max as f64, contig_end_exclusion);

        // adds contig info to variant struct
        variant_struct.add_contig(
            variant_matrix.variants_of_contig(ref_idx, last_tid),
            last_tid.clone(),
            total_indels_in_current_contig,
            contig_name.clone(),
            contig_len,
            sample_idx,
            coverages,
            ups_and_downs,
        );

        match mode {
            "polymorph" => {
                // calculates minimum number of genotypes possible for each variant location
                //                variant_struct.generate_minimum_genotypes();
                if variant_struct.len() > 0 {

                    //                    variant_struct.print_variants(&ref_sequence, stoit_name);
                }
            }
            "summarize" | "genotype" | "evolve" | "polish" => {
                // Add samples contig information to main struct
                debug!("Adding in new info for contig...");
                variant_matrix.add_contig(variant_struct, sample_count, sample_idx, ref_idx);
            }
            // "polish" => {
            //     let stoit_name = stoit_name
            //         .split("..").last().unwrap()
            //         .split("/").last().unwrap();
            //     let output_prefix = format!("{}/{}", output_prefix.to_string(), stoit_name);
            //     variant_struct.polish_contig(&ref_sequence,
            //                                  &output_prefix);
            // }
            _ => {
                panic!("unknown mode {}", mode);
            }
        }
    }
}
