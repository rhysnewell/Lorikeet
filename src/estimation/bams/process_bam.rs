use std;
use std::collections::{HashMap, HashSet};
use rust_htslib::bam::{self, record::Cigar};

use estimation::contig_variants::*;
use estimation::variant_matrix::*;
use estimation::codon_structs::*;
use coverm::bam_generator::*;
use rayon::prelude::*;
use model::variants::*;
use utils::*;

use std::str;
use std::path::Path;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::genomes_and_contigs::*;
use coverm::FlagFilter;
use bio::io::gff::Record;


/// Process all reads in a BAM file
#[allow(unused)]
pub fn process_bam<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_generator: G,
    sample_count: usize,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    variant_matrix_map: &mut HashMap<usize, VariantMatrix>,
    gff_map: &mut HashMap<usize, HashMap<String, Vec<Record>>>,
    split_threads: usize,
    m: &clap::ArgMatches,
    output_prefix: &str,
    coverage_fold: f32,
    codon_table: &CodonTable,
    min_var_depth: usize,
    contig_end_exclusion: u64,
    min: f32, max: f32,
    ani: f32,
    mode: &str,
    include_soft_clipping: bool,
    include_indels: bool,
    flag_filters: &FlagFilter,
    mapq_threshold: u8,
    method: &str,
    sample_groups: &HashMap<&str, HashSet<String>>,
    genomes_and_contigs: &GenomesAndContigs,
    reference_map: &HashMap<usize, String>) {

    let mut bam_generated = bam_generator.start();

    // Adjust the sample index if the bam is from long reads
    let longread;


//    let mut bam_properties =
//        AlignmentProperties::default(InsertSize::default());

    let stoit_name = bam_generated.name().to_string().replace("/", ".");
    // adjust sample index for longread bams


    bam_generated.set_threads(split_threads);
    debug!("Managed to set threads.");

    let header = bam_generated.header().clone(); // bam header
    let target_names = header.target_names(); // contig names
    let mut record: bam::record::Record = bam::Record::new();
    let mut ups_and_downs: Vec<i32> = Vec::new();
    let mut num_mapped_reads_total: u64 = 0;
    let mut num_mapped_reads_in_current_contig: u64 = 0;
    let mut total_edit_distance_in_current_contig: u64 = 0;

    let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig
    let mut last_tid: i32 = -2; // no such tid in a real BAM file
    let mut total_indels_in_current_contig = 0;

    let reference_stem = genomes_and_contigs.genome_of_contig(
        &str::from_utf8(&target_names[0]).unwrap().to_string()).unwrap();
    let ref_idx = genomes_and_contigs.genome_index(&reference_stem).unwrap();
    let reference_path = reference_map.get(&ref_idx).expect("Unable to retrieve reference path");

    let gff_map = gff_map.entry(ref_idx).or_insert(HashMap::new());
    let mut variant_matrix = variant_matrix_map.entry(ref_idx)
        .or_insert(VariantMatrix::new_matrix(sample_count));
    let mut reference = match bio::io::fasta::IndexedReader::from_file(&Path::new(&reference_path)) {
        Ok(reader) => reader,
        Err(_e) => generate_faidx(&reference_path),
    };

    let sample_idx = match variant_matrix {
        VariantMatrix::VariantContigMatrix {
            sample_names,
            ..
        } => {
            sample_names.iter().position(|p| {p == &stoit_name}).unwrap()
        }
    };

    if sample_groups.contains_key("long") {
        if sample_groups["long"].contains(&stoit_name) {
            longread = true;
            debug!("Longread {} {}", stoit_name, sample_idx)
        } else {
            longread = false;
            debug!(" SOME Longread {} {}", stoit_name, sample_idx)

        }
    } else {
        longread = false;
        debug!(" NO Longread {} {}", stoit_name, sample_idx)

    }


    // for record in records
    let mut skipped_reads = 0;

    while bam_generated.read(&mut record)
        .expect("Error while reading BAM record") == true {
        if (!flag_filters.include_supplementary && record.is_supplementary() && !longread) ||
            (!flag_filters.include_secondary && record.is_secondary() && !longread) ||
            (!flag_filters.include_improper_pairs && !record.is_proper_pair() && !longread) {
            skipped_reads += 1;
            continue;
        } else if !flag_filters.include_secondary && record.is_secondary() && longread {
            skipped_reads += 1;
            continue;
        }
        // if reference has changed, print the last record
        let tid = record.tid();
        if !record.is_unmapped() { // if mapped
            if record.seq().len() == 0 {
                continue
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
                    let contig_len = header.target_len(last_tid as u32)
                        .expect("Corrupt BAM file?") as usize;
                    let contig_name = target_names[last_tid as usize].to_vec();
                    let total_mismatches = total_edit_distance_in_current_contig -
                        total_indels_in_current_contig;

                    match reference.fetch_all(std::str::from_utf8(target_names[last_tid as usize]).unwrap()) {
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

                    process_previous_contigs_var(
                        mode,
                        ani,
                        last_tid,
                        ups_and_downs,
                        coverage_estimators,
                        min, max,
                        total_indels_in_current_contig as usize,
                        contig_end_exclusion,
                        min_var_depth,
                        contig_len,
                        contig_name,
                        variant_matrix,
                        ref_seq,
                        sample_idx,
                        method,
                        total_mismatches,
                        gff_map,
                        &codon_table,
                        coverage_fold,
                        num_mapped_reads_in_current_contig,
                        sample_count,
                        output_prefix,
                        &stoit_name,
                        longread);
                }

                ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                debug!("Working on new reference {}",
                       std::str::from_utf8(target_names[tid as usize]).unwrap());
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
                        let mut mnv = vec!();
                        let mut mnv_cursor = 0;
                        for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                            // See if read is match MNV
                            if potential_mnv {
                                let read_char = record.seq()[qpos];
                                debug!("MNV searching {} {:?} {}", &mnv_pos, &mnv, &read_char);
                                if mnv[mnv_pos] == read_char {
                                    mnv_pos += 1;
                                    debug!("pos {} length {}", &mnv_pos, &mnv.len());
                                    if mnv_pos == mnv.len() {
                                        match variant_matrix.variants(tid, mnv_cursor) {

                                            Some(current_variants) => {
                                                current_variants.iter_mut().for_each(|(variant, base)| {
                                                    match variant {
                                                        Variant::MNV(alt) => {
                                                            debug!("alt {:?} found {:?}", &alt, &mnv);

                                                            if *alt == mnv {
                                                                base.assign_read(record.qname().to_vec());
                                                                base.truedepth[sample_idx] += 1;
                                                                mnv = vec!();
                                                                mnv_pos = 0;
                                                                potential_mnv = false;
                                                            } else {
                                                                mnv = vec!();
                                                                mnv_pos = 0;
                                                                potential_mnv = false
                                                            }
                                                        },
                                                        _ => {
                                                            debug!("found {:?} but not confirmed wrongly", &mnv);
                                                        },
                                                    }
                                                })
                                            },
                                            None => {
                                                mnv = vec!();
                                                mnv_pos = 0;
                                                potential_mnv = false
                                            },
                                        }
                                    }
                                } else {
                                    mnv = vec!();
                                    mnv_pos = 0;
                                    potential_mnv = false
                                }
                            }
                            match variant_matrix.variants(tid, cursor as i64) {
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
                                                if alt[mnv_pos] == read_char {

                                                    mnv = alt.clone();
                                                    debug!("Potential MNV {} {:?} {}", &mnv_pos, &mnv, &read_char);

                                                    mnv_pos += 1;
                                                    potential_mnv = true;
                                                    mnv_cursor = cursor as i64;
                                                }
                                            },
                                            Variant::None => {
                                                if base.refr[0] == read_char {
                                                    base.assign_read(record.qname().to_vec());
                                                    base.truedepth[sample_idx] += 1;
                                                } else {
                                                    mnv = vec!();
                                                    mnv_pos = 0;
                                                    potential_mnv = false
                                                }
                                            },
                                            _ => {}
                                        }
                                    });
                                },
                                _ => {},
                            }
                            cursor += 1;
                        }
                        if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        read_cursor += cig.len() as usize;
                    },
                    Cigar::Del(del) => {
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Deletion(alt) => {
                                            if alt == del {
                                                base.assign_read(record.qname().to_vec());
                                                base.truedepth[sample_idx] += 1;
                                            }
                                        },
                                        _ => {}
                                    }
                                });
                            },
                            _ => {},
                        }

                        cursor += cig.len() as usize;
                    },
                    Cigar::RefSkip(_) => {
                        // if D or N, move the cursor
                        cursor += cig.len() as usize;
                    },
                    Cigar::Ins(ins) => {
                        let insertion = &record.seq().as_bytes()[read_cursor..read_cursor+cig.len() as usize];
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Insertion(alt) => {
                                            if String::from_utf8(alt.to_vec()).expect("Unable to convert to string")
                                                .contains(&String::from_utf8(insertion.to_vec()).expect("Unable to convert to string")) {
                                                base.assign_read(record.qname().to_vec());
                                                base.truedepth[sample_idx] += 1;
                                            }
                                        },
                                        _ => {}
                                    }
                                });
                            },
                            _ => {},
                        }
                        read_cursor += cig.len() as usize;
                        total_indels_in_current_contig += cig.len() as u64;
                    },
                    Cigar::SoftClip(_) => {
                        // soft clipped portions of reads can be included as structural variants
                        // not sure if this correct protocol or not
//                        if longreads {
//                            let refr = (ref_seq[cursor as usize] as char).to_string();
//                            let insert = match str::from_utf8(&record.seq().as_bytes()[
//                                read_cursor..read_cursor + cig.len() as usize]) {
//                                Ok(ins) => {ins.to_string()},
//                                Err(_e) => {"".to_string()},
//                            };
//                            let indel_map = indels.entry(cursor as i32)
//                                .or_insert(BTreeMap::new());
//
//                            let id = indel_map.entry(refr + &insert)
//                                .or_insert(BTreeSet::new());
//                            let mut read_to_id = read_to_id.lock().unwrap();
//                            id.insert(read_to_id[&record.qname().to_vec()]);
//                            total_indels_in_current_contig += cig.len() as u64;
//                        }
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
                    aux.integer() as u64
                },
                None => {
                    panic!("Mapping record encountered that does not have an 'NM' \
                            auxiliary tag in the SAM/BAM format. This is required \
                            to work out some coverage statistics");
                }
            };
        }
    }
    if last_tid != -2 {
        let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
        let contig_name = target_names[last_tid as usize].to_vec();
        let total_mismatches = total_edit_distance_in_current_contig -
            total_indels_in_current_contig;
        match reference.fetch_all(std::str::from_utf8(target_names[last_tid as usize]).unwrap()) {
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

        process_previous_contigs_var(
            mode,
            ani,
            last_tid,
            ups_and_downs,
            coverage_estimators,
            min, max,
            total_indels_in_current_contig as usize,
            contig_end_exclusion,
            min_var_depth,
            contig_len,
            contig_name,
            variant_matrix,
            ref_seq,
            sample_idx,
            method,
            total_mismatches,
            gff_map,
            &codon_table,
            coverage_fold,
            num_mapped_reads_in_current_contig,
            sample_count,
            output_prefix,
            &stoit_name,
            longread);

        num_mapped_reads_total += num_mapped_reads_in_current_contig;
    }


    info!("In sample '{}', found {} reads mapped out of {} total ({:.*}%) and filtered {}",
          stoit_name, num_mapped_reads_total,
          bam_generated.num_detected_primary_alignments(), 2,
          (num_mapped_reads_total * 100) as f64 /
              bam_generated.num_detected_primary_alignments() as f64, skipped_reads);


    if bam_generated.num_detected_primary_alignments() == 0 {
        warn!("No primary alignments were observed for sample {} \
           - perhaps something went wrong in the mapping?",
              stoit_name);
    }
    bam_generated.finish();
}

#[allow(unused)]
pub fn process_previous_contigs_var(
    mode: &str,
    ani: f32,
    last_tid: i32,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    min: f32, max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u64,
    min_var_depth: usize,
    contig_len: usize,
    contig_name: Vec<u8>,
    variant_matrix: &mut VariantMatrix,
    ref_sequence: Vec<u8>,
    sample_idx: usize,
    method: &str,
    total_mismatches: u64,
    gff_map: &mut HashMap<String, Vec<Record>>,
    codon_table: &CodonTable,
    coverage_fold: f32,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
    output_prefix: &str,
    stoit_name: &str,
    longread: bool,) {

    if last_tid != -2 {

        coverage_estimators.par_iter_mut().for_each(|estimator|{
            estimator.setup()
        });

        coverage_estimators.par_iter_mut().for_each(|estimator|{
            estimator.add_contig(
                &ups_and_downs,
                num_mapped_reads_in_current_contig,
                total_mismatches)
        });

        let coverages: Vec<f64> = coverage_estimators.iter_mut()
            .map(|estimator| estimator.calculate_coverage(&vec![0]) as f64).collect();

        let mut variant_struct = VariantStats::new_contig_stats(min as f64,
                                                                max as f64,
                                                                contig_end_exclusion);


        // adds contig info to variant struct
        variant_struct.add_contig(variant_matrix.variants_of_contig(last_tid),
                                  last_tid.clone(),
                                  total_indels_in_current_contig,
                                  contig_name.clone(),
                                  contig_len,
                                  sample_idx,
                                  coverages,
                                  ups_and_downs);

        if ani == 0. {
//            variant_struct.calc_error(ani);


            // filters variants across contig
//            variant_struct.calc_variants(
//                min_var_depth,
//                coverage_fold as f64);
        } else {
//            let min_var_depth = variant_struct.calc_error(ani);

            info!("Minimum Variant Depth set to {} for strain ANI of {}", min_var_depth, ani);

            // filters variants across contig
//            variant_struct.calc_variants(
//                min_var_depth,
//                coverage_fold as f64);
        }


        match mode {
            "polymorph" => {

                // calculates minimum number of genotypes possible for each variant location
//                variant_struct.generate_minimum_genotypes();
                if variant_struct.len() > 0 {

//                    variant_struct.print_variants(&ref_sequence, stoit_name);
                }

            },
            "summarize" | "genotype" | "evolve" => {
                // Add samples contig information to main struct
                debug!("Adding in new info for contig...");
                variant_matrix.add_contig(variant_struct,
                                          sample_count,
                                          sample_idx);
            },
            "polish" => {
                let stoit_name = stoit_name
                    .split("..").last().unwrap()
                    .split("/").last().unwrap();
                let output_prefix = format!("{}/{}", output_prefix.to_string(), stoit_name);
                variant_struct.polish_contig(&ref_sequence,
                                             &output_prefix);
            }
            _ => {panic!("unknown mode {}", mode);},
        }
    }
}
