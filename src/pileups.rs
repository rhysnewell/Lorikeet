use std;
use std::collections::{HashMap, BTreeMap, BTreeSet};
use rust_htslib::bam::{self, record::Cigar};

use external_command_checker;
use pileup_structs::*;
use pileup_matrix::*;
use codon_structs::*;
use coverm::bam_generator::*;
use rayon::prelude::*;
use crate::estimation::alignment_properties::{InsertSize, AlignmentProperties};

use crate::*;
use std::str;
use std::fs::File;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use bio::io::gff;
use bio::io::gff::Record;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use tempfile;
use std::sync::{Arc, Mutex};


pub fn pileup_variants<R: NamedBamReader + Send,
    G: NamedBamReaderGenerator<R> + Send>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    _print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter,
    mapq_threshold: u8,
    min_var_depth: usize,
    min: f32, max: f32,
    contig_end_exclusion: u64,
    output_prefix: &str,
    n_threads: usize,
    method: &str,
    coverage_fold: f32,
    include_indels: bool,
    include_soft_clipping: bool) {

    let sample_count = bam_readers.len();
    let reference = Arc::new(Mutex::new(reference));
    let coverage_estimators = Arc::new(Mutex::new(coverage_estimators));
    let mut ani = 0.;
    // Print file header
    let mut pileup_matrix = Arc::new(Mutex::new(PileupMatrix::new_matrix(sample_count)));
    let mut gff_map = Arc::new(Mutex::new(HashMap::new()));
    let mut codon_table = CodonTable::setup();
    let nanopore = false;


    match mode {
        "evolve" => {

            codon_table.get_codon_table(11);
            ani = 0.;

            let mut gff_reader;
            if m.is_present("gff") {
                let gff_file = m.value_of("gff").unwrap();
                gff_reader = gff::Reader::from_file(gff_file,
                                                    bio::io::gff::GffType::GFF3)
                    .expect("GFF File not found");
            } else {
                external_command_checker::check_for_prodigal();
                let tmp_dir = TempDir::new("lorikeet_fifo")
                    .expect("Unable to create temporary directory");
                let fifo_path = tmp_dir.path().join("foo.pipe");

                // create new fifo and give read, write and execute rights to the owner.
                // This is required because we cannot open a Rust stream as a BAM file with
                // rust-htslib.
                unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
                    .expect(&format!("Error creating named pipe {:?}", fifo_path));

                let gff_file = tempfile::Builder::new()
                    .prefix("lorikeet-prodigal-gff")
                    .tempfile_in(tmp_dir.path())
                    .expect(&format!("Failed to create distances tempfile"));
                let cmd_string = format!(
                    "set -e -o pipefail; \
                     prodigal -f gff -i {} -o {} {}",
                    // prodigal
                    m.value_of("reference").unwrap(),
                    gff_file.path().to_str()
                        .expect("Failed to convert tempfile path to str"),
                    m.value_of("prodigal-params").unwrap_or(""));
                info!("Queuing cmd_string: {}", cmd_string);
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&cmd_string)
                    .output()
                    .expect("Unable to execute bash");

                gff_reader = gff::Reader::from_file(gff_file.path().to_str()
                                                        .expect("Failed to convert tempfile path to str"),
                                                    bio::io::gff::GffType::GFF3)
                    .expect("Failed to read prodigal output");

                tmp_dir.close().expect("Failed to close temp directory");

            }
            gff_reader.records().into_iter().for_each(|record| {
                let rec = record.unwrap();
                let mut gff_map = gff_map.lock().unwrap();
                let contig_genes = gff_map.entry(rec.seqname().to_owned())
                    .or_insert(Vec::new());
                contig_genes.push(rec);
            });
        },
        "genotype" | "summarize" => {
            if m.is_present("strain-ani") {
                ani = parse_percentage(m, "strain-ani");
            }
        },
        _ => {
//            min_cluster_size = m.value_of("min-cluster-size").unwrap().parse().unwrap();
//            epsilon = m.value_of("epsilon").unwrap().parse().unwrap();
        }
    }
    let mut read_cnt_id: Arc<Mutex<i64>> = Arc::new(Mutex::new(0));
    let mut read_to_id = Arc::new(Mutex::new(HashMap::new()));
    // Loop through bam generators in parallel
    bam_readers.into_par_iter().enumerate().for_each(|(sample_idx, bam_generator)|{
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(n_threads / sample_count);

        let bam_properties =
            AlignmentProperties::default(InsertSize::default());

        let stoit_name = bam_generated.name().to_string();

        let header = bam_generated.header().clone(); // bam header
        let target_names = header.target_names(); // contig names
        let mut record: bam::record::Record = bam::Record::new();
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u64 = 0;

        let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        let mut nuc_freq = HashMap::new();
        let mut indels = HashMap::new();

        let mut last_tid: i32 = -2; // no such tid in a real BAM file
        let mut total_indels_in_current_contig = 0;


        // for record in records
        let mut skipped_reads = 0;

        while bam_generated.read(&mut record)
            .expect("Error while reading BAM record") == true {

            if (!flag_filters.include_supplementary && record.is_supplementary()) ||
                (!flag_filters.include_secondary && record.is_secondary()) ||
                (!flag_filters.include_improper_pairs && !record.is_proper_pair()){
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
                // Check if new read to id
                if !read_to_id.lock().unwrap().contains_key(&record.qname().to_vec()) {
                    let mut read_to_id = read_to_id.lock().unwrap();
                    let mut read_cnt_id = read_cnt_id.lock().unwrap();
                    read_to_id.entry(record.qname().to_vec())
                              .or_insert(*read_cnt_id);
                    *read_cnt_id += 1;
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
                            mode,
                            ani,
                            last_tid,
                            nuc_freq,
                            indels,
                            ups_and_downs,
                            &coverage_estimators,
                            min, max,
                            total_indels_in_current_contig as usize,
                            contig_end_exclusion,
                            min_var_depth,
                            contig_len,
                            contig_name,
                            &pileup_matrix,
                            ref_seq,
                            sample_idx,
                            method,
                            total_mismatches,
                            &gff_map,
                            &codon_table,
                            coverage_fold,
                            num_mapped_reads_in_current_contig,
                            sample_count,
                            output_prefix,
                            &stoit_name);
                    }
                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    debug!("Working on new reference {}",
                           std::str::from_utf8(target_names[tid as usize]).unwrap());
                    last_tid = tid;
                    num_mapped_reads_total += num_mapped_reads_in_current_contig;
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    nuc_freq = HashMap::new();
                    indels = HashMap::new();

                    let mut reference = reference.lock().unwrap();
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
                let mut cursor: usize = record.pos() as usize;
                let _quals = record.qual();
                let mut read_cursor: usize = 0;
                for cig in record.cigar().iter() {
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or = increment start and decrement end index
                            ups_and_downs[cursor] += 1;
                            let final_pos = cursor + cig.len() as usize;
                            if !nanopore {
                                for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                                    let base = record.seq()[qpos] as char;
                                    let refr = ref_seq[cursor as usize] as char;
//                                let mut nuc_freq = nuc_freq.lock().unwrap();
                                    let nuc_map = nuc_freq
                                        .entry(cursor as i32).or_insert(BTreeMap::new());

                                    if base != refr {
//                                    let nuc_freq = Arc::clone(nuc_freq.lock().unwrap());
                                        let id = nuc_map.entry(base).or_insert(BTreeSet::new());
                                        let mut read_to_id = read_to_id.lock().unwrap();
                                        id.insert(read_to_id[&record.qname().to_vec()]);
                                    } else {
                                        let id = nuc_map.entry("R".as_bytes()[0] as char).or_insert(BTreeSet::new());
                                        let mut read_to_id = read_to_id.lock().unwrap();
                                        id.insert(read_to_id[&record.qname().to_vec()]);
                                    }
                                    cursor += 1;
                                }
                            } else {
                                cursor = final_pos;
                            }


                            if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                                ups_and_downs[final_pos] -= 1;
                            }
                            read_cursor += cig.len() as usize;
                        },
                        Cigar::Del(_) => {
                            if include_indels {
                                let refr = (ref_seq[cursor as usize] as char).to_string();
                                let insert = refr +
                                    &std::iter::repeat("N").take(cig.len() as usize).collect::<String>();
                                let refr = str::from_utf8(&ref_seq[cursor as usize..
                                    cursor as usize + cig.len() as usize]).unwrap().to_string();
                                if refr != insert {
                                    let indel_map = indels
                                        .entry(cursor as i32).or_insert(BTreeMap::new());
                                    let id = indel_map.entry(insert)
                                        .or_insert(BTreeSet::new());
                                    let mut read_to_id = read_to_id.lock().unwrap();
                                    id.insert(read_to_id[&record.qname().to_vec()]);
                                    total_indels_in_current_contig += cig.len() as u64;
                                }
                            }

                            cursor += cig.len() as usize;

                        },
                        Cigar::RefSkip(_) => {
                            // if D or N, move the cursor
                            cursor += cig.len() as usize;
                        },
                        Cigar::Ins(_) => {
                            if include_indels {
                                let refr = (ref_seq[cursor as usize] as char).to_string();
                                let insert = match str::from_utf8(&record.seq().as_bytes()[
                                    read_cursor..read_cursor + cig.len() as usize]) {
                                    Ok(ins) => { ins.to_string() },
                                    Err(_e) => { "".to_string() },
                                };

                                let indel_map = indels.entry(cursor as i32)
                                    .or_insert(BTreeMap::new());

                                let id = indel_map.entry(refr + &insert)
                                    .or_insert(BTreeSet::new());
                                let mut read_to_id = read_to_id.lock().unwrap();
                                id.insert(read_to_id[&record.qname().to_vec()]);
                            }
                            read_cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len() as u64;
                        },
                        Cigar::SoftClip(_) => {
                            // soft clipped portions of reads can be included as insertions
                            // not sure if this correct protocol or not
                            if include_soft_clipping {
                                let refr = (ref_seq[cursor as usize] as char).to_string();
                                let insert = match str::from_utf8(&record.seq().as_bytes()[
                                    read_cursor..read_cursor + cig.len() as usize]) {
                                    Ok(ins) => {ins.to_string()},
                                    Err(_e) => {"".to_string()},
                                };
                                let indel_map = indels.entry(cursor as i32)
                                    .or_insert(BTreeMap::new());

                                let id = indel_map.entry(refr + &insert)
                                    .or_insert(BTreeSet::new());
                                let mut read_to_id = read_to_id.lock().unwrap();
                                id.insert(read_to_id[&record.qname().to_vec()]);
                                total_indels_in_current_contig += cig.len() as u64;
                            }
                            read_cursor += cig.len() as usize;
                        },
                        Cigar::HardClip(_) | Cigar::Pad(_) => {
                        }
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
        } if last_tid != -2 {
            let contig_len = header.target_len(last_tid as u32).expect("Corrupt BAM file?") as usize;
            let contig_name = target_names[last_tid as usize].to_vec();
            let total_mismatches = total_edit_distance_in_current_contig -
                total_indels_in_current_contig;

            process_previous_contigs_var(
                mode,
                ani,
                last_tid,
                nuc_freq,
                indels,
                ups_and_downs,
                &coverage_estimators,
                min, max,
                total_indels_in_current_contig as usize,
                contig_end_exclusion,
                min_var_depth,
                contig_len,
                contig_name,
                &pileup_matrix,
                ref_seq,
                sample_idx,
                method,
                total_mismatches,
                &gff_map,
                &codon_table,
                coverage_fold,
                num_mapped_reads_in_current_contig,
                sample_count,
                output_prefix,
                &stoit_name);

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
        let mut pileup_matrix = pileup_matrix.lock().unwrap();
        pileup_matrix.add_sample(stoit_name, sample_idx);
    });
    if mode=="genotype" {
        let mut pileup_matrix = pileup_matrix.lock().unwrap();
        pileup_matrix.generate_distances(n_threads, output_prefix);
        let e_min: f64 = m.value_of("e-min").unwrap().parse().unwrap();
        let e_max: f64 = m.value_of("e-max").unwrap().parse().unwrap();
        let pts_min: f64 = m.value_of("pts-min").unwrap().parse().unwrap();
        let pts_max: f64 = m.value_of("pts-max").unwrap().parse().unwrap();
        let phi: f64 = m.value_of("phi").unwrap().parse().unwrap();

        pileup_matrix.run_fuzzy_scan(e_min, e_max, pts_min, pts_max, phi);
        pileup_matrix.generate_genotypes(output_prefix);
    } else if mode=="summarize" {
        let mut pileup_matrix = pileup_matrix.lock().unwrap();
        pileup_matrix.print_variant_stats(output_prefix);
    }
}

fn process_previous_contigs_var(
    mode: &str,
    ani: f32,
    last_tid: i32,
    nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i64>>>,
    indels: HashMap<i32, BTreeMap<String, BTreeSet<i64>>>,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &Arc<Mutex<&mut Vec<CoverageEstimator>>>,
    min: f32, max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u64,
    min_var_depth: usize,
    contig_len: usize,
    contig_name: Vec<u8>,
    pileup_matrix: &Arc<Mutex<PileupMatrix>>,
    ref_sequence: Vec<u8>,
    sample_idx: usize,
    method: &str,
    total_mismatches: u64,
    gff_map: &Arc<Mutex<HashMap<String, Vec<Record>>>>,
    codon_table: &CodonTable,
    coverage_fold: f32,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
    output_prefix: &str,
    stoit_name: &str) {

    if last_tid != -2 {
        let mut coverage_estimators = coverage_estimators.lock().unwrap();
        coverage_estimators.par_iter_mut().for_each(|estimator|{
            estimator.setup()
        });

        coverage_estimators.par_iter_mut().for_each(|estimator|{
            estimator.add_contig(
                &ups_and_downs,
                num_mapped_reads_in_current_contig,
                total_mismatches)
        });

        let coverages: Vec<f64> = coverage_estimators.par_iter_mut()
            .map(|estimator| estimator.calculate_coverage(&vec![0]) as f64).collect();

        let mut pileup_struct = PileupStats::new_contig_stats(min as f64,
                                                              max as f64,
                                                              contig_end_exclusion);

        // adds contig info to pileup struct
        pileup_struct.add_contig(nuc_freq,
                                 indels,
                                 last_tid.clone(),
                                 total_indels_in_current_contig,
                                 contig_name.clone(),
                                 contig_len,
                                 method,
                                 coverages,
                                 ups_and_downs);

        if ani == 0. {
            pileup_struct.calc_error(ani);


            // filters variants across contig
            pileup_struct.calc_variants(
                min_var_depth,
                coverage_fold as f64);
        } else {
            let min_var_depth = pileup_struct.calc_error(ani);

            info!("Minimum Variant Depth set to {} for strain ANI of {}", min_var_depth, ani);

            // filters variants across contig
            pileup_struct.calc_variants(
                min_var_depth,
                coverage_fold as f64);
        }


        match mode {
            "polymorph" => {

                // calculates minimum number of genotypes possible for each variant location
//                pileup_struct.generate_minimum_genotypes();
                if pileup_struct.len() > 0 {

                    pileup_struct.print_variants(&ref_sequence, stoit_name);
                }

            },
            "summarize" | "genotype" => {
                let mut pileup_matrix = pileup_matrix.lock().unwrap();
                // calculates minimum number of genotypes possible for each variant location
                pileup_matrix.add_contig(pileup_struct,
                                         sample_count,
                                         sample_idx,
                                        ref_sequence);
            },
            "evolve" => {
                let gff_map = gff_map.lock().unwrap();
                pileup_struct.calc_gene_mutations(&*gff_map, &ref_sequence, codon_table);
            },
            "polish" => {
                let stoit_name = stoit_name
                    .split("..").last().unwrap()
                    .split("/").last().unwrap();
                let output_prefix = output_prefix.to_string() + "_" + stoit_name;
                pileup_struct.polish_contig(&ref_sequence,
                                            &output_prefix);
            }
            _ => {panic!("unknown mode {}", mode);},
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::str;
    use std::fs::File;
    use coverm::mapping_parameters::*;
    use coverm::shard_bam_reader::*;
    use coverm::genome_exclusion::*;
    use coverm::bam_generator::*;
//
//    fn test_with_stream<R: NamedBamReader + Send,
//        G: NamedBamReaderGenerator<R> + Send>(
//        expected: &str,
//        bam_readers: Vec<G>,
//        mut reference: bio::io::fasta::IndexedReader<File>,
//        proper_pairs_only: bool,
//        n_threads: usize,
//        coverage_fold: f32,
//        min_var_depth: usize,
//        min: f32,
//        max: f32,
//        mode: &str,
//        include_indels: bool,
//        include_soft_clipping: bool) {
////        let mut stream = Cursor::new(Vec::new());
//        {
//            reads_mapped_vec = pileup_variants(
//                bam_readers,
//                &mut coverage_taker,
//                coverage_estimators,
//                print_zero_coverage_contigs,
//                flag_filters,
//                false,
//            );
//        }
////        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
//    }
}