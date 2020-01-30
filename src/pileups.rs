use std;
use std::collections::{HashMap, BTreeMap, BTreeSet};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use external_command_checker;
use pileup_structs::*;
use pileup_matrix::*;
use codon_structs::*;
use bam_generator::*;
use FlagFilter;
use rayon::prelude::*;

use std::str;
use std::fs::File;
use mosdepth_genome_coverage_estimators::*;
use bio::io::gff;
use bio::io::gff::Record;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use tempfile;
use std::sync::{Arc, Mutex};



pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    mut reference: bio::io::fasta::IndexedReader<File>,
    print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter,
    mapq_threshold: u8,
    min_var_depth: usize,
    min: f32, max: f32,
    contig_end_exclusion: u32,
    output_prefix: &str,
    n_threads: usize,
    method: &str,
    coverage_fold: f32,
    include_indels: bool,
    include_soft_clipping: bool) {

    let sample_idx = 0;
    let sample_count = bam_readers.len();
    let mut sample_idx = 0;
    // Print file header
    let mut pileup_matrix = PileupMatrix::new_matrix();

    let mut gff_map = HashMap::new();
    let mut epsilon = 0.05;
    let mut min_cluster_size= 10;

    match mode {
        "evolve" => {
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

                let mut gff_file = tempfile::Builder::new()
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

                tmp_dir.close().expect("Failed to close temo directory");

            }
            for record in gff_reader.records() {
                let rec = record.unwrap();
                let contig_genes = gff_map.entry(rec.seqname().to_owned())
                    .or_insert(Vec::new());
                contig_genes.push(rec);
            }
        },
        _ => {
            min_cluster_size = m.value_of("min-cluster-size").unwrap().parse().unwrap();
            epsilon = m.value_of("epsilon").unwrap().parse().unwrap();
        }
    }
    let mut read_cnt_id: i64 = 0;
    let mut read_to_id = HashMap::new();
    let mut codon_table = CodonTable::setup();
    codon_table.get_codon_table(11);
    // Loop through bam generators in parallel
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(n_threads);
        let stoit_name = bam_generated.name().to_string();

        let header = bam_generated.header().clone(); // bam header
        let target_names = header.target_names(); // contig names
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u32 = 0;

        let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig

        // for each genomic position, only has hashmap when variants are present. Includes read ids
        let mut nuc_freq: Arc<Mutex<HashMap<i32, BTreeMap<char, BTreeSet<i64>>>>> = Arc::new(Mutex::new(HashMap::new()));
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
                            mode,
                            last_tid,
                            nuc_freq.lock().unwrap().clone(),
                            indels,
                            ups_and_downs,
                            coverage_estimators,
                            min, max,
                            total_indels_in_current_contig as usize,
                            contig_end_exclusion,
                            min_var_depth,
                            contig_len,
                            contig_name,
                            &mut pileup_matrix,
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
                            epsilon,
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
                    nuc_freq = Arc::new(Mutex::new(HashMap::new()));
                    indels = HashMap::new();

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
                let quals = record.qual();
                let mut read_cursor: usize = 0;
                for cig in record.cigar().iter() {
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or = increment start and decrement end index
                            ups_and_downs[cursor] += 1;
                            let final_pos = cursor + cig.len() as usize;

                            for qpos in read_cursor..(read_cursor+cig.len() as usize) {
                                let base = record.seq()[qpos] as char;
                                let refr = ref_seq[cursor as usize] as char;
                                let mut nuc_freq = nuc_freq.lock().unwrap();
                                let nuc_map = nuc_freq
                                    .entry(cursor as i32).or_insert(BTreeMap::new());

                                if base != refr {
//                                    let nuc_freq = Arc::clone(nuc_freq.lock().unwrap());

                                    let id = nuc_map.entry(base)
                                        .or_insert(BTreeSet::new());
                                    let id = nuc_map.entry(base).or_insert(BTreeSet::new());
                                    id.insert(read_to_id[&record.qname().to_vec()]);
                                }
                                cursor += 1;
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
                                    id.insert(read_to_id[&record.qname().to_vec()]);
                                    total_indels_in_current_contig += cig.len();
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
                                id.insert(read_to_id[&record.qname().to_vec()]);
                            }
                            read_cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len();
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

                                id.insert(read_to_id[&record.qname().to_vec()]);
                                total_indels_in_current_contig += cig.len();
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
                        aux.integer() as u32
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
                last_tid,
                nuc_freq.lock().unwrap().clone(),
                indels,
                ups_and_downs,
                coverage_estimators,
                min, max,
                total_indels_in_current_contig as usize,
                contig_end_exclusion,
                min_var_depth,
                contig_len,
                contig_name,
                &mut pileup_matrix,
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
                epsilon,
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
        pileup_matrix.add_sample(stoit_name);

        sample_idx += 1;
    };
    if mode=="genotype" {
//        pileup_matrix.dbscan_cluster(epsilon, min_cluster_size);
//        pileup_matrix.generate_genotypes(output_prefix);
        pileup_matrix.generate_distances(n_threads, output_prefix);
//        pileup_matrix.print_matrix();
    } else if mode=="summarize" {
        pileup_matrix.print_variant_stats(output_prefix);
    }
}

fn process_previous_contigs_var(
    mode: &str,
    last_tid: i32,
    nuc_freq: HashMap<i32, BTreeMap<char, BTreeSet<i64>>>,
    indels: HashMap<i32, BTreeMap<String, BTreeSet<i64>>>,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    min: f32, max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u32,
    min_var_depth: usize,
    contig_len: usize,
    contig_name: Vec<u8>,
    pileup_matrix: &mut PileupMatrix,
    ref_sequence: Vec<u8>,
    sample_idx: i32,
    method: &str,
    total_mismatches: u32,
    gff_map: &HashMap<String, Vec<Record>>,
    codon_table: &CodonTable,
    coverage_fold: f32,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
    output_prefix: &str,
    epsilon: f64,
    stoit_name: &str) {

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

        let coverages: Vec<f32> = coverage_estimators.par_iter_mut()
            .map(|estimator| estimator.calculate_coverage(&vec![0])).collect();

        let mut pileup_struct = PileupStats::new_contig_stats(min,
                                                              max,
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

        pileup_struct.calc_error();


        // filters variants across contig
        pileup_struct.calc_variants(
            min_var_depth,
            coverage_fold);

        match mode {
            "polymorph" => {
                // calculates minimum number of genotypes possible for each variant location
//                pileup_struct.generate_minimum_genotypes();
                if pileup_struct.len() > 0 {
//                    pileup_struct.cluster_variants(epsilon);


                    // prints results of variants calling
//                    pileup_struct.generate_variant_contig(&ref_sequence, output_prefix);

                    pileup_struct.print_variants(&ref_sequence, stoit_name);
                }
//                pileup_struct.generate_svd();

            },
            "summarize" | "genotype" => {
                // calculates minimum number of genotypes possible for each variant location
//                pileup_struct.generate_minimum_genotypes();
                pileup_matrix.add_contig(pileup_struct,
                                         sample_count,
                                         sample_idx as usize,
                                        ref_sequence);
            },
            "evolve" => {
                pileup_struct.calc_gene_mutations(gff_map, &ref_sequence, codon_table);
            },
            "polish" => {
                pileup_struct.polish_contig(&ref_sequence,
                                            output_prefix);
            }
            _ => {panic!("unknown mode {}", mode);},
        }
    }
}