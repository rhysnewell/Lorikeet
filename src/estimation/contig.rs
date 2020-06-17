use std;
use std::collections::{HashMap};
use rust_htslib::bam::{self, record::Cigar};
use rust_htslib::{bcf::Reader, bcf::*};
use bird_tool_utils::command;
use glob::glob;

use external_command_checker;
use estimation::contig_variants::*;
use estimation::variant_matrix::*;
use estimation::codon_structs::*;
use coverm::bam_generator::*;
use rayon::prelude::*;
use estimation::alignment_properties::{InsertSize, AlignmentProperties};
use model::variants::*;

use crate::*;
use std::str;
use std::fs::File;
use std::path::Path;
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
    G: NamedBamReaderGenerator<R> + Send,
    S: NamedBamReader + Send,
    U: NamedBamReaderGenerator<S> + Send,>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    long_readers: Option<Vec<U>>,
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
    include_soft_clipping: bool,
    is_long_read: bool) {

    let mut sample_count = bam_readers.len();
    let reference = Arc::new(Mutex::new(reference));
    let coverage_estimators = Arc::new(Mutex::new(coverage_estimators));
    let mut ani = 0.;


    let gff_map = Arc::new(Mutex::new(HashMap::new()));
    let mut codon_table = CodonTable::setup();

    // Get long reads bams if they exist
    let mut longreads = match long_readers {
        Some(vec) => {
            // update sample count
            sample_count += vec.len();
            vec
        },
        _ => vec!(),
    };

    let variant_matrix = Arc::new(Mutex::new(VariantMatrix::new_matrix(sample_count)));

    info!("{} Longread BAM files and {} Shortread BAM files {} Total BAMs",
          longreads.len(), bam_readers.len(), sample_count);

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
                external_command_checker::check_for_prokka();
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
                     prokka -f gff -i {} -o {} {}",
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
//    let read_cnt_id: Arc<Mutex<i64>> = Arc::new(Mutex::new(0));
//    let read_to_id = Arc::new(Mutex::new(HashMap::new()));
    // Loop through bam generators in parallel
    let split_threads = std::cmp::max(n_threads / sample_count, 1);
    let short_threads = std::cmp::max(n_threads / bam_readers.len(), 1);
    let mut long_threads = 1;

    // Get the total amound of bases in reference using the index
    let reference_length = bio::io::fasta::Index::with_fasta_file(
        &Path::new(m.value_of("reference").unwrap())).unwrap()
        .sequences().iter().fold(0, |acc, seq| acc + seq.len);

    info!("Running SNP calling on {} shortread samples", bam_readers.len());
    bam_readers.into_par_iter().enumerate().for_each(|(sample_idx, bam_generator)|{
        process_vcf(bam_generator,
                    short_threads,
                    sample_idx,
                    sample_count,
                    &variant_matrix,
                    false,
                    m,
                    reference_length)
    });

    if longreads.len() > 0 && m.is_present("include-longread-svs"){
        long_threads = std::cmp::max(n_threads / longreads.len(), 1);
        info!("Running structural variant detection on {} longread samples", longreads.len());
        longreads.into_par_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            process_vcf(bam_generator,
                        long_threads,
                        sample_idx,
                        sample_count,
                        &variant_matrix,
                        true,
                        m,
                        reference_length)
        });
    }


    // Annoyingly read in bam file again
    let mut bam_path = "".to_string();
    let mut longreads = vec!();
    let mut bam_readers = vec!();

    if m.is_present("bam-files") {
        let bam_paths = m.values_of("bam-files").unwrap().collect::<Vec<&str>>();
        bam_readers = generate_named_bam_readers_from_bam_files(bam_paths);

    } else {
        let cache = m.value_of("outdir").unwrap().to_string() + "/*.bam";
        let bam_paths = glob(&cache).expect("Failed to read cache")
            .map(|p| p.expect("Failed to read cached bam path")
                .to_str().unwrap().to_string()).collect::<Vec<String>>();
        let bam_paths = bam_paths.iter().map(|p| &**p).collect::<Vec<&str>>();
        bam_readers = generate_named_bam_readers_from_bam_files(bam_paths);
    }

    // Process Short Read BAMs
    info!("Performing guided variant calling...");
    bam_readers.into_par_iter().enumerate().for_each(|(sample_idx, bam_generator)|{
        process_bam(bam_generator,
                    sample_idx,
                    sample_count,
                    &reference,
                    &coverage_estimators,
                    &variant_matrix,
                    &gff_map,
                    split_threads,
                    m,
                    output_prefix,
                    coverage_fold,
                    &codon_table,
                    min_var_depth,
                    contig_end_exclusion,
                    min, max, ani,
                    mode,
                    include_soft_clipping,
                    include_indels,
                    &flag_filters,
                    mapq_threshold, method, false)
    });

    // Process Long Read BAMs if they are present
    if m.is_present("longread-bam-files") {
        let longreads_path = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>();
        longreads = generate_named_bam_readers_from_bam_files(longreads_path);
    }
    if longreads.len() > 0 {
        long_threads = std::cmp::max(n_threads / longreads.len(), 1);
        longreads.into_par_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            process_bam(bam_generator,
                        sample_idx,
                        sample_count,
                        &reference,
                        &coverage_estimators,
                        &variant_matrix,
                        &gff_map,
                        split_threads,
                        m,
                        output_prefix,
                        coverage_fold,
                        &codon_table,
                        min_var_depth,
                        contig_end_exclusion,
                        min, max, ani,
                        mode,
                        include_soft_clipping,
                        include_indels,
                        &flag_filters,
                        mapq_threshold, method, true)
        });
    }

    if mode=="genotype" {
        let mut variant_matrix = variant_matrix.lock().unwrap();
        variant_matrix.generate_distances(n_threads, output_prefix);
        let e_min: f64 = m.value_of("e-min").unwrap().parse().unwrap();
        let e_max: f64 = m.value_of("e-max").unwrap().parse().unwrap();
        let pts_min: f64 = m.value_of("pts-min").unwrap().parse().unwrap();
        let pts_max: f64 = m.value_of("pts-max").unwrap().parse().unwrap();
        let phi: f64 = m.value_of("phi").unwrap().parse().unwrap();

        variant_matrix.run_fuzzy_scan(e_min, e_max, pts_min, pts_max, phi);
        variant_matrix.generate_genotypes(output_prefix);
    } else if mode=="summarize" {
        let mut variant_matrix = variant_matrix.lock().unwrap();
        variant_matrix.print_variant_stats(output_prefix);
    }
}

fn process_vcf<R: NamedBamReader + Send,
    G: NamedBamReaderGenerator<R> + Send>(
    bam_generator: G,
    split_threads: usize,
    sample_idx: usize,
    sample_count: usize,
    variant_matrix: &Arc<Mutex<VariantMatrix>>,
    longread: bool,
    m: &clap::ArgMatches,
    reference_length: u64) {
    let mut bam_generated = bam_generator.start();

    let mut bam_properties =
        AlignmentProperties::default(InsertSize::default());

    let stoit_name = bam_generated.name().to_string();


    bam_generated.set_threads(split_threads);
    let header = bam_generated.header().clone(); // bam header
    let target_names = header.target_names(); // contig names

    // for each genomic position, only has hashmap when variants are present. Includes read ids
    let mut variant_map = HashMap::new();

    // Write Bam Early and then reread in later
    bam_generated.finish();

    // Get VCF file from BAM using freebayes of SVIM
    let mut vcf_reader = get_vcf(&stoit_name,
                                 &m,
                                 sample_idx,
                                 split_threads,
                                 longread,
                                 reference_length);
    vcf_reader.set_threads(split_threads);
    let mut sample_idx = sample_idx;
    if longread {
        sample_idx = sample_count - sample_idx - 1;
    }

    info!("Collecting VCF records for sample {}", sample_idx);
    vcf_reader.records().into_iter().for_each(|vcf_record| {
        let mut vcf_record = vcf_record.unwrap();
        let header = vcf_record.header();
        let variant_rid = vcf_record.rid().unwrap();
        // Check bam header names and vcf header names are in same order
        // Sanity check
        if target_names[variant_rid as usize]
            == header.rid2name(variant_rid).unwrap() {
            let base_option = Base::from_vcf_record(&mut vcf_record,
                                                    sample_count,
                                                    sample_idx,
                                                    longread);
            match base_option {

                Some(bases) => {
                    for base in bases {
                        let variant_con = variant_map.entry(variant_rid as i32).or_insert(HashMap::new());
                        let variant_pos = variant_con.entry(base.pos).or_insert(HashMap::new());
                        variant_pos.entry(base.variant.to_owned()).or_insert(base);
                    }
                },
                None => {},
            }
        } else {
            panic!("Bug: VCF record reference ids do not match BAM reference ids. Perhaps BAM is unsorted?")
        }
    });

    let mut variant_matrix = variant_matrix.lock().unwrap();

    variant_matrix.add_sample(stoit_name.clone(), sample_idx, &variant_map, &header);
}

/// Process all reads in a BAM file
fn process_bam<R: NamedBamReader + Send,
                G: NamedBamReaderGenerator<R> + Send>(
    bam_generator: G,
    sample_idx: usize,
    sample_count: usize,
    reference: &Arc<Mutex<bio::io::fasta::IndexedReader<File>>>,
    coverage_estimators: &Arc<Mutex<&mut Vec<CoverageEstimator>>>,
    variant_matrix: &Arc<Mutex<VariantMatrix>>,
    gff_map: &Arc<Mutex<HashMap<String, Vec<Record>>>>,
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
    longread: bool) {

    let mut bam_generated = bam_generator.start();

    // Adjust the sample index if the bam is from long reads
    let mut sample_idx = sample_idx;
    if longread {
        sample_idx = sample_count - sample_idx - 1;
    }

    let mut bam_properties =
        AlignmentProperties::default(InsertSize::default());

    let stoit_name = bam_generated.name().to_string();

    bam_generated.set_threads(split_threads);

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


    // for record in records
    let mut skipped_reads = 0;

    while bam_generated.read(&mut record)
        .expect("Error while reading BAM record") == true {
        if (!flag_filters.include_supplementary && record.is_supplementary() && !longread) ||
            (!flag_filters.include_secondary && record.is_secondary() && !longread) ||
            (!flag_filters.include_improper_pairs && !record.is_proper_pair() && !longread) {
            skipped_reads += 1;
            continue;
        } else if (!flag_filters.include_secondary && record.is_secondary() && longread) {
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
                        ups_and_downs,
                        &coverage_estimators,
                        min, max,
                        total_indels_in_current_contig as usize,
                        contig_end_exclusion,
                        min_var_depth,
                        contig_len,
                        contig_name,
                        &variant_matrix,
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

                let mut reference = reference.lock().unwrap();
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
            }

            if !record.is_supplementary() {
                num_mapped_reads_in_current_contig += 1;
            }

            // Lock variant matrix to collect variants
            let mut variant_matrix = variant_matrix.lock().unwrap();
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
                        for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                            match variant_matrix.variants(tid, cursor as i64) {
                                Some(current_variants) => {
                                    let read_char = record.seq()[qpos];
                                    let refr_char = ref_seq[cursor as usize];
                                    current_variants.par_iter_mut().for_each(|(variant, base)| {
                                        match variant {
                                            Variant::SNV(alt) => {
                                                if *alt != refr_char && *alt == read_char {
                                                    base.assign_read(record.qname().to_vec());
                                                    base.truedepth[sample_idx] += 1;
                                                }
                                            },
                                            Variant::None => {
                                                if refr_char == read_char {
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
                            cursor += 1;
                        }
                        if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        read_cursor += cig.len() as usize;
                    },
                    Cigar::Del(_) => {
                        if longread {
//                            let refr = (ref_seq[cursor as usize] as char).to_string();
//                            let insert = refr +
//                                &std::iter::repeat("N").take(cig.len() as usize).collect::<String>();
//                            let refr = str::from_utf8(&ref_seq[cursor as usize..
//                                cursor as usize + cig.len() as usize]).unwrap().to_string();
//                            if refr != insert {
//                                let indel_map = indels
//                                    .entry(cursor as i32).or_insert(BTreeMap::new());
//                                let id = indel_map.entry(insert)
//                                    .or_insert(BTreeSet::new());
//                                let mut read_to_id = read_to_id.lock().unwrap();
//                                id.insert(read_to_id[&record.qname().to_vec()]);
//                                total_indels_in_current_contig += cig.len() as u64;
//                            }
                        }

                        cursor += cig.len() as usize;
                    },
                    Cigar::RefSkip(_) => {
                        // if D or N, move the cursor
                        cursor += cig.len() as usize;
                    },
                    Cigar::Ins(_) => {
//                        if longreads {
//                            let refr = (ref_seq[cursor as usize] as char).to_string();
//                            let insert = match str::from_utf8(&record.seq().as_bytes()[
//                                read_cursor..read_cursor + cig.len() as usize]) {
//                                Ok(ins) => { ins.to_string() },
//                                Err(_e) => { "".to_string() },
//                            };
//
//                            let indel_map = indels.entry(cursor as i32)
//                                .or_insert(BTreeMap::new());
//
//                            let id = indel_map.entry(refr + &insert)
//                                .or_insert(BTreeSet::new());
//                            let mut read_to_id = read_to_id.lock().unwrap();
//                            id.insert(read_to_id[&record.qname().to_vec()]);
//                        }
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

        process_previous_contigs_var(
            mode,
            ani,
            last_tid,
            ups_and_downs,
            &coverage_estimators,
            min, max,
            total_indels_in_current_contig as usize,
            contig_end_exclusion,
            min_var_depth,
            contig_len,
            contig_name,
            &variant_matrix,
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

fn process_previous_contigs_var(
    mode: &str,
    ani: f32,
    last_tid: i32,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &Arc<Mutex<&mut Vec<CoverageEstimator>>>,
    min: f32, max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u64,
    min_var_depth: usize,
    contig_len: usize,
    contig_name: Vec<u8>,
    variant_matrix: &Arc<Mutex<VariantMatrix>>,
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
    stoit_name: &str,
    longread: bool,) {

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

        let mut variant_struct = VariantStats::new_contig_stats(min as f64,
                                                              max as f64,
                                                              contig_end_exclusion);

        let mut variant_matrix = variant_matrix.lock().unwrap();

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
            "summarize" | "genotype" => {
                // calculates minimum number of genotypes possible for each variant location
                variant_matrix.add_contig(variant_struct,
                                          sample_count,
                                          sample_idx,
                                          ref_sequence);
            },
            "evolve" => {
                let gff_map = gff_map.lock().unwrap();
                variant_struct.calc_gene_mutations(&*gff_map, &ref_sequence, codon_table);
            },
            "polish" => {
                let stoit_name = stoit_name
                    .split("..").last().unwrap()
                    .split("/").last().unwrap();
                let output_prefix = output_prefix.to_string() + "_" + stoit_name;
                variant_struct.polish_contig(&ref_sequence,
                                            &output_prefix);
            }
            _ => {panic!("unknown mode {}", mode);},
        }
    }
}

/// Get or generate vcf file
pub fn get_vcf(stoit_name: &str, m: &clap::ArgMatches,
               sample_idx: usize, threads: usize, longread: bool, reference_length: u64) -> Reader {
    // if vcfs are already provided find correct one first
    if m.is_present("vcfs") {
        let vcf_paths: Vec<&str> = m.values_of("vcfs").unwrap().collect();
        // Filter out bams that don't have stoit name and get their sample idx
        let vcf_path: Vec<&str> = vcf_paths.iter().cloned()
            .filter(|x| x.contains(stoit_name)).collect();

        if vcf_path.len() > 1 || vcf_path.len() == 0 {
            info!("Could not associate VCF file with current BAM file. Re-running variant calling");
            if longread {
                let bam_path: &str = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
                return generate_vcf(bam_path, m, threads, longread, reference_length)
            } else {
                let bam_path: &str = m.values_of("bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
                return generate_vcf(bam_path, m, threads, longread, reference_length)
            }
        } else {
            let vcf_path = vcf_path[0];
            let vcf = Reader::from_path(vcf_path).unwrap();
            return vcf
        }
    } else if longread {
        let bam_path: &str = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
        return generate_vcf(bam_path, m, threads, longread, reference_length)
    } else if m.is_present("bam-files") {
        let bam_path: &str = m.values_of("bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
        return generate_vcf(bam_path, m, threads, longread, reference_length)
    } else {
        // We are streaming a generated bam file, so we have had to cache the bam for this to work
        let cache = m.value_of("outdir").unwrap().to_string() + "/";
        let stoit_name: Vec<&str> = stoit_name.split("/").collect();
        let stoit_name = stoit_name.join(".");
        let stoit_name = stoit_name.replace(|c: char| !c.is_ascii(), "");
        let bam_path = cache + &(stoit_name + ".bam");
        info!("Cached bam path {} ", bam_path);

        return generate_vcf(&bam_path, m, threads, longread, reference_length)
    }

}

/// Makes direct call to freebayes or SVIM
pub fn generate_vcf(bam_path: &str, m: &clap::ArgMatches, threads: usize, longread: bool, reference_length: u64) -> Reader {

    // setup temp directory
    let tmp_dir = TempDir::new("lorikeet_fifo")
        .expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    if !longread {
        external_command_checker::check_for_freebayes();
        external_command_checker::check_for_freebayes_parallel();
        external_command_checker::check_for_fasta_generate_regions();
        external_command_checker::check_for_samclip();
        external_command_checker::check_for_samtools();
        external_command_checker::check_for_vt();
        external_command_checker::check_for_bcftools();

        let region_size = reference_length / threads as u64;

        let index_path = m.value_of("reference").unwrap().to_string() + ".fai";

        let freebayes_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/freebayes.vcf");
//        let freebayes_path = &("freebayes.vcf");
        let tmp_bam_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/tmp.bam");

        // Generate uncompressed filtered SAM file
        let sam_cmd_string = format!(
            "samtools view -h -O SAM {} -@ {} | \
            samclip --max 10 --ref {} | \
            samtools sort -@ {} -n -l 0 -T /tmp | \
            samtools fixmate -@ {} -m - - | \
            samtools sort -@ {} -l 0 -T /tmp | \
            samtools markdup -@ {} -T /tmp -r -s - - > {}",
            bam_path,
            threads-1,
            index_path,
            threads-1,
            threads-1,
            threads-1,
            threads-1,
            tmp_bam_path);
        debug!("Queuing cmd_string: {}", sam_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&sam_cmd_string)
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .spawn()
                .expect("Unable to execute bash"), "samtools");

        // check and build bam index if it doesn't exist
        if !Path::new(&(tmp_bam_path.to_string() + ".bai")).exists() {
            bam::index::build(tmp_bam_path, Some(&(tmp_bam_path.to_string() + ".bai")),
                              bam::index::Type::BAI, threads as u32).expect(
                &format!("Unable to index bam at {}", &tmp_bam_path));
        }

        // Variant calling pipeline adapted from Snippy but without all of the rewriting of BAM files
        let vcf_cmd_string = format!(
            "set -e -o pipefail;  \
            freebayes-parallel <(fasta_generate_regions.py {} {}) {} -f {} -C {} -q {} \
            --min-repeat-entropy {} --strict-vcf -m {} {} | \
            vt normalize -n -r {} - | \
            bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > {}",
            index_path,
            region_size,
            threads,
            m.value_of("reference").unwrap(),
            m.value_of("min-variant-depth").unwrap(),
            m.value_of("base-quality-threshold").unwrap(),
            m.value_of("min-repeat-entropy").unwrap(),
            m.value_of("mapq-threshold").unwrap(),
            tmp_bam_path,
            m.value_of("reference").unwrap(),
            freebayes_path);
        debug!("Queuing cmd_string: {}", vcf_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&vcf_cmd_string)
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .spawn()
                .expect("Unable to execute bash"), "freebayes");
        debug!("VCF Path {:?}", freebayes_path);
        let vcf_reader = Reader::from_path(freebayes_path)
            .expect("Failed to read pilon vcf output");

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader
    } else {
        external_command_checker::check_for_svim();
        let svim_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/svim");

        // check and build bam index if it doesn't exist
        if !Path::new(&(bam_path.to_string() + ".bai")).exists() {
            bam::index::build(bam_path, Some(&(bam_path.to_string() + ".bai")),
                              bam::index::Type::BAI, threads as u32).expect(
                &format!("Unable to index bam at {}", &bam_path));
        }

        let cmd_string = format!(
            "set -e -o pipefail; svim alignment --read_names --skip_genotyping \
            --min_mapq {} --minimum_depth {} --sequence_alleles {} {} {}",
            m.value_of("mapq-threshold").unwrap(),
            m.value_of("min-variant-depth").unwrap(),
            svim_path,
            bam_path,
            m.value_of("reference").unwrap());
        debug!("Queuing cmd_string: {}", cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&cmd_string)
                .stderr(std::process::Stdio::null())
//                .stdout(std::process::Stdio::null())
                .spawn()
                .expect("Unable to execute bash"), "svim");
        let vcf_path = &(svim_path.to_string() + "/variants.vcf");
        debug!("VCF Path {:?}", vcf_path);
        let vcf_reader = Reader::from_path(vcf_path)
            .expect("Failed to read SVIM vcf output");

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader
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
//            reads_mapped_vec = variant_variants(
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