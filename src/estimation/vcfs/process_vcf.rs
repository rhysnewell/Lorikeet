use bird_tool_utils::command;
use rust_htslib::{bcf, bcf::Read};
use std;
use std::collections::HashMap;

use coverm::bam_generator::*;
use estimation::variant_matrix::*;
use external_command_checker;
use model::variants::*;
use utils::*;

use crate::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use scoped_threadpool::Pool;
use std::io::Write;
use std::path::Path;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::Builder;

#[allow(unused)]
pub fn process_vcf<'b, R: IndexedNamedBamReader + Send, G: NamedBamReaderGenerator<R> + Send>(
    bam_generator: G,
    split_threads: usize,
    ref_idx: usize,
    sample_idx: usize,
    mut sample_count: usize,
    variant_matrix: &'b mut VariantMatrix,
    longread: bool,
    m: &'b clap::ArgMatches,
    // sample_groups: &mut HashMap<&str, HashSet<String>>,
    genomes_and_contigs: &'b GenomesAndContigs,
    reference_map: &'b HashMap<usize, String>,
    mut short_sample_count: usize,
    concatenated_genomes: &'b Option<String>,
) {
    let mut bam_generated = bam_generator.start();
    let mut stoit_name = bam_generated.name().to_string();

    // if longread {
    //     let group = sample_groups.entry("long").or_insert(HashSet::new());
    //     group.insert(stoit_name.clone());
    // } else {
    //     let group = sample_groups.entry("short").or_insert(HashSet::new());
    //     group.insert(stoit_name.clone());
    // }

    let reference = &genomes_and_contigs.genomes[ref_idx];
    let mut reference_file = retrieve_reference(concatenated_genomes);

    bam_generated.set_threads(split_threads);

    let header = bam_generated.header().clone(); // bam header
    let target_lens: Vec<u64> = (0..header.target_count())
        .into_iter()
        .map(|tid| header.target_len(tid).unwrap())
        .collect();
    let target_names = header.target_names();
    //     .into_iter()
    //     .map(|target| String::from_utf8(target.to_vec()).unwrap())
    //     .collect(); // contig names

    let bam_path = bam_generated.path().to_string();

    // minimum PHRED base quality
    let bq = m
        .value_of("base-quality-threshold")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    // Minimum MAPQ value
    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());
    // Minimum count for a SNP to be considered
    let min_variant_depth: i32 = m.value_of("min-variant-depth").unwrap().parse().unwrap();
    let min_variant_quality = m
        .value_of("min-variant-quality")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let mut total_records = 0;
    variant_matrix.add_sample_name(stoit_name.to_string(), sample_idx);
    let mut ref_target_names = Vec::new();
    let mut ref_target_lengths = Vec::new();

    // for each genomic position, only has hashmap when variants are present. Includes read ids
    target_names
        .iter()
        .enumerate()
        .for_each(|(tid, contig_name)| {
            let target_name = String::from_utf8(contig_name.to_vec()).unwrap();
            if target_name.contains(reference)
                || match genomes_and_contigs.contig_to_genome.get(&target_name) {
                    Some(ref_id) => *ref_id == ref_idx,
                    None => false,
                }
            {
                // use pileups to call SNPs for low quality variants
                // That are usually skipped by GATK
                ref_target_names.push(target_name.clone());
                let target_len = target_lens[tid];
                ref_target_lengths.push(target_len);
                variant_matrix.add_info(ref_idx, tid, target_name.as_bytes().to_vec(), target_len);
                {
                    bam_generated.fetch(tid as u32, 0, target_len);
                    match bam_generated.pileup() {
                        Some(pileups) => {
                            let mut ref_seq = Vec::new();
                            // Update all contig information
                            fetch_contig_from_reference(
                                &mut reference_file,
                                &contig_name.to_vec(),
                                genomes_and_contigs,
                                ref_idx as usize,
                            );
                            ref_seq = Vec::new();
                            read_sequence_to_vec(
                                &mut ref_seq,
                                &mut reference_file,
                                &contig_name.to_vec(),
                            );

                            for p in pileups {
                                // if {
                                let pileup = p.unwrap();
                                let tid = pileup.tid() as i32;
                                let pos = pileup.pos() as usize;
                                let depth = pileup.depth();

                                let refr_base = ref_seq[pos];
                                // let mut base = Base::new(tid, pos, );
                                // info!("Base dict {:?} {}", &pos_dict, pos);

                                let mut base_dict = HashMap::new();

                                let mut refr_depth = 0;
                                let mut refr_qual = 0.;

                                for alignment in pileup.alignments() {
                                    let record = alignment.record();
                                    if record.mapq() >= mapq_thresh {
                                        if !alignment.is_del() && !alignment.is_refskip() {
                                            // query position in read
                                            let qpos = alignment.qpos().unwrap();
                                            let record_qual = record.qual()[qpos] as f64;
                                            if record_qual >= bq {
                                                let read_base = alignment.record().seq()[qpos];
                                                if read_base != refr_base {
                                                    let mut base = base_dict
                                                        .entry(read_base)
                                                        .or_insert(Base::new(
                                                            tid as u32,
                                                            pos as i64,
                                                            sample_count,
                                                            vec![refr_base],
                                                        ));

                                                    base.depth[sample_idx] += 1;
                                                    base.quals[sample_idx] += record_qual;
                                                    if base.variant == Variant::None {
                                                        base.variant = Variant::SNV(read_base);
                                                    }
                                                }
                                            } else {
                                                refr_depth += 1;
                                                refr_qual += record_qual;
                                            }
                                        }
                                    }
                                }

                                // Collect refr base information
                                {
                                    let mut base = base_dict.entry(refr_base).or_insert(Base::new(
                                        tid as u32,
                                        pos as i64,
                                        sample_count,
                                        vec![refr_base],
                                    ));

                                    base.depth[sample_idx] = refr_depth;
                                    base.quals[sample_idx] = refr_qual;
                                }

                                // If more than one variant at location (including reference)
                                // Collect the variants
                                if base_dict.keys().len() > 1 {
                                    let mut variant_found = false;
                                    let mut refr_base = 0;
                                    for (var_char, base) in base_dict.iter() {
                                        if base.depth[sample_idx] >= min_variant_depth &&
                                        // base.quals[sample_idx] >= min_variant_quality &&
                                        base.variant != Variant::None
                                        {
                                            total_records += 1;
                                            variant_found = true;
                                            refr_base = base.refr[0];
                                            variant_matrix.add_variant_to_matrix(
                                                sample_idx,
                                                base,
                                                tid as usize,
                                                ref_idx,
                                            );
                                        }
                                    }
                                    // Add reference
                                    match base_dict.get(&refr_base) {
                                        Some(base) => {
                                            variant_matrix.add_variant_to_matrix(
                                                sample_idx,
                                                base,
                                                tid as usize,
                                                ref_idx,
                                            );
                                        }
                                        None => {}
                                    }
                                }
                            }
                            // }
                        }
                        None => println!("no bam for pileups"),
                    };
                    // bam_generated.set_threads(1);
                };
            }
        });

    bam_generated.finish();
    let mut variant_matrix_sync = Arc::new(Mutex::new(variant_matrix.clone()));
    // let mut thread_locker = Mutex::new(true);
    let freebayes_threads = std::cmp::max(split_threads, 1);
    // target_names.par_iter().enumerate().for_each(|(tid, target_name)| {
    //         if target_name.contains(reference)
    //             || match genomes_and_contigs.contig_to_genome.get(target_name) {
    //             Some(ref_id) => *ref_id == ref_idx,
    //             None => false,
    //         }
    //         {

    // use pileups to call SNPs for low quality variants
    // That are usually skipped by GATK
    // let target_len = target_lens[tid];

    // Get VCF file from BAM using freebayes of SVIM
    let mut vcf_reader = get_vcf(
        &stoit_name,
        &m,
        freebayes_threads,
        longread,
        ref_target_lengths,
        &concatenated_genomes.as_ref().unwrap(),
        &bam_path,
        &ref_target_names,
    );

    match vcf_reader {
        Ok(ref mut reader) => {
            // let vcf_header = reader.header();
            // let vcf_target_names: Vec<String> = vcf_header.samples()
            //     .into_iter()
            //     .map(|target| String::from_utf8(target.to_vec()).unwrap())
            //     .collect();
            // println!("{:?} {:?}", &vcf_target_names, &target_names);

            let min_qual = m.value_of("min-variant-quality").unwrap().parse().unwrap();

            let mut pool = Pool::new(freebayes_threads as u32);
            // let total_records = Arc::new(Mutex::new(total_records));

            pool.scoped(|scope| {
                for vcf_record in reader.records().into_iter() {
                    scope.execute(|| {
                        let mut vcf_record = vcf_record.unwrap();
                        // let vcf_header = vcf_record.header();
                        let variant_rid = vcf_record.rid().unwrap();
                        // let vcf_target = String::from_utf8(vcf_header.rid2name(variant_rid).unwrap().to_vec()).unwrap();
                        // let variant_rid: &usize = &target_names.iter().position(|t| t==&vcf_target).unwrap();
                        // Check bam header names and vcf header names are in same order
                        // Sanity check
                        // total_records += 1;
                        // if target_name.as_bytes()
                        //     == vcf_header.rid2name(variant_rid).unwrap() {
                        let base_option = Base::from_vcf_record(
                            &mut vcf_record,
                            sample_count,
                            sample_idx,
                            longread,
                            min_qual,
                        );
                        match base_option {
                            Some(bases) => {
                                for base in bases {
                                    let mut variant_matrix_sync =
                                        variant_matrix_sync.lock().unwrap();
                                    variant_matrix_sync.add_variant_to_matrix(
                                        sample_idx,
                                        &base,
                                        variant_rid as usize,
                                        ref_idx,
                                    );
                                }
                            }
                            None => {}
                        }
                        // } else {
                        //     panic!("Bug: VCF record reference ids do not match BAM reference ids. Perhaps BAM is unsorted?")
                        // }
                    });
                }
            });
        }
        Err(_) => {
            info!("No VCF records found for sample {}", &stoit_name);
        }
    }

    //     }
    // });
    *variant_matrix = variant_matrix_sync.lock().unwrap().clone();
}

/// Get or generate vcf file
#[allow(unused)]
pub fn get_vcf(
    stoit_name: &str,
    m: &clap::ArgMatches,
    threads: usize,
    longread: bool,
    target_lengths: Vec<u64>,
    reference: &str,
    bam_path: &str,
    target_names: &Vec<String>,
) -> std::result::Result<bcf::Reader, rust_htslib::bcf::Error> {
    // if vcfs are already provided find correct one first
    if m.is_present("vcfs") {
        let vcf_paths: Vec<&str> = m.values_of("vcfs").unwrap().collect();
        // Filter out bams that don't have stoit name and get their sample idx
        let vcf_path: Vec<&str> = vcf_paths
            .iter()
            .cloned()
            .filter(|x| x.contains(stoit_name))
            .collect();

        if vcf_path.len() > 1 || vcf_path.len() == 0 {
            info!("Could not associate VCF file with current BAM file. Re-running variant calling");
            if longread {
                // let bam_path: &str = *m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()
                //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
                return generate_vcf(
                    bam_path,
                    m,
                    threads,
                    longread,
                    target_lengths,
                    reference,
                    target_names,
                );
            } else {
                // let bam_path: &str = *m.values_of("bam-files").unwrap().collect::<Vec<&str>>()
                //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
                return generate_vcf(
                    bam_path,
                    m,
                    threads,
                    longread,
                    target_lengths,
                    reference,
                    target_names,
                );
            }
        } else {
            let vcf_path = vcf_path[0];
            let vcf = bcf::Reader::from_path(vcf_path);
            return vcf;
        }
    } else if longread && m.is_present("longread-bam-files") {
        // let bam_path: &str = *m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()
        //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
        return generate_vcf(
            bam_path,
            m,
            threads,
            longread,
            target_lengths,
            reference,
            target_names,
        );
    } else if m.is_present("bam-files") {
        // let bam_path: &str = *m.values_of("bam-files").unwrap().collect::<Vec<&str>>()
        //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
        return generate_vcf(
            bam_path,
            m,
            threads,
            longread,
            target_lengths,
            reference,
            target_names,
        );
    } else {
        // We are streaming a generated bam file, so we have had to cache the bam for this to work
        // let cache = m.value_of("bam-file-cache-directory").unwrap().to_string() + "/";
        //
        // let bam_path = cache + &(stoit_name.to_string() + ".bam");

        return generate_vcf(
            &bam_path,
            m,
            threads,
            longread,
            target_lengths,
            reference,
            target_names,
        );
    }
}

/// Makes direct call to freebayes or SVIM
#[allow(unused)]
pub fn generate_vcf(
    bam_path: &str,
    m: &clap::ArgMatches,
    threads: usize,
    longread: bool,
    target_lengths: Vec<u64>,
    reference: &str,
    target_names: &Vec<String>,
) -> std::result::Result<bcf::Reader, rust_htslib::bcf::Error> {
    // setup temp directory
    let tmp_dir = TempDir::new("lorikeet_fifo").expect("Unable to create temporary directory");

    let tmp_bam_path1 = Builder::new()
        .prefix(&(tmp_dir.path().to_str().unwrap().to_string() + "/"))
        .suffix(".bam")
        .tempfile()
        .unwrap();

    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());

    // if !longread {
    external_command_checker::check_for_freebayes();
    external_command_checker::check_for_freebayes_parallel();
    external_command_checker::check_for_samtools();
    external_command_checker::check_for_vt();

    // Old way of calulating region size, not a good method
    // let region_size = reference_length / threads as u64;
    // Now we just set it to be 100000, doesn't seem necessary to make this user defined?
    // let mut region_size = target_length / threads as u64
    //     + if target_length % threads as u64 != 0 {
    //         1
    //     } else {
    //         0
    //     };
    // if target_length <= 100000 {
    //     region_size = 100000
    // }

    let region_size = 100000;

    let index_path = format!("{}.fai", reference);

    let vcf_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/output.vcf");
    let vcf_path_prenormalization =
        &(tmp_dir.path().to_str().unwrap().to_string() + "/output_prenormalization.vcf");

    let mut region_tmp_file = Builder::new()
        .prefix(&(tmp_dir.path().to_str().unwrap().to_string() + "/"))
        .suffix(".txt")
        .tempfile()
        .unwrap();

    for (target_name, target_length) in target_names.iter().zip(target_lengths) {
        let mut total_region_covered = 0;
        let total_region_chunks = target_length / region_size
            + if target_length % region_size != 0 {
                1
            } else {
                0
            };
        for idx in (0..total_region_chunks).into_iter() {
            total_region_covered += region_size;
            if total_region_covered > target_length {
                writeln!(
                    region_tmp_file,
                    "{}:{}-{}",
                    &target_name,
                    total_region_covered - region_size,
                    target_length,
                )
                .unwrap();
            } else {
                writeln!(
                    region_tmp_file,
                    "{}:{}-{}",
                    &target_name,
                    total_region_covered - region_size,
                    total_region_covered,
                )
                .unwrap();
            }
        }
    }

    // Variant calling pipeline adapted from Snippy but without all of the rewriting of BAM files
    let vcf_cmd_string = format!(
        "set -e -o pipefail;  \
            freebayes-parallel {:?} {} -f {} -C {} -q {} \
            -p {} --strict-vcf -m {} {} > {}",
        region_tmp_file.path(),
        threads,
        &reference,
        m.value_of("min-variant-depth").unwrap(),
        m.value_of("base-quality-threshold").unwrap(),
        // m.value_of("min-repeat-entropy").unwrap(),
        m.value_of("ploidy").unwrap(),
        mapq_thresh,
        &bam_path,
        &vcf_path_prenormalization,
    );
    let vt_cmd_string = format!(
        "vt normalize {} -n -r {} -o {}",
        &vcf_path_prenormalization, &reference, vcf_path,
    );
    debug!("Queuing cmd_string: {}", vcf_cmd_string);
    command::finish_command_safely(
        std::process::Command::new("bash")
            .arg("-c")
            .arg(&vcf_cmd_string)
            .stderr(std::process::Stdio::piped())
            // .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Unable to execute bash"),
        "freebayes",
    );
    debug!("Queuing cmd_string: {}", vt_cmd_string);
    command::finish_command_safely(
        std::process::Command::new("bash")
            .arg("-c")
            .arg(&vt_cmd_string)
            .stderr(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .expect("Unable to execute bash"),
        "vt",
    );
    debug!("VCF Path {:?}", vcf_path);
    let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

    tmp_dir.close().expect("Failed to close temp directory");
    return vcf_reader;
    // } else {
    //     external_command_checker::check_for_svim();
    // let svim_path = tmp_dir.path().to_str().unwrap().to_string();
    //
    // // check and build bam index if it doesn't exist
    // if !Path::new(&(bam_path.to_string() + ".bai")).exists() {
    //     bam::index::build(
    //         bam_path,
    //         Some(&(bam_path.to_string() + ".bai")),
    //         bam::index::Type::BAI,
    //         threads as u32,
    //     )
    //     .expect(&format!("Unable to index bam at {}", &bam_path));
    // }
    //
    // let cmd_string = format!(
    //     "set -e -o pipefail; samtools view -bh {} {} > {} && \
    //     svim alignment --read_names --skip_genotyping \
    //     --tandem_duplications_as_insertions --interspersed_duplications_as_insertions \
    //     --min_mapq {} --sequence_alleles {} {} {}",
    //     bam_path,
    //     &target_name,
    //     tmp_bam_path1.path().to_str().unwrap().to_string(),
    //     mapq_thresh,
    //     &svim_path,
    //     tmp_bam_path1.path().to_str().unwrap().to_string(),
    //     &reference
    // );
    // debug!("Queuing cmd_string: {}", cmd_string);
    // command::finish_command_safely(
    //     std::process::Command::new("bash")
    //         .arg("-c")
    //         .arg(&cmd_string)
    //         .stderr(std::process::Stdio::piped())
    //         //                .stdout(std::process::Stdio::null())
    //         .spawn()
    //         .expect("Unable to execute bash"),
    //     "svim",
    // );
    // let vcf_path = &(svim_path + "/variants.vcf");
    // debug!("VCF Path {:?}", vcf_path);
    // let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));
    //
    // tmp_dir.close().expect("Failed to close temp directory");
    // return vcf_reader;
    // }
}
