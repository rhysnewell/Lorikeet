use bird_tool_utils::command;
use rust_htslib::{bam, bcf, bcf::Read};
use std;
use std::collections::{HashMap, HashSet};

use coverm::bam_generator::*;
use estimation::variant_matrix::*;
use external_command_checker;
use model::variants::*;
use utils::*;

use crate::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use nix::sys::stat;
use nix::unistd;
use std::path::Path;
use std::str;
use tempdir::TempDir;
use tempfile::NamedTempFile;

#[allow(unused)]
pub fn process_vcf<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bam_generator: G,
    split_threads: usize,
    mut prev_ref_idx: &mut i32,
    mut per_ref_sample_idx: &mut i32,
    mut sample_count: usize,
    variant_matrix: &mut VariantMatrix,
    longread: bool,
    m: &clap::ArgMatches,
    sample_groups: &mut HashMap<&str, HashSet<String>>,
    genomes_and_contigs: &GenomesAndContigs,
    reference_map: &HashMap<usize, String>,
    mut short_sample_count: usize,
    concatenated_genomes: &Option<NamedTempFile>,
) {
    let mut bam_generated = bam_generator.start();
    let mut stoit_name = bam_generated.name().to_string();

    if longread {
        let group = sample_groups.entry("long").or_insert(HashSet::new());
        group.insert(stoit_name.clone());
    } else {
        let group = sample_groups.entry("short").or_insert(HashSet::new());
        group.insert(stoit_name.clone());
    }

    bam_generated.set_threads(split_threads);
    let header = bam_generated.header().clone(); // bam header
    let target_names = header.target_names(); // contig names

    let bam_path = bam_generated.path().to_string();

    // Set AUX tags used by GATK
    // Also saves BAM file to disk
    let read_group = stoit_name.as_bytes();
    let mut record: bam::record::Record = bam::record::Record::new();
    while bam_generated
        .read(&mut record)
        .expect("Failure to read BAM record")
        == true
    {
        record.push_aux("SM".as_bytes(), &bam::record::Aux::String(read_group));
        record.push_aux("LB".as_bytes(), &bam::record::Aux::String("N".as_bytes()));
        record.push_aux("PL".as_bytes(), &bam::record::Aux::String("N".as_bytes()));
        record.push_aux("PU".as_bytes(), &bam::record::Aux::String("N".as_bytes()));
    }
    bam_generated.finish();

    // Adjust indices based on whether or not we are using a concatenated reference or not
    let (reference, ref_idx) = if stoit_name.contains(".fna") && reference_map.len() > 1 {
        debug!("Stoit_name {:?} {:?}", &stoit_name, &reference_map);
        let reference_stem = stoit_name.split(".fna").next().unwrap();
        debug!("possible reference stem {:?}", reference_stem);
        let ref_idx = match genomes_and_contigs.genome_index(&reference_stem.to_string()) {
            Some(idx) => idx,
            None => panic!("Unable to retrieve reference index"),
        };
        debug!("Actual reference idx {:?}", ref_idx);

        let reference = reference_map
            .get(&ref_idx)
            .expect("Unable to retrieve reference path")
            .clone();
        (reference, ref_idx)
    } else {
        match concatenated_genomes {
            Some(ref temp_file) => {
                stoit_name = format!(
                    "{}.{}",
                    temp_file.path().file_name().unwrap().to_str().unwrap(),
                    stoit_name,
                );
                debug!("Renamed Stoit_name {:?}", &stoit_name);

                (temp_file.path().to_str().unwrap().to_string(), 0)
            }
            None => {
                retrieve_genome_from_contig(target_names[0], genomes_and_contigs, reference_map)
            }
        }
    };

    debug!(
        "retrieving genome id with contig {:?} from {} for sample {}",
        str::from_utf8(&target_names[0]),
        &reference,
        &stoit_name
    );

    let reference_length = match bio::io::fasta::Index::with_fasta_file(&Path::new(&reference)) {
        Ok(index) => index.sequences().iter().fold(0, |acc, seq| acc + seq.len),
        Err(_e) => {
            generate_faidx(&reference);
            bio::io::fasta::Index::with_fasta_file(&Path::new(&reference))
                .unwrap()
                .sequences()
                .iter()
                .fold(0, |acc, seq| acc + seq.len)
        }
    };

    // Get VCF file from BAM using freebayes of SVIM
    let mut vcf_reader = get_vcf(
        &stoit_name,
        &m,
        split_threads,
        longread,
        reference_length,
        &reference,
        &bam_path,
    );

    // for each genomic position, only has hashmap when variants are present. Includes read ids
    let mut variant_map = HashMap::new();

    match vcf_reader {
        Ok(ref mut reader) => {
            reader
                .set_threads(split_threads)
                .expect("Unable to set threads on VCF reader");

            let min_qual = m.value_of("min-variant-quality").unwrap().parse().unwrap();

            let mut total_records = 0;
            reader.records().into_iter().for_each(|vcf_record| {
                let mut vcf_record = vcf_record.unwrap();
                let vcf_header = vcf_record.header();
                let variant_rid = vcf_record.rid().unwrap();
                // Check bam header names and vcf header names are in same order
                // Sanity check
                total_records += 1;
                if target_names[variant_rid as usize]
                    == vcf_header.rid2name(variant_rid).unwrap() {
                    let base_option =
                        Base::from_vcf_record(
                            &mut vcf_record,
                            sample_count,
                            *per_ref_sample_idx as usize,
                            longread,
                            min_qual,
                        );
                    match base_option {
                        Some(bases) => {
                            for base in bases {
                                let variant_con = variant_map
                                    .entry(variant_rid as i32).or_insert(HashMap::new());
                                let variant_pos = variant_con
                                    .entry(base.pos).or_insert(HashMap::new());
                                variant_pos.entry(base.variant.to_owned()).or_insert(base);
                            }
                        },
                        None => {},
                    }
                } else {
                    panic!("Bug: VCF record reference ids do not match BAM reference ids. Perhaps BAM is unsorted?")
                }
            });
            info!(
                "Collected {} VCF records for sample {}",
                total_records, &stoit_name
            );
            if longread {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    (short_sample_count as i32 + *per_ref_sample_idx) as usize,
                    &variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            } else {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    *per_ref_sample_idx as usize,
                    &variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            }
        }
        Err(_) => {
            info!("No VCF records found for sample {}", &stoit_name);
            if longread {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    (short_sample_count as i32 + *per_ref_sample_idx) as usize,
                    &variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            } else {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    *per_ref_sample_idx as usize,
                    &variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            }
        }
    }
}

/// Get or generate vcf file
#[allow(unused)]
pub fn get_vcf(
    stoit_name: &str,
    m: &clap::ArgMatches,
    threads: usize,
    longread: bool,
    reference_length: u64,
    reference: &String,
    bam_path: &str,
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
                return generate_vcf(bam_path, m, threads, longread, reference_length, reference);
            } else {
                // let bam_path: &str = *m.values_of("bam-files").unwrap().collect::<Vec<&str>>()
                //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
                return generate_vcf(bam_path, m, threads, longread, reference_length, reference);
            }
        } else {
            let vcf_path = vcf_path[0];
            let vcf = bcf::Reader::from_path(vcf_path);
            return vcf;
        }
    } else if longread && m.is_present("longread-bam-files") {
        // let bam_path: &str = *m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()
        //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
        return generate_vcf(bam_path, m, threads, longread, reference_length, reference);
    } else if m.is_present("bam-files") {
        // let bam_path: &str = *m.values_of("bam-files").unwrap().collect::<Vec<&str>>()
        //     .iter().filter(|bam| bam.contains(&stoit_name)).collect::<Vec<&&str>>()[0];
        return generate_vcf(bam_path, m, threads, longread, reference_length, reference);
    } else {
        // We are streaming a generated bam file, so we have had to cache the bam for this to work
        // let cache = m.value_of("bam-file-cache-directory").unwrap().to_string() + "/";
        //
        // let bam_path = cache + &(stoit_name.to_string() + ".bam");

        return generate_vcf(&bam_path, m, threads, longread, reference_length, reference);
    }
}

/// Makes direct call to freebayes or SVIM
#[allow(unused)]
pub fn generate_vcf(
    bam_path: &str,
    m: &clap::ArgMatches,
    threads: usize,
    longread: bool,
    reference_length: u64,
    reference: &String,
) -> std::result::Result<bcf::Reader, rust_htslib::bcf::Error> {
    // setup temp directory
    let tmp_dir = TempDir::new("lorikeet_fifo").expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());

    if !longread {
        external_command_checker::check_for_gatk();
        external_command_checker::check_for_samtools();
        external_command_checker::check_for_vt();

        // Old way of calulating region size, not a good method
        // let region_size = reference_length / threads as u64;
        // Now we just set it to be 100000, doesn't seem necessary to make this user defined?
        let region_size = 10000;
        let index_path = reference.clone() + ".fai";

        let vcf_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/output.vcf");
        let vcf_path_prenormalization =
            &(tmp_dir.path().to_str().unwrap().to_string() + "/output_prenormalization.vcf");

        //        let freebayes_path = &("freebayes.vcf");
        let tmp_bam_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/tmp.bam");

        // Generate uncompressed filtered SAM file
        let sam_cmd_string = format!(
            "samtools sort -@ {} -n -l 0 -T /tmp {} | \
            samtools fixmate -@ {} -m - - | \
            gatk AddOrReplaceReadGroups -I - -O {} -SM 1 -LB N -PL N -PU N",
            threads - 1,
            bam_path,
            threads - 1,
            tmp_bam_path,
        );
        debug!("Queuing cmd_string: {}", sam_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&sam_cmd_string)
                .stderr(std::process::Stdio::piped())
                // .stdout(std::process::Stdio::piped())
                .spawn()
                .expect("Unable to execute bash"),
            "samtools",
        );

        // check and build bam index if it doesn't exist
        bam::index::build(
            bam_path,
            Some(&(bam_path.to_string() + ".bai")),
            bam::index::Type::BAI,
            threads as u32,
        )
        .expect(&format!("Unable to index bam at {}", &bam_path));

        // Variant calling pipeline adapted from Snippy but without all of the rewriting of BAM files
        let vcf_cmd_string = format!(
            "set -e -o pipefail;  \
            gatk HaplotypeCaller -I {} -R {} -O {} --native-pair-hmm-threads {} --sample-ploidy {} -mbq {} \
            --annotation AlleleFraction --annotation DepthPerAlleleBySample --minimum-mapping-quality {}",
            tmp_bam_path,
            &reference,
            &vcf_path_prenormalization,
            threads,
            m.value_of("ploidy").unwrap(),
            m.value_of("base-quality-threshold").unwrap(),
            mapq_thresh,
        );
        let vt_cmd_string = format!(
            "vt normalize -n -r {} {} > {}",
            &reference, &vcf_path_prenormalization, vcf_path,
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
            "gatk",
        );
        debug!("Queuing cmd_string: {}", vt_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&vt_cmd_string)
                .stderr(std::process::Stdio::piped())
                // .stdout(std::process::Stdio::piped())
                .spawn()
                .expect("Unable to execute bash"),
            "vt",
        );
        debug!("VCF Path {:?}", vcf_path);
        let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader;
    } else {
        external_command_checker::check_for_svim();
        let svim_path = tmp_dir.path().to_str().unwrap().to_string();

        // check and build bam index if it doesn't exist
        if !Path::new(&(bam_path.to_string() + ".bai")).exists() {
            bam::index::build(
                bam_path,
                Some(&(bam_path.to_string() + ".bai")),
                bam::index::Type::BAI,
                threads as u32,
            )
            .expect(&format!("Unable to index bam at {}", &bam_path));
        }

        let cmd_string = format!(
            "set -e -o pipefail; svim alignment --read_names --skip_genotyping \
            --tandem_duplications_as_insertions --interspersed_duplications_as_insertions \
            --min_mapq {} --sequence_alleles {} {} {}",
            mapq_thresh, &svim_path, &bam_path, &reference
        );
        debug!("Queuing cmd_string: {}", cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&cmd_string)
                .stderr(std::process::Stdio::piped())
                //                .stdout(std::process::Stdio::null())
                .spawn()
                .expect("Unable to execute bash"),
            "svim",
        );
        let vcf_path = &(svim_path + "/variants.vcf");
        debug!("VCF Path {:?}", vcf_path);
        let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader;
    }
}
