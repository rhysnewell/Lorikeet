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
pub fn process_vcf<R: IndexedNamedBamReader>(
    mut bam_generated: R,
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
    // let mut bam_generated = bam_generator.start();
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

    // Adjust indices based on whether or not we are using a concatenated reference or not
    let (reference, ref_idx) = if stoit_name.contains(".fna") && reference_map.len() > 1 {
        debug!("Stoit_name {:?} {:?}", &stoit_name, &reference_map);
        let reference_stem = stoit_name.split(".fna").next().unwrap();
        debug!("possible reference stem {:?}", reference_stem);

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
            None => match genomes_and_contigs.genome_index(&reference_stem.to_string()) {
                Some(ref_idx) => {
                    debug!("Actual reference idx {:?}", ref_idx);

                    let reference = reference_map
                        .get(&ref_idx)
                        .expect("Unable to retrieve reference path")
                        .clone();
                    (reference, ref_idx)
                }
                None => {
                    retrieve_genome_from_contig(target_names[0], genomes_and_contigs, reference_map)
                }
            },
        }
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

    let mut reference_file = retrieve_reference(concatenated_genomes);

    // minimum PHRED base quality
    let bq = m
        .value_of("base-quality-threshold")
        .unwrap()
        .parse()
        .unwrap();
    // Minimum MAPQ value
    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());
    // Minimum count for a SNP to be considered
    let min_variant_depth: i32 = m.value_of("min-variant-depth").unwrap().parse().unwrap();
    let min_variant_quality: f32 = m.value_of("min-variant-quality").unwrap().parse().unwrap();

    // for each genomic position, only has hashmap when variants are present. Includes read ids
    let mut variant_map = HashMap::new();

    // use pileups to call SNPs for low quality variants
    // That are usually skipped by GATK
    match bam_generated.pileup() {
        Some(pileups) => {
            let mut last_tid = -2;
            let mut contig_name = Vec::new();
            let mut ref_idx = 0;
            let mut ref_seq = Vec::new();

            for p in pileups {
                if !longread {
                    let pileup = p.unwrap();
                    let tid = pileup.tid() as i32;
                    let pos = pileup.pos() as usize;
                    let depth = pileup.depth();
                    if tid != last_tid {
                        if tid < last_tid {
                            error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                            panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                        }
                        // Update all contig information
                        contig_name = target_names[tid as usize].to_vec();
                        fetch_contig_from_reference(
                            &mut reference_file,
                            &contig_name,
                            genomes_and_contigs,
                            ref_idx,
                        );
                        ref_seq = Vec::new();
                        read_sequence_to_vec(&mut ref_seq, &mut reference_file, &contig_name);
                        last_tid = tid;
                    }

                    let refr_base = ref_seq[pos];
                    // let mut base = Base::new(tid, pos, );

                    let mut base_dict = HashMap::new();

                    for alignment in pileup.alignments() {
                        let record = alignment.record();
                        if record.mapq() >= mapq_thresh {
                            if !alignment.is_del() && !alignment.is_refskip() {
                                // query position in read
                                let qpos = alignment.qpos().unwrap();
                                let record_qual = record.qual()[qpos];
                                if record_qual >= bq {
                                    let read_base = alignment.record().seq()[qpos];
                                    if read_base != refr_base {
                                        let mut base =
                                            base_dict.entry(read_base).or_insert(Base::new(
                                                tid as u32,
                                                pos as i64,
                                                sample_count,
                                                vec![refr_base],
                                            ));

                                        if base.variant == Variant::None {
                                            base.variant = Variant::SNV(read_base);
                                        }

                                        base.depth[*per_ref_sample_idx as usize] += 1;
                                        base.quals[*per_ref_sample_idx as usize] +=
                                            record_qual as f32;
                                    }
                                }
                            }
                        }
                    }

                    for (var_char, base) in base_dict {
                        if base.depth[*per_ref_sample_idx as usize] >= min_variant_depth
                            && base.quals[*per_ref_sample_idx as usize] >= min_variant_quality
                        {
                            let variant_con =
                                variant_map.entry(tid as i32).or_insert(HashMap::new());
                            let variant_pos = variant_con.entry(base.pos).or_insert(HashMap::new());

                            // Overwrite any existing variants called by mpileup
                            variant_pos.insert(base.variant.to_owned(), base);
                        }
                    }
                }
            }
        }
        None => println!("no bam for pileups"),
    }

    bam_generated.finish();

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

                                // Overwrite any existing variants called by mpileup
                                variant_pos.insert(base.variant.to_owned(), base);
                            }
                        },
                        None => {},
                    }
                } else {
                    panic!("Bug: VCF record reference ids do not match BAM reference ids. Perhaps BAM is unsorted?")
                }
            });
            info!(
                "Collected {} variant positions for sample {}",
                total_records, // remove tmp file name from sample id
                match &stoit_name[..4] {
                    ".tmp" => &stoit_name[15..],
                    _ => &stoit_name,
                },
            );
            if longread {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    (short_sample_count as i32 + *per_ref_sample_idx) as usize,
                    &mut variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            } else {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    *per_ref_sample_idx as usize,
                    &mut variant_map,
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
                    &mut variant_map,
                    &header,
                    &genomes_and_contigs,
                );
            } else {
                variant_matrix.add_sample(
                    stoit_name.to_string(),
                    *per_ref_sample_idx as usize,
                    &mut variant_map,
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
        let tmp_bam_path1 = &(tmp_dir.path().to_str().unwrap().to_string() + "/tmp1.bam");
        let tmp_bam_path2 = &(tmp_dir.path().to_str().unwrap().to_string() + "/tmp2.bam");

        // Generate uncompressed filtered SAM file
        let sam_cmd_string = format!(
            "samtools sort -@ {} {} > {} && \
            gatk AddOrReplaceReadGroups -I {} -O {} -SM 1 -LB N -PL N -PU N && \
            samtools index -@ {} {}",
            threads - 1,
            bam_path,
            tmp_bam_path1,
            tmp_bam_path1,
            tmp_bam_path2,
            threads - 1,
            tmp_bam_path2,
        );
        debug!("Queuing cmd_string: {}", sam_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&sam_cmd_string)
                .stderr(std::process::Stdio::piped())
                .stdout(std::process::Stdio::piped())
                .spawn()
                .expect("Unable to execute bash"),
            "samtools",
        );

        // Variant calling pipeline adapted from Snippy but without all of the rewriting of BAM files
        let vcf_cmd_string = format!(
            "set -e -o pipefail;  \
            gatk HaplotypeCaller -I {} -R {} -O {} --native-pair-hmm-threads {} --sample-ploidy {} -mbq {} \
            --annotation AlleleFraction --annotation DepthPerAlleleBySample --minimum-mapping-quality {} \
            --heterozygosity {} --indel-heterozygosity {} \
            --pcr-indel-model CONSERVATIVE \
            --base-quality-score-threshold 6 --max-reads-per-alignment-start 0 --force-call-filtered-alleles false",
            tmp_bam_path2,
            &reference,
            &vcf_path_prenormalization,
            threads,
            m.value_of("ploidy").unwrap(),
            m.value_of("base-quality-threshold").unwrap(),
            mapq_thresh,
            m.value_of("heterozygosity").unwrap(),
            m.value_of("indel-heterozygosity").unwrap(),
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
