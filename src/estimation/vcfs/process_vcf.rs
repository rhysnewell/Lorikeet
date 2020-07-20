use std;
use std::collections::{HashMap, HashSet};
use rust_htslib::{bam, bcf, bcf::Read};
use bird_tool_utils::command;

use external_command_checker;
use estimation::variant_matrix::*;
use coverm::bam_generator::*;
use model::variants::*;
use utils::*;

use crate::*;
use std::str;
use std::path::Path;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use coverm::genomes_and_contigs::GenomesAndContigs;


#[allow(unused)]
pub fn process_vcf<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>>(
    bam_generator: G,
    split_threads: usize,
    sample_idx: usize,
    sample_count: usize,
    variant_matrix_map: &mut HashMap<usize, VariantMatrix>,
    longread: bool,
    m: &clap::ArgMatches,
    sample_groups: &mut HashMap<&str, HashSet<String>>,
    genomes_and_contigs: &GenomesAndContigs,
    reference_map: &HashMap<usize, String>) {

    let mut bam_generated = bam_generator.start();

    let stoit_name = bam_generated.name().to_string().replace("/", ".");

    if longread {
        let group = sample_groups.entry("long").or_insert(HashSet::new());
        group.insert(stoit_name.clone());
    } else {
        let group = sample_groups.entry("short").or_insert(HashSet::new());
        group.insert(stoit_name.clone());
    }

    debug!("Setting threads...");
    bam_generated.set_threads(split_threads);
    debug!("Managed to set threads.");
    let header = bam_generated.header().clone(); // bam header
    let target_names = header.target_names(); // contig names

    debug!("retrieving genome id with contig {:?}", str::from_utf8(&target_names[0]));

    // Get the total amound of bases in reference using the index
    let reference_stem = genomes_and_contigs.genome_of_contig(
        &str::from_utf8(&target_names[0]).unwrap().to_string())
        .expect(&format!("Found invalid contig in bam, {:?}. Please provide corresponding reference genomes", str::from_utf8(&target_names[0]).unwrap()));
    let ref_idx = genomes_and_contigs.genome_index(&reference_stem).unwrap();
    let reference = reference_map.get(&ref_idx).expect("Unable to retrieve reference path");

    let reference_length = match bio::io::fasta::Index::with_fasta_file(&Path::new(&reference)) {
        Ok(index) => index.sequences().iter().fold(0, |acc, seq| acc + seq.len),
        Err(_e) => {
            generate_faidx(&reference);
            bio::io::fasta::Index::with_fasta_file(&Path::new(&reference)).unwrap().sequences().iter().fold(0, |acc, seq| acc + seq.len)
        }
    };

    // Get the variant matrix for current reference
    let mut variant_matrix = variant_matrix_map.entry(ref_idx).or_insert(
        VariantMatrix::new_matrix(sample_count));

    // for each genomic position, only has hashmap when variants are present. Includes read ids
    let mut variant_map = HashMap::new();


    // Get VCF file from BAM using freebayes of SVIM
    let mut vcf_reader = get_vcf(&stoit_name,
                                 &m,
                                 sample_idx,
                                 split_threads,
                                 longread,
                                 reference_length,
                                 reference);
    let mut sample_idx = sample_idx;
    if longread {
        sample_idx = sample_count - sample_idx - 1;
    }
    match vcf_reader {
        Ok(ref mut reader) => {
            reader.set_threads(split_threads).expect("Unable to set threads on VCF reader");

            let min_qual = m.value_of("min-variant-quality").unwrap().parse().unwrap();
            info!("Collecting VCF records for sample {} against {}", sample_idx, &reference_stem);
            reader.records().into_iter().for_each(|vcf_record| {
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
                                                            longread,
                                                            min_qual);
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
            variant_matrix.add_sample(stoit_name.to_string(), sample_idx, &variant_map, &header);
        },
        Err(_) => {
            info!("No VCF records found for sample {} against {}", sample_idx, &reference_stem);
            variant_matrix.add_sample(stoit_name.to_string(), sample_idx, &variant_map, &header);
        }
    }

}



/// Get or generate vcf file
#[allow(unused)]
pub fn get_vcf(stoit_name: &str, m: &clap::ArgMatches,
               sample_idx: usize, threads: usize, longread: bool,
               reference_length: u64, reference: &String) -> std::result::Result<bcf::Reader, rust_htslib::bcf::Error> {
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
                return generate_vcf(bam_path, m, threads, longread, reference_length, reference)
            } else {
                let bam_path: &str = m.values_of("bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
                return generate_vcf(bam_path, m, threads, longread, reference_length, reference)
            }
        } else {
            let vcf_path = vcf_path[0];
            let vcf = bcf::Reader::from_path(vcf_path);
            return vcf
        }
    } else if longread && m.is_present("longread-bam-files") {
        let bam_path: &str = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
        return generate_vcf(bam_path, m, threads, longread, reference_length, reference)
    } else if m.is_present("bam-files") {
        let bam_path: &str = m.values_of("bam-files").unwrap().collect::<Vec<&str>>()[sample_idx];
        return generate_vcf(bam_path, m, threads, longread, reference_length, reference)
    } else {
        // We are streaming a generated bam file, so we have had to cache the bam for this to work
        let cache = m.value_of("bam-file-cache-directory").unwrap().to_string() + "/";
        let stoit_name: Vec<&str> = stoit_name.split("/").collect();
        let stoit_name = stoit_name.join(".");
        let stoit_name = stoit_name.replace(|c: char| !c.is_ascii(), "");
        let bam_path = cache + &(stoit_name + ".bam");

        return generate_vcf(&bam_path, m, threads, longread, reference_length, reference)
    }

}

/// Makes direct call to freebayes or SVIM
#[allow(unused)]
pub fn generate_vcf(bam_path: &str, m: &clap::ArgMatches,
                    threads: usize, longread: bool, reference_length: u64, reference: &String) -> std::result::Result<bcf::Reader, rust_htslib::bcf::Error> {

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
        external_command_checker::check_for_samtools();
        external_command_checker::check_for_vt();
        external_command_checker::check_for_bcftools();

        let region_size = reference_length / threads as u64;

        let index_path = reference.clone() + ".fai";

        let freebayes_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/freebayes.vcf");
//        let freebayes_path = &("freebayes.vcf");
        let tmp_bam_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/tmp.bam");

        // Generate uncompressed filtered SAM file
        let sam_cmd_string = format!(
            "samtools sort -@ {} -n -l 0 -T /tmp {} | \
            samtools fixmate -@ {} -m - - | \
            samtools sort -@ {} -l 0 -T /tmp | \
            samtools markdup -@ {} -T /tmp -r -s - - > {}",
            threads-1,
            bam_path,
            threads-1,
            threads-1,
            threads-1,
            tmp_bam_path);
        debug!("Queuing cmd_string: {}", sam_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&sam_cmd_string)
                .stderr(std::process::Stdio::piped())
                .stdout(std::process::Stdio::piped())
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
            &reference,
            m.value_of("min-variant-depth").unwrap(),
            m.value_of("base-quality-threshold").unwrap(),
            m.value_of("min-repeat-entropy").unwrap(),
            m.value_of("mapq-threshold").unwrap(),
            tmp_bam_path,
            &reference,
            freebayes_path);
        debug!("Queuing cmd_string: {}", vcf_cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&vcf_cmd_string)
                .stderr(std::process::Stdio::piped())
                .stdout(std::process::Stdio::piped())
                .spawn()
                .expect("Unable to execute bash"), "freebayes");
        debug!("VCF Path {:?}", freebayes_path);
        let vcf_reader = bcf::Reader::from_path(&Path::new(freebayes_path));

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader
    } else {
        external_command_checker::check_for_svim();
        let svim_path = tmp_dir.path().to_str().unwrap().to_string();

        // check and build bam index if it doesn't exist
        if !Path::new(&(bam_path.to_string() + ".bai")).exists() {
            bam::index::build(bam_path, Some(&(bam_path.to_string() + ".bai")),
                              bam::index::Type::BAI, threads as u32).expect(
                &format!("Unable to index bam at {}", &bam_path));
        }

        let cmd_string = format!(
            "set -e -o pipefail; svim alignment --read_names --skip_genotyping \
            --tandem_duplications_as_insertions --interspersed_duplications_as_insertions \
            --min_mapq {} --sequence_alleles {} {} {}",
            m.value_of("mapq-threshold").unwrap(),
            &svim_path,
            &bam_path,
            &reference);
        debug!("Queuing cmd_string: {}", cmd_string);
        command::finish_command_safely(
            std::process::Command::new("bash")
                .arg("-c")
                .arg(&cmd_string)
                .stderr(std::process::Stdio::piped())
//                .stdout(std::process::Stdio::null())
                .spawn()
                .expect("Unable to execute bash"), "svim");
        let vcf_path = &(svim_path + "/variants.vcf");
        debug!("VCF Path {:?}", vcf_path);
        let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

        tmp_dir.close().expect("Failed to close temp directory");
        return vcf_reader
    }
}
