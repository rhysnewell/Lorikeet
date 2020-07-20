use std;
use std::collections::{HashMap, HashSet};
use glob::glob;

use external_command_checker;
use estimation::variant_matrix::*;
use estimation::vcfs::process_vcf::*;
use estimation::bams::process_bam::*;
use estimation::codon_structs::*;
use coverm::bam_generator::*;
use bird_tool_utils::{command};
use utils::*;

use crate::*;
use std::str;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use bio::io::gff;
use tempdir::TempDir;
use std::process::Stdio;
use coverm::genome_parsing::read_genome_fasta_files;
use rayon::prelude::*;
use std::sync::Mutex;
use std::path::Path;

#[allow(unused)]
pub fn pileup_variants<R: NamedBamReader,
    G: NamedBamReaderGenerator<R>,
    S: NamedBamReader,
    U: NamedBamReaderGenerator<S>>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    longreads: Option<Vec<U>>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
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

    let references = match m.values_of("genome-fasta-files") {
        Some(vec) => {
            let reference_paths = vec.map(|p| p.to_string()).collect::<Vec<String>>();
            debug!("Reference files {:?}", reference_paths);
            reference_paths
        },
        None => {
            match m.value_of("genome-fasta-directory") {
                Some(path) => {
                    let ext = m.value_of("genome-fasta-extension").unwrap();
                    let reference_glob = format!("{}/*.{}", path, ext);
                    let reference_paths = glob(&reference_glob).expect("Failed to read cache")
                        .map(|p| p.expect("Failed to read cached bam path")
                            .to_str().unwrap().to_string()).collect::<Vec<String>>();
                    debug!("Reference files {:?}", reference_paths);
                    reference_paths

                }
                None => panic!("Can't find suitable references for variant calling")
            }
        }
    };
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();
    let genomes_and_contigs = read_genome_fasta_files(&references);

    let mut sample_count = bam_readers.len();
    let mut reference_count = references.len();
    let mut coverage_estimators = coverage_estimators;
    let mut ani = 0.;

    // Gff map lists coding regions
    let mut gff_map = HashMap::new();
    let mut codon_table = CodonTable::setup();

    let longreads = match longreads {
        Some(vec) => {
            sample_count += vec.len();
            vec
        },
        None => {
            vec!()
        }
    };

    let mut variant_matrix_map = HashMap::new();
    // Put reference index in the variant map and initialize matrix
    let mut reference_map = HashMap::new();
    for reference in references.iter() {
        debug!("Genomes {:?} contigs {:?}", &genomes_and_contigs.genomes, &genomes_and_contigs.contig_to_genome);
        let ref_idx = genomes_and_contigs.genome_index(&Path::new(reference)
            .file_stem().expect("problem determining file stem").to_str().unwrap().to_string()).unwrap();
        reference_map.entry(ref_idx).or_insert(reference.to_string());
        variant_matrix_map.entry(ref_idx).or_insert(
            VariantMatrix::new_matrix(sample_count / references.len()));
    }
    info!("{} Longread BAM files and {} Shortread BAM files {} Total BAMs over {} genome(s)",
          longreads.len(),
          bam_readers.len(), sample_count, reference_count);

    match mode {
        "evolve" => {

            codon_table.get_codon_table(11);
            ani = 0.;

            let mut gff_reader;
            if m.is_present("gff") || !m.is_present("gff") {
                external_command_checker::check_for_prokka();

                // create new fifo and give read, write and execute rights to the owner.
                let gff_dir = TempDir::new("lorikeet-prokka")
                    .expect("unable to create prokka directory");
                for reference in references.iter() {
                    let cmd_string = format!(
                    "set -e -o pipefail; \
                     prokka --outdir {} --prefix {} --force {} {}",
                    // prodigal
                    gff_dir.path().to_str()
                        .expect("Failed to convert tempfile path to str"),
                    &reference,
                    m.value_of("prokka-params").unwrap_or(""),
                    &reference);
                    info!("Queuing cmd_string: {}", cmd_string);
                    command::finish_command_safely(
                        std::process::Command::new("bash")
                            .arg("-c")
                            .arg(&cmd_string)
                            .stdout(Stdio::piped())
                            .stderr(Stdio::piped())
                            .spawn()
                            .expect("Unable to execute bash"), "prokka");

                    // Read in newly created gff
                    gff_reader = gff::Reader::from_file(format!("{}/{}.gff", gff_dir.path().to_str()
                        .expect("Failed to convert tempfile path to str"), &reference),
                                                        bio::io::gff::GffType::GFF3)
                        .expect("Failed to read prodigal output");

                    // Map to reference id
                    gff_reader.records().into_iter().for_each(|record| {
                        match record {
                            Ok(rec) => {
                                let gff_ref = gff_map.entry(
                                    genomes_and_contigs.genome_index(
                                        &reference.to_string()).unwrap()).or_insert(HashMap::new());
                                let contig_genes = gff_ref.entry(rec.seqname().to_owned())
                                    .or_insert(Vec::new());
                                contig_genes.push(rec);
                            },
                            _ => {},
                        };
                    });
                }

                gff_dir.close().expect("Failed to close temp directory");

            }

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



    let mut sample_groups = HashMap::new();

    info!("Running SNP calling on {} shortread samples", bam_readers.len());
    bam_readers.into_iter().enumerate().for_each(|(sample_idx, bam_generator)|{
        // Get the appropriate sample index based on how many references we are using
        let sample_idx = if references.len() == 1 {
            sample_idx
        } else {
            sample_idx % references.len()
        };
        process_vcf(bam_generator,
                    n_threads,
                    sample_idx,
                    sample_count / references.len(),
                    &mut variant_matrix_map,
                    false,
                    m,
                    &mut sample_groups,
                    &genomes_and_contigs,
                    &reference_map)
    });

    if m.is_present("include-longread-svs") && (m.is_present("longreads") | m.is_present("longread-bam-files")){
//        long_threads = std::cmp::max(n_threads / longreads.len(), 1);
        info!("Running structural variant detection...");
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            // Get the appropriate sample index based on how many references we are using
            let sample_idx = if references.len() == 1 {
                sample_idx
            } else {
                sample_idx % references.len()
            };
            process_vcf(bam_generator,
                        n_threads,
                        sample_idx,
                        sample_count / references.len(),
                        &mut variant_matrix_map,
                        true,
                        m,
                        &mut sample_groups, &genomes_and_contigs, &reference_map)
        });
    } else if m.is_present("longreads") | m.is_present("longread-bam-files") {
        // We need update the variant matrix anyway
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            let bam_generated = bam_generator.start();
            let header = bam_generated.header().clone(); // bam header
            let target_names = header.target_names(); // contig names
            let mut variant_map = HashMap::new();

            let stoit_name = bam_generated.name().to_string().replace("/", ".");
            let group = sample_groups.entry("long").or_insert(HashSet::new());
            group.insert(stoit_name.clone());

            let reference_stem = genomes_and_contigs.genome_of_contig(
                &str::from_utf8(&target_names[0]).unwrap().to_string()).unwrap();
            let ref_idx = genomes_and_contigs.genome_index(&reference_stem).unwrap();

            // Adjust sample index based on number of references
            let sample_idx = sample_idx % references.len();
            let sample_count = sample_count / references.len();
            // Longread adjusted sample index
            let sample_idx_l = sample_count - sample_idx - 1;

            let mut variant_matrix = variant_matrix_map.entry(ref_idx).or_insert(VariantMatrix::new_matrix(sample_count));
            variant_matrix.
                add_sample(stoit_name, sample_idx_l, &variant_map, &header);

        });

    }


    // Annoyingly read in bam file again
    let mut bam_readers = vec!();

    // This is going to catch cached longread bam files from mapping
    if m.is_present("bam-files") {
        let bam_paths = m.values_of("bam-files").unwrap().collect::<Vec<&str>>();
        bam_readers = generate_named_bam_readers_from_bam_files(bam_paths);

    } else {

        let mut all_bam_paths = vec!();

        for ref_name in references.iter() {
            let cache = format!("{}/{}*.bam",
                                m.value_of("bam-file-cache-directory").unwrap().to_string(),
                                ref_name);
            let bam_paths = glob(&cache).expect("Failed to read cache")
                .map(|p| p.expect("Failed to read cached bam path")
                    .to_str().unwrap().to_string()).collect::<Vec<String>>();
            all_bam_paths.extend(bam_paths);
        }
        let all_bam_paths = all_bam_paths.iter().map(|p| p.as_str()).collect::<Vec<&str>>();
        let bam_cnts = all_bam_paths.len();
        bam_readers = generate_named_bam_readers_from_bam_files(all_bam_paths);
        if m.value_of("bam-file-cache-directory").unwrap() == "/tmp/lorikeet-bam-cache" {
            info!("Removing {} cached BAM files...", bam_cnts);
            let remove_directory_cmd = format!("rm -r {}", m.value_of("bam-file-cache-directory").unwrap());
            command::finish_command_safely(
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(remove_directory_cmd)
                    .spawn()
                    .expect("Unable to remove cached BAM files..."),
                "rm"
            );
        }
    }

    // Process Short Read BAMs
    if bam_readers.len() > 0 {
        info!("Performing guided variant calling...");
        bam_readers.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {

            // Get the appropriate sample index based on how many references we are using
            let sample_idx = if references.len() == 1 {
                sample_idx
            } else {
                sample_idx % references.len()
            };
            process_bam(bam_generator,
                        sample_idx,
                        sample_count / references.len(),
                        &mut coverage_estimators,
                        &mut variant_matrix_map,
                        &mut gff_map,
                        n_threads,
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
                        mapq_threshold, method, &sample_groups, &genomes_and_contigs, &reference_map)
        });
    }

    // Process Long Read BAMs if they are present
    if m.is_present("longread-bam-files") {

        let longreads_path = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>();
        let longreads = generate_named_bam_readers_from_bam_files(longreads_path);
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            // Get the appropriate sample index based on how many references we are using
            let sample_idx = if references.len() == 1 {
                sample_idx
            } else {
                sample_idx % references.len()
            };
            process_bam(bam_generator,
                        sample_idx,
                        sample_count / references.len(),
                        &mut coverage_estimators,
                        &mut variant_matrix_map,
                        &mut gff_map,
                        n_threads,
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
                        mapq_threshold, method, &sample_groups, &genomes_and_contigs, &reference_map)
        });
    }

    let gff_map = Mutex::new(gff_map);
    variant_matrix_map.par_iter_mut().for_each(|(ref_idx, variant_matrix)|{
        if mode == "genotype" {
            variant_matrix.generate_distances(n_threads, output_prefix);
            let e_min: f64 = m.value_of("e-min").unwrap().parse().unwrap();
            let e_max: f64 = m.value_of("e-max").unwrap().parse().unwrap();
            let pts_min: f64 = m.value_of("pts-min").unwrap().parse().unwrap();
            let pts_max: f64 = m.value_of("pts-max").unwrap().parse().unwrap();
            let phi: f64 = m.value_of("phi").unwrap().parse().unwrap();
            let anchor_size: usize = m.value_of("minimum-seed-size").unwrap().parse().unwrap();
            let anchor_similarity: f64 = m.value_of("maximum-seed-similarity").unwrap().parse().unwrap();
            let minimum_reads_in_link: usize = m.value_of("minimum-reads-in-link").unwrap().parse().unwrap();
            let reference_path = reference_map.get(&ref_idx).expect("Unable to retrieve reference path");
            let mut reference = match bio::io::fasta::IndexedReader::from_file(&Path::new(&reference_path)) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(&reference_path),
            };

            variant_matrix.run_fuzzy_scan(e_min, e_max, pts_min, pts_max, phi,
                                          anchor_size, anchor_similarity, minimum_reads_in_link);
            variant_matrix.generate_genotypes(&genomes_and_contigs.genomes[*ref_idx], &mut reference);
            variant_matrix.write_vcf(&genomes_and_contigs.genomes[*ref_idx]);
            if m.is_present("plot") {
                let window_size = m.value_of("window-size").unwrap().parse().unwrap();
                variant_matrix.print_variant_stats(&genomes_and_contigs.genomes[*ref_idx], window_size);
            }
        } else if mode == "summarize" {
            let window_size = m.value_of("window-size").unwrap().parse().unwrap();
            variant_matrix.write_vcf(&genomes_and_contigs.genomes[*ref_idx]);
            variant_matrix.print_variant_stats(&genomes_and_contigs.genomes[*ref_idx], window_size);
        } else if mode == "evolve" {
            let reference_path = reference_map.get(&ref_idx).expect("Unable to retrieve reference path");
            let mut reference = match bio::io::fasta::IndexedReader::from_file(&Path::new(&reference_path)) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(&reference_path),
            };
            let mut gff_map = gff_map.lock().unwrap();
            let mut gff_ref = gff_map.get_mut(ref_idx).expect(&format!("No GFF records for reference {:?}", reference_path));
            variant_matrix.calc_gene_mutation(&mut gff_ref, &mut reference, &codon_table,
                                              Path::new(&reference_path).file_stem().unwrap().to_str().unwrap(), &genomes_and_contigs.genomes[*ref_idx])
        }
    });
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