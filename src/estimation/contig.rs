use std;
use std::collections::HashMap;
use glob::glob;

use external_command_checker;
use estimation::variant_matrix::*;
use estimation::vcfs::process_vcf::*;
use estimation::bams::process_bam::*;
use estimation::codon_structs::*;
use coverm::bam_generator::*;
use bird_tool_utils::command;

use crate::*;
use std::str;
use std::fs::File;
use std::path::Path;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use bio::io::gff;
use tempdir::TempDir;
use std::process::Stdio;

#[allow(unused)]
pub fn pileup_variants<R: NamedBamReader + Send,
    G: NamedBamReaderGenerator<R> + Send,
    S: NamedBamReader + Send,
    U: NamedBamReaderGenerator<S> + Send,>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    long_readers: Option<Vec<U>>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    reference: bio::io::fasta::IndexedReader<File>,
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
    let mut reference = reference;
    let mut coverage_estimators = coverage_estimators;
    let mut ani = 0.;

    // Gff map lists coding regions
    let mut gff_map = HashMap::new();
    let mut codon_table = CodonTable::setup();

    // Get long reads bams if they exist
    let longreads = match long_readers {
        Some(vec) => {
            // update sample count
            debug!("Longread bams {:?}", vec.len());

            sample_count += vec.len();
            vec
        },
        _ => vec!(),
    };


    let mut variant_matrix = VariantMatrix::new_matrix(sample_count);

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

                // create new fifo and give read, write and execute rights to the owner.
                let gff_dir = TempDir::new("lorikeet-prokka")
                    .expect("unable to create prokka directory");

                let cmd_string = format!(
                    "set -e -o pipefail; \
                     prokka --outdir {} --prefix prokka_out --force {} {}",
                    // prodigal
                    gff_dir.path().to_str()
                        .expect("Failed to convert tempfile path to str"),
                    m.value_of("prokka-params").unwrap_or(""),
                    m.value_of("reference").unwrap());
                info!("Queuing cmd_string: {}", cmd_string);
                command::finish_command_safely(
                    std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&cmd_string)
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .expect("Unable to execute bash"), "prokka");

                gff_reader = gff::Reader::from_file(format!("{}/prokka_out.gff", gff_dir.path().to_str()
                                                        .expect("Failed to convert tempfile path to str")),
                                                    bio::io::gff::GffType::GFF3)
                    .expect("Failed to read prodigal output");

                gff_dir.close().expect("Failed to close temp directory");

            }
            gff_reader.records().into_iter().for_each(|record| {
                match record {
                    Ok(rec) => {
                        let contig_genes = gff_map.entry(rec.seqname().to_owned())
                            .or_insert(Vec::new());
                        contig_genes.push(rec);
                    },
                    _ => {},
                };
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

    // Get the total amound of bases in reference using the index
    let reference_length = bio::io::fasta::Index::with_fasta_file(
        &Path::new(m.value_of("reference").unwrap())).unwrap()
        .sequences().iter().fold(0, |acc, seq| acc + seq.len);

    info!("Running SNP calling on {} shortread samples", bam_readers.len());
    bam_readers.into_iter().enumerate().for_each(|(sample_idx, bam_generator)|{
        process_vcf(bam_generator,
                    n_threads,
                    sample_idx,
                    sample_count,
                    &mut variant_matrix,
                    false,
                    m,
                    reference_length)
    });

    if longreads.len() > 0 && m.is_present("include-longread-svs"){
//        long_threads = std::cmp::max(n_threads / longreads.len(), 1);
        info!("Running structural variant detection on {} longread samples", longreads.len());
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            process_vcf(bam_generator,
                        n_threads,
                        sample_idx,
                        sample_count,
                        &mut variant_matrix,
                        true,
                        m,
                        reference_length)
        });
    } else if longreads.len() > 0 {
        // We need update the variant matrix anyway
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            let bam_generated = bam_generator.start();
            let header = bam_generated.header().clone(); // bam header
            let mut variant_map = HashMap::new();

            let stoit_name = bam_generated.name().to_string();
            // Longread adjusted sample index
            let sample_idx_l = sample_count - sample_idx - 1;
            variant_matrix.
                add_sample(stoit_name.clone(), sample_idx_l, &variant_map, &header);

        });

    }


    // Annoyingly read in bam file again
    let mut longreads = vec!();
    let mut bam_readers = vec!();

    if m.is_present("bam-files") {
        let bam_paths = m.values_of("bam-files").unwrap().collect::<Vec<&str>>();
        bam_readers = generate_named_bam_readers_from_bam_files(bam_paths);

    } else {
        let ref_name = Path::new(m.value_of("reference").unwrap())
            .file_name().unwrap()
            .to_str().unwrap();
        let cache = format!("{}/{}*.bam",
                            m.value_of("bam-file-cache-directory").unwrap().to_string(),
                            ref_name);
        let bam_paths = glob(&cache).expect("Failed to read cache")
            .map(|p| p.expect("Failed to read cached bam path")
                .to_str().unwrap().to_string()).collect::<Vec<String>>();
        let bam_paths = bam_paths.iter().map(|p| &**p).collect::<Vec<&str>>();
        let bam_cnts = bam_paths.len();
        bam_readers = generate_named_bam_readers_from_bam_files(bam_paths);
        if m.value_of("bam-file-cache-directory").unwrap() == "/tmp" {
            info!("Removing {} cached BAM files...", bam_cnts);
            let remove_directory_cmd = format!("rm {}", cache);
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
            process_bam(bam_generator,
                        sample_idx,
                        sample_count,
                        &mut reference,
                        &mut coverage_estimators,
                        &mut variant_matrix,
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
                        mapq_threshold, method, false)
        });
    }

    // Process Long Read BAMs if they are present
    if m.is_present("longread-bam-files") {
        let longreads_path = m.values_of("longread-bam-files").unwrap().collect::<Vec<&str>>();
        longreads = generate_named_bam_readers_from_bam_files(longreads_path);
    }
    if longreads.len() > 0 {
//        long_threads = std::cmp::max(n_threads / longreads.len(), 1);
        longreads.into_iter().enumerate().for_each(|(sample_idx, bam_generator)| {
            process_bam(bam_generator,
                        sample_idx,
                        sample_count,
                        &mut reference,
                        &mut coverage_estimators,
                        &mut variant_matrix,
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
                        mapq_threshold, method, true)
        });
    }

    if mode=="genotype" {
        variant_matrix.generate_distances(n_threads, output_prefix);
        let e_min: f64 = m.value_of("e-min").unwrap().parse().unwrap();
        let e_max: f64 = m.value_of("e-max").unwrap().parse().unwrap();
        let pts_min: f64 = m.value_of("pts-min").unwrap().parse().unwrap();
        let pts_max: f64 = m.value_of("pts-max").unwrap().parse().unwrap();
        let phi: f64 = m.value_of("phi").unwrap().parse().unwrap();
        let anchor_size: usize = m.value_of("minimum-seed-size").unwrap().parse().unwrap();
        let anchor_similarity: f64 = m.value_of("maximum-seed-similarity").unwrap().parse().unwrap();
        let minimum_reads_in_link: usize = m.value_of("minimum-reads-in-link").unwrap().parse().unwrap();


        variant_matrix.run_fuzzy_scan(e_min, e_max, pts_min, pts_max, phi,
                                      anchor_size, anchor_similarity, minimum_reads_in_link);
        variant_matrix.generate_genotypes(output_prefix, &mut reference);
        variant_matrix.write_vcf(output_prefix);
        if m.is_present("plot") {
            let window_size = m.value_of("window-size").unwrap().parse().unwrap();
            variant_matrix.print_variant_stats(output_prefix, window_size);
        }
    } else if mode=="summarize" {
        let window_size = m.value_of("window-size").unwrap().parse().unwrap();
        variant_matrix.write_vcf(output_prefix);
        variant_matrix.print_variant_stats(output_prefix, window_size);
    } else if mode=="evolve" {
        let ref_name = Path::new(m.value_of("reference").unwrap())
            .file_name().unwrap()
            .to_str().unwrap();
        variant_matrix.calc_gene_mutation(&mut gff_map, &mut reference, &codon_table,
                                          ref_name, output_prefix)
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