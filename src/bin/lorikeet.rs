extern crate lorikeet_genome;

use lorikeet_genome::cli::*;
use lorikeet_genome::estimation::contig;
use lorikeet_genome::external_command_checker;
use lorikeet_genome::utils::*;
use lorikeet_genome::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate bio;
use bio::alignment::sparse::*;

extern crate coverm;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use coverm::*;
use coverm::{bam_generator::*, filter};

use std::collections::BTreeMap;
use std::env;
use std::path::Path;
use std::process;
use std::str;

extern crate tempdir;
extern crate tempfile;
use tempfile::NamedTempFile;
extern crate clap;
use clap::*;

extern crate itertools;
use itertools::Itertools;

#[macro_use]
extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            if m.is_present("full-help") {
                println!("{}", filter_full_help());
                process::exit(1);
            }
            set_log_level(m, true);

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let output_bam_files: Vec<&str> = m.values_of("output-bam-files").unwrap().collect();
            if bam_files.len() != output_bam_files.len() {
                panic!("The number of input BAM files must be the same as the number output")
            }

            let filter_params = FilterParameters::generate_from_clap(m);

            let num_threads = value_t!(m.value_of("threads"), u16).unwrap();

            for (bam, output) in bam_files.iter().zip(output_bam_files.iter()) {
                let reader =
                    bam::Reader::from_path(bam).expect(&format!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer =
                    bam::Writer::from_path(output, &header, rust_htslib::bam::Format::BAM)
                        .expect(&format!("Failed to write BAM file {}", output));
                writer
                    .set_threads(num_threads as usize)
                    .expect("Failed to set num threads in writer");
                let mut filtered = filter::ReferenceSortedBamFilter::new(
                    reader,
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    !m.is_present("inverse"),
                );

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        }
        Some("summarize") => {
            let m = matches.subcommand_matches("summarize").unwrap();
            let mode = "summarize";
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        }
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            let mode = "genotype";
            if m.is_present("full-help") {
                println!("{}", genotype_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        }
        Some("evolve") => {
            let m = matches.subcommand_matches("evolve").unwrap();
            let mode = "evolve";
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        }
        Some("polish") => {
            let m = matches.subcommand_matches("polish").unwrap();
            let mode = "polish";
            if m.is_present("full-help") {
                println!("{}", polish_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        }
        Some("kmer") => {
            let m = matches.subcommand_matches("kmer").unwrap();
            if m.is_present("full-help") {
                //                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let reference_path = Path::new(m.value_of("reference").unwrap());
            let fasta_reader = match bio::io::fasta::Reader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(e) => {
                    eprintln!("Missing or corrupt fasta file {}", e);
                    process::exit(1);
                }
            };
            let contigs = fasta_reader.records().collect_vec();
            // Initialize bound contig variable
            let mut tet_freq = BTreeMap::new();
            let contig_count = contigs.len();
            let mut contig_idx = 0 as usize;
            let mut contig_names = vec![String::new(); contig_count];
            for contig in contigs {
                let contig = contig.unwrap();
                contig_names[contig_idx] = contig.id().to_string();
                debug!("Parsing contig: {}", contig.id());
                let kmers = hash_kmers(contig.seq(), 4);
                // Get kmer counts in a contig
                for (kmer, pos) in kmers {
                    let k = tet_freq
                        .entry(kmer.to_vec())
                        .or_insert(vec![0; contig_count]);
                    k[contig_idx] = pos.len();
                }
                contig_idx += 1;
            }
            for (idx, contig) in contig_names.iter().enumerate() {
                print!("{}\t", contig);
                for (_kmer, counts) in &tet_freq {
                    print!("{}\t", counts[idx])
                }
                print!("\n");
            }
        }
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn prepare_pileup(m: &clap::ArgMatches, mode: &str) {
    // This function is amazingly painful. It handles every combination of longread and short read
    // mapping or bam file reading. Could not make it smaller using dynamic or static dispatch
    set_log_level(m, true);
    let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
    let filter_params = FilterParameters::generate_from_clap(m);
    let threads = m.value_of("threads").unwrap().parse().unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let references = parse_references(m);
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

    // Temp directory that will house all cached bams for variant calling
    let tmp_dir = match m.is_present("bam-file-cache-directory") {
        false => {
            let tmp_direct = tempdir::TempDir::new("lorikeet_fifo")
                .expect("Unable to create temporary directory");
            debug!("Temp directory {}", tmp_direct.as_ref().to_str().unwrap());
            std::fs::create_dir(format!("{}/long", &tmp_direct.as_ref().to_str().unwrap()))
                .unwrap();
            std::fs::create_dir(format!("{}/short", &tmp_direct.as_ref().to_str().unwrap()))
                .unwrap();
            std::fs::create_dir(format!(
                "{}/assembly",
                &tmp_direct.as_ref().to_str().unwrap()
            ))
            .unwrap();

            Some(tmp_direct)
        }
        true => None,
    };

    let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
    debug!("Found genomes_and_contigs {:?}", genomes_and_contigs_option);
    if m.is_present("bam-files") {
        let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();

        // Associate genomes and contig names, if required
        if filter_params.doing_filtering() {
            let bam_readers = bam_generator::generate_filtered_bam_readers_from_bam_files(
                bam_files,
                filter_params.flag_filters.clone(),
                filter_params.min_aligned_length_single,
                filter_params.min_percent_identity_single,
                filter_params.min_aligned_percent_single,
                filter_params.min_aligned_length_pair,
                filter_params.min_percent_identity_pair,
                filter_params.min_aligned_percent_pair,
            );

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else if m.is_present("longreads") {
                // Perform mapping
                let mut long_generators = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    bam_readers,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            }
        } else {
            let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_readers),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else if m.is_present("longreads") {
                // Perform mapping
                let mut long_generators = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        bam_readers,
                        filter_params.flag_filters,
                        Some(long_generators),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    bam_readers,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            }
        }
    } else {
        let mapping_program = parse_mapping_program(m.value_of("mapper"));
        external_command_checker::check_for_samtools();

        if filter_params.doing_filtering() {
            debug!("Filtering..");
            let readtype = ReadType::Short;
            let generator_sets = get_streamed_filtered_bam_readers(
                m,
                mapping_program,
                &concatenated_genomes,
                &filter_params,
                &readtype,
                &Some(references.clone()),
                &tmp_dir,
            );
            let mut all_generators = vec![];
            let mut indices = vec![]; // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            debug!("Finished collecting generators.");
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else if m.is_present("longreads") {
                // Perform mapping
                let mut long_generators = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    all_generators,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            }
        } else {
            debug!("Not filtering..");
            let readtype = ReadType::Short;
            let generator_sets = get_streamed_bam_readers(
                m,
                mapping_program,
                &concatenated_genomes,
                &readtype,
                &Some(references.clone()),
                &tmp_dir,
            );
            let mut all_generators = vec![];
            let mut indices = vec![]; // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_readers),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else if m.is_present("longreads") {
                // Perform mapping
                let mut long_generators = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                if m.is_present("assembly-bam-files") {
                    let assembly_bam_files = m.values_of("assembly-bam-files").unwrap().collect();
                    let assembly_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        assembly_bam_files,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_readers),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else if m.is_present("assembly") {
                    // Perform mapping
                    let mut assembly_generators = assembly_generator_setup(
                        &m,
                        &concatenated_genomes,
                        &Some(references.clone()),
                        &tmp_dir,
                    );
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        Some(assembly_generators),
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                } else {
                    run_pileup(
                        m,
                        mode,
                        &mut estimators,
                        all_generators,
                        filter_params.flag_filters,
                        Some(long_generators),
                        None::<Vec<PlaceholderBamFileReader>>,
                        genomes_and_contigs_option,
                        tmp_dir,
                        concatenated_genomes,
                    )
                }
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    all_generators,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            }
        }
    }
}

struct EstimatorsAndTaker {
    estimators: Vec<CoverageEstimator>,
}

impl EstimatorsAndTaker {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> EstimatorsAndTaker {
        let mut estimators = vec![];
        let min_fraction_covered = parse_percentage(&m, "min-covered-fraction");
        let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();

        let methods: Vec<&str> = m.values_of("method").unwrap().collect();

        if doing_metabat(&m) {
            estimators.push(CoverageEstimator::new_estimator_length());
            estimators.push(CoverageEstimator::new_estimator_mean(
                min_fraction_covered,
                contig_end_exclusion,
                false,
            ));
            estimators.push(CoverageEstimator::new_estimator_variance(
                min_fraction_covered,
                contig_end_exclusion,
            ));

            debug!("Cached regular coverage taker for metabat mode being used");
        } else {
            for (_i, method) in methods.iter().enumerate() {
                match method {
                    &"mean" => {
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        )); // TODO: Parameterise exclude_mismatches

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"trimmed_mean" => {
                        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                            error!(
                                "error: Trim bounds must be between 0 and 1, and \
                                 min must be less than max, found {} and {}",
                                min, max
                            );
                            process::exit(1);
                        }
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min,
                            max,
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    _ => unreachable!(),
                };
            }
        }

        // Check that min-covered-fraction is being used as expected
        if min_fraction_covered != 0.0 {
            let die = |estimator_name| {
                error!(
                    "The '{}' coverage estimator cannot be used when \
                     --min-covered-fraction is > 0 as it does not calculate \
                     the covered fraction. You may wish to set the \
                     --min-covered-fraction to 0 and/or run this estimator \
                     separately.",
                    estimator_name
                );
                process::exit(1)
            };
            for e in &estimators {
                match e {
                    CoverageEstimator::ReadCountCalculator { .. } => die("counts"),
                    CoverageEstimator::ReferenceLengthCalculator { .. } => die("length"),
                    CoverageEstimator::ReadsPerBaseCalculator { .. } => die("reads_per_base"),
                    _ => {}
                }
            }
        }

        return EstimatorsAndTaker {
            estimators: estimators,
        };
    }
}

fn run_pileup<
    'a,
    R: NamedBamReader,
    S: NamedBamReaderGenerator<R>,
    T: NamedBamReader,
    U: NamedBamReaderGenerator<T>,
    V: NamedBamReader,
    W: NamedBamReaderGenerator<V>,
>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
    bam_readers: Vec<S>,
    flag_filters: FlagFilter,
    long_readers: Option<Vec<U>>,
    assembly_readers: Option<Vec<W>>,
    genomes_and_contigs_option: Option<GenomesAndContigs>,
    tmp_bam_file_cache: Option<tempdir::TempDir>,
    concatenated_genomes: Option<NamedTempFile>,
) {
    match mode {
        "genotype" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();

            let output_prefix = match m.is_present("output-directory") {
                true => {
                    match std::fs::create_dir_all(
                        m.value_of("output-directory").unwrap().to_string(),
                    ) {
                        Ok(_) => {}
                        Err(err) => panic!(format!("Unable to create output directory {:?}", err)),
                    };
                    m.value_of("output-directory").unwrap()
                }
                false => "./",
            };
            let include_indels = m.is_present("include-indels");
            let include_soft_clipping = m.is_present("include-soft-clipping");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!(
                    "error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}",
                    min, max
                );
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                assembly_readers,
                mode,
                &mut estimators.estimators,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                include_soft_clipping,
                m.is_present("longread-bam-files"),
                genomes_and_contigs,
                tmp_bam_file_cache,
                concatenated_genomes,
            );
        }
        "summarize" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();

            let output_prefix = match m.is_present("output-directory") {
                true => {
                    match std::fs::create_dir_all(
                        m.value_of("output-directory").unwrap().to_string(),
                    ) {
                        Ok(_) => {}
                        Err(err) => panic!(format!("Unable to create output directory {:?}", err)),
                    };
                    m.value_of("output-directory").unwrap()
                }
                false => "./",
            };
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!(
                    "error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}",
                    min, max
                );
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                assembly_readers,
                mode,
                &mut estimators.estimators,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                m.is_present("include-soft-clipping"),
                m.is_present("longread-bam-files"),
                genomes_and_contigs,
                tmp_bam_file_cache,
                concatenated_genomes,
            );
        }
        "evolve" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let output_prefix = match m.is_present("output-directory") {
                true => {
                    match std::fs::create_dir_all(
                        m.value_of("output-directory").unwrap().to_string(),
                    ) {
                        Ok(_) => {}
                        Err(err) => panic!(format!("Unable to create output directory {:?}", err)),
                    };
                    m.value_of("output-directory").unwrap()
                }
                false => "./",
            };
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!(
                    "error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}",
                    min, max
                );
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                assembly_readers,
                mode,
                &mut estimators.estimators,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                m.is_present("include-soft-clipping"),
                m.is_present("longread-bam-files"),
                genomes_and_contigs,
                tmp_bam_file_cache,
                concatenated_genomes,
            );
        }
        "polish" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let output_prefix = match m.is_present("output-directory") {
                true => {
                    match std::fs::create_dir_all(
                        m.value_of("output-directory").unwrap().to_string(),
                    ) {
                        Ok(_) => {}
                        Err(err) => panic!(format!("Unable to create output directory {:?}", err)),
                    };
                    m.value_of("output-directory").unwrap()
                }
                false => "./",
            };

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();
            let include_indels = m.is_present("include-indels");

            //            let index_path = reference_path.clone().to_owned() + ".fai";
            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!(
                    "error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}",
                    min, max
                );
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                assembly_readers,
                mode,
                &mut estimators.estimators,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                &output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                false,
                m.is_present("longread-bam-files"),
                genomes_and_contigs,
                tmp_bam_file_cache,
                concatenated_genomes,
            );
        }
        _ => panic!("Unknown lorikeet mode"),
    }
}

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("lorikeet version {}", crate_version!());
    }
}
