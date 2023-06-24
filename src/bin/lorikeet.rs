use clap::crate_version;
use lorikeet_genome::cli::*;
use lorikeet_genome::external_command_checker;
use lorikeet_genome::utils::utils::*;
use lorikeet_genome::bam_parsing::bam_generator::*;
use lorikeet_genome::processing::lorikeet_engine::{
    run_summarize, start_lorikeet_engine, ReadType
};
use lorikeet_genome::reference::reference_reader_utils::{ReferenceReaderUtils, GenomesAndContigs};
use lorikeet_genome::utils::errors::BirdToolError;
use lorikeet_genome::bam_parsing::FlagFilter;

use log::{info, warn};
use std::env;
use tempfile::NamedTempFile;
use clap_complete::{generate, Shell};
use log::LevelFilter;
use env_logger::Builder;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("summarise") => {
            let m = matches.subcommand_matches("summarise").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, summarise_full_help());
            rayon::ThreadPoolBuilder::new()
                .num_threads(*m.get_one::<usize>("threads").unwrap())
                .build_global()
                .unwrap();
            run_summarize(m);
        }
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, genotype_full_help());
            let mode = "genotype";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Genotype complete."),
                Err(e) => warn!("Genotype failed with error: {:?}", e),
            };
        }
        Some("call") => {
            let m = matches.subcommand_matches("call").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, call_full_help());
            let mode = "call";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Call complete."),
                Err(e) => warn!("Call failed with error: {:?}", e),
            };
        }
        Some("consensus") => {
            let m = matches.subcommand_matches("consensus").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, consensus_full_help());
            let mode = "consensus";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Consensus complete."),
                Err(e) => warn!("Consensus failed with error: {:?}", e),
            };
        }
        Some("shell-completion") => {
            let m = matches.subcommand_matches("shell-completion").unwrap();
            set_log_level(m, true);
            let mut file = std::fs::File::create(m.get_one::<String>("output-file").unwrap())
                .expect("failed to open output file");

            if let Some(generator) = m.get_one::<Shell>("shell").copied() {
                let mut cmd = build_cli();
                info!("Generating completion script for shell {}", generator);
                let name = cmd.get_name().to_string();
                generate(generator, &mut cmd, name, &mut file);
            }
        }
        _ => {
            app.print_help().unwrap();
        }
    }
}

fn prepare_pileup(m: &clap::ArgMatches, mode: &str) -> Result<(), BirdToolError> {
    // This function is amazingly painful. It handles every combination of longread and short read
    // mapping or bam file reading. Could not make it smaller using dynamic or static dispatch
    set_log_level(m, true);
    let filter_params = FilterParameters::generate_from_clap(m);
    let threads = *m.get_one::<usize>("threads").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let references = ReferenceReaderUtils::parse_references(m);
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

    // Temp directory that will house all cached bams for variant calling
    let tmp_dir = match m.contains_id("bam-file-cache-directory") {
        false => {
            let tmp_direct = tempdir::TempDir::new("lorikeet_fifo")
                .expect("Unable to create temporary directory");
            // debug!("Temp directory {}", tmp_direct.as_ref().to_str().unwrap());
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

    let (concatenated_genomes, genomes_and_contigs_option) =
        ReferenceReaderUtils::setup_genome_fasta_files(m);
    // debug!("Found genomes_and_contigs {:?}", genomes_and_contigs_option);
    if m.contains_id("bam-files") {
        let bam_files: Vec<&str> = m.get_many::<String>("bam-files").unwrap().map(|s| &**s).collect();

        // Associate genomes and contig names, if required
        if filter_params.doing_filtering() {
            let bam_readers = generate_filtered_bam_readers_from_bam_files(
                bam_files,
                filter_params.flag_filters.clone(),
                filter_params.min_aligned_length_single,
                filter_params.min_percent_identity_single,
                filter_params.min_aligned_percent_single,
                filter_params.min_aligned_length_pair,
                filter_params.min_percent_identity_pair,
                filter_params.min_aligned_percent_pair,
            );

            if m.contains_id("longread-bam-files") {
                let bam_files = m.get_many::<String>("longread-bam-files").unwrap().map(|s| &**s).collect();
                let long_readers =
                    generate_named_bam_readers_from_bam_files(bam_files);
                run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            } else if m.contains_id("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_generators),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else {
                return run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            }
        } else {
            let bam_readers = generate_named_bam_readers_from_bam_files(bam_files);

            if m.contains_id("longread-bam-files") {
                let bam_files = m.get_many::<String>("longread-bam-files").unwrap().map(|s| &**s).collect();
                let long_readers =
                    generate_named_bam_readers_from_bam_files(bam_files);
                run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            } else if m.contains_id("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_generators),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else {
                return run_pileup(
                    m,
                    mode,
                    bam_readers,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            }
        }
    } else {
        let mapping_program = parse_mapping_program(m.get_one::<String>("mapper").map(|s| &**s));
        external_command_checker::check_for_samtools();

        if filter_params.doing_filtering() {
            // debug!("Filtering..");
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
            // debug!("Finished collecting generators.");
            if m.contains_id("longread-bam-files") {
                let bam_files = m.get_many::<String>("longread-bam-files").unwrap().map(|s| &**s).collect();
                let long_readers =
                    generate_named_bam_readers_from_bam_files(bam_files);
                run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            } else if m.contains_id("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_generators),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else {
                return run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            }
        } else {
            // debug!("Not filtering..");
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

            if m.contains_id("longread-bam-files") {
                let bam_files = m.get_many::<String>("longread-bam-files").unwrap().map(|s| &**s).collect();
                let long_readers =
                    generate_named_bam_readers_from_bam_files(bam_files);

                run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                )
            } else if m.contains_id("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_generators),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else {
                return run_pileup(
                    m,
                    mode,
                    all_generators,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            }
        }
    }
}

fn run_pileup<
    'a,
    R: NamedBamReader,
    S: NamedBamReaderGenerator<R>,
    T: NamedBamReader,
    U: NamedBamReaderGenerator<T>,
>(
    m: &clap::ArgMatches,
    mode: &str,
    bam_readers: Vec<S>,
    flag_filters: FlagFilter,
    long_readers: Option<Vec<U>>,
    genomes_and_contigs_option: Option<GenomesAndContigs>,
    tmp_bam_file_cache: Option<tempdir::TempDir>,
    concatenated_genomes: Option<NamedTempFile>,
) -> Result<(), BirdToolError> {
    let genomes_and_contigs = genomes_and_contigs_option.unwrap();

    start_lorikeet_engine(
        m,
        bam_readers,
        long_readers,
        mode,
        flag_filters,
        genomes_and_contigs,
        tmp_bam_file_cache,
        concatenated_genomes,
    )?;
    Ok(())
}

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.get_flag("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.get_flag("quiet") {
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
