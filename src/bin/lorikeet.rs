extern crate lorikeet_genome;

use lorikeet_genome::*;
use lorikeet_genome::estimation::contig;
use lorikeet_genome::external_command_checker;
use lorikeet_genome::cli::*;
use lorikeet_genome::utils::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate bio;
use bio::alignment::sparse::*;

extern crate coverm;
use coverm::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::shard_bam_reader::*;
use coverm::FlagFilter;
use coverm::genome_exclusion::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::{bam_generator::*, mapping_parameters::*, filter, mapping_index_maintenance};

use std::env;
use std::str;
use std::process;
use std::collections::{BTreeMap, HashSet};
use std::process::Stdio;
use std::fs::File;
use std::io::Write;
use std::path::Path;

extern crate clap;
use clap::*;

extern crate itertools;
use itertools::Itertools;

#[macro_use]
extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;

extern crate tempfile;
use tempfile::NamedTempFile;

fn galah_command_line_definition(
) -> galah::cluster_argument_parsing::GalahClustererCommandDefinition {
    galah::cluster_argument_parsing::GalahClustererCommandDefinition {
        dereplication_ani_argument: "dereplication-ani".to_string(),
        dereplication_prethreshold_ani_argument: "dereplication-prethreshold-ani".to_string(),
        dereplication_quality_formula_argument: "dereplication-quality-formula".to_string(),
        dereplication_precluster_method_argument: "dereplication-precluster-method".to_string(),
    }
}

fn main(){
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
                let reader = bam::Reader::from_path(bam).expect(
                    &format!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer = bam::Writer::from_path(
                    output,
                    &header,
                    rust_htslib::bam::Format::BAM)
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
                    !m.is_present("inverse"));

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        },
        Some("polymorph") => {
            let m = matches.subcommand_matches("polymorph").unwrap();
            let mode = "polymorph";
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        },
        Some("summarize") => {
            let m = matches.subcommand_matches("summarize").unwrap();
            let mode = "summarize";
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        },
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            let mode = "genotype";
            if m.is_present("full-help") {
                println!("{}", genotype_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);

        },
        Some("evolve") => {
            let m = matches.subcommand_matches("evolve").unwrap();
            let mode = "evolve";
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        },
        Some("polish") => {
            let m = matches.subcommand_matches("polish").unwrap();
            let mode = "polish";
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        },
        Some("kmer") => {
            let m = matches.subcommand_matches("kmer").unwrap();
            if m.is_present("full-help") {
//                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let reference_path = Path::new(m.value_of("reference").unwrap());
            let fasta_reader = match bio::io::fasta::Reader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(e) => {
                    eprintln!("Missing or corrupt fasta file {}", e);
                    process::exit(1);
                },
            };
            let contigs = fasta_reader.records().collect_vec();
            // Initialize bound contig variable
            let mut tet_freq = BTreeMap::new();
            let contig_count = contigs.len();
            let mut contig_idx = 0 as usize;
            let mut contig_names = vec![String::new(); contig_count];
            for contig in contigs{
                let contig = contig.unwrap();
                contig_names[contig_idx] = contig.id().to_string();
                debug!("Parsing contig: {}", contig.id());
                let kmers = hash_kmers(contig.seq(), 4);
                // Get kmer counts in a contig
                for (kmer, pos) in kmers {
                    let k = tet_freq.entry(kmer.to_vec()).or_insert(vec![0; contig_count]);
                    k[contig_idx] = pos.len();
                }
                contig_idx += 1;
            }
            for (idx, contig) in contig_names.iter().enumerate() {
                print!("{}\t", contig);
                for (_kmer, counts) in &tet_freq{
                    print!("{}\t", counts[idx])
                }
                print!("\n");
            }

        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

// Enum for exclusion out here so long read can find it
enum GenomeExclusionTypes {
    SeparatorType,
    NoneType,
    GenomesAndContigsType,
}

fn prepare_pileup
    (m: &clap::ArgMatches, mode: &str) {

    // This function is amazingly painful. It handles every combination of longread and short read
    // mapping or bam file reading. Could not make it smaller using dynamic or static dispatch
    set_log_level(m, true);
    let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
    let filter_params = FilterParameters::generate_from_clap(m);
    let threads = m.value_of("threads").unwrap().parse().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    let genome_names_content: Vec<u8>;

    let separator = parse_separator(m);

    let genomes_and_contigs_option_predereplication = if !m.is_present("separator")
        && !m.is_present("dereplicate")
        && !m.is_present("single-genome")
    {
        parse_all_genome_definitions(&m)
    } else {
        None
    };

    // This would be better as a separate function to make this function
    // smaller, but I find this hard because functions cannot return a
    // trait.
    let mut genome_exclusion_filter_separator_type: Option<SeparatorGenomeExclusionFilter> =
        None;
    let mut genome_exclusion_filter_non_type: Option<NoExclusionGenomeFilter> = None;
    let mut genome_exclusion_genomes_and_contigs: Option<GenomesAndContigsExclusionFilter> =
        None;


    let genome_exclusion_type = {
        if m.is_present("sharded") {
            if m.is_present("exclude-genomes-from-deshard") {
                let filename = m.value_of("exclude-genomes-from-deshard").unwrap();
                genome_names_content = std::fs::read(filename).expect(&format!(
                    "Failed to open file '{}' containing list of excluded genomes",
                    filename
                ));
                let mut genome_names_hash: HashSet<&[u8]> = HashSet::new();
                for n in genome_names_content.split(|s| *s == b"\n"[0]) {
                    if n != b"" {
                        genome_names_hash.insert(n);
                    }
                }
                if genome_names_hash.is_empty() {
                    warn!("No genomes read in that are to be excluded from desharding process");
                    genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                    GenomeExclusionTypes::NoneType
                } else {
                    info!("Read in {} distinct genomes to exclude from desharding process e.g. '{}'",
                          genome_names_hash.len(),
                          std::str::from_utf8(genome_names_hash.iter().next().unwrap())
                              .unwrap());
                    if separator.is_some() {
                        genome_exclusion_filter_separator_type =
                            Some(SeparatorGenomeExclusionFilter {
                                split_char: separator.unwrap(),
                                excluded_genomes: genome_names_hash,
                            });
                        GenomeExclusionTypes::SeparatorType
                    } else {
                        match genomes_and_contigs_option_predereplication {
                            Some(ref gc) => {
                                genome_exclusion_genomes_and_contigs =
                                    Some(GenomesAndContigsExclusionFilter {
                                        genomes_and_contigs: gc,
                                        excluded_genomes: genome_names_hash,
                                    });
                                GenomeExclusionTypes::GenomesAndContigsType
                            }
                            None => unreachable!(),
                        }
                    }
                }
            } else {
                debug!("Not excluding any genomes during the deshard process");
                genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                GenomeExclusionTypes::NoneType
            }
        } else {
            genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
            GenomeExclusionTypes::NoneType
        }
    };

    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();

    if m.is_present("bam-files") {
        let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();

        // Associate genomes and contig names, if required
        let genomes_and_contigs_option = parse_all_genome_definitions(&m);

        if filter_params.doing_filtering() {
            let bam_readers = bam_generator::generate_filtered_bam_readers_from_bam_files(
                bam_files,
                filter_params.flag_filters.clone(),
                filter_params.min_aligned_length_single,
                filter_params.min_percent_identity_single,
                filter_params.min_aligned_percent_single,
                filter_params.min_aligned_length_pair,
                filter_params.min_percent_identity_pair,
                filter_params.min_aligned_percent_pair);

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           Some(long_readers),
                           separator,
                           &genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &concatenated_genomes, true);
                let mut all_generators = vec!();
                let mut indices = vec!(); // Prevent indices from being dropped
                for set in generator_sets {
                    indices.push(set.index);
                    for g in set.generators {
                        all_generators.push(g)
                    }
                }
                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           Some(all_generators),
                           separator,
                           &genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           separator,
                           &genomes_and_contigs_option)
            }


        } else if m.is_present("sharded") {
            external_command_checker::check_for_samtools();

            match genome_exclusion_type {
                GenomeExclusionTypes::NoneType => {
                    let genome_exclusion_non_type = genome_exclusion_filter_non_type.unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads,
                        &genome_exclusion_non_type);
                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_non_type);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_non_type,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }
                },
                GenomeExclusionTypes::SeparatorType => {
                    let genome_exclusion_separator = genome_exclusion_filter_separator_type.unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads,
                        &genome_exclusion_separator);

                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_separator);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_separator,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }

                }
                GenomeExclusionTypes::GenomesAndContigsType => {
                    let genome_exclusion_genomes = genome_exclusion_genomes_and_contigs.unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads,
                        &genome_exclusion_genomes);
                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_genomes);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_genomes,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   bam_readers,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }
                }

            }
        } else {
            let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                bam_files);

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    bam_files);
                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           Some(long_readers),
                           separator,
                           &genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &concatenated_genomes, true);
                let mut all_generators = vec!();
                let mut indices = vec!(); // Prevent indices from being dropped
                for set in generator_sets {
                    indices.push(set.index);
                    for g in set.generators {
                        all_generators.push(g)
                    }
                }
                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           Some(all_generators),
                           separator,
                           &genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           separator,
                           &genomes_and_contigs_option)
            }
        }

    } else {
        let mapping_program = parse_mapping_program(m.value_of("mapper"));
        external_command_checker::check_for_samtools();
        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
        if filter_params.doing_filtering() {
            debug!("Filtering..");
            let generator_sets = get_streamed_filtered_bam_readers(
                m,
                mapping_program,
                &concatenated_genomes,
                &filter_params,
                false);
            let mut all_generators = vec!();
            let mut indices = vec!(); // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            debug!("Finished collecting generators.");
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    bam_files);
                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           Some(long_readers),
                           separator,
                           &genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &concatenated_genomes, true);
                let mut long_generators = vec!();
                let mut indices = vec!(); // Prevent indices from being dropped
                for set in generator_sets {
                    indices.push(set.index);
                    for g in set.generators {
                        long_generators.push(g)
                    }
                }
                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           Some(long_generators),
                           separator,
                           &genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           separator,
                           &genomes_and_contigs_option)
            }
        } else if m.is_present("sharded") {
            match genome_exclusion_type {
                GenomeExclusionTypes::NoneType => {
                    let genome_exclusion_filter_none = genome_exclusion_filter_non_type.unwrap();
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &concatenated_genomes,
                        &genome_exclusion_filter_none, false);
                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_filter_none);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_filter_none,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }
                }
                GenomeExclusionTypes::SeparatorType => {
                    let genome_exclusion_separator = genome_exclusion_filter_separator_type.unwrap();
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &concatenated_genomes,
                        &genome_exclusion_separator, false);
                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_separator);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_separator,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }
                }
                GenomeExclusionTypes::GenomesAndContigsType => {
                    let genome_exclusion_genomes = genome_exclusion_genomes_and_contigs.unwrap();
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &concatenated_genomes,
                        &genome_exclusion_genomes, false);
                    if m.is_present("longread-bam-files") {
                        let long_files = m.values_of("longread-bam-files").unwrap().collect();
                        let long_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            long_files, sort_threads,
                            &genome_exclusion_genomes);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);

                    } else if m.is_present("longreads") {
                        let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                        external_command_checker::check_for_samtools();
                        let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                        let long_readers = get_sharded_bam_readers(
                            m,
                            mapping_program,
                            &concatenated_genomes,
                            &genome_exclusion_genomes,
                            true);
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   Some(long_readers),
                                   separator,
                                   &genomes_and_contigs_option);
                    } else {
                        run_pileup(m,
                                   mode,
                                   &mut estimators,
                                   generator_sets,
                                   filter_params.flag_filters,
                                   None::<Vec<PlaceholderBamFileReader>>,
                                   separator,
                                   &genomes_and_contigs_option);
                    }
                }
            }
        } else {
            debug!("Not filtering..");
            let generator_sets = get_streamed_bam_readers(m, mapping_program, &concatenated_genomes, false);
            let mut all_generators = vec!();
            let mut indices = vec!(); // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    bam_files);
                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           Some(long_readers),
                           separator,
                           &genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &concatenated_genomes, true);
                let mut long_generators = vec!();
                let mut indices = vec!(); // Prevent indices from being dropped
                for set in generator_sets {
                    indices.push(set.index);
                    for g in set.generators {
                        long_generators.push(g)
                    }
                }
                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           Some(long_generators),
                           separator,
                           &genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           separator,
                           &genomes_and_contigs_option)
            }
        }
    }
}

fn extract_genomes_and_contigs_option(
    m: &clap::ArgMatches,
    genome_fasta_files: &Vec<&str>,
) -> Option<GenomesAndContigs> {
    match m.is_present("genome-definition") {
        true => Some(coverm::genome_parsing::read_genome_definition_file(
            m.value_of("genome-definition").unwrap(),
        )),
        false => Some(coverm::genome_parsing::read_genome_fasta_files(
            &genome_fasta_files,
        )),
    }
}

fn parse_all_genome_definitions(m: &clap::ArgMatches) -> Option<GenomesAndContigs> {
    if m.is_present("single-genome") || m.is_present("separator") {
        None
    } else if m.is_present("genome-definition") {
        Some(coverm::genome_parsing::read_genome_definition_file(
            m.value_of("genome-definition").unwrap(),
        ))
    } else {
        extract_genomes_and_contigs_option(
            m,
            &bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(&m, true)
                .expect("Failed to parse genome paths")
                .iter()
                .map(|s| s.as_str())
                .collect(),
        )
    }
}

fn parse_separator(m: &clap::ArgMatches) -> Option<u8> {
    let single_genome = m.is_present("single-genome");
    if single_genome {
        Some("0".as_bytes()[0])
    } else if m.is_present("separator") {
        let separator_str = m.value_of("separator").unwrap().as_bytes();
        if separator_str.len() != 1 {
            eprintln!(
                "error: Separator can only be a single character, found {} ({}).",
                separator_str.len(),
                str::from_utf8(separator_str).unwrap()
            );
            process::exit(1);
        }
        Some(separator_str[0])
    } else if m.is_present("bam-files") || m.is_present("reference") {
        // Argument parsing enforces that genomes have been specified as FASTA
        // files.
        None
    } else {
        // Separator is set by CoverM and written into the generated reference
        // fasta file.
        Some(CONCATENATED_FASTA_FILE_SEPARATOR.as_bytes()[0])
    }
}

fn setup_genome_fasta_files(m: &clap::ArgMatches) -> (Option<NamedTempFile>, Option<GenomesAndContigs>){
    let genome_fasta_files_opt = {
        match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(&m, false) {
            Ok(paths) => {
                if paths.len() == 0 {
                    error!(
                        "Genome paths were described, but ultimately none were found"
                    );
                    process::exit(1);
                }
                if m.is_present("checkm-tab-table") || m.is_present("genome-info") {
                    let genomes_after_filtering =
                        galah::cluster_argument_parsing::filter_genomes_through_checkm(
                            &paths,
                            &m,
                            &galah_command_line_definition(),
                        )
                            .expect("Error parsing CheckM-related options");
                    info!(
                        "After filtering by CheckM, {} genomes remained",
                        genomes_after_filtering.len()
                    );
                    if genomes_after_filtering.len() == 0 {
                        error!("All genomes were filtered out, so none remain to be mapped to");
                        process::exit(1);
                    }
                    Some(
                        genomes_after_filtering
                            .iter()
                            .map(|s| s.to_string())
                            .collect(),
                    )
                } else {
                    Some(paths)
                }
            }
            Err(_) => None,
        }
    };

    let (concatenated_genomes, genomes_and_contigs_option) =
        match m.is_present("reference") {
            true => match genome_fasta_files_opt {
                Some(genome_paths) => (
                    None,
                    extract_genomes_and_contigs_option(
                        &m,
                        &genome_paths.iter().map(|s| s.as_str()).collect(),
                    ),
                ),
                None => (None, None),
            },
            false => {
                // Dereplicate if required
                let dereplicated_genomes: Vec<String> = if m.is_present("dereplicate") {
                    dereplicate(&m, &genome_fasta_files_opt.unwrap())
                } else {
                    genome_fasta_files_opt.unwrap()
                };
                info!("Profiling {} genomes", dereplicated_genomes.len());

                let list_of_genome_fasta_files = &dereplicated_genomes;
                info!(
                    "Generating concatenated reference FASTA file of {} genomes ..",
                    list_of_genome_fasta_files.len()
                );

                (
                    Some(
                        coverm::mapping_index_maintenance::generate_concatenated_fasta_file(
                            list_of_genome_fasta_files,
                        ),
                    ),
                    extract_genomes_and_contigs_option(
                        &m,
                        &dereplicated_genomes
                            .clone()
                            .iter()
                            .map(|s| s.as_str())
                            .collect(),
                    ),
                )
            }
        };
    return (concatenated_genomes, genomes_and_contigs_option)
}

fn dereplicate(m: &clap::ArgMatches, genome_fasta_files: &Vec<String>) -> Vec<String> {
    info!(
        "Found {} genomes specified before dereplication",
        genome_fasta_files.len()
    );

    // Generate clusterer and check for dependencies
    let clusterer = galah::cluster_argument_parsing::generate_galah_clusterer(
        genome_fasta_files,
        &m,
        &galah_command_line_definition(),
    ).expect("Failed to parse galah clustering arguments correctly");
    galah::external_command_checker::check_for_dependencies();
    info!("Dereplicating genome at {}% ANI ..", clusterer.ani*100.);

    let cluster_indices = clusterer.cluster();
    info!(
        "Finished dereplication, finding {} representative genomes.",
        cluster_indices.len()
    );
    debug!("Found cluster indices: {:?}", cluster_indices);
    let reps = cluster_indices
        .iter()
        .map(|cluster| genome_fasta_files[cluster[0]].clone())
        .collect::<Vec<_>>();
    debug!("Found cluster representatives: {:?}", reps);

    if m.is_present("output-dereplication-clusters") {
        let path = m.value_of("output-dereplication-clusters").unwrap();
        info!("Writing dereplication cluster memberships to {}", path);
        let mut f =
            std::fs::File::create(path).expect("Error creating dereplication cluster output file");
        for cluster in cluster_indices.iter() {
            let rep = cluster[0];
            for member in cluster {
                writeln!(
                    f,
                    "{}\t{}",
                    genome_fasta_files[rep], genome_fasta_files[*member]
                )
                    .expect("Failed to write a specific line to dereplication cluster file");
            }
        }
    }
    reps
}

fn parse_mapping_program(mapper: Option<&str>) -> MappingProgram {
    let mapping_program = match mapper {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        Some("ngmlr-ont") => MappingProgram::NGMLR_ONT,
        Some("ngmlr-pb") => MappingProgram::NGMLR_PB,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!(
            "Unexpected definition for --mapper: {:?}",
            mapper
        ),
    };
    match mapping_program {
        MappingProgram::BWA_MEM => {
            external_command_checker::check_for_bwa();
        }
        MappingProgram::MINIMAP2_SR |
        MappingProgram::MINIMAP2_ONT |
        MappingProgram::MINIMAP2_PB |
        MappingProgram::MINIMAP2_NO_PRESET => {
            external_command_checker::check_for_minimap2();
        }
        MappingProgram::NGMLR_ONT | MappingProgram::NGMLR_PB => {
            external_command_checker::check_for_ngmlr();
        }
    }
    return mapping_program;
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


fn generate_faidx(m: &clap::ArgMatches) -> bio::io::fasta::IndexedReader<File> {
    external_command_checker::check_for_samtools();
    info!("Generating reference index");
    let cmd_string = format!(
        "set -e -o pipefail; \
                     samtools faidx {}",
        m.value_of("reference").unwrap());
    debug!("Queuing cmd_string: {}", cmd_string);

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd_string)
        .stdout(Stdio::piped())
        .output()
        .expect("Unable to execute bash");

    let reference_path = Path::new(m.value_of("reference").unwrap());
    return bio::io::fasta::IndexedReader::from_file(&reference_path).expect(
        "Unable to generate index")
}

fn run_pileup<'a,
    R: NamedBamReader,
    T: NamedBamReaderGenerator<R>,
    S: NamedBamReader,
    U: NamedBamReaderGenerator<S>>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
    bam_readers: Vec<T>,
    flag_filters: FlagFilter,
    long_readers: Option<Vec<U>>,
    separator: Option<u8>,
    genomes_and_contigs_option: &Option<GenomesAndContigs>) {

    match mode {
        "polymorph" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let output_prefix = m.value_of("output-prefix").unwrap().to_string();
            let include_indels = m.is_present("include-indels");

            let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let threads = m.value_of("threads").unwrap().parse().unwrap();


            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }


            info!("Beginning polymorph with {} bam readers and {} threads", bam_readers.len(), threads);
            println!("sample\ttid\tpos\tvariant\treference\tvariant_depth\tdepth\tgenotypes\tvaf_cluster\tgenotype");
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
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
                false);
        },
        "genotype" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let reference_path = Path::new(m.value_of("reference").unwrap());

            let output_prefix = m.value_of("output-prefix").unwrap();
            let include_indels = m.is_present("include-indels");
            let include_soft_clipping = m.is_present("include-soft-clipping");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };

            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
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
                m.is_present("longread-bam-files"));
        },
        "summarize" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let reference_path = Path::new(m.value_of("reference").unwrap());

            let output_prefix = m.value_of("output-prefix").unwrap();
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };

            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
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
                m.is_present("longread-bam-files"));
        },
        "evolve" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let reference_path = Path::new(m.value_of("reference").unwrap());
            let output_prefix = m.value_of("output-prefix").unwrap();
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            info!("Beginning evolve with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
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
                m.is_present("longread-bam-files"));
        },
        "polish" => {
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let output_prefix = m.value_of("reference").unwrap().to_string();
            let output_prefix = output_prefix
                .split("..").last().unwrap()
                .split("/").last().unwrap().split(".").next().unwrap();


            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();
            let include_indels = m.is_present("include-indels");

            let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let threads = m.value_of("threads").unwrap().parse().unwrap();


            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }


            info!("Beginning polishing with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
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
                m.is_present("longread-bam-files"));
        },
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

