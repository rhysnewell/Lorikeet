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
use coverm::FlagFilter;
use coverm::genome_exclusion::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::{bam_generator::*, filter};

use std::env;
use std::str;
use std::process;
use std::collections::{BTreeMap, HashSet};
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

fn prepare_pileup
    (m: &clap::ArgMatches, mode: &str) {

    // This function is amazingly painful. It handles every combination of longread and short read
    // mapping or bam file reading. Could not make it smaller using dynamic or static dispatch
    set_log_level(m, true);
    let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
    let filter_params = FilterParameters::generate_from_clap(m);
    let threads = m.value_of("threads").unwrap().parse().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    let separator = parse_separator(m);

    let references = parse_references(m);
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

    let genomes_and_contigs_option = extract_genomes_and_contigs_option(&m, &references);
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
                           genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &None, true);
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
                           genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           genomes_and_contigs_option)
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
                           genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &None, true);
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
                           genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           bam_readers,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           genomes_and_contigs_option)
            }
        }

    } else {
        let mapping_program = parse_mapping_program(m.value_of("mapper"));
        external_command_checker::check_for_samtools();
        if filter_params.doing_filtering() {
            debug!("Filtering..");
            let generator_sets = get_streamed_filtered_bam_readers(
                m,
                mapping_program,
                &None,
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
                           genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &None, true);
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
                           genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           genomes_and_contigs_option)
            }
        } else {
            debug!("Not filtering..");
            let generator_sets = get_streamed_bam_readers(m, mapping_program, &None, false);
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
                           genomes_and_contigs_option)
            } else if m.is_present("longreads") {
                // Perform mapping
                let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
                external_command_checker::check_for_samtools();
                let generator_sets = get_streamed_bam_readers(m, mapping_program, &None, true);
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
                           genomes_and_contigs_option)
            } else {

                run_pileup(m,
                           mode,
                           &mut estimators,
                           all_generators,
                           filter_params.flag_filters,
                           None::<Vec<PlaceholderBamFileReader>>,
                           genomes_and_contigs_option)
            }
        }
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
    genomes_and_contigs_option: Option<GenomesAndContigs>) {

    match mode {
        "polymorph" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();
            let output_prefix = m.value_of("output-prefix").unwrap().to_string();
            let include_indels = m.is_present("include-indels");
            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            println!("sample\ttid\tpos\tvariant\treference\tvariant_depth\tdepth\tgenotypes\tvaf_cluster\tgenotype");
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
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
                false,
                genomes_and_contigs);
        },
        "genotype" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();

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
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();

            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
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
                genomes_and_contigs);
        },
        "summarize" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();

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
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();


            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
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
                genomes_and_contigs);
        },
        "evolve" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let output_prefix = m.value_of("output-prefix").unwrap();
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();


            info!("Beginning evolve with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
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
                genomes_and_contigs);
        },
        "polish" => {
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let output_prefix = m.value_of("reference").unwrap().to_string();
            let output_prefix = output_prefix
                .split("..").last().unwrap()
                .split("/").last().unwrap().split(".").next().unwrap();


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
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }
            let genomes_and_contigs = genomes_and_contigs_option.unwrap();



            info!("Beginning polishing with {} bam readers and {} threads", bam_readers.len(), threads);
            contig::pileup_variants(
                m,
                bam_readers,
                long_readers,
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
                genomes_and_contigs);
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

