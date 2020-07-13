extern crate lorikeet_genome;

use lorikeet_genome::*;
use lorikeet_genome::estimation::contig;
use lorikeet_genome::external_command_checker;
use lorikeet_genome::cli::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate bio;
use bio::alignment::sparse::*;

extern crate coverm;
use coverm::*;
use coverm::shard_bam_reader::*;
use coverm::FlagFilter;
use coverm::genome_exclusion::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::{bam_generator::*, mapping_parameters::*, filter, mapping_index_maintenance};

use std::env;
use std::str;
use std::process;
use std::collections::BTreeMap;
use std::process::Stdio;
use std::fs::File;
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

const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";

const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;


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
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
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
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters,
                               None);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters,
                               None);
                } else {
                    let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters,
                               None);
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters,
                               None);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters,
                               None);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
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
                               all_generators,
                               filter_params.flag_filters.clone(),
                               None);
                }
            }
        },
        Some("summarize") => {
            let m = matches.subcommand_matches("summarize").unwrap();
            let mode = "summarize";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            let mut long_readers = vec!();
            if m.is_present("longread-bam-files") {
                let longreads = m.values_of("longread-bam-files").unwrap().collect();
                long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    longreads);
            };
            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
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
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters, Some(long_readers));
                } else {
                    let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters, Some(long_readers));
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       all_generators,
                                       filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       generator_sets,
                                       filter_params.flag_filters, Some(long_readers));
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
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
                                       all_generators,
                                       filter_params.flag_filters.clone(),
                               Some(long_readers));
                }
            }
        },
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            let mode = "genotype";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", genotype_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
            let mut long_readers = vec!();
            if m.is_present("longread-bam-files") {
                let longreads = m.values_of("longread-bam-files").unwrap().collect();
                long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    longreads);
            };

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
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
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers))
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers))

                } else {
                    let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers))
                }

            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters, Some(long_readers));
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
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
                               all_generators,
                               filter_params.flag_filters.clone(), Some(long_readers));
                }
            }
        },
        Some("evolve") => {
            let m = matches.subcommand_matches("evolve").unwrap();
            let mode = "evolve";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            let mut long_readers = vec!();
            if m.is_present("longread-bam-files") {
                let longreads = m.values_of("longread-bam-files").unwrap().collect();
                long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    longreads);
            };
            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();

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

                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m, mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers));
                } else {
                    let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers));
                }
            } else if m.is_present("read1") |
                m.is_present("interleaved") |
                m.is_present("coupled") |
                m.is_present("single"){
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters, Some(long_readers));
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
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
                               all_generators,
                               filter_params.flag_filters.clone(), Some(long_readers));
                }
            }
        },
        Some("polish") => {
            let m = matches.subcommand_matches("polish").unwrap();
            let mode = "polish";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            let mut long_readers = vec!();
            if m.is_present("longread-bam-files") {
                let longreads = m.values_of("longread-bam-files").unwrap().collect();
                long_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                    longreads);
            };

            if m.is_present("bam-file") {
                let bam_files: Vec<&str> = m.values_of("bam-file").unwrap().collect();
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
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers));
                } else {
                    let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters, Some(long_readers));
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters, Some(long_readers));
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters, Some(long_readers));
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
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
                               all_generators,
                               filter_params.flag_filters.clone(), Some(long_readers));
                }
            }
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

fn setup_mapping_index(
    reference_wise_params: &SingleReferenceMappingParameters,
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
) -> Option<Box<dyn mapping_index_maintenance::MappingIndex>> {
    match mapping_program {
        MappingProgram::BWA_MEM => Some(mapping_index_maintenance::generate_bwa_index(
            reference_wise_params.reference,
            None
        )),
        MappingProgram::MINIMAP2_SR |
        MappingProgram::MINIMAP2_ONT |
        MappingProgram::MINIMAP2_PB |
        MappingProgram::MINIMAP2_NO_PRESET => {
            if m.is_present("minimap2-reference-is-index") || reference_wise_params.len() == 1 {
                info!("Not pre-generating minimap2 index");
                if m.is_present("minimap2-reference-is-index") {
                    warn!("Minimap2 uses mapping parameters defined when the index was created, \
                    not parameters defined when mapping. Proceeding on the assumption that you \
                    passed the correct parameters when creating the minimap2 index.");
                }
                None
            } else {
                Some(mapping_index_maintenance::generate_minimap2_index(
                    reference_wise_params.reference,
                    Some(m.value_of("threads").unwrap().parse::<usize>().unwrap()),
                    Some(m.value_of("minimap2-params").unwrap_or("")),
                    mapping_program,
                ))
            }
        }
    }
}

fn parse_mapping_program(m: &clap::ArgMatches) -> MappingProgram {
    let mapping_program = match m.value_of("mapper") {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!(
            "Unexpected definition for --mapper: {:?}",
            m.value_of("mapper")
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
    }
    return mapping_program;
}


fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("method") { return false; }
            let methods: &str = m.value_of("method").unwrap();
            if methods.contains(&"metabat") {
                return true
            }
            return false
        },
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false
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

#[derive(Debug)]
struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        let mut f = FilterParameters {
            flag_filters: FlagFilter {
                include_improper_pairs: !m.is_present("proper-pairs-only"),
                include_secondary: m.is_present("include-secondary"),
                include_supplementary: m.is_present("include-supplementary"),
            },
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => value_t!(m.value_of("min-read-aligned-length"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => value_t!(m.value_of("min-read-aligned-length-pair"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_pair: parse_percentage(&m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(&m, "min-read-aligned-percent-pair"),

        };
        if m.is_present("nanopore") {
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
        }
        if doing_metabat(&m) {
            debug!(
                "Setting single read percent identity threshold at 0.97 for \
                 MetaBAT adjusted coverage."
            );
            // we use >= where metabat uses >. Gah.
            f.min_percent_identity_single = 0.97001;
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
            f.flag_filters.include_secondary = true;
        }
        debug!("Filter parameters set as {:?}", f);
        return f;
    }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0
            || self.min_percent_identity_pair > 0.0
            || self.min_aligned_percent_single > 0.0
            || self.min_aligned_percent_pair > 0.0
            || self.min_aligned_length_single > 0
            || self.min_aligned_length_pair > 0;
    }
}

fn get_sharded_bam_readers<'a, 'b, T>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    genome_exclusion: &'b T,
) -> Vec<ShardedBamReaderGenerator<'b, T>>
    where
        T: GenomeExclusion,
{
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");
    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut bam_readers = vec![];
    let mut concatenated_reference_name: Option<String> = None;
    let mut concatenated_read_names: Option<String> = None;

    for reference_wise_params in params {
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let reference_name = std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string();
        concatenated_reference_name = match concatenated_reference_name {
            Some(prev) => Some(format!("{}|{}", prev, reference_name)),
            None => Some(reference_name),
        };
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };
        let n_samples = reference_wise_params.len() as u16;

        for p in reference_wise_params {
            bam_readers.push(
                shard_bam_reader::generate_named_sharded_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads / n_samples,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                ),
            );
            let name = &std::path::Path::new(p.read1)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string();
            concatenated_read_names = match concatenated_read_names {
                Some(prev) => Some(format!("{}|{}", prev, name)),
                None => Some(name.to_string()),
            };
        }

        debug!("Finished BAM setup");
    }
    let gen = ShardedBamReaderGenerator {
        stoit_name: format!(
            "{}/{}",
            concatenated_reference_name.unwrap(),
            concatenated_read_names.unwrap()
        ),
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec![gen];
}

fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        let n_samples = reference_wise_params.len() as u16;

        for p in reference_wise_params {
            bam_readers.push(
                bam_generator::generate_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

fn generate_cached_bam_file_name(directory: &str, reference: &str, read1_path: &str) -> String {
    debug!("Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
           directory, reference, read1_path);
    std::path::Path::new(directory).to_str()
        .expect("Unable to covert bam-file-cache-directory name into str").to_string()+"/"+
        &std::path::Path::new(reference).file_name()
        .expect("Unable to convert reference to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+"."+
        &std::path::Path::new(read1_path).file_name()
        .expect("Unable to convert read1 name to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+".bam"
}

fn setup_bam_cache_directory(cache_directory: &str) {
    let path = std::path::Path::new(cache_directory);
    if path.is_dir() {
        if path.metadata().expect("Unable to read metadata for cache directory")
            .permissions().readonly() {
                panic!("Cache directory {} does not appear to be writeable, not continuing",
                       cache_directory);
            } else {
                info!("Writing BAM files to already existing directory {}", cache_directory)
            }
    } else {
        match path.parent() {
            Some(parent) => {
                let parent2 = match parent == std::path::Path::new("") {
                    true => std::path::Path::new("."),
                    false => parent
                };
                if parent2.canonicalize().expect(
                    &format!("Unable to canonicalize parent of cache directory {}", cache_directory)).is_dir() {
                    if parent2.metadata().expect(
                        &format!("Unable to get metadata for parent of cache directory {}",
                                 cache_directory))
                        .permissions().readonly() {
                            panic!(
                                "The parent directory of the (currently non-existent) \
                                 cache directory {} is not writeable, not continuing",
                                cache_directory);
                        } else {
                            info!("Creating cache directory {}", cache_directory);
                            std::fs::create_dir(path).expect("Unable to create cache directory");
                        }
                } else {
                    panic!("The parent directory of the cache directory {} does not \
                            yet exist, so not creating that cache directory, and not continuing.",
                           cache_directory)
                }
            },
            None => {
                panic!("Cannot create root directory {}", cache_directory)
            }
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        panic!("Failed to create test file in bam cache directory: {}", tf_result.err().unwrap())
    }
}

fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };
        let n_samples = reference_wise_params.len() as u16;

        for p in reference_wise_params {
            bam_readers.push(
                bam_generator::generate_filtered_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    std::cmp::max(p.threads / n_samples, 1),
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    p.mapping_options,
                    discard_unmapped,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
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
    R: bam_generator::NamedBamReader + Send,
    T: bam_generator::NamedBamReaderGenerator<R> + Send>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
    bam_readers: Vec<T>,
    flag_filters: FlagFilter,
    long_readers: Option<Vec<bam_generator::BamFileNamedReader>>) {
    match mode {
        "polymorph" => {
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let output_prefix = m.value_of("output-prefix").unwrap().to_string();
            let include_indels = m.is_present("include-indels");

            // Make sure we are dealing with a fresh file
            let file_name = output_prefix.to_string()
                + &".fna".to_owned();

            let file_path = Path::new(&file_name);

            File::create(file_path)
                .expect("No Read or Write Permission in current directory");
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
                print_zeros,
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
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
//            let kmer_size = m.value_of("kmer-size").unwrap().parse().unwrap();
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
                print_zeros,
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
            let print_zeros = !m.is_present("no-zeros");
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
                print_zeros,
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
            let print_zeros = !m.is_present("no-zeros");
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
                print_zeros,
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
                print_zeros,
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

