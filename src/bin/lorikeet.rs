extern crate openssl;
extern crate openssl_sys;

extern crate galah;
use galah::PreclusterDistanceFinder;
use galah::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use galah::{dashing::DashingPreclusterer, finch::FinchPreclusterer};

extern crate lorikeet_genome;
use lorikeet_genome::cli::*;
use lorikeet_genome::external_command_checker;
use lorikeet_genome::pair_hmm::pair_hmm_likelihood_calculation_engine::AVXMode;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::SmithWatermanAligner;
use lorikeet_genome::smith_waterman::smith_waterman_aligner::NEW_SW_PARAMETERS;
use lorikeet_genome::utils::utils::*;
use lorikeet_genome::*;

extern crate gkl;
use gkl::smithwaterman::OverhangStrategy;

extern crate needletail;
use needletail::parse_fastx_file;

extern crate coverm;
use coverm::mapping_index_maintenance::generate_concatenated_fasta_file;

use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use coverm::*;

extern crate bird_tool_utils;

use std::env;
use std::path;
use std::process;
use std::process::Stdio;
use std::sync::{Mutex, Arc};
use std::time::Duration;

extern crate tempfile;
use tempfile::NamedTempFile;

extern crate clap;
use clap::*;

extern crate clap_complete;
use clap_complete::{generate, Shell};

#[macro_use]
extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;
use lorikeet_genome::processing::lorikeet_engine::{
    run_summarize, start_lorikeet_engine, ReadType, Elem
};
use lorikeet_genome::reference::reference_reader_utils::ReferenceReaderUtils;
use lorikeet_genome::utils::errors::BirdToolError;

extern crate rayon;
use rayon::prelude::*;

extern crate indicatif;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};

extern crate partitions;
use partitions::partition_vec::PartitionVec;

// contant variables
// max i16 size minus 1
const MAX_I16: i16 = 32767;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("summarise") => {
            let m = matches.subcommand_matches("summarise").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(&m, summarise_full_help());
            rayon::ThreadPoolBuilder::new()
                .num_threads(m.value_of("threads").unwrap().parse().unwrap())
                .build_global()
                .unwrap();
            run_summarize(m);
        }
        Some("pangenome") => {
            let m = matches.subcommand_matches("pangenome").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(&m, summarise_full_help());

            match run_pangenome(m) {
                Ok(_) => info!("Pangenome created."),
                Err(e) => panic!("Pangenome creation failed with error: {:?}", e),
            };
        }
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(&m, genotype_full_help());
            let mode = "genotype";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Genotype complete."),
                Err(e) => warn!("Genotype failed with error: {:?}", e),
            };
        }
        Some("call") => {
            let m = matches.subcommand_matches("call").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(&m, call_full_help());
            let mode = "call";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Call complete."),
                Err(e) => warn!("Call failed with error: {:?}", e),
            };
        }
        Some("consensus") => {
            let m = matches.subcommand_matches("consensus").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(&m, consensus_full_help());
            let mode = "consensus";

            match prepare_pileup(m, mode) {
                Ok(_) => info!("Consensus complete."),
                Err(e) => warn!("Consensus failed with error: {:?}", e),
            };
        }
        Some("shell-completion") => {
            let m = matches.subcommand_matches("shell-completion").unwrap();
            set_log_level(m, true);
            let mut file = std::fs::File::create(m.value_of("output-file").unwrap())
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
            println!();
        }
    }
}

fn run_pangenome(m: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    if !m.is_present("use-avx") {
        external_command_checker::check_for_pggb();
        set_log_level(m, true);
        let threads = m.value_of("threads").unwrap().parse()?;
        let percent_identity =
            galah::cluster_argument_parsing::parse_percentage(m, "percent-identity")?
                .unwrap_or(0.95);
        let kmer_size: usize = m.value_of("kmer-size").unwrap().parse()?;
        let segment_length: usize = m.value_of("segment-length").unwrap().parse()?;

        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
        debug!("Parsing references");
        let references = ReferenceReaderUtils::parse_references(m);
        // let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

        debug!("Building concatenated reference");
        let concatenated_genomes = generate_concatenated_fasta_file(&references);
        let concat_path = concatenated_genomes.path().to_str().unwrap();
        ReferenceReaderUtils::generate_faidx(concat_path);

        let pggb_params = m.value_of("pggb-params").unwrap();
        let output_dir = m.value_of("output").unwrap();
        let cmd_string = format!(
            "pggb -i {concat_path} -n {} -t {threads} -o {output_dir} -p {percent_identity} -k {kmer_size} -s {segment_length} {pggb_params}", 
            references.len()
        );

        info!("Running pggb.");
        std::process::Command::new("bash")
            .arg("-c")
            .arg(&cmd_string)
            .stdout(Stdio::piped())
            .output()?;
    } else {
        // read in references and loop through entries and align entries pairwise
        // using SmithWatermanAligner
        set_log_level(m, true);
        let threads = m.value_of("threads").unwrap().parse()?;
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
        debug!("Parsing references");
        let references = ReferenceReaderUtils::parse_references(m);
        let references = references.iter().map(|s| s.as_str()).collect::<Vec<&str>>();
        debug!("Building concatenated reference");
        // let concatenated_genomes  = generate_concatenated_fasta_file(&references);

        // create GalahClusterer
        // match preclusterer, dashing or finch
        let preclusterer: Box<dyn PreclusterDistanceFinder> = 
            match m.value_of("precluster-method").unwrap() {
                "dashing" => {
                    external_command_checker::check_for_dashing();
                    Box::new(
                        DashingPreclusterer {
                            min_ani: galah::cluster_argument_parsing::parse_percentage(
                                m,
                                "precluster-identity",
                            )?
                            .unwrap_or(0.95),
                            threads,
                        })
                }
                "finch" => Box::new(
                    FinchPreclusterer {
                        min_ani: galah::cluster_argument_parsing::parse_percentage(
                            m,
                            "precluster-identity",
                        )?
                        .unwrap_or(0.95),
                        num_kmers: 1000,
                        kmer_length: m.value_of("kmer-size").unwrap().parse()?,
                    }),
                _ => unreachable!(),
        };
        
        let dashing_cache = preclusterer.distances(&references);
        info!("Preclustering with {}", m.value_of("precluster-method").unwrap());
        let minhash_preclusters = partition_sketches(
            &references, 
            &dashing_cache
        );
        trace!("Found preclusters: {:?}", minhash_preclusters);
        // Convert single linkage data structure into just a list of list of indices
        let mut clusters: Vec<Vec<usize>> = minhash_preclusters
            .all_sets()
            .map(|cluster| {
                let mut indices: Vec<_> = cluster.map(|cluster_genome| *cluster_genome.1).collect();
                indices.sort_unstable();
                indices
            })
            .collect();

        // Sort preclusters so bigger clusters are started before smaller
        clusters.sort_unstable_by(|c1, c2| c2.len().cmp(&c1.len()));
        debug!("After sorting, found preclusters {:?}", clusters);
        info!(
            "Found {} preclusters. The largest contained {} genomes",
            clusters.len(),
            clusters[0].len()
        );

        // set up ProgressBar the length of num_contigs
        let sty_eta = ProgressStyle::default_bar().template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]",
        )?;

        // Set up multi progress bars
        let multi_inner = Arc::new(MultiProgress::new());
        let mut progress_bars = vec![
            Elem {
                key: "Genome clusters to align".to_string(),
                index: 0,
                progress_bar: ProgressBar::new(clusters.len() as u64),
            };
            clusters.len() + 1
        ];

        let progress_bar = multi_inner.add(progress_bars[0].progress_bar.clone());
        progress_bar.set_style(sty_eta.clone());
        progress_bar.set_message("Genome clusters to align");
        progress_bar.enable_steady_tick(Duration::from_millis(200));
        progress_bars[0] = Elem {
            key: format!("Genome clusters to align"),
            index: 0,
            progress_bar,
        };


        debug!("Setting progress bar styles");
        clusters.iter().enumerate().for_each(|(i, cluster)| {
            
            let progress_bar = multi_inner.add(ProgressBar::new(cluster.len() as u64));
            progress_bar.set_style(sty_eta.clone());
            progress_bar.set_message(format!("Aligning cluster {}", i));
            progress_bars[i + 1] = Elem {
                key: format!("Aligning cluster {}", i),
                index: i + 1,
                progress_bar,
            };
        });
        
        debug!("Inserting bars into tree");
        let tree: Arc<Mutex<Vec<&Elem>>> =
            Arc::new(Mutex::new(Vec::with_capacity(progress_bars.len())));
        {
            let mut tree = tree.lock().unwrap();
            for pb in progress_bars.iter() {
                tree.push(pb)
            }
        }
        
        // begin ticks on Elem in progress_bars
        // debug!("Beginning ticks on progress bars");
        // (0..clusters.len()+1).into_iter().for_each(|i| {
        //     let elem = &progress_bars[i];
        //     let pb = multi_inner.insert(i, elem.progress_bar.clone());
        //     pb.enable_steady_tick(Duration::from_millis(200));
        // });

        // take indices from clusters and index into references
        // align pairwise refrences in each cluster
        // using SmithWatermanAligner
        info!("Aligning genomes");
        clusters.par_iter().enumerate().for_each(|(cluster_index, cluster)| {
            // let aligner = SmithWatermanAligner::new();
            // let mut alignments = Vec::new();
            {
                let tree = tree.lock().unwrap();
                let elem = tree[cluster_index + 1];
                elem.progress_bar.enable_steady_tick(Duration::from_millis(200));
            }
            // loop through cluster indices in parallel
            cluster.par_iter().for_each(|i| {
                cluster[*i+1..].par_iter().for_each(|j| {
                    let mut ref1 = parse_fastx_file(path::Path::new(&references[*i]))
                        .expect("Unable to open reference file");
                    let mut ref2 = parse_fastx_file(path::Path::new(&references[*j]))
                        .expect("Unable to open reference file");
                    while let Some(record1) = ref1.next() {
                        let seqrec1 = record1.expect("invalid record");
                        let seq1 = seqrec1.seq();
                        while let Some(record2) = ref2.next() {
                            let seqrec2 = record2.expect("invalid record");
                            
                            let seq2 = seqrec2.seq();
                            // iterate seq1 and seq2 in chunks of size MAX_I16 and align
                            // each chunk
                            seq1.chunks(MAX_I16 as usize).for_each(|chunk1| {
                                seq2.chunks(MAX_I16 as usize).for_each(|chunk2| {
                                    let _ = SmithWatermanAligner::align(
                                        &chunk1,
                                        &chunk2,
                                        &NEW_SW_PARAMETERS,
                                        OverhangStrategy::SoftClip,
                                        AVXMode::detect_mode(),
                                    );
                                });
                            });
                            // alignments.push(alignment);
                        }
                    }
                });
                let tree = tree.lock().unwrap();
                tree[cluster_index + 1].progress_bar.inc(1);
            });
            // increment top level of multi progress bar in tree
            let tree = tree.lock().unwrap();
            tree[0].progress_bar.inc(1);
            tree[cluster_index + 1].progress_bar.finish();
        });
        // finish all progress bars in tree
        let tree = tree.lock().unwrap();
        tree[0].progress_bar.finish();
        // pb.finish();
    }
    Ok(())
}

/// Create sub-sets by single linkage clustering
fn partition_sketches(
    genomes: &[&str],
    dashing_cache: &SortedPairGenomeDistanceCache,
) -> PartitionVec<usize> {
    let mut to_return: PartitionVec<usize> = PartitionVec::with_capacity(genomes.len());
    for (i, _) in genomes.iter().enumerate() {
        to_return.push(i);
    }

    genomes.iter().enumerate().for_each(|(i, _)| {
        genomes[0..i].iter().enumerate().for_each(|(j, _)| {
            trace!("Testing precluster between {} and {}", i, j);
            if dashing_cache.contains_key(&(i, j)) {
                debug!(
                    "During preclustering, found a match between genomes {} and {}",
                    i, j
                );
                to_return.union(i, j)
            }
        });
    });
    return to_return;
}

fn prepare_pileup(m: &clap::ArgMatches, mode: &str) -> Result<(), BirdToolError> {
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

    let references = ReferenceReaderUtils::parse_references(m);
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

    let (concatenated_genomes, genomes_and_contigs_option) =
        ReferenceReaderUtils::setup_genome_fasta_files(&m);
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
                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else if m.is_present("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
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
                    &mut estimators,
                    bam_readers,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            }
        } else {
            let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    bam_readers,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else if m.is_present("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
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
                    &mut estimators,
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
                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else if m.is_present("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
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
                    &mut estimators,
                    all_generators,
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
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

                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    all_generators,
                    filter_params.flag_filters,
                    Some(long_readers),
                    genomes_and_contigs_option,
                    tmp_dir,
                    concatenated_genomes,
                );
            } else if m.is_present("longreads") {
                // Perform mapping
                let (long_generators, _indices) = long_generator_setup(
                    &m,
                    &concatenated_genomes,
                    &Some(references.clone()),
                    &tmp_dir,
                );

                return run_pileup(
                    m,
                    mode,
                    &mut estimators,
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
                    &mut estimators,
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
>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
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
        estimators.estimators.clone(),
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
