use crate::*;

use coverm::bam_generator::*;
use coverm::genomes_and_contigs::*;
use coverm::mapping_index_maintenance;
use coverm::mapping_parameters::*;
use coverm::FlagFilter;

use bio::io::fasta::IndexedReader;
use glob::glob;
use nix::{sys::stat, unistd};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::process::Stdio;
use std::str;
use tempdir::TempDir;
use tempfile::NamedTempFile;

pub const NUMERICAL_EPSILON: f64 = 1e-3;
pub const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";
pub const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

pub fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    longread: bool,
    references: &'a Option<Vec<&'a str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        if longread {
            setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
            setup_bam_cache_directory(&format!(
                "{}/long/",
                m.value_of("bam-file-cache-directory").unwrap()
            ));
        } else {
            setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
            setup_bam_cache_directory(&format!(
                "{}/short/",
                m.value_of("bam-file-cache-directory").unwrap()
            ));
        }
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = if !longread {
        MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile, &references)
    } else {
        MappingParameters::generate_longread_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        )
    };
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        debug!("Reference file {:?}", &reference);
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        tmp_bam_file_cache
                            .as_ref()
                            .unwrap()
                            .path()
                            .to_str()
                            .unwrap(),
                        reference,
                        naming_readset,
                        longread,
                    );
                    Some(bam_file_cache_path)
                }
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                        longread,
                    );
                    debug!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        let _n_samples = reference_wise_params.len() as u16;

        for p in reference_wise_params {
            bam_readers.push(generate_named_bam_readers_from_reads(
                mapping_program,
                match index {
                    Some(ref index) => index.index_path(),
                    None => {
                        warn!("Not using reference index...");
                        reference
                    }
                },
                p.read1,
                p.read2,
                p.read_format.clone(),
                p.threads,
                bam_file_cache(p.read1).as_ref().map(String::as_ref),
                discard_unmapped,
                p.mapping_options,
                reference_tempfile.is_none(),
            ));
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

pub fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
    longread: bool,
    references: &Option<Vec<&str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        if longread {
            setup_bam_cache_directory(&format!(
                "{}/long/",
                m.value_of("bam-file-cache-directory").unwrap()
            ));
        } else {
            setup_bam_cache_directory(&format!(
                "{}/short/",
                m.value_of("bam-file-cache-directory").unwrap()
            ));
        }
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = if !longread {
        MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile, &references)
    } else {
        MappingParameters::generate_longread_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        )
    };
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        debug!("Reference file {:?}", &reference);

        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        tmp_bam_file_cache
                            .as_ref()
                            .unwrap()
                            .path()
                            .to_str()
                            .unwrap(),
                        reference,
                        naming_readset,
                        longread,
                    );
                    Some(bam_file_cache_path)
                }
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                        longread,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(generate_filtered_named_bam_readers_from_reads(
                mapping_program,
                match index {
                    Some(ref index) => index.index_path(),
                    None => {
                        warn!("Not using reference index...");
                        reference
                    }
                },
                p.read1,
                p.read2,
                p.read_format.clone(),
                p.threads,
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
            ));
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

pub fn setup_mapping_index(
    reference_wise_params: &SingleReferenceMappingParameters,
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
) -> Option<Box<dyn mapping_index_maintenance::MappingIndex>> {
    match mapping_program {
        MappingProgram::BWA_MEM => Some(mapping_index_maintenance::generate_bwa_index(
            reference_wise_params.reference,
            None,
        )),
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => {
            if m.is_present("minimap2-reference-is-index") || reference_wise_params.len() == 1 {
                info!("Not pre-generating minimap2 index");
                if m.is_present("minimap2-reference-is-index") {
                    warn!(
                        "Minimap2 uses mapping parameters defined when the index was created, \
                    not parameters defined when mapping. Proceeding on the assumption that you \
                    passed the correct parameters when creating the minimap2 index."
                    );
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
        MappingProgram::NGMLR_ONT | MappingProgram::NGMLR_PB => {
            // NGMLR won't let us create a mapping index, we use --skip-write to avoid index being written
            // to disk
            None
        }
    }
}

pub fn setup_bam_cache_directory(cache_directory: &str) {
    let path = std::path::Path::new(cache_directory);
    if path.is_dir() {
        if path
            .metadata()
            .expect("Unable to read metadata for cache directory")
            .permissions()
            .readonly()
        {
            panic!(
                "Cache directory {} does not appear to be writeable, not continuing",
                cache_directory
            );
        } else {
            debug!(
                "Writing BAM files to already existing directory {}",
                cache_directory
            )
        }
    } else {
        match path.parent() {
            Some(parent) => {
                let parent2 = match parent == std::path::Path::new("") {
                    true => std::path::Path::new("."),
                    false => parent,
                };
                if parent2
                    .canonicalize()
                    .expect(&format!(
                        "Unable to canonicalize parent of cache directory {}",
                        cache_directory
                    ))
                    .is_dir()
                {
                    if parent2
                        .metadata()
                        .expect(&format!(
                            "Unable to get metadata for parent of cache directory {}",
                            cache_directory
                        ))
                        .permissions()
                        .readonly()
                    {
                        panic!(
                            "The parent directory of the (currently non-existent) \
                                 cache directory {} is not writeable, not continuing",
                            cache_directory
                        );
                    } else {
                        info!("Creating cache directory {}", cache_directory);
                        std::fs::create_dir(path).expect("Unable to create cache directory");
                    }
                } else {
                    panic!(
                        "The parent directory of the cache directory {} does not \
                            yet exist, so not creating that cache directory, and not continuing.",
                        cache_directory
                    )
                }
            }
            None => panic!("Cannot create root directory {}", cache_directory),
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        panic!(
            "Failed to create test file in bam cache directory: {}",
            tf_result.err().unwrap()
        )
    }
}

pub fn generate_cached_bam_file_name(
    directory: &str,
    reference: &str,
    read1_path: &str,
    longread: bool,
) -> String {
    debug!(
        "Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
        directory, reference, read1_path
    );
    if longread {
        std::path::Path::new(&format!("{}/long/", directory))
            .to_str()
            .expect("Unable to covert bam-file-cache-directory name into str")
            .to_string()
            + "/"
            + &std::path::Path::new(reference)
                .file_name()
                .expect("Unable to convert reference to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
            + "."
            + &std::path::Path::new(read1_path)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
            + ".bam"
    } else {
        std::path::Path::new(&format!("{}/short/", directory))
            .to_str()
            .expect("Unable to covert bam-file-cache-directory name into str")
            .to_string()
            + "/"
            + &std::path::Path::new(reference)
                .file_name()
                .expect("Unable to convert reference to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
            + "."
            + &std::path::Path::new(read1_path)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
            + ".bam"
    }
}

#[derive(Debug)]
pub struct FilterParameters {
    pub flag_filters: FlagFilter,
    pub min_aligned_length_single: u32,
    pub min_percent_identity_single: f32,
    pub min_aligned_percent_single: f32,
    pub min_aligned_length_pair: u32,
    pub min_percent_identity_pair: f32,
    pub min_aligned_percent_pair: f32,
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

pub fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("method") {
                return false;
            }
            let methods: &str = m.value_of("method").unwrap();
            if methods.contains(&"metabat") {
                return true;
            }
            return false;
        }
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false;
        }
    }
}

pub fn setup_genome_fasta_files(
    m: &clap::ArgMatches,
) -> (Option<NamedTempFile>, Option<GenomesAndContigs>) {
    let genome_fasta_files_opt = {
        match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(&m, false) {
            Ok(paths) => {
                if paths.len() == 0 {
                    error!("Genome paths were described, but ultimately none were found");
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

    debug!("Found paths {:?}", &genome_fasta_files_opt);

    let (concatenated_genomes, genomes_and_contigs_option) = match m.is_present("reference") {
        true => match genome_fasta_files_opt {
            Some(genome_paths) => (
                Some(
                    coverm::mapping_index_maintenance::generate_concatenated_fasta_file(
                        &genome_paths,
                    ),
                ),
                extract_genomes_and_contigs_option(
                    &m,
                    &genome_paths.iter().map(|s| s.as_str()).collect(),
                ),
            ),
            None => (None, None),
        },
        false => {
            // Dereplicate if required
            // TODO: Properly implement dereplication in cli.rs and make sure function works
            let dereplicated_genomes: Vec<String> = if m.is_present("dereplicate") {
                dereplicate(&m, &genome_fasta_files_opt.unwrap())
            } else {
                genome_fasta_files_opt.unwrap()
            };
            debug!("Profiling {} genomes", dereplicated_genomes.len());

            let list_of_genome_fasta_files = &dereplicated_genomes;

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

    debug!("Found genome_and_contigs {:?}", &genomes_and_contigs_option);
    return (concatenated_genomes, genomes_and_contigs_option);
}

pub fn parse_references(m: &clap::ArgMatches) -> Vec<String> {
    let references = match m.values_of("genome-fasta-files") {
        Some(vec) => {
            let reference_paths = vec.map(|p| p.to_string()).collect::<Vec<String>>();
            debug!("Reference files {:?}", reference_paths);
            reference_paths
        }
        None => match m.value_of("genome-fasta-directory") {
            Some(path) => {
                let ext = m.value_of("genome-fasta-extension").unwrap();
                let reference_glob = format!("{}/*.{}", path, ext);
                let reference_paths = glob(&reference_glob)
                    .expect("Failed to read cache")
                    .map(|p| {
                        p.expect("Failed to read cached bam path")
                            .to_str()
                            .unwrap()
                            .to_string()
                    })
                    .collect::<Vec<String>>();
                debug!("Reference files {:?}", reference_paths);
                reference_paths
            }
            None => panic!("Can't find suitable references for variant calling"),
        },
    };
    return references;
}

pub fn extract_genomes_and_contigs_option(
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

pub fn generate_faidx(reference_path: &str) -> bio::io::fasta::IndexedReader<File> {
    external_command_checker::check_for_samtools();
    info!("Generating reference index");
    let cmd_string = format!(
        "set -e -o pipefail; \
                     samtools faidx {}",
        &reference_path
    );
    debug!("Queuing cmd_string: {}", cmd_string);

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd_string)
        .stdout(Stdio::piped())
        .output()
        .expect("Unable to execute bash");

    return bio::io::fasta::IndexedReader::from_file(&reference_path)
        .expect("Unable to generate index");
}

pub fn galah_command_line_definition(
) -> galah::cluster_argument_parsing::GalahClustererCommandDefinition {
    galah::cluster_argument_parsing::GalahClustererCommandDefinition {
        dereplication_ani_argument: "dereplication-ani".to_string(),
        dereplication_prethreshold_ani_argument: "dereplication-prethreshold-ani".to_string(),
        dereplication_quality_formula_argument: "dereplication-quality-formula".to_string(),
        dereplication_precluster_method_argument: "dereplication-precluster-method".to_string(),
    }
}

pub fn dereplicate(m: &clap::ArgMatches, genome_fasta_files: &Vec<String>) -> Vec<String> {
    info!(
        "Found {} genomes specified before dereplication",
        genome_fasta_files.len()
    );

    // Generate clusterer and check for dependencies
    let clusterer = galah::cluster_argument_parsing::generate_galah_clusterer(
        genome_fasta_files,
        &m,
        &galah_command_line_definition(),
    )
    .expect("Failed to parse galah clustering arguments correctly");
    galah::external_command_checker::check_for_dependencies();
    info!("Dereplicating genome at {}% ANI ..", clusterer.ani * 100.);

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

pub fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    trace!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 std::str::from_utf8(target_name).unwrap()));
    return &target_name[(0..offset)];
}

pub fn retrieve_genome_from_contig<'a>(
    target_name: &'a [u8],
    genomes_and_contigs: &'a GenomesAndContigs,
    reference_map: &'a HashMap<usize, String>,
) -> (String, usize) {
    let genome_from_contig = || -> &'a String {
        genomes_and_contigs
            .genome_of_contig(&str::from_utf8(&target_name).unwrap().to_string())
            .expect(&format!(
                "Found invalid contig in bam, {:?}. \
                Please provide corresponding reference genomes",
                str::from_utf8(&target_name).unwrap()
            ))
    };

    // Concatenated references have the reference file name in front of the contig name
    // separated by the "~" symbol by default.
    // TODO: Parse as a separator value to this function
    let reference_stem = match str::from_utf8(&target_name).unwrap().splitn(2, "~").next() {
        Some(ref_stem) => ref_stem,
        None => genome_from_contig(),
    };

    debug!("possible reference stem {:?}", reference_stem);
    let ref_idx = match genomes_and_contigs.genome_index(&reference_stem.to_string()) {
        Some(idx) => idx,
        None => genomes_and_contigs
            .genome_index(genome_from_contig())
            .expect("Unable to parse genome name"),
    };
    debug!("Actual reference idx {:?}", ref_idx);

    let reference = reference_map
        .get(&ref_idx)
        .expect("Unable to retrieve reference path")
        .clone();
    (reference, ref_idx)
}

// Splits a contig name based on the ~
pub fn split_contig_name(target_name: &Vec<u8>) -> String {
    String::from_utf8(target_name.clone())
        .unwrap()
        .splitn(2, "~")
        .skip(1)
        .next()
        .unwrap_or(std::str::from_utf8(&target_name).unwrap())
        .to_string()
}

pub fn retrieve_reference_index_from_contig(
    target_name: &Vec<u8>,
    genomes_and_contigs: &GenomesAndContigs,
) -> usize {
    let target_name_str = String::from_utf8(target_name.clone()).unwrap();

    match genomes_and_contigs.genome_index_of_contig(&target_name_str) {
        Some(idx) => idx,
        None => {
            let split_name = split_contig_name(target_name);
            genomes_and_contigs
                .genome_index_of_contig(&split_name)
                .unwrap()
        }
    }
}

pub fn generate_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    discard_unmapped: bool,
    mapping_options: Option<&str>,
    include_reference_in_stoit_name: bool,
) -> StreamingNamedBamReaderGenerator {
    let tmp_dir = TempDir::new("lorikeet_fifo").expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let mapping_log = tempfile::NamedTempFile::new().expect(&format!(
        "Failed to create {:?} log tempfile",
        mapping_program
    ));
    let samtools2_log =
        tempfile::NamedTempFile::new().expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log =
        tempfile::NamedTempFile::new().expect("Failed to create cache samtools view log tempfile");

    let cached_bam_file_args = match cached_bam_file {
        Some(path) => {
            format!(
                "|tee {:?} |samtools view {} -@ {} -b -o '{}' 2>{}",
                // tee
                fifo_path,
                // samtools view
                match discard_unmapped {
                    true => "-F4",
                    false => "",
                },
                threads - 1,
                path,
                samtools_view_cache_log
                    .path()
                    .to_str()
                    .expect("Failed to convert tempfile path to str")
            )
        }
        None => format!("> {:?}", fifo_path),
    };

    let mapping_command = build_mapping_command(
        mapping_program,
        read_format,
        threads,
        read1_path,
        reference,
        read2_path,
        mapping_options,
    );
    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("lorikeet-samtools-sort")
        .tempfile_in(tmp_dir.path())
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} {}\
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         {}",
        // Mapping program
        mapping_command,
        mapping_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // remove extraneous @SQ lines
        match mapping_program {
            MappingProgram::BWA_MEM | MappingProgram::NGMLR_ONT | MappingProgram::NGMLR_PB => {
                ""
            }
            // Required because of https://github.com/lh3/minimap2/issues/527
            _ => " | remove_minimap2_duplicated_headers",
        },
        // samtools
        bwa_sort_prefix
            .path()
            .to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // Caching (or not)
        cached_bam_file_args
    );
    debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    let mut log_descriptions = vec![
        format!("{:?}", mapping_program).to_string(),
        "samtools sort".to_string(),
    ];
    let mut log_files = vec![mapping_log, samtools2_log];
    if cached_bam_file.is_some() {
        log_descriptions.push("samtools view for cache".to_string());
        log_files.push(samtools_view_cache_log);
    }

    let stoit_name = match include_reference_in_stoit_name {
        true => {
            std::path::Path::new(reference)
                .file_name()
                .expect("Unable to convert reference to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
                + "/"
        }
        false => "".to_string(),
    } + &std::path::Path::new(read1_path)
        .file_name()
        .expect("Unable to convert read1 name to file name")
        .to_str()
        .expect("Unable to covert file name into str")
        .to_string();

    return StreamingNamedBamReaderGenerator {
        stoit_name: stoit_name,
        tempdir: tmp_dir,
        fifo_path: fifo_path,
        cache_path: match cached_bam_file {
            Some(cache_path) => cache_path.to_string(),
            None => "".to_string(),
        },
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c \"{}\"", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files: log_files,
    };
}

pub fn generate_filtered_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    bwa_options: Option<&str>,
    discard_unmapped: bool,
    include_reference_in_stoit_name: bool,
) -> StreamingFilteredNamedBamReaderGenerator {
    let streaming = generate_named_bam_readers_from_reads(
        mapping_program,
        reference,
        read1_path,
        read2_path,
        read_format,
        threads,
        cached_bam_file,
        discard_unmapped,
        bwa_options,
        include_reference_in_stoit_name,
    );
    return StreamingFilteredNamedBamReaderGenerator {
        stoit_name: streaming.stoit_name,
        tempdir: streaming.tempdir,
        fifo_path: streaming.fifo_path,
        cache_path: streaming.cache_path,
        pre_processes: streaming.pre_processes,
        command_strings: streaming.command_strings,
        log_file_descriptions: streaming.log_file_descriptions,
        log_files: streaming.log_files,
        flag_filters: flag_filters,
        min_aligned_length_single: min_aligned_length_single,
        min_percent_identity_single: min_percent_identity_single,
        min_aligned_percent_single: min_aligned_percent_single,
        min_aligned_length_pair: min_aligned_length_pair,
        min_percent_identity_pair: min_percent_identity_pair,
        min_aligned_percent_pair: min_aligned_percent_pair,
    };
}

pub fn retrieve_reference(concatenated_genomes: &Option<String>) -> IndexedReader<File> {
    let reference = match concatenated_genomes {
        Some(reference_path) => match bio::io::fasta::IndexedReader::from_file(&reference_path) {
            Ok(reader) => reader,
            Err(_e) => generate_faidx(&reference_path),
        },
        None => panic!("Concatenated reference file does not exist"),
    };

    reference
}

pub fn fetch_contig_from_reference(
    reference: &mut IndexedReader<File>,
    contig_name: &Vec<u8>,
    genomes_and_contigs: &GenomesAndContigs,
    ref_idx: usize,
) {
    match reference.fetch_all(std::str::from_utf8(&contig_name[..]).unwrap()) {
        Ok(reference) => reference,
        Err(_e) => match reference.fetch_all(&format!(
            "{}~{}",
            &genomes_and_contigs.genomes[ref_idx],
            std::str::from_utf8(&contig_name[..]).unwrap()
        )) {
            Ok(reference) => reference,
            Err(e) => {
                println!(
                    "Cannot read sequence from reference {} {:?}",
                    format!(
                        "{}~{}",
                        &genomes_and_contigs.genomes[ref_idx],
                        std::str::from_utf8(&contig_name[..]).unwrap()
                    ),
                    e,
                );
                std::process::exit(1);
            }
        },
    };
}

pub fn read_sequence_to_vec(
    ref_seq: &mut Vec<u8>,
    reference: &mut IndexedReader<File>,
    contig_name: &Vec<u8>,
) {
    match reference.read(ref_seq) {
        Ok(reference) => reference,
        Err(e) => {
            println!(
                "Cannot read sequence from reference {} {:?}",
                std::str::from_utf8(&contig_name[..]).unwrap(),
                e,
            );
            std::process::exit(1)
        }
    };
}

pub fn mean(data: &[i32]) -> Option<f32> {
    let sum = data.par_iter().sum::<i32>() as f32;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f32),
        _ => None,
    }
}

pub fn std_deviation(data: &[i32]) -> Option<f32> {
    match (mean(data), data.len()) {
        (Some(data_mean), count) if count > 0 => {
            let variance = data
                .par_iter()
                .map(|value| {
                    let diff = data_mean - (*value as f32);

                    diff * diff
                })
                .sum::<f32>()
                / count as f32;

            Some(variance.sqrt())
        }
        _ => None,
    }
}
