use crate::*;

use coverm::bam_generator::{self, *};
use coverm::shard_bam_reader::{self, ShardedBamReaderGenerator};
use coverm::mapping_parameters::*;
use coverm::genome_exclusion::*;
use coverm::genomes_and_contigs::*;
use coverm::mapping_index_maintenance;
use coverm::FlagFilter;

use std::fs::File;
use std::process::Stdio;
use std::io::Write;
use tempfile::NamedTempFile;
use glob::glob;

pub const NUMERICAL_EPSILON: f64 = 1e-3;
pub const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";
pub const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

pub fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    longread: bool,
    references: &'a Option<Vec<&'a str>>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params =
        if !longread {
            MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile, &references)
        } else {
            MappingParameters::generate_longread_from_clap(&m, mapping_program, &reference_tempfile, &references)
        };    let mut generator_set = vec![];
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

        let _n_samples = reference_wise_params.len() as u16;

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


pub fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
    longread: bool,
    references: &Option<Vec<&str>>,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params =
        if !longread {
            MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile, &references)
        } else {
            MappingParameters::generate_longread_from_clap(&m, mapping_program, &reference_tempfile, &references)
        };
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

pub fn setup_mapping_index(
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
        },
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

pub fn generate_cached_bam_file_name(directory: &str, reference: &str, read1_path: &str) -> String {
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

pub fn doing_metabat(m: &clap::ArgMatches) -> bool {
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

pub fn setup_genome_fasta_files(m: &clap::ArgMatches) -> (Option<NamedTempFile>, Option<GenomesAndContigs>){
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

    debug!("Found paths {:?}", &genome_fasta_files_opt);


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
    return (concatenated_genomes, genomes_and_contigs_option)
}

pub fn parse_references(m: &clap::ArgMatches) -> Vec<String> {
    let references = match m.values_of("reference") {
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
    return references
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
        &reference_path);
    debug!("Queuing cmd_string: {}", cmd_string);

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd_string)
        .stdout(Stdio::piped())
        .output()
        .expect("Unable to execute bash");

    return bio::io::fasta::IndexedReader::from_file(&reference_path).expect(
        "Unable to generate index")
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

pub fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    trace!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 std::str::from_utf8(target_name).unwrap()));
    return &target_name[(0..offset)];
}