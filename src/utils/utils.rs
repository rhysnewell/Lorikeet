use crate::*;

use coverm::bam_generator::*;
use coverm::genomes_and_contigs::*;
use coverm::mapping_index_maintenance;
use coverm::mapping_parameters::*;
use coverm::FlagFilter;

use nix::{sys::stat, unistd};
use rayon::prelude::*;
use std::collections::HashMap;
use std::str;
use tempdir::TempDir;
use tempfile::NamedTempFile;
use std::path::Path;
use std::sync::Arc;

pub const NUMERICAL_EPSILON: f64 = 1e-3;
pub const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";
pub const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

// pub fn log10_binomial_coefficient(n: i64, k: i64) -> i64 {
//
// }

// pub fn factorial<T: Sized + Send + Add + Div + Mul + PartialEq + PartialOrd>(x: T) -> T {
//     if let Some(mut factorial) =
// }

// pub fn finish_and_clear()

pub fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    readtype: &ReadType,
    references: &'a Option<Vec<&'a str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        match readtype {
            &ReadType::Long => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/long/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
            &ReadType::Short => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/short/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
            &ReadType::Assembly => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/assembly/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
        }
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = match readtype {
        &ReadType::Short => MappingParameters::generate_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Long => MappingParameters::generate_longread_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Assembly => MappingParameters::generate_assembly_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
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
                        readtype,
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
                        readtype,
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

pub fn long_generator_setup(
    m: &clap::ArgMatches,
    reference_tempfile: &Option<NamedTempFile>,
    references: &Option<Vec<&str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<coverm::bam_generator::StreamingNamedBamReaderGenerator> {
    // Perform mapping
    let mapping_program = parse_mapping_program(m.value_of("longread-mapper"));
    let readtype = ReadType::Long;
    external_command_checker::check_for_samtools();
    let generator_sets = get_streamed_bam_readers(
        m,
        mapping_program,
        reference_tempfile,
        &readtype,
        references,
        tmp_bam_file_cache,
    );
    let mut long_generators = vec![];
    let mut indices = vec![]; // Prevent indices from being dropped
    for set in generator_sets {
        indices.push(set.index);
        for g in set.generators {
            long_generators.push(g)
        }
    }

    return long_generators;
}

pub fn assembly_generator_setup(
    m: &clap::ArgMatches,
    reference_tempfile: &Option<NamedTempFile>,
    references: &Option<Vec<&str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<coverm::bam_generator::StreamingNamedBamReaderGenerator> {
    // Perform mapping
    let mapping_program = parse_mapping_program(Some("minimap2-assembly"));
    let readtype = ReadType::Assembly;
    external_command_checker::check_for_samtools();
    let generator_sets = get_streamed_bam_readers(
        m,
        mapping_program,
        reference_tempfile,
        &readtype,
        references,
        tmp_bam_file_cache,
    );
    let mut long_generators = vec![];
    let mut indices = vec![]; // Prevent indices from being dropped
    for set in generator_sets {
        indices.push(set.index);
        for g in set.generators {
            long_generators.push(g)
        }
    }

    return long_generators;
}

pub fn parse_mapping_program(mapper: Option<&str>) -> MappingProgram {
    let mapping_program = match mapper {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-assembly") => MappingProgram::MINIMAP2_ASS,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        Some("ngmlr-ont") => MappingProgram::NGMLR_ONT,
        Some("ngmlr-pb") => MappingProgram::NGMLR_PB,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!("Unexpected definition for --mapper: {:?}", mapper),
    };
    match mapping_program {
        MappingProgram::BWA_MEM => {
            external_command_checker::check_for_bwa();
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_ASS
        | MappingProgram::MINIMAP2_NO_PRESET => {
            external_command_checker::check_for_minimap2();
        }
        MappingProgram::NGMLR_ONT | MappingProgram::NGMLR_PB => {
            external_command_checker::check_for_ngmlr();
        }
    }
    return mapping_program;
}

pub fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
    readtype: &ReadType,
    references: &Option<Vec<&str>>,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        match readtype {
            &ReadType::Long => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/long/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
            &ReadType::Short => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/short/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
            &ReadType::Assembly => {
                setup_bam_cache_directory(&m.value_of("bam-file-cache-directory").unwrap());
                setup_bam_cache_directory(&format!(
                    "{}/assembly/",
                    m.value_of("bam-file-cache-directory").unwrap()
                ));
            }
        }
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = match readtype {
        &ReadType::Short => MappingParameters::generate_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Long => MappingParameters::generate_longread_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Assembly => MappingParameters::generate_assembly_from_clap(
            &m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
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
                        readtype,
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
                        readtype,
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
        | MappingProgram::MINIMAP2_ASS
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
    readtype: &ReadType,
) -> String {
    debug!(
        "Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
        directory, reference, read1_path
    );
    match readtype {
        &ReadType::Long => {
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
        }
        &ReadType::Short => {
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
        &ReadType::Assembly => {
            std::path::Path::new(&format!("{}/assembly/", directory))
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
