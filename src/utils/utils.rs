use crate::*;
use coverm::bam_generator::*;
use coverm::mapping_index_maintenance;
use coverm::mapping_parameters::*;
use coverm::FlagFilter;
use processing::lorikeet_engine::ReadType;

use coverm::mapping_index_maintenance::MappingIndex;
use rayon::prelude::*;
use std::str;
use tempdir::TempDir;
use tempfile::NamedTempFile;

pub const NUMERICAL_EPSILON: f64 = 1e-3;
pub const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";
pub const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

// pub fn log10_binomial_coefficient(n: i64, k: i64) -> i64 {
//
// }

// pub fn factorial<T: Sized + Send + Add + Div + Mul + PartialEq + PartialOrd>(x: T) -> T {
//     if let Some(mut factorial) =
// }

/// Finds the first occurence of element in a slice
pub fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
where
    T: std::cmp::PartialEq<T>,
{
    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index);
        }
        index += 1;
    }
    return Err("Element not found in slice");
}

// pub fn finish_and_clear()

/**
 * Make all combinations of N size of objects
 *
 * if objects = [A, B, C]
 * if N = 1 => [[A], [B], [C]]
 * if N = 2 => [[A, A], [B, A], [C, A], [A, B], [B, B], [C, B], [A, C], [B, C], [C, C]]
 *
 * @param objects         list of objects
 * @param n               size of each combination
 * @param withReplacement if false, the resulting permutations will only contain unique objects from objects
 * @return a list with all combinations with size n of objects.
 */
pub fn make_permutations<T: Clone + PartialEq>(
    objects: &[T],
    n: usize,
    with_replacement: bool,
) -> Vec<Vec<T>> {
    let mut combinations = Vec::new();

    if n == 1 {
        for o in objects {
            combinations.push(vec![o.clone()]);
        }
    } else if n > 1 {
        let sub = make_permutations(objects, n - 1, with_replacement);
        for sub_i in sub {
            for a in objects {
                if with_replacement || !sub_i.contains(a) {
                    combinations.push(cons(a.clone(), sub_i.clone()));
                }
            }
        }
    }

    return combinations;
}

pub fn cons<T>(elt: T, l: Vec<T>) -> Vec<T> {
    let mut l2 = vec![elt];
    l2.extend(l);
    return l2;
}

pub fn clean_sample_name(sample_idx: usize, samples: &[String]) -> &str {
    if samples[sample_idx].contains(".tmp") {
        samples[sample_idx]
            .split("/.tmp")
            .skip(1)
            .next()
            .unwrap()
            .split(".fna.")
            .skip(1)
            .next()
            .unwrap()
    } else {
        &samples[sample_idx]
    }
}

pub fn get_cleaned_sample_names(samples: &[String]) -> Vec<&str> {
    (0..samples.len())
        .into_iter()
        .map(|i| clean_sample_name(i, samples))
        .collect()
}

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
    debug!("Reference Tempfile: {:?}", &reference_tempfile);
    let params = match readtype {
        &ReadType::Short => MappingParameters::generate_from_clap(
            m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Long => MappingParameters::generate_longread_from_clap(
            m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Assembly => MappingParameters::generate_assembly_from_clap(
            m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
    };

    let mut generator_set = vec![];
    for reference_wise_params in params {
        debug!("Ref Wise Params: {:?}", &reference_wise_params.len());
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
            index,
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
) -> (
    Vec<coverm::bam_generator::StreamingNamedBamReaderGenerator>,
    Vec<Option<Box<dyn MappingIndex>>>,
) {
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

    return (long_generators, indices);
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
        Some("bwa-mem2") => MappingProgram::BWA_MEM2,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-hifi") => MappingProgram::MINIMAP2_HIFI,
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
        MappingProgram::BWA_MEM2 => {
            external_command_checker::check_for_bwa_mem2();
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
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
            m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Long => MappingParameters::generate_longread_from_clap(
            m,
            mapping_program,
            &reference_tempfile,
            &references,
        ),
        &ReadType::Assembly => MappingParameters::generate_assembly_from_clap(
            m,
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
            index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

pub fn setup_mapping_index(
    reference_wise_params: &SingleReferenceMappingParameters,
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
) -> Option<Box<dyn coverm::mapping_index_maintenance::MappingIndex>> {
    match mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
            Some(coverm::mapping_index_maintenance::generate_bwa_index(
                reference_wise_params.reference,
                None,
                mapping_program,
            ))
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
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
                Some(coverm::mapping_index_maintenance::generate_minimap2_index(
                    reference_wise_params.reference,
                    Some(m.value_of("threads").unwrap().parse::<usize>().unwrap()),
                    Some(m.value_of("minimap2-params").unwrap_or("")),
                    mapping_program,
                ))
            }
        }
        MappingProgram::NGMLR_PB | MappingProgram::NGMLR_ONT => {
            // NGMLR won't let us just create mapping index, so we will use the --skip-write parameter
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
            error!(
                "Cache directory {} does not appear to be writeable, not continuing",
                cache_directory
            );
            process::exit(1);
        } else {
            info!(
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
                        error!(
                            "The parent directory of the (currently non-existent) \
                             cache directory {} is not writeable, not continuing",
                            cache_directory
                        );
                        process::exit(1);
                    } else {
                        info!("Creating cache directory {}", cache_directory);
                        std::fs::create_dir(path).expect("Unable to create cache directory");
                    }
                } else {
                    error!(
                        "The parent directory of the cache directory {} does not \
                         yet exist, so not creating that cache directory, and not continuing.",
                        cache_directory
                    );
                    process::exit(1);
                }
            }
            None => {
                error!("Cannot create root directory {}", cache_directory);
                process::exit(1);
            }
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        error!(
            "Failed to create test file in bam cache directory: {}",
            tf_result.err().unwrap()
        );
        process::exit(1);
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
                include_supplementary: !m.is_present("exclude-supplementary"),
            },
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => m
                    .value_of("min-read-aligned-length")
                    .unwrap()
                    .parse()
                    .unwrap(),
                false => 0,
            },
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => m
                    .value_of("min-read-aligned-length-pair")
                    .unwrap()
                    .parse()
                    .unwrap(),
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

pub fn table_roff(strings: &[&[&str]]) -> String {
    //start with a new line so the first .IP starts at the first char of the row
    let mut s: String = "\n.TS\n\
        tab(@);\n"
        .to_string();
    for row in strings {
        for _ in *row {
            s.push_str("l ");
        }
        break;
    }
    s.push_str(".\n");

    let mut first_row = true;
    for e in strings {
        let mut first_column = true;
        for cell in *e {
            if !first_column {
                s.push_str("@");
            }
            s.push_str("T{\n");
            if !first_column {
                s.push_str("|")
            } else {
                first_column = false
            }
            s.push_str(cell.clone());
            s.push_str("\nT}");
        }
        s.push_str("\n");
        if first_row {
            first_row = false;
            s.push_str("T{\n----------------------\nT}@T{\n|----------------------------------------\nT}\n");
        }
    }
    s.push_str(".TE\n");
    s
}
