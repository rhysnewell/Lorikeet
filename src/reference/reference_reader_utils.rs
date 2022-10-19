use bio::io::fasta::IndexedReader;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::genomes_and_contigs::*;
use external_command_checker;
use galah::cluster_argument_parsing::GalahClustererCommandDefinition;
use glob::glob;
use process::Stdio;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process;
use tempfile::NamedTempFile;

lazy_static! {
    static ref GALAH_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
        GalahClustererCommandDefinition {
            dereplication_ani_argument: "ani".to_string(),
            dereplication_prethreshold_ani_argument: "precluster-ani".to_string(),
            dereplication_quality_formula_argument: "quality-formula".to_string(),
            dereplication_precluster_method_argument: "precluster-method".to_string(),
            dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
            dereplication_fraglen_argument: "fragment-length".to_string(),
            dereplication_output_cluster_definition_file: "output-cluster-definition".to_string(),
            dereplication_output_representative_fasta_directory:
                "output-representative-fasta-directory".to_string(),
            dereplication_output_representative_fasta_directory_copy:
                "output-representative-fasta-directory-copy".to_string(),
            dereplication_output_representative_list: "output-representative-list".to_string(),
        }
    };
}

pub struct ReferenceReaderUtils {}

impl ReferenceReaderUtils {
    pub fn retrieve_reference(concatenated_genomes: &Option<String>) -> IndexedReader<File> {
        let reference = match concatenated_genomes {
            Some(reference_path) => match IndexedReader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(_e) => Self::generate_faidx(reference_path.as_str()),
            },
            None => panic!("Concatenated reference file does not exist"),
        };

        reference
    }

    pub fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
        let target_name = target_names[tid as usize];
        trace!("target name {:?}, separator {:?}", target_name, split_char);
        let offset = find_first(target_name, split_char).unwrap_or_else(|_|
            panic!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                     std::str::from_utf8(target_name).unwrap()));
        return &target_name[(0..offset)];
    }

    // pub fn get_reference_path(args: &Vec<&str>, genomes_and_contigs: &GenomesAndContigs, ref_idx: usize) -> String {
    //     let ref_stub = &genomes_and_contigs.genomes[ref_idx];
    //
    // }

    pub fn retrieve_genome_from_contig<'a>(
        target_name: &'a [u8],
        genomes_and_contigs: &'a GenomesAndContigs,
        reference_map: &'a HashMap<usize, String>,
    ) -> (String, usize) {
        let genome_from_contig = || -> &'a String {
            genomes_and_contigs
                .genome_of_contig(&std::str::from_utf8(&target_name).unwrap().to_string())
                .unwrap_or_else(|| {
                    panic!(
                        "Found invalid contig in bam, {:?}. \
                Please provide corresponding reference genomes",
                        std::str::from_utf8(&target_name).unwrap()
                    )
                })
        };

        // Concatenated references have the reference file name in front of the contig name
        // separated by the "~" symbol by default.
        // TODO: Parse as a separator value to this function
        let reference_stem = match std::str::from_utf8(&target_name)
            .unwrap()
            .splitn(2, '~')
            .next()
        {
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
            .splitn(2, '~')
            .nth(1)
            .unwrap_or_else(|| std::str::from_utf8(&target_name).unwrap())
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
                let split_name = Self::split_contig_name(target_name);
                genomes_and_contigs
                    .genome_index_of_contig(&split_name)
                    .unwrap()
            }
        }
    }

    pub fn setup_genome_fasta_files(
        m: &clap::ArgMatches,
    ) -> (Option<NamedTempFile>, Option<GenomesAndContigs>) {
        let genome_fasta_files_opt = {
            match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(m, false) {
                Ok(paths) => {
                    if paths.len() == 0 {
                        error!("Genome paths were described, but ultimately none were found");
                        process::exit(1);
                    }
                    if m.is_present("checkm-tab-table") || m.is_present("genome-info") {
                        let genomes_after_filtering =
                            galah::cluster_argument_parsing::filter_genomes_through_checkm(
                                &paths,
                                m,
                                &Self::galah_command_line_definition(),
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
                    Self::extract_genomes_and_contigs_option(
                        m,
                        &genome_paths.iter().map(|s| s.as_str()).collect(),
                    ),
                ),
                None => (None, None),
            },
            false => {
                // Dereplicate if required
                // TODO: Properly implement dereplication in cli.rs and make sure function works
                let dereplicated_genomes: Vec<String> = if m.is_present("dereplicate") {
                    // Self::dereplicate(m, &genome_fasta_files_opt.unwrap())
                    genome_fasta_files_opt.unwrap()
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
                    Self::extract_genomes_and_contigs_option(
                        m,
                        &dereplicated_genomes
                            // .clone()
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
                false,
            )),
        }
    }

    pub fn generate_faidx(reference_path: &str) -> IndexedReader<File> {
        external_command_checker::check_for_samtools();
        debug!("Generating reference index");
        let cmd_string = format!(
            "set -e -o pipefail; \
                     samtools faidx {}",
            reference_path
        );
        debug!("Queuing cmd_string: {}", cmd_string);

        std::process::Command::new("bash")
            .arg("-c")
            .arg(&cmd_string)
            .stdout(Stdio::piped())
            .output()
            .expect("Unable to execute bash");

        return IndexedReader::from_file(&reference_path).expect("Unable to generate index");
    }

    pub fn galah_command_line_definition() -> &'static GalahClustererCommandDefinition {
        &*GALAH_COMMAND_DEFINITION
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
            Self::galah_command_line_definition(),
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
            let mut f = std::fs::File::create(path)
                .expect("Error creating dereplication cluster output file");
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
}
