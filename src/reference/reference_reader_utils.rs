use bio::io::fasta::IndexedReader;
use glob::glob;
use needletail::parse_fastx_file;
use std::process::{Stdio, exit, self};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;
use tempfile::NamedTempFile;

use crate::external_command_checker;
use crate::bam_parsing::mapping_index_maintenance::generate_concatenated_fasta_file;
use crate::utils::utils::find_first;

// lazy_static! {
//     static ref GALAH_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
//         GalahClustererCommandDefinition {
//             dereplication_ani_argument: "ani".to_string(),
//             dereplication_prethreshold_ani_argument: "precluster-ani".to_string(),
//             dereplication_quality_formula_argument: "quality-formula".to_string(),
//             dereplication_precluster_method_argument: "precluster-method".to_string(),
//             dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
//             dereplication_fraglen_argument: "fragment-length".to_string(),
//             dereplication_output_cluster_definition_file: "output-cluster-definition".to_string(),
//             dereplication_output_representative_fasta_directory:
//                 "output-representative-fasta-directory".to_string(),
//             dereplication_output_representative_fasta_directory_copy:
//                 "output-representative-fasta-directory-copy".to_string(),
//             dereplication_output_representative_list: "output-representative-list".to_string(),
//         }
//     };
// }

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
        return &target_name[0..offset];
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

    pub fn setup_genome_fasta_files(
        m: &clap::ArgMatches,
    ) -> (Option<NamedTempFile>, Option<GenomesAndContigs>) {
        let genome_fasta_files_opt = {
            match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(m, false) {
                Ok(paths) => {
                    if paths.len() == 0 {
                        error!("Genome paths were described, but ultimately none were found");
                        exit(1);
                    }
                    // if m.is_present("checkm-tab-table") || m.is_present("genome-info") {
                    //     let genomes_after_filtering =
                    //         galah::cluster_argument_parsing::filter_genomes_through_checkm(
                    //             &paths,
                    //             m,
                    //             &Self::galah_command_line_definition(),
                    //         )
                    //         .expect("Error parsing CheckM-related options");
                    //     info!(
                    //         "After filtering by CheckM, {} genomes remained",
                    //         genomes_after_filtering.len()
                    //     );
                    //     if genomes_after_filtering.len() == 0 {
                    //         error!("All genomes were filtered out, so none remain to be mapped to");
                    //         exit(1);
                    //     }
                    //     Some(
                    //         genomes_after_filtering
                    //             .iter()
                    //             .map(|s| s.to_string())
                    //             .collect(),
                    //     )
                    // } else {
                    Some(paths)
                    // }
                }
                Err(_) => None,
            }
        };

        // debug!("Found paths {:?}", &genome_fasta_files_opt);

        let (concatenated_genomes, genomes_and_contigs_option) = match m.contains_id("genome-fasta-files") {
            true => match genome_fasta_files_opt {
                Some(genome_paths) => (
                    Some(
                        generate_concatenated_fasta_file(
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
                let dereplicated_genomes: Vec<String> = if m.contains_id("dereplicate") {
                    // Self::dereplicate(m, &genome_fasta_files_opt.unwrap())
                    genome_fasta_files_opt.unwrap()
                } else {
                    genome_fasta_files_opt.unwrap()
                };
                // debug!("Profiling {} genomes", dereplicated_genomes.len());

                let list_of_genome_fasta_files = &dereplicated_genomes;

                (
                    Some(
                        generate_concatenated_fasta_file(
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

        // debug!("Found genome_and_contigs {:?}", &genomes_and_contigs_option);
        return (concatenated_genomes, genomes_and_contigs_option);
    }

    pub fn parse_references(m: &clap::ArgMatches) -> Vec<String> {
        let references = match m.get_many::<String>("genome-fasta-files") {
            Some(vec) => {
                let reference_paths = vec.map(|p| p.to_string()).collect::<Vec<String>>();
                // debug!("Reference files {:?}", reference_paths);
                reference_paths
            }
            None => match m.get_one::<String>("genome-fasta-directory") {
                Some(path) => {
                    let ext = m.get_one::<String>("genome-fasta-extension").unwrap();
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
                    // debug!("Reference files {:?}", reference_paths);
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
        match m.contains_id("genome-definition") {
            true => Some(read_genome_definition_file(
                m.get_one::<String>("genome-definition").unwrap(),
            )),
            false => Some(read_genome_fasta_files(
                &genome_fasta_files,
                false,
            )),
        }
    }

    pub fn generate_faidx(reference_path: &str) -> IndexedReader<File> {
        // debug!("Generating reference index");
        let cmd_string = format!(
            "set -e -o pipefail; \
            samtools faidx {}",
            reference_path
        );
        // debug!("Queuing cmd_string: {}", cmd_string);
        // check if {reference_path}.fai exists 
        let fai_path = format!("{}.fai", reference_path);
        if Path::new(&fai_path).exists() {
            // debug!("Found existing index file at {}", fai_path);
            return IndexedReader::from_file(&reference_path).expect("Unable to generate index");
        }
                    
        external_command_checker::check_for_samtools();
        std::process::Command::new("bash")
            .arg("-c")
            .arg(&cmd_string)
            .stdout(Stdio::piped())
            .output()
            .expect("Unable to execute bash");

        IndexedReader::from_file(&reference_path).expect("Unable to generate index")
    }

    // pub fn galah_command_line_definition() -> &'static GalahClustererCommandDefinition {
    //     &*GALAH_COMMAND_DEFINITION
    // }

    // pub fn dereplicate(m: &clap::ArgMatches, genome_fasta_files: &Vec<String>) -> Vec<String> {
    //     info!(
    //         "Found {} genomes specified before dereplication",
    //         genome_fasta_files.len()
    //     );

    //     // Generate clusterer and check for dependencies
    //     let clusterer = galah::cluster_argument_parsing::generate_galah_clusterer(
    //         genome_fasta_files,
    //         &m,
    //         Self::galah_command_line_definition(),
    //     )
    //     .expect("Failed to parse galah clustering arguments correctly");
    //     galah::external_command_checker::check_for_dependencies();
    //     info!("Dereplicating genome at {}% ANI ..", clusterer.ani * 100.);

    //     let cluster_indices = clusterer.cluster();
    //     info!(
    //         "Finished dereplication, finding {} representative genomes.",
    //         cluster_indices.len()
    //     );
    //     debug!("Found cluster indices: {:?}", cluster_indices);
    //     let reps = cluster_indices
    //         .iter()
    //         .map(|cluster| genome_fasta_files[cluster[0]].clone())
    //         .collect::<Vec<_>>();
    //     debug!("Found cluster representatives: {:?}", reps);

    //     if m.is_present("output-dereplication-clusters") {
    //         let path = m.value_of("output-dereplication-clusters").unwrap();
    //         info!("Writing dereplication cluster memberships to {}", path);
    //         let mut f = std::fs::File::create(path)
    //             .expect("Error creating dereplication cluster output file");
    //         for cluster in cluster_indices.iter() {
    //             let rep = cluster[0];
    //             for member in cluster {
    //                 writeln!(
    //                     f,
    //                     "{}\t{}",
    //                     genome_fasta_files[rep], genome_fasta_files[*member]
    //                 )
    //                 .expect("Failed to write a specific line to dereplication cluster file");
    //             }
    //         }
    //     }
    //     reps
    // }
}

#[derive(Debug, Clone)]
pub struct GenomesAndContigs {
    pub genomes: Vec<String>,
    pub contigs: usize,
}

impl GenomesAndContigs {
    pub fn new() -> Self {
        GenomesAndContigs {
            genomes: Vec::new(),
            contigs: 0,
        }
    }

    pub fn establish_genome(&mut self, genome_name: String) -> usize {
        let index = self.genomes.len();
        self.genomes.push(genome_name);
        return index;
    }

    pub fn genome_index(&self, genome_name: &str) -> Option<usize> {
        match find_first(self.genomes.as_slice(), genome_name.to_string()) {
            Ok(index) => Some(index),
            Err(_) => {
                debug!("Genome {} not found in {:?}", genome_name, self.genomes);
                None
            },
        }
    }
}

pub fn read_genome_fasta_files(
    fasta_file_paths: &Vec<&str>,
    _use_full_sequence_name: bool,
) -> GenomesAndContigs {
    let mut contig_to_genome = GenomesAndContigs::new();

    // NOTE: A lot of this code is shared with mapping_index_maintenance.rs#generate_concatenated_fasta_file
    for file in fasta_file_paths {
        let path = Path::new(file);
        let mut reader =
            parse_fastx_file(path).expect(&format!("Unable to read fasta file {}", file));

        // Remove .gz .bz .xz from file names if present
        let genome_name1 =
            String::from(path.to_str().expect("File name string conversion problem"));
        // if genome_name1 ends with .gz, .bz, or .xz, the Rust HTSLIB indexed reader
        // will not be able to parse it. So we need to gracefully exit and inform the user
        // that they need to decompress the file first. Additionally, if they have provided BAM
        // files the BAM files will need to be regenerated with decompressed input or reheadered
        // to remove the .gz, .bz, or .xz extension.
        if genome_name1.ends_with(".gz") || genome_name1.ends_with(".bz") || genome_name1.ends_with(".xz") {
            error!("The genome file {} is compressed. Please decompress the file before running lorikeet.", genome_name1);
            error!("If you are using BAM files, please regenerate the BAM files with decompressed input");
            error!("Or reheader the BAM files to remove the .gz, .bz, or .xz extension from genome names");
            process::exit(1);
        }
        let path1 = Path::new(&genome_name1);

        let genome_name = String::from(
            path1
                .file_stem()
                .expect("Problem while determining file stem")
                .to_str()
                .expect("File name string conversion problem"),
        );
        if contig_to_genome.genome_index(&genome_name).is_some() {
            error!("The genome name {} was derived from >1 file", genome_name);
            exit(1);
        }
        let _genome_index = contig_to_genome.establish_genome(genome_name);
        while let Some(record) = reader.next() {
            let record_expected =
                record.expect(&format!("Failed to parse record in fasta file {:?}", path));

            if record_expected.format() != needletail::parser::Format::Fasta {
                panic!(
                    "File {:?} is not a fasta file, but a {:?}",
                    path,
                    record_expected.format()
                );
            }

            contig_to_genome.contigs += 1;
        }
    }
    return contig_to_genome;
}

pub fn read_genome_definition_file(definition_file_path: &str) -> GenomesAndContigs {
    let f = std::fs::File::open(definition_file_path).expect(&format!(
        "Unable to find/read genome definition file {}",
        definition_file_path
    ));
    let file = std::io::BufReader::new(&f);
    let mut genome_to_contig: HashMap<String, Vec<String>> = HashMap::new();
    // Maintain the same order as the input file.
    let mut genome_order: Vec<String> = vec![];
    let mut contig_to_genome = GenomesAndContigs::new();
    for line_res in file.lines() {
        let line = line_res.expect("Read error on genome definition file");
        let v: Vec<&str> = line.split("\t").collect();
        if v.len() == 2 {
            let genome = v[0].trim();
            let contig = v[1]
                .split_ascii_whitespace()
                .next()
                .expect("Failed to split contig name by whitespace in genome definition file");
            contig_to_genome.contigs += 1;

            if genome_to_contig.contains_key(genome) {
                genome_to_contig
                    .get_mut(genome)
                    .unwrap()
                    .push(contig.to_string());
            } else {
                genome_to_contig.insert(genome.to_string(), vec![contig.to_string()]);
                genome_order.push(genome.to_string());
            }
        } else if v.len() == 0 {
            continue;
        } else {
            error!(
                "The line \"{}\" in the genome definition file is not a \
                    genome name and contig name separated by a tab",
                line
            );
            exit(1);
        }
    }

    info!(
        "Found {} contigs assigned to {} different genomes from \
           the genome definition file",
        contig_to_genome.contigs,
        genome_to_contig.len()
    );

    for genome in genome_order {
        let _contigs = &genome_to_contig[&genome];
        let _genome_index = contig_to_genome.establish_genome(genome);
    }
    return contig_to_genome;
}