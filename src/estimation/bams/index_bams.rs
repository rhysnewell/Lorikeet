use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use external_command_checker;
use glob::glob;
use rust_htslib::bam;
use tempdir::TempDir;
use tempfile::NamedTempFile;

pub fn finish_bams<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bams: Vec<G>,
    n_threads: usize,
) {
    external_command_checker::check_for_gatk();
    let mut record: bam::record::Record = bam::Record::new();
    for bam_generator in bams {
        let mut bam = bam_generator.start();
        bam.set_threads(n_threads);

        let path = bam.path().to_string();
        let stoit_name = bam.name().to_string().replace("/", ".");

        let header = bam.header();
        let mut tmp_header = bam::header::Header::from_template(&header);

        // Check if bam already has read group, if it doesn't then add one by writing to a new bam
        if !tmp_header.to_hashmap().contains_key("RG") {
            let tmp_bam = tempfile::NamedTempFile::new().unwrap();
            {
                // Setup tmp bam writer
                let mut tmp_header_record = bam::header::HeaderRecord::new("RG".as_bytes());
                // Push tags
                tmp_header_record.push_tag("ID".as_bytes(), &1);
                tmp_header_record.push_tag("SM".as_bytes(), &1);
                tmp_header_record.push_tag("LB".as_bytes(), &"N");
                tmp_header_record.push_tag("PL".as_bytes(), &"N");
                tmp_header_record.push_tag("PU".as_bytes(), &"N");

                tmp_header.push_record(&tmp_header_record);

                let mut bam_writer =
                    bam::Writer::from_path(tmp_bam.path(), &tmp_header, bam::Format::BAM)
                        .expect("Unable to create bam");

                bam_writer
                    .set_threads(n_threads)
                    .expect("Unable to set threads for BAM writer");
                bam_writer
                    .set_compression_level(bam::CompressionLevel::Uncompressed)
                    .expect("Unexpected compression level");

                let rg = bam::record::Aux::String("1".as_bytes());

                while bam
                    .read(&mut record)
                    .expect("Error while reading BAM record")
                    == true
                {
                    // push aux flags
                    record.push_aux("RG".as_bytes(), &rg);

                    // Write to bam
                    bam_writer.write(&record).expect("Unable to write to BAM");
                }
            }

            info!(
                "Finished Mapping Sample {}",
                match &stoit_name[..4] {
                    ".tmp" => &stoit_name[15..],
                    _ => &stoit_name,
                }
            );
            bam.finish();

            std::fs::copy(&tmp_bam.path(), &path).expect("Unable to move BAM");

            // index the bam file
            bam::index::build(&path, None, bam::index::Type::BAI, n_threads as u32)
                .expect("Unable to index BAM");
        } else {
            while bam
                .read(&mut record)
                .expect("Error while reading BAM record")
                == true
            {
                // do nothing
            }

            info!(
                "Finished Mapping Sample {}",
                match &stoit_name[..4] {
                    ".tmp" => &stoit_name[15..],
                    _ => &stoit_name,
                }
            );
            bam.finish();
        }
    }
}

pub fn recover_bams(
    m: &clap::ArgMatches,
    concatenated_genomes: &Option<NamedTempFile>,
    short_sample_count: usize,
    long_sample_count: usize,
    genomes_and_contigs: &GenomesAndContigs,
    n_threads: u32,
    tmp_bam_file_cache: &Option<TempDir>,
) -> Vec<IndexedBamFileNamedReader> {
    // Annoyingly read in bam file again
    let mut bam_readers = vec![];

    // This is going to catch cached longread bam files from mapping
    if m.is_present("bam-files") {
        let bam_paths = m.values_of("bam-files").unwrap().collect::<Vec<&str>>();
        bam_readers.extend(generate_indexed_named_bam_readers_from_bam_files(
            bam_paths, n_threads,
        ));
    } else if m.is_present("read1")
        | m.is_present("single")
        | m.is_present("coupled")
        | m.is_present("interleaved")
    {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/short/*.bam",
                    match m.is_present("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache
                                .as_ref()
                                .unwrap()
                                .path()
                                .to_str()
                                .unwrap()
                        }
                        true => {
                            m.value_of("bam-file-cache-directory").unwrap()
                        }
                    },
                );
                let bam_paths = glob(&cache)
                    .expect("Failed to read cache")
                    .map(|p| {
                        p.expect("Failed to read cached bam path")
                            .to_str()
                            .unwrap()
                            .to_string()
                    })
                    .collect::<Vec<String>>();
                all_bam_paths.extend(bam_paths);
            }
            None => {
                for ref_name in genomes_and_contigs.genomes.iter() {
                    let cache = format!(
                        "{}/short/{}*.bam",
                        match m.is_present("bam-file-cache-directory") {
                            false => {
                                tmp_bam_file_cache
                                    .as_ref()
                                    .unwrap()
                                    .path()
                                    .to_str()
                                    .unwrap()
                            }
                            true => {
                                m.value_of("bam-file-cache-directory").unwrap()
                            }
                        },
                        ref_name
                    );
                    let bam_paths = glob(&cache)
                        .expect("Failed to read cache")
                        .map(|p| {
                            p.expect("Failed to read cached bam path")
                                .to_str()
                                .unwrap()
                                .to_string()
                        })
                        .collect::<Vec<String>>();
                    all_bam_paths.extend(bam_paths);
                }
            }
        }

        debug!("Rereading in {}", all_bam_paths.len());
        let all_bam_paths = all_bam_paths
            .iter()
            .map(|p| p.as_str())
            .collect::<Vec<&str>>();

        let _bam_cnts = all_bam_paths.len();
        bam_readers.extend(generate_indexed_named_bam_readers_from_bam_files(
            all_bam_paths,
            n_threads,
        ));
    }

    if m.is_present("longread-bam-files") {
        let bam_paths = m
            .values_of("longread-bam-files")
            .unwrap()
            .collect::<Vec<&str>>();
        bam_readers.extend(generate_indexed_named_bam_readers_from_bam_files(
            bam_paths, n_threads,
        ));
    } else if m.is_present("longread") {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/long/*.bam",
                    match m.is_present("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache
                                .as_ref()
                                .unwrap()
                                .path()
                                .to_str()
                                .unwrap()
                        }
                        true => {
                            m.value_of("bam-file-cache-directory").unwrap()
                        }
                    },
                );
                let bam_paths = glob(&cache)
                    .expect("Failed to read cache")
                    .map(|p| {
                        p.expect("Failed to read cached bam path")
                            .to_str()
                            .unwrap()
                            .to_string()
                    })
                    .collect::<Vec<String>>();
                all_bam_paths.extend(bam_paths);
            }
            None => {
                for ref_name in genomes_and_contigs.genomes.iter() {
                    let cache = format!(
                        "{}/long/{}*.bam",
                        match m.is_present("bam-file-cache-directory") {
                            false => {
                                tmp_bam_file_cache
                                    .as_ref()
                                    .unwrap()
                                    .path()
                                    .to_str()
                                    .unwrap()
                            }
                            true => {
                                m.value_of("bam-file-cache-directory").unwrap()
                            }
                        },
                        ref_name
                    );
                    let bam_paths = glob(&cache)
                        .expect("Failed to read cache")
                        .map(|p| {
                            p.expect("Failed to read cached bam path")
                                .to_str()
                                .unwrap()
                                .to_string()
                        })
                        .collect::<Vec<String>>();
                    all_bam_paths.extend(bam_paths);
                }
            }
        }

        debug!("Rereading in {}", all_bam_paths.len());
        let all_bam_paths = all_bam_paths
            .iter()
            .map(|p| p.as_str())
            .collect::<Vec<&str>>();

        let _bam_cnts = all_bam_paths.len();
        bam_readers.extend(generate_indexed_named_bam_readers_from_bam_files(
            all_bam_paths,
            n_threads,
        ));
    }

    if bam_readers.len() == (short_sample_count + long_sample_count) {
        debug!("Reread in correct bam count")
    } else {
        panic!(format!(
            "Original sample count {} does not match new sample count {}, \
                please clear bam cache directory or ask for support at: github.com/rhysnewell/Lorikeet",
            (short_sample_count + long_sample_count),
            bam_readers.len(),
        ))
    }
    return bam_readers;
}
