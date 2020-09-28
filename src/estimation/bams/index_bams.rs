use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use external_command_checker;
use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam;
use std::path::Path;

pub fn finish_bams<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bams: Vec<G>,
    n_threads: usize,
) {
    external_command_checker::check_for_gatk();
    let mut record: bam::record::Record = bam::Record::new();

    // progress bar
    let sty = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}");
    let pb1 = ProgressBar::new(bams.len() as u64);
    pb1.set_style(sty.clone());
    for bam_generator in bams {
        let mut bam = bam_generator.start();
        bam.set_threads(n_threads / 2);

        let path = bam.path().to_string();
        let stoit_name = bam.name().to_string().replace("/", ".");

        pb1.set_message(&format!(
            "Processing sample: {}",
            match &stoit_name[..4] {
                ".tmp" => &stoit_name[15..],
                _ => &stoit_name,
            },
        ));

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
                    .set_threads(n_threads / 2)
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

            debug!(
                "Finished mapping sample {}",
                match &stoit_name[..4] {
                    ".tmp" => &stoit_name[15..],
                    _ => &stoit_name,
                }
            );
            bam.finish();

            std::fs::copy(&tmp_bam.path(), &path).expect("Unable to move BAM");

            // index the bam file
            bam::index::build(
                &path,
                Some(&format!("{}.bai", path)),
                bam::index::Type::BAI,
                n_threads as u32,
            )
            .expect("Unable to index BAM");
        } else {
            // while bam
            //     .read(&mut record)
            //     .expect("Error while reading BAM record")
            //     == true
            // {
            //     do nothing
            // }

            debug!(
                "Read groups already present for sample {}",
                match &stoit_name[..4] {
                    "lori" => &stoit_name[16..],
                    _ => &stoit_name,
                }
            );

            if !Path::new(&format!("{}.bai", path)).exists() {
                bam::index::build(
                    &path,
                    Some(&format!("{}.bai", path)),
                    bam::index::Type::BAI,
                    n_threads as u32,
                )
                .expect(&format!("Unable to index bam at {}", &path));
            }
        }
        pb1.inc(1);
    }
    pb1.finish_with_message(&format!("Reads and BAM files processed..."));
}

pub fn recover_bams(
    m: &clap::ArgMatches,
    concatenated_genomes: &Option<String>,
    short_sample_count: usize,
    long_sample_count: usize,
    genomes_and_contigs: &GenomesAndContigs,
    _n_threads: u32,
    tmp_bam_file_cache: &Option<String>,
) -> Vec<String> {
    // Annoyingly read in bam file again
    let mut bam_readers = vec![];

    // This is going to catch cached longread bam files from mapping
    if m.is_present("bam-files") {
        let bam_paths = m
            .values_of("bam-files")
            .unwrap()
            .map(|b| b.to_string())
            .collect::<Vec<String>>();
        bam_readers.extend(bam_paths);
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
                            tmp_bam_file_cache.as_ref().unwrap()
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
                                tmp_bam_file_cache.as_ref().unwrap()
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

        let _bam_cnts = all_bam_paths.len();
        bam_readers.extend(all_bam_paths);
    }

    if m.is_present("longread-bam-files") {
        let bam_paths = m
            .values_of("longread-bam-files")
            .unwrap()
            .map(|b| b.to_string())
            .collect::<Vec<String>>();
        bam_readers.extend(bam_paths);
    } else if m.is_present("longreads") {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/long/*.bam",
                    match m.is_present("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache.as_ref().unwrap()
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
                                tmp_bam_file_cache.as_ref().unwrap()
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

        let _bam_cnts = all_bam_paths.len();
        bam_readers.extend(all_bam_paths);
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
    let bam_readers = bam_readers.into_iter().map(|b| b.to_string()).collect();
    return bam_readers;
}
