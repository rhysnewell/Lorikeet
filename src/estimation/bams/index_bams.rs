use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use glob::glob;
use tempdir::TempDir;
use tempfile::NamedTempFile;
use rust_htslib::bam;

pub fn finish_bams<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bams: Vec<G>,
    n_threads: usize,
) {
    let mut record: bam::record::Record = bam::Record::new();
    for bam_generator in bams {
        let mut bam = bam_generator.start();
        bam.set_threads(n_threads);

        while bam
            .read(&mut record)
            .expect("Error while reading BAM record")
            == true
        {
            // do nothing
        }
        bam.finish();
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
