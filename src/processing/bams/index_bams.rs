use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam;
use std::path::Path;
use std::result::Result::Err;
use std::time::Duration;

use crate::bam_parsing::bam_generator::*;
use crate::reference::reference_reader_utils::GenomesAndContigs;
use crate::utils::errors::BirdToolError;

/// Ensures mapping is completed for provided bams. If multiple references are provided and
/// the user has asked to run genomes in parallel, then the bams are split per reference to
/// avoid file locking when reading bams in parallel.
pub fn finish_bams<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bams: Vec<G>,
    n_threads: usize,
    _references: &GenomesAndContigs,
    _parallel_genomes: bool,
    mapping: bool,
) -> Result<(), BirdToolError> {
    let mut record: bam::Record = bam::Record::new();

    // progress bar
    let sty = match ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]")
    {
        Ok(s) => s,
        Err(e) => return Err(BirdToolError::DebugError(e.to_string())),
    };
    let pb1 = ProgressBar::new(bams.len() as u64);
    pb1.set_style(sty);
    pb1.enable_steady_tick(Duration::from_millis(200));
    for bam_generator in bams {
        let mut bam = bam_generator.start();
        // bam.set_threads(std::cmp::max(n_threads / 2, 1));

        let path = bam.path().to_string();
        let stoit_name = bam.name().to_string().replace("/", ".");

        pb1.set_message(format!(
            "Processing sample: {}",
            match &stoit_name[..4] {
                ".tmp" => &stoit_name[15..],
                _ => &stoit_name,
            },
        ));

        if mapping {
            while bam.read(&mut record).is_some() {
                continue;
            }

            bam.finish();
        }

        if !Path::new(&format!("{}.bai", path)).exists() || mapping {
            bam::index::build(
                &path,
                Some(&format!("{}.bai", path)),
                bam::index::Type::Bai,
                n_threads as u32,
            )
            .unwrap_or_else(|_| panic!("Unable to index bam at {}", &path));
        }
        // }
        pb1.inc(1);
    }
    pb1.finish_with_message(format!("Reads and BAM files processed..."));
    Ok(())
}

pub fn recover_bams(
    m: &clap::ArgMatches,
    concatenated_genomes: &Option<String>,
    short_sample_count: usize,
    long_sample_count: usize,
    genomes_and_contigs: &GenomesAndContigs,
    _n_threads: u32,
    tmp_bam_file_cache: &Option<String>,
    bams_are_split: bool,
    ref_idx: usize,
) -> Vec<String> {
    // Annoyingly read in bam file again
    let mut bam_readers = vec![];

    debug!("Bams are split {}", bams_are_split);

    // This is going to catch cached longread bam files from mapping
    if m.get_flag("bam-files") {
        let bam_paths = m
            .get_many::<String>("bam-files")
            .unwrap()
            .map(|b| b.to_string())
            .collect::<Vec<String>>();
        if !bams_are_split {
            bam_readers.extend(bam_paths);
        } else {
            for path in bam_paths {
                let mut path_split = path.rsplitn(2, '/');
                let path_suffix = path_split.next().unwrap();
                let path_prefix = path_split.next().unwrap();

                bam_readers.push(format!(
                    "{}/{}/{}",
                    path_prefix, &genomes_and_contigs.genomes[ref_idx], path_suffix
                ))
            }
        }
    } else if m.get_flag("read1")
        | m.get_flag("single")
        | m.get_flag("coupled")
        | m.get_flag("interleaved")
    {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/short/{}/*.bam",
                    match m.get_flag("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache.as_ref().unwrap()
                        }
                        true => {
                            m.get_one::<String>("bam-file-cache-directory").unwrap()
                        }
                    },
                    if bams_are_split {
                        &genomes_and_contigs.genomes[ref_idx]
                    } else {
                        ""
                    }
                );
                debug!("Cache: {}", &cache);
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
                        match m.get_flag("bam-file-cache-directory") {
                            false => {
                                tmp_bam_file_cache.as_ref().unwrap()
                            }
                            true => {
                                m.get_one::<String>("bam-file-cache-directory").unwrap()
                            }
                        },
                        ref_name
                    );
                    debug!("Cache: {}", &cache);

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

    if m.get_flag("longread-bam-files") {
        let bam_paths = m
            .get_many::<String>("longread-bam-files")
            .unwrap()
            .map(|b| b.to_string())
            .collect::<Vec<String>>();
        if !bams_are_split {
            bam_readers.extend(bam_paths);
        } else {
            for path in bam_paths {
                let mut path_split = path.rsplitn(2, '/');
                let path_suffix = path_split.next().unwrap();
                let path_prefix = path_split.next().unwrap();

                bam_readers.push(format!(
                    "{}/{}/{}",
                    path_prefix, &genomes_and_contigs.genomes[ref_idx], path_suffix
                ))
            }
        }
    } else if m.get_flag("longreads") {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/long/{}/*.bam",
                    match m.get_flag("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache.as_ref().unwrap()
                        }
                        true => {
                            m.get_one::<String>("bam-file-cache-directory").unwrap()
                        }
                    },
                    if bams_are_split {
                        &genomes_and_contigs.genomes[ref_idx]
                    } else {
                        ""
                    }
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
                        match m.get_flag("bam-file-cache-directory") {
                            false => {
                                tmp_bam_file_cache.as_ref().unwrap()
                            }
                            true => {
                                m.get_one::<String>("bam-file-cache-directory").unwrap()
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
        panic!(
            "Original sample count {} does not match new sample count {}, \
                please clear bam cache directory or ask for support at: github.com/rhysnewell/Lorikeet",
            (short_sample_count + long_sample_count),
            bam_readers.len(),
        )
    }

    return bam_readers;
}
