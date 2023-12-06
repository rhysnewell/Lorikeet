use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};
// use rayon::prelude::*;
use rust_htslib::bam;
use std::collections::HashMap;
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
    references: &GenomesAndContigs,
    split_bams: bool,
    mapping: bool,
) -> Result<(), BirdToolError> {
    let mut record: bam::Record = bam::Record::new();

    // progress bar
    let sty = match ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
    {
        Ok(s) => s,
        Err(e) => return Err(BirdToolError::DebugError(e.to_string())),
    };
    let pb1 = ProgressBar::new(bams.len() as u64);
    pb1.set_style(sty);
    pb1.enable_steady_tick(Duration::from_millis(200));

    // bams.into_iter().par_bridge(|bam_generator| {

    // });

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

        if split_bams {
            split_bams_to_references(bam, references, &path, n_threads)?;
        } else if mapping {
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

/// Splits bams by reference if the user has requested split bams.
/// This is done to avoid file locking when reading bams in parallel.
fn split_bams_to_references<R: NamedBamReader>(
    mut bam_generator: R,
    references: &GenomesAndContigs,
    bam_path: &str,
    n_threads: usize
) -> Result<(), BirdToolError> {
    let mut bam_writer_map: HashMap<String, bam::Writer> = HashMap::with_capacity(references.genomes.len());
    let bam_header = bam_generator.header();
    let new_header = bam::Header::from_template(bam_header);
    debug!("Beginning bam splitting...");
    // ensure output file exists for each reference
    for ref_name in references.genomes.iter() {
        let mut path_split = bam_path.rsplitn(2, '/');
        let path_suffix = path_split.next().unwrap();
        let path_prefix = path_split.next().unwrap();

        let path = format!(
            "{}/{}/{}",
            path_prefix, ref_name, path_suffix
        );
        // ensure path exists
        std::fs::create_dir_all(Path::new(&path).parent().unwrap())
            .map_err(|_| BirdToolError::IOError(format!("Unable to create path {}", &path)))?;
        
        let writer = bam::Writer::from_path(&path, &new_header, bam::Format::Bam)
            .map_err(|_| BirdToolError::IOError(format!("Unable to write bam at {}", &path)))?;

        bam_writer_map.insert(ref_name.to_string(), writer);
    }

    debug!("Splitting bam to {} references", bam_writer_map.len());

    let mut record: bam::Record = bam::Record::new();

    while bam_generator.read(&mut record).is_some() {
        let record_tid = record.tid();
        if record_tid < 0 {
            continue;
        }
        let ref_name = std::str::from_utf8(bam_generator.header().tid2name(record_tid as u32)).expect("Cannot read reference name from bam file").split('~').next().unwrap();

        let writer = bam_writer_map.get_mut(ref_name).unwrap();
        writer.write(&record).unwrap();
    }

    bam_generator.finish();
    
    let mut paths_to_index = Vec::with_capacity(bam_writer_map.len());
    // build indices for writers in two loops
    // we use two loops as the destructor for the writer needs to run
    // otherwise there is no EOF marker in the bam file
    for (ref_name, _) in bam_writer_map.into_iter() {
        let mut path_split = bam_path.rsplitn(2, '/');
        let path_suffix = path_split.next().unwrap();
        let path_prefix = path_split.next().unwrap();

        let path = format!(
            "{}/{}/{}",
            path_prefix, ref_name, path_suffix
        );
        paths_to_index.push(path);
        
    }

    // build indices for writers
    for path in paths_to_index.into_iter() {
        bam::index::build(
            &path,
            Some(&format!("{}.bai", path)),
            bam::index::Type::Bai,
            n_threads as u32,
        )
        .unwrap_or_else(|_| panic!("Unable to index bam at {}", &path));
    }

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
    if m.contains_id("bam-files") {
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
    } else if m.contains_id("read1")
        | m.contains_id("single")
        | m.contains_id("coupled")
        | m.contains_id("interleaved")
    {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/short/{}/*.bam",
                    match m.contains_id("bam-file-cache-directory") {
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
                        match m.contains_id("bam-file-cache-directory") {
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

    if m.contains_id("longread-bam-files") {
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
    } else if m.contains_id("longreads") {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/long/{}/*.bam",
                    match m.contains_id("bam-file-cache-directory") {
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
                        match m.contains_id("bam-file-cache-directory") {
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
