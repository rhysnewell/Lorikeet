use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam;
use std::fs::create_dir_all;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Ensures mapping is completed for provided bams. If multiple references are provided and
/// the user has asked to run genomes in parallel, then the bams are split per reference to
/// avoid file locking when reading bams in parallel.
pub fn finish_bams<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
    bams: Vec<G>,
    n_threads: usize,
    references: &GenomesAndContigs,
    parallel_genomes: bool,
    mapping: bool,
) {
    if !parallel_genomes {
        let mut record: bam::record::Record = bam::Record::new();

        // progress bar
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");
        let pb1 = ProgressBar::new(bams.len() as u64);
        pb1.set_style(sty);
        pb1.enable_steady_tick(500);
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
    } else {
        // we need to subset the bams
        let mut record: bam::record::Record = bam::Record::new();

        // progress bar
        let sty = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");
        let pb1 = Arc::new(Mutex::new(ProgressBar::new(bams.len() as u64)));
        {
            let mut pb1 = pb1.lock().unwrap();
            pb1.set_style(sty);
            pb1.enable_steady_tick(500);
            pb1.set_message("Splitting bam files per reference...");
        }
        let mut paths = Vec::with_capacity(bams.len());

        for bam_generator in bams {
            let mut bam = bam_generator.start();
            // bam.set_threads(std::cmp::max(n_threads / 2, 1));

            let path = bam.path().to_string();

            paths.push(path);
            if mapping {
                while bam.read(&mut record).is_some() {
                    continue;
                }
                bam.finish();
            }
        }
        paths.into_par_iter().for_each(|path| {
            let mut record: bam::record::Record = bam::Record::new();

            let mut bam = generate_named_bam_readers_from_bam_files(vec![&path])
                .into_iter()
                .next()
                .unwrap();
            let header = bam.header();
            let stoit_name = bam.name().to_string().replace('/', ".");
            let path_prefix = path.rsplitn(2, '/').nth(1).unwrap();
            let mut tmp_header = bam::header::Header::from_template(&header);
            let mut need_to_write = false;
            let mut writers = references
                .genomes
                .iter()
                .map(|ref_stub| {
                    let output_path = format!("{}/{}", &path_prefix, ref_stub);
                    create_dir_all(&output_path);
                    let bam_path = format!("{}/{}.bam", output_path, &stoit_name);
                    let mut out = None;
                    if !Path::new(&bam_path).exists() {
                        out = Some(
                            bam::Writer::from_path(&bam_path, &tmp_header, bam::Format::Bam)
                                .unwrap(),
                        );
                        need_to_write = true;
                    }
                    out
                })
                .collect::<Vec<Option<bam::Writer>>>();

            if need_to_write {
                while bam.read(&mut record).is_some() {
                    if !record.is_unmapped() {
                        let header = bam.header();
                        let contig_name = std::str::from_utf8(header.tid2name(record.tid() as u32))
                            .unwrap()
                            .to_string()
                            .split_once("~")
                            .unwrap()
                            .1
                            .to_string();
                        let ref_idx = references.genome_index_of_contig(&contig_name);
                        match ref_idx {
                            Some(ref_idx) => {
                                match &mut writers[ref_idx] {
                                    None => {
                                        // bam file already written
                                    }
                                    Some(writer) => writer.write(&record).unwrap(),
                                }
                            }
                            None => {} // shall pass
                        }
                    }
                }
                bam.finish();
            }
            {
                let mut pb1 = pb1.lock().unwrap();
                pb1.inc(1);
            }
        });
        {
            let mut pb1 = pb1.lock().unwrap();
            pb1.finish_with_message(format!("Reads and BAM files processed..."));
        }
    }
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
    if m.is_present("bam-files") {
        let bam_paths = m
            .values_of("bam-files")
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
    } else if m.is_present("read1")
        | m.is_present("single")
        | m.is_present("coupled")
        | m.is_present("interleaved")
    {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/short/{}/*.bam",
                    match m.is_present("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache.as_ref().unwrap()
                        }
                        true => {
                            m.value_of("bam-file-cache-directory").unwrap()
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

    if m.is_present("longread-bam-files") {
        let bam_paths = m
            .values_of("longread-bam-files")
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
    } else if m.is_present("longreads") {
        let mut all_bam_paths = vec![];

        match concatenated_genomes {
            Some(ref _tmp_file) => {
                let cache = format!(
                    "{}/long/{}/*.bam",
                    match m.is_present("bam-file-cache-directory") {
                        false => {
                            tmp_bam_file_cache.as_ref().unwrap()
                        }
                        true => {
                            m.value_of("bam-file-cache-directory").unwrap()
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
        panic!(
            "Original sample count {} does not match new sample count {}, \
                please clear bam cache directory or ask for support at: github.com/rhysnewell/Lorikeet",
            (short_sample_count + long_sample_count),
            bam_readers.len(),
        )
    }

    return bam_readers;
}
