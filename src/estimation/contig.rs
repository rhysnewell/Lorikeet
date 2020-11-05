use std;
use std::collections::HashMap;

use bird_tool_utils::command;
use coverm::bam_generator::*;
use estimation::bams::{index_bams::*, process_bam::*};
use estimation::codon_structs::*;
use estimation::variant_matrix::*;
use estimation::vcfs::process_vcf::*;
use external_command_checker;
use utils::*;

use crate::*;
use bio::io::gff;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use glob::glob;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use scoped_threadpool::Pool;
use std::path::Path;
use std::process::Stdio;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::NamedTempFile;

#[derive(Clone, Debug)]
struct Elem {
    key: String,
    index: usize,
    progress_bar: ProgressBar,
}

#[allow(unused)]
pub fn pileup_variants<
    'a,
    R: NamedBamReader,
    G: NamedBamReaderGenerator<R>,
    S: NamedBamReader,
    U: NamedBamReaderGenerator<S>,
>(
    m: &clap::ArgMatches,
    bam_readers: Vec<G>,
    longreads: Option<Vec<U>>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    flag_filters: FlagFilter,
    mapq_threshold: u8,
    min_var_depth: usize,
    min: f32,
    max: f32,
    contig_end_exclusion: u64,
    output_prefix: &str,
    n_threads: usize,
    method: &str,
    coverage_fold: f32,
    include_indels: bool,
    include_soft_clipping: bool,
    is_long_read: bool,
    genomes_and_contigs: GenomesAndContigs,
    tmp_bam_file_cache: Option<TempDir>,
    concatenated_genomes: Option<NamedTempFile>,
) {
    // TODO: We need to split up analyses per reference to help contain memory issues
    let references = parse_references(&m);
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

    retrieve_reference(&Some(
        concatenated_genomes
            .as_ref()
            .unwrap()
            .path()
            .to_str()
            .unwrap()
            .to_string(),
    ));
    // let tmp_dir = TempDir::new("lorikeet_references").expect("Unable to create temporary directory");

    // All different counts of samples I need. Changes depends on when using concatenated genomes or not
    let mut short_sample_count = bam_readers.len();
    let mut long_sample_count = 0;
    let mut reference_count = references.len();

    let mut ani = 0.;

    let longreads = match longreads {
        Some(vec) => {
            long_sample_count += vec.len();
            vec
        }
        None => vec![],
    };

    let alpha: f64 = m.value_of("fdr-threshold").unwrap().parse().unwrap();

    // Finish each BAM source
    if m.is_present("longreads") || m.is_present("longread-bam-files") {
        info!("Processing long reads...");
        finish_bams(longreads, n_threads);
    }
    // if !m.is_present("bam-files") {
    info!("Processing short reads...");
    finish_bams(bam_readers, n_threads);
    // }

    // Put reference index in the variant map and initialize matrix
    let mut progress_bars = vec![
        Elem {
            key: "Genomes complete".to_string(),
            index: 1,
            progress_bar: ProgressBar::new(references.len() as u64),
        };
        references.len() + 2
    ];

    let mut reference_map = HashMap::new();
    for reference in references.iter() {
        debug!(
            "Genomes {:?} contigs {:?}",
            &genomes_and_contigs.genomes, &genomes_and_contigs.contig_to_genome,
        );

        let ref_idx = genomes_and_contigs
            .genome_index(
                &Path::new(reference)
                    .file_stem()
                    .expect("problem determining file stem")
                    .to_str()
                    .unwrap()
                    .to_string(),
            )
            .unwrap();

        progress_bars[ref_idx + 2] = Elem {
            key: genomes_and_contigs.genomes[ref_idx].clone(),
            index: ref_idx,
            progress_bar: ProgressBar::new((short_sample_count + long_sample_count) as u64),
        };
        debug!("Reference {}", reference,);
        reference_map
            .entry(ref_idx)
            .or_insert(reference.to_string());
    }

    progress_bars[0] = Elem {
        key: "Operations remaining".to_string(),
        index: 0,
        progress_bar: ProgressBar::new(
            ((references.len() * (short_sample_count + long_sample_count)) * 2 + references.len())
                as u64,
        ),
    };

    info!(
        "{} Longread BAM files and {} Shortread BAM files {} Total BAMs over {} genome(s)",
        long_sample_count,
        short_sample_count,
        (short_sample_count + long_sample_count),
        reference_count
    );

    let parallel_genomes = m.value_of("parallel-genomes").unwrap().parse().unwrap();
    let mut pool = Pool::new(parallel_genomes);
    let n_threads = std::cmp::max(n_threads / parallel_genomes as usize, 2);
    // Set up multi progress bars
    let multi = Arc::new(MultiProgress::new());
    let sty_eta = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");

    let sty_aux = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");
    progress_bars
        .par_iter()
        .for_each(|pb| pb.progress_bar.set_style(sty_aux.clone()));
    progress_bars[0].progress_bar.set_style(sty_eta.clone());

    // let pb_main = multi.add(ProgressBar::new(reference_map.keys().len() as u64));
    // pb_main.set_style(sty_eta.clone());

    let tree: Arc<Mutex<Vec<&Elem>>> =
        Arc::new(Mutex::new(Vec::with_capacity(progress_bars.len())));
    {
        let mut tree = tree.lock().unwrap();
        for pb in progress_bars.iter() {
            tree.push(pb)
        }
    }
    // let tree2 = Arc::clone(&tree);

    let multi_inner = Arc::clone(&multi);

    pool.scoped(|scope| {
        {
            // Total steps eta progress bar
            let elem = &progress_bars[0];
            let pb = multi_inner.insert(0, elem.progress_bar.clone());

            pb.enable_steady_tick(500);

            pb.set_message(&format!("{}...", &elem.key,));
        }
        {
            // completed genomes progress bar
            let elem = &progress_bars[1];
            let pb = multi_inner.insert(1, elem.progress_bar.clone());

            pb.enable_steady_tick(500);

            pb.set_message(&format!("{}...", &elem.key,));
        }
        for (ref_idx, reference_stem) in reference_map.clone().into_iter() {
            // let _ = std::thread::spawn(move || {
            //     multi.join().unwrap();
            // });

            // let ref_idx = Arc::new(ref_idx);

            // let ref_idx = Arc::new(ref_idx);
            let multi_inner = &multi_inner;
            let tree = &tree;
            let progress_bars = &progress_bars;
            let flag_filters = &flag_filters;
            let reference_map = &reference_map;
            let references = references.clone();
            let tmp_bam_file_cache = match tmp_bam_file_cache.as_ref() {
                Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
                None => None,
            };
            let concatenated_genomes = match concatenated_genomes.as_ref() {
                Some(file) => Some(file.path().to_str().unwrap().to_string()),
                None => None,
            };
            let mut coverage_estimators = coverage_estimators.clone();
            let genomes_and_contigs = &genomes_and_contigs;

            let output_prefix = format!(
                "{}/{}",
                &output_prefix,
                Path::new(&reference_stem)
                    .file_stem()
                    .unwrap()
                    .to_str()
                    .unwrap(),
            );

            if Path::new(&output_prefix).exists() && !m.is_present("force") {
                let cache = glob(&format!("{}/*.vcf", &output_prefix))
                    .expect("failed to interpret glob")
                    .map(|p| {
                        p.expect("Failed to read cached bam path")
                            .to_str()
                            .unwrap()
                            .to_string()
                    })
                    .collect::<Vec<String>>();
                if cache.len() > 0 {
                    {
                        let elem = &progress_bars[ref_idx + 2];
                        let pb = multi_inner.insert(ref_idx + 2, elem.progress_bar.clone());
                    }
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];

                        pb.progress_bar.set_message(&format!(
                            "{}: Output already present. Run with --force to overwrite",
                            &genomes_and_contigs.genomes[ref_idx]
                        ));
                        pb.progress_bar.finish_and_clear();
                    }
                    {
                        let pb = &tree.lock().unwrap()[1];
                        pb.progress_bar.inc(1);
                        let pos = pb.progress_bar.position();
                        let len = pb.progress_bar.length();
                        if pos >= len {
                            pb.progress_bar
                                .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                        }
                    }
                    {
                        let pb = &tree.lock().unwrap()[0];
                        pb.progress_bar
                            .inc(((short_sample_count + long_sample_count) as u64) * 2 + 1);
                        let pos = pb.progress_bar.position();
                        let len = pb.progress_bar.length();
                        if pos >= len {
                            pb.progress_bar
                                .finish_with_message(&format!("All steps completed {}", "✔",));
                        }
                    }
                    continue;
                }
            }
            // pb_main.tick();
            //
            // pb_main.inc(1);
            // pb_main.set_message(&format!(
            //     "Staging reference: {}",
            //     &genomes_and_contigs.genomes[ref_idx],
            // ));

            scope.execute(move || {
                let reference = &genomes_and_contigs.genomes[ref_idx];

                {
                    let elem = &progress_bars[ref_idx + 2];
                    let pb = multi_inner.insert(ref_idx + 2, elem.progress_bar.clone());

                    pb.enable_steady_tick(500);

                    pb.set_message(&format!("{}: Preparing variants...", &elem.key,));
                    // multi.join().unwrap();

                    // tree.lock().unwrap().insert(elem.index, &elem);
                }
                let mut codon_table = CodonTable::setup();

                // Read BAMs back in as indexed
                let mut indexed_bam_readers = recover_bams(
                    m,
                    &concatenated_genomes,
                    short_sample_count,
                    long_sample_count,
                    &genomes_and_contigs,
                    n_threads as u32,
                    &tmp_bam_file_cache,
                );
                let mut per_reference_samples = 0;
                let mut per_reference_short_samples = 0;
                let mut variant_matrix = match concatenated_genomes {
                    Some(ref concat) => {
                        per_reference_samples = short_sample_count + long_sample_count;
                        per_reference_short_samples = short_sample_count;
                        debug!(
                            "Per reference samples concatenated {}",
                            &per_reference_samples
                        );

                        VariantMatrix::new_matrix(per_reference_samples)
                    }
                    None => {
                        per_reference_samples =
                            (short_sample_count + long_sample_count) / references.len();
                        per_reference_short_samples = short_sample_count / references.len();
                        debug!(
                            "Per reference samples not concatenated {}",
                            &per_reference_samples
                        );
                        VariantMatrix::new_matrix(per_reference_samples)
                    }
                };

                // Gff map lists coding regions
                let mut gff_map = HashMap::new();
                match mode {
                    "evolve" => {
                        codon_table.get_codon_table(11);
                        // ani = 0.;

                        let mut gff_reader;
                        if m.is_present("gff") || !m.is_present("gff") {
                            external_command_checker::check_for_prodigal();

                            // create new fifo and give read, write and execute rights to the owner.
                            let gff_dir = TempDir::new("lorikeet-prokka")
                                .expect("unable to create prokka directory");
                            for reference in references.iter() {
                                let cmd_string = format!(
                                    "set -e -o pipefail; \
                                    prodigal -o {}/{}.gff -i {} -f gff -p meta",
                                    // prodigal
                                    gff_dir
                                        .path()
                                        .to_str()
                                        .expect("Failed to convert tempfile path to str"),
                                    Path::new(&reference).file_stem().unwrap().to_str().unwrap(),
                                    &reference
                                );
                                debug!("Queuing cmd_string: {}", cmd_string);
                                command::finish_command_safely(
                                    std::process::Command::new("bash")
                                        .arg("-c")
                                        .arg(&cmd_string)
                                        .stdout(Stdio::piped())
                                        .stderr(Stdio::piped())
                                        .spawn()
                                        .expect("Unable to execute bash"),
                                    "prodigal",
                                );

                                // Read in newly created gff
                                gff_reader = gff::Reader::from_file(
                                    format!(
                                        "{}/{}.gff",
                                        gff_dir
                                            .path()
                                            .to_str()
                                            .expect("Failed to convert tempfile path to str"),
                                        Path::new(&reference)
                                            .file_stem()
                                            .unwrap()
                                            .to_str()
                                            .unwrap()
                                    ),
                                    bio::io::gff::GffType::GFF3,
                                )
                                .expect("Failed to read GFF file");

                                // Map to reference id
                                gff_reader.records().into_iter().for_each(|record| {
                                    match record {
                                        Ok(rec) => {
                                            let gff_ref = gff_map
                                                .entry(
                                                    genomes_and_contigs
                                                        .genome_index(
                                                            &Path::new(reference)
                                                                .file_stem()
                                                                .expect(
                                                                    "problem determining file stem",
                                                                )
                                                                .to_str()
                                                                .unwrap()
                                                                .to_string(),
                                                        )
                                                        .unwrap(),
                                                )
                                                .or_insert(HashMap::new());
                                            let contig_genes = gff_ref
                                                .entry(rec.seqname().to_owned())
                                                .or_insert(Vec::new());
                                            contig_genes.push(rec);
                                        }
                                        _ => {}
                                    };
                                });
                            }

                            gff_dir.close().expect("Failed to close temp directory");
                        }
                    }
                    "genotype" | "summarize" => {
                        // if m.is_present("strain-ani") {
                        //     ani = parse_percentage(m, "strain-ani");
                        // }
                    }
                    _ => {
                        //            min_cluster_size = m.value_of("min-cluster-size").unwrap().parse().unwrap();
                        //            epsilon = m.value_of("epsilon").unwrap().parse().unwrap();
                    }
                }

                // let mut sample_groups = HashMap::new();

                debug!(
                    "Running SNP calling on {} shortread samples",
                    indexed_bam_readers.len()
                );

                // let mut prev_ref_idx = -1;
                // let mut per_ref_sample_idx = 0;

                // let threads = std::cmp::max(n_threads / indexed_bam_readers.len(), 1);

                indexed_bam_readers.into_iter().enumerate().for_each(
                    |(sample_idx, bam_generator)| {
                        // Get the appropriate sample index based on how many references we are using
                        let mut bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                            vec![&bam_generator],
                            n_threads as u32,
                        )
                        .into_iter()
                        .next()
                        .unwrap();
                        if sample_idx < short_sample_count {
                            process_vcf(
                                bam_generator,
                                n_threads,
                                ref_idx,
                                sample_idx,
                                per_reference_samples,
                                &mut variant_matrix,
                                false,
                                m,
                                // &mut sample_groups,
                                &genomes_and_contigs,
                                &reference_map,
                                per_reference_short_samples,
                                &concatenated_genomes,
                                &flag_filters,
                            );
                        } else if m.is_present("include-longread-svs")
                            && (m.is_present("longreads") | m.is_present("longread-bam-files"))
                        {
                            debug!("Running structural variant detection...");
                            // Get the appropriate sample index based on how many references we are using by tracking
                            // changes in references
                            process_vcf(
                                bam_generator,
                                n_threads,
                                ref_idx,
                                sample_idx,
                                per_reference_samples,
                                &mut variant_matrix,
                                true,
                                m,
                                // &mut sample_groups,
                                &genomes_and_contigs,
                                &reference_map,
                                per_reference_short_samples,
                                &concatenated_genomes,
                                &flag_filters,
                            );
                        } else if m.is_present("longreads") | m.is_present("longread-bam-files") {
                            // We need update the variant matrix anyway
                            let bam_generated = bam_generator.start();
                            let header = bam_generated.header().clone(); // bam header
                            let target_names = header.target_names(); // contig names
                            let reference = &genomes_and_contigs.genomes[ref_idx];

                            let mut stoit_name = bam_generated.name().to_string().replace("/", ".");
                            debug!("Stoit_name {:?}", &stoit_name);
                            variant_matrix.add_sample_name(stoit_name.to_string(), sample_idx);

                            // let group = sample_groups.entry("long").or_insert(HashSet::new());
                            // group.insert(stoit_name.clone());
                            // for each genomic position, only has hashmap when variants are present. Includes read ids
                            for (tid, target) in target_names.iter().enumerate() {
                                let target_name = String::from_utf8(target.to_vec()).unwrap();
                                if target_name.contains(reference) {
                                    let target_len = header.target_len(tid as u32).unwrap();
                                    variant_matrix.add_info(
                                        ref_idx,
                                        tid,
                                        target_name.as_bytes().to_vec(),
                                        target_len,
                                    );
                                }
                            }
                        }
                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];

                            pb.progress_bar.set_message(&format!(
                                "{}: Variant calling on sample: {}",
                                pb.key,
                                variant_matrix.get_sample_name(sample_idx),
                            ));
                            pb.progress_bar.inc(1);
                        }
                        {
                            let pb = &tree.lock().unwrap()[0];
                            pb.progress_bar.inc(1);
                        }
                    },
                );
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];
                    pb.progress_bar
                        .set_message(&format!("{}: Initial variant calling complete...", pb.key));
                }
                // // Read BAMs back in as indexed
                let mut indexed_bam_readers = recover_bams(
                    m,
                    &concatenated_genomes,
                    short_sample_count,
                    long_sample_count,
                    &genomes_and_contigs,
                    n_threads as u32,
                    &tmp_bam_file_cache,
                );

                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];
                    pb.progress_bar.reset();
                    pb.progress_bar.enable_steady_tick(1000);
                    pb.progress_bar
                        .set_message(&format!("{}: Performing guided variant calling...", pb.key));
                }
                // let mut variant_matrix = Mutex::new(variant_matrix);
                if variant_matrix.get_variant_count(ref_idx) > 0 {
                    // if there are variants, perform guided variant calling
                    std::fs::create_dir_all(&output_prefix).unwrap();

                    indexed_bam_readers.into_iter().enumerate().for_each(
                        |(sample_idx, bam_generator)| {
                            let mut bam_generator =
                                generate_indexed_named_bam_readers_from_bam_files(
                                    vec![&bam_generator],
                                    n_threads as u32,
                                )
                                .into_iter()
                                .next()
                                .unwrap();
                            if sample_idx < short_sample_count {
                                process_bam(
                                    bam_generator,
                                    sample_idx,
                                    per_reference_samples,
                                    &mut coverage_estimators,
                                    &mut variant_matrix,
                                    n_threads,
                                    m,
                                    &output_prefix,
                                    coverage_fold,
                                    &codon_table,
                                    min_var_depth,
                                    contig_end_exclusion,
                                    min,
                                    max,
                                    ref_idx,
                                    mode,
                                    include_soft_clipping,
                                    include_indels,
                                    &flag_filters,
                                    mapq_threshold,
                                    method,
                                    false,
                                    &genomes_and_contigs,
                                    &reference_map,
                                    &concatenated_genomes,
                                )
                            } else {
                                process_bam(
                                    bam_generator,
                                    sample_idx,
                                    per_reference_samples,
                                    &mut coverage_estimators,
                                    &mut variant_matrix,
                                    n_threads,
                                    m,
                                    &output_prefix,
                                    coverage_fold,
                                    &codon_table,
                                    min_var_depth,
                                    contig_end_exclusion,
                                    min,
                                    max,
                                    ref_idx,
                                    mode,
                                    include_soft_clipping,
                                    include_indels,
                                    &flag_filters,
                                    mapq_threshold,
                                    method,
                                    true,
                                    &genomes_and_contigs,
                                    &reference_map,
                                    &concatenated_genomes,
                                )
                            }

                            {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar.set_message(&format!(
                                    "{}: Guided variant calling on sample: {}",
                                    pb.key,
                                    variant_matrix.get_sample_name(sample_idx),
                                ));
                                pb.progress_bar.inc(1);
                            }
                            {
                                let pb = &tree.lock().unwrap()[0];
                                pb.progress_bar.inc(1);
                            }
                        },
                    );
                } else {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(indexed_bam_readers.len() as u64);
                    pb.progress_bar.reset_eta();
                }
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];
                    pb.progress_bar
                        .set_message(&format!("{}: Guided variant calling complete...", pb.key));
                }

                // Collects info about variants across samples to check whether they are genuine or not
                // using FDR

                // TODO: Make sure that this is fixed. It seems to work appropriately now
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];
                    pb.progress_bar
                        .set_message(&format!("{}: Setting FDR threshold...", pb.key));
                }
                variant_matrix
                    .remove_false_discoveries(alpha, &genomes_and_contigs.genomes[ref_idx]);

                if mode == "genotype" {
                    let e_min: f64 = m.value_of("e-min").unwrap().parse().unwrap();
                    let e_max: f64 = m.value_of("e-max").unwrap().parse().unwrap();
                    let pts_min: f64 = m.value_of("pts-min").unwrap().parse().unwrap();
                    let pts_max: f64 = m.value_of("pts-max").unwrap().parse().unwrap();
                    let phi: f64 = m.value_of("phi").unwrap().parse().unwrap();
                    let anchor_size: usize = m.value_of("n-neighbors").unwrap().parse().unwrap();
                    let anchor_similarity: f64 = m
                        .value_of("maximum-seed-similarity")
                        .unwrap()
                        .parse()
                        .unwrap();
                    let minimum_reads_in_link: usize = m
                        .value_of("minimum-reads-in-link")
                        .unwrap()
                        .parse()
                        .unwrap();

                    // Calculate the geometric mean values and CLR for each variant, reference specific
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating variant distances...",
                            &reference,
                        ));
                    }
                    variant_matrix.generate_distances();

                    // Generate initial read linked clusters
                    // Cluster each variant using phi-D and fuzzy DBSCAN, reference specific
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Running UMAP and HDBSCAN...", &reference,));
                    }
                    variant_matrix.run_fuzzy_scan(
                        e_min,
                        e_max,
                        pts_min,
                        pts_max,
                        phi,
                        anchor_size,
                        anchor_similarity,
                        minimum_reads_in_link,
                        &reference_map,
                        &multi_inner,
                        &output_prefix,
                        anchor_size,
                    );

                    // Write genotypes to disk, reference specific
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating genotypes...", &reference,));
                    }
                    variant_matrix.generate_genotypes(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                        &concatenated_genomes,
                    );

                    // Get strain abundances
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Calculating genotype abundances...",
                            &reference,
                        ));
                    }
                    variant_matrix.calculate_strain_abundances(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                    );

                    // Get sample distances
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating adjacency matrix...",
                            &reference,
                        ));
                    }
                    variant_matrix.calculate_sample_distances(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                    );

                    // If flagged, then create plots using CMplot
                    if m.is_present("plot") {
                        let window_size = m.value_of("window-size").unwrap().parse().unwrap();
                        variant_matrix.print_variant_stats(
                            window_size,
                            &output_prefix,
                            &genomes_and_contigs,
                        );
                    };

                    // Write variants in VCF format, reference specific
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating VCF file...", &reference,));
                    }
                    variant_matrix.write_vcf(&output_prefix, &genomes_and_contigs);
                } else if mode == "summarize" {
                    variant_matrix.generate_distances();

                    let window_size = m.value_of("window-size").unwrap().parse().unwrap();

                    // Get sample distances
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating adjacency matrix...",
                            &reference,
                        ));
                    }
                    variant_matrix.calculate_sample_distances(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                    );

                    variant_matrix.print_variant_stats(
                        window_size,
                        &output_prefix,
                        &genomes_and_contigs,
                    );

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating VCF file...", &reference,));
                    }
                    variant_matrix.write_vcf(&output_prefix, &genomes_and_contigs);
                } else if mode == "evolve" {
                    variant_matrix.generate_distances();

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Calculating dN/dS values...", &reference,));
                    }
                    variant_matrix.calc_gene_mutation(
                        &mut gff_map,
                        &genomes_and_contigs,
                        &reference_map,
                        &codon_table,
                        &output_prefix,
                        &concatenated_genomes,
                    );

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating VCF file...", &reference,));
                    }
                    variant_matrix.write_vcf(&output_prefix, &genomes_and_contigs);
                } else if mode == "polish" {
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating consensus genomes...",
                            &reference,
                        ));
                    }
                    variant_matrix.generate_distances();
                    variant_matrix.polish_genomes(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                        &concatenated_genomes,
                    );

                    // Get sample distances
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating adjacency matrix...",
                            &reference,
                        ));
                    }
                    variant_matrix.calculate_sample_distances(
                        &output_prefix,
                        &reference_map,
                        &genomes_and_contigs,
                    );
                    // If flagged, then create plots using CMplot
                    if m.is_present("plot") {
                        let window_size = m.value_of("window-size").unwrap().parse().unwrap();
                        variant_matrix.print_variant_stats(
                            window_size,
                            &output_prefix,
                            &genomes_and_contigs,
                        );
                    };

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating VCF file...", &reference,));
                    }
                    variant_matrix.write_vcf(&output_prefix, &genomes_and_contigs);
                };
                {
                    let pb = &tree.lock().unwrap()[ref_idx + 2];
                    pb.progress_bar
                        .set_message(&format!("{}: All steps completed {}", &reference, "✔",));
                    pb.progress_bar.finish_and_clear();
                }
                {
                    let pb = &tree.lock().unwrap()[1];
                    pb.progress_bar.inc(1);
                    let pos = pb.progress_bar.position();
                    let len = pb.progress_bar.length();
                    if pos >= len {
                        pb.progress_bar
                            .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                    }
                }
                {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(1);
                    let pos = pb.progress_bar.position();
                    let len = pb.progress_bar.length();
                    if pos >= len {
                        pb.progress_bar
                            .finish_with_message(&format!("All steps completed {}", "✔",));
                    }
                }
            });
        }

        // pb_main.finish_with_message("All genomes staged...");
        multi.join().unwrap();
    });
    // reference_map.iter().for_each(|(ref_idx, reference_stem)| {
    //     pb2.reset();
    //     pb3.reset();
    //     pb3.set_message("Waiting for variant matrix...");
    //     pb4.reset();
    //
    // });

    info!("Analysis finished!");
    // multi.join_and_clear().unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use coverm::bam_generator::*;
    use coverm::genome_exclusion::*;
    use coverm::mapping_parameters::*;
    use coverm::shard_bam_reader::*;
    use std::fs::File;
    use std::str;
    //
    //    fn test_with_stream<R: NamedBamReader + Send,
    //        G: NamedBamReaderGenerator<R> + Send>(
    //        expected: &str,
    //        bam_readers: Vec<G>,
    //        mut reference: bio::io::fasta::IndexedReader<File>,
    //        proper_pairs_only: bool,
    //        n_threads: usize,
    //        coverage_fold: f32,
    //        min_var_depth: usize,
    //        min: f32,
    //        max: f32,
    //        mode: &str,
    //        include_indels: bool,
    //        include_soft_clipping: bool) {
    ////        let mut stream = Cursor::new(Vec::new());
    //        {
    //            reads_mapped_vec = variant_variants(
    //                bam_readers,
    //                &mut coverage_taker,
    //                coverage_estimators,
    //                print_zero_coverage_contigs,
    //                flag_filters,
    //                false,
    //            );
    //        }
    ////        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
    //    }
}
