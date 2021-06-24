use std;
use std::collections::HashMap;

use bird_tool_utils::command;
use coverm::bam_generator::*;
use estimation::bams::{index_bams::*, process_bam::*};
use estimation::codon_structs::*;
use estimation::variant_matrix::*;
use external_command_checker;
use utils::utils::{*, Elem};
use utils::reference_reader_utils::ReferenceReaderUtils;

use crate::*;
use bio::io::gff;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use glob::glob;
use indicatif::MultiProgress;
use rayon::prelude::*;
use scoped_threadpool::Pool;
use std::path::Path;
use std::process::Stdio;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::NamedTempFile;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;

#[allow(unused)]
pub fn pileup_variants<
    'a,
    R: NamedBamReader,
    S: NamedBamReaderGenerator<R>,
    T: NamedBamReader,
    U: NamedBamReaderGenerator<T>,
>(
    m: &clap::ArgMatches,
    bam_readers: Vec<S>,
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
    threads: usize,
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
    let alpha: f64 = m.value_of("fdr-threshold").unwrap().parse().unwrap();
    let parallel_genomes = m.value_of("parallel-genomes").unwrap().parse().unwrap();
    let mut pool = Pool::new(parallel_genomes);
    let n_threads = std::cmp::max(threads / parallel_genomes as usize, 2);

    let references = ReferenceReaderUtils::parse_references(&m);
    let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();

    ReferenceReaderUtils::retrieve_reference(&Some(
        concatenated_genomes
            .as_ref()
            .unwrap()
            .path()
            .to_str()
            .unwrap()
            .to_string(),
    ));

    // All different counts of samples I need. Changes depends on when using concatenated genomes or not
    let mut short_sample_count = bam_readers.len();
    let mut long_sample_count = 0;
    let mut reference_count = references.len();

    let longreads = match longreads {
        Some(vec) => {
            long_sample_count += vec.len();
            vec
        }
        None => vec![],
    };



    // Finish each BAM source
    if m.is_present("longreads") || m.is_present("longread-bam-files") {
        info!("Processing long reads...");
        finish_bams(longreads, threads);
    }

    info!("Processing short reads...");
    finish_bams(bam_readers, threads);


    let mut reference_map = HashMap::new();


    // Set up multi progress bars
    let multi = Arc::new(MultiProgress::new());

    let multi_inner = Arc::clone(&multi);
    let mut progress_bars = setup_progress_bars(
        &references,
        &mut reference_map,
        &genomes_and_contigs,
        short_sample_count,
        long_sample_count,
    );

    let tree: Arc<Mutex<Vec<&Elem>>> =
        Arc::new(Mutex::new(Vec::with_capacity(progress_bars.len())));
    {
        let mut tree = tree.lock().unwrap();
        for pb in progress_bars.iter() {
            tree.push(pb)
        }
    }

    debug!(
        "{} Longread BAM files, {} Shortread BAM files {} Total BAMs over {} genome(s)",
        long_sample_count,
        short_sample_count,
        (short_sample_count + long_sample_count),
        reference_count
    );

    pool.scoped(|scope| {

        begin_tick(0, &progress_bars, &multi_inner, "");
        begin_tick(1, &progress_bars, &multi_inner, "");

        for (ref_idx, reference_stem) in reference_map.clone().into_iter() {

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
                        pb.progress_bar.reset_eta();
                        let pos = pb.progress_bar.position();
                        let len = pb.progress_bar.length();
                        if pos >= len {
                            pb.progress_bar
                                .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                        }
                    }
                    {
                        let pb = &tree.lock().unwrap()[0];
                        pb.progress_bar.inc(
                            ((short_sample_count + long_sample_count)
                                as u64)
                                * 2
                                + 1,
                        );
                        pb.progress_bar.reset_eta();
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

            scope.execute(move || {
                let reference = &genomes_and_contigs.genomes[ref_idx];
                begin_tick(ref_idx + 2, &progress_bars, &multi_inner, "Preparing variants");
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
                        per_reference_samples =
                            short_sample_count + long_sample_count;
                        per_reference_short_samples = short_sample_count;
                        debug!(
                            "Per reference samples concatenated {}",
                            &per_reference_samples
                        );

                        VariantMatrix::new_matrix(per_reference_samples)
                    }
                    None => {
                        per_reference_samples =
                            (short_sample_count + long_sample_count)
                                / references.len();
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
                            let gff_dir = TempDir::new("lorikeet-prodigal")
                                .expect("unable to create prodigal directory");
                            for reference in references.iter() {
                                let cmd_string = format!(
                                    "set -e -o pipefail; \
                                    prodigal -o {}/{}.gff -i {} -f gff {}",
                                    // prodigal
                                    gff_dir
                                        .path()
                                        .to_str()
                                        .expect("Failed to convert tempfile path to str"),
                                    Path::new(&reference).file_stem().unwrap().to_str().unwrap(),
                                    &reference,
                                    m.value_of("prodigal-params").unwrap(),
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
                                let ref_stem = Path::new(reference)
                                    .file_stem()
                                    .expect("problem determining file stem")
                                    .to_str()
                                    .unwrap()
                                    .to_string();

                                let ref_idx = genomes_and_contigs.genome_index(&ref_stem).unwrap();
                                gff_reader.records().into_iter().for_each(|record| {
                                    match record {
                                        Ok(rec) => {
                                            let gff_ref =
                                                gff_map.entry(ref_idx).or_insert(HashMap::new());
                                            let contig_genes = gff_ref
                                                .entry(format!("{}~{}", &ref_stem, rec.seqname()))
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
                    _ => {}
                }

                // let mut sample_groups = HashMap::new();

                debug!(
                    "Running SNP calling on {} samples",
                    indexed_bam_readers.len()
                );

                let mut hc_engine = HaplotypeCallerEngine::new(
                    m,
                    ref_idx,
                    indexed_bam_readers.clone(),
                    false,
                    m.value_of("ploidy").unwrap().parse().unwrap()
                );

                hc_engine.collect_activity_profile(
                    &indexed_bam_readers,
                    short_sample_count,
                    long_sample_count,
                    n_threads,
                    ref_idx,
                    per_reference_samples,
                    m,
                    genomes_and_contigs,
                    &concatenated_genomes,
                    flag_filters,
                    tree
                );

                // if genotype_likelihoods.len() == 1 {
                //     // Faster implementation for single sample analysis
                // } else {
                //
                // }

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
                                    ReadType::Short,
                                    &genomes_and_contigs,
                                    &reference_map,
                                    &concatenated_genomes,
                                )
                            } else if sample_idx >= short_sample_count
                                && sample_idx < (short_sample_count + long_sample_count)
                            {
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
                                    ReadType::Long,
                                    &genomes_and_contigs,
                                    &reference_map,
                                    &concatenated_genomes,
                                )
                            }

                            if mode == "evolve" || mode == "full" {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar.set_message(&format!(
                                    "{}: Calculating dN/dS values...",
                                    &reference,
                                ));
                                variant_matrix.add_gene_info(
                                    &mut gff_map,
                                    &genomes_and_contigs,
                                    &reference_map,
                                    &codon_table,
                                    &concatenated_genomes,
                                    sample_idx,
                                );
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
                    let anchor_similarity: f64 =
                        m.value_of("cluster-distance").unwrap().parse().unwrap();
                    let minimum_reads_in_link: usize = m
                        .value_of("minimum-reads-in-link")
                        .unwrap()
                        .parse()
                        .unwrap();
                    let n_components: usize = m.value_of("n-components").unwrap().parse().unwrap();

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
                        n_components,
                        n_threads,
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
