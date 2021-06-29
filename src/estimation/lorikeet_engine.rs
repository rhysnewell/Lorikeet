use coverm::mosdepth_genome_coverage_estimators::CoverageEstimator;
use coverm::bam_generator::*;
use coverm::FlagFilter;
use utils::reference_reader_utils::{ReferenceReaderUtils, ReferenceReader};
use utils::utils::*;
use estimation::bams::index_bams::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use std::sync::{Arc, Mutex};
use std::collections::HashMap;
use scoped_threadpool::Pool;
use std::path::Path;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use assembly::assembly_region_walker::AssemblyRegionWalker;


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReadType {
    Short,
    Long,
    Assembly,
}

#[derive(Clone, Debug)]
pub struct Elem {
    pub key: String,
    pub index: usize,
    pub progress_bar: ProgressBar,
}

/**
The main lorikeet engine, takes input files and performs activity profiling and local reassembly
per reference and per contig
*/
pub struct LorikeetEngine<'a> {
    args: &'a clap::ArgMatches,
    short_read_bam_count: usize,
    long_read_bam_count: usize,
    coverage_estimators: Vec<CoverageEstimator>,
    flag_filters: FlagFilter,
    genomes_and_contigs: GenomesAndContigs,
    concatenated_genomes: Option<NamedTempFile>,
    tmp_bam_file_cache: Option<TempDir>,
    reference_map: HashMap<usize, String>,
    references: Vec<&str>,
    multi: Arc<MultiProgress>,
    multi_inner: Arc<MultiProgress>,
    progress_bars: Vec<Elem>,
    tree: Arc<Mutex<&Elem>>,
    reference_reader: ReferenceReader,
}

impl LorikeetEngine {
    pub fn start<
        R: NamedBamReader,
        S: NamedBamReaderGenerator<R>,
        T: NamedBamReader,
        U: NamedBamReaderGenerator<T>,
    >(
        m: &clap::ArgMatches,
        bam_readers: Vec<S>,
        longreads: Option<Vec<U>>,
        mode: &str,
        coverage_estimators: Vec<CoverageEstimator>,
        flag_filters: FlagFilter,
        genomes_and_contigs: GenomesAndContigs,
        tmp_bam_file_cache: Option<TempDir>,
        concatenated_genomes: Option<NamedTempFile>,
    ) {

        let references = ReferenceReaderUtils::parse_references(&m);
        let references = references.par_iter().map(|p| &**p).collect::<Vec<&str>>();

        ReferenceReaderUtils::retrieve_reference(&Some(
            concatenated_genomes
                .as_ref()
                .unwrap()
                .path()
                .to_str()
                .unwrap()
                .to_string(),
        ));

        let reference_reader = ReferenceReader::new(
            concatenated_genomes,
            genomes_and_contigs.clone(),
            genomes_and_contigs.contig_to_genome.len()
        );

        // All different counts of samples I need. Changes depends on when using concatenated genomes or not
        let short_read_bam_count = bam_readers.len();
        let mut long_read_bam_count = 0;
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
        let mut progress_bars = Self::setup_progress_bars(
            &references,
            &mut reference_map,
            &genomes_and_contigs,
            short_read_bam_count,
            long_read_bam_count,
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
            long_read_bam_count,
            short_read_bam_count,
            (short_read_bam_count + long_read_bam_count),
            reference_count
        );


        let mut lorikeet_engine = LorikeetEngine {
            args: m,
            short_read_bam_count,
            long_read_bam_count,
            coverage_estimators,
            flag_filters,
            genomes_and_contigs,
            concatenated_genomes,
            tmp_bam_file_cache,
            reference_map,
            references,
            multi,
            multi_inner,
            progress_bars,
            tree,
            reference_reader,
        };

        lorikeet_engine.apply_per_reference()
    }

    fn apply_per_reference(&self) {
        let alpha: f64 = self.args.value_of("fdr-threshold").unwrap().parse().unwrap();
        let parallel_genomes = self.args.value_of("parallel-genomes").unwrap().parse().unwrap();
        let mut pool = Pool::new(parallel_genomes);
        let n_threads = std::cmp::max(threads / parallel_genomes as usize, 2);
        let output_prefix = match self.args.is_present("output-directory") {
            true => {
                match std::fs::create_dir_all(
                    self.args.value_of("output-directory").unwrap().to_string(),
                ) {
                    Ok(_) => {}
                    Err(err) => panic!(format!("Unable to create output directory {:?}", err)),
                };
                self.args.value_of("output-directory").unwrap()
            }
            false => "./",
        };

        pool.scoped(|scope| {

            Self::begin_tick(0, &progress_bars, &multi_inner, "");
            Self::begin_tick(1, &progress_bars, &multi_inner, "");

            for (ref_idx, reference_stem) in self.reference_map.clone().into_iter() {

                let multi_inner = &self.multi_inner;
                let tree = &self.tree;
                let progress_bars = &self.progress_bars;
                let flag_filters = &self.flag_filters;
                let reference_map = &self.reference_map;
                let references = self.references.clone();
                let tmp_bam_file_cache = match self.tmp_bam_file_cache.as_ref() {
                    Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
                    None => None,
                };
                let concatenated_genomes = match self.concatenated_genomes.as_ref() {
                    Some(file) => Some(file.path().to_str().unwrap().to_string()),
                    None => None,
                };
                let mut coverage_estimators = self.coverage_estimators.clone();
                let mut reference_reader = self.reference_reader.clone();
                let genomes_and_contigs = &self.genomes_and_contigs;

                let output_prefix = format!(
                    "{}/{}",
                    &output_prefix,
                    Path::new(&reference_stem)
                        .file_stem()
                        .unwrap()
                        .to_str()
                        .unwrap(),
                );


                if Path::new(&output_prefix).exists() && !self.args.is_present("force") {
                    let cache = glob::glob(&format!("{}/*.vcf", &output_prefix))
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
                    Self::begin_tick(ref_idx + 2, &progress_bars, &multi_inner, "Preparing variants");

                    // Read BAMs back in as indexed
                    let mut indexed_bam_readers = recover_bams(
                        &self.args,
                        &concatenated_genomes,
                        self.short_read_bam_count,
                        self.long_read_bam_count,
                        &genomes_and_contigs,
                        n_threads as u32,
                        &tmp_bam_file_cache,
                    );
                    let mut per_reference_samples = 0;
                    let mut per_reference_short_samples = 0;

                    debug!(
                        "Running SNP calling on {} samples",
                        indexed_bam_readers.len()
                    );


                    let mut assembly_engine = AssemblyRegionWalker::start(
                        self.args,
                        ref_idx,
                        self.short_read_bam_count,
                        self.long_read_bam_count,
                        &indexed_bam_readers,
                        genomes_and_contigs,
                        &concatenated_genomes,
                        flag_filters,
                        tree,
                        reference_reader,
                        n_threads,
                    );

                    assembly_engine.traverse();
                    // let mut hc_engine = HaplotypeCallerEngine::new(
                    //     &self.args,
                    //     ref_idx,
                    //     indexed_bam_readers.clone(),
                    //     false,
                    //     self.args.value_of("ploidy").unwrap().parse().unwrap()
                    // );
                    //
                    // hc_engine.apply(
                    //     &indexed_bam_readers,
                    //     self.short_read_bam_count,
                    //     self.long_read_bam_count,
                    //     n_threads,
                    //     ref_idx,
                    //     per_reference_samples,
                    //     &self.args,
                    //     genomes_and_contigs,
                    //     &concatenated_genomes,
                    //     flag_filters,
                    //     tree
                    // );

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(&format!("{}: Initial variant calling complete...", pb.key));
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
    }

    pub fn setup_progress_bars(
        references: &Vec<&str>,
        reference_map: &mut HashMap<usize, String>,
        genomes_and_contigs: &GenomesAndContigs,
        short_sample_count: usize,
        long_sample_count: usize,
    ) -> Vec<Elem> {
        // Put reference index in the variant map and initialize matrix
        let mut progress_bars = vec![
            Elem {
                key: "Genomes complete".to_string(),
                index: 1,
                progress_bar: ProgressBar::new(references.len() as u64),
            };
            references.len() + 2
        ];

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
                progress_bar: ProgressBar::new(
                    (short_sample_count + long_sample_count) as u64,
                ),
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
                ((references.len() * (short_sample_count + long_sample_count))
                    * 2
                    + references.len()) as u64,
            ),
        };

        let sty_eta = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");

        let sty_aux = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");
        progress_bars
            .par_iter()
            .for_each(|pb| pb.progress_bar.set_style(sty_aux.clone()));
        progress_bars[0].progress_bar.set_style(sty_eta.clone());

        return progress_bars
    }

    pub fn begin_tick(
        index: usize,
        progress_bars: &Vec<Elem>,
        multi_inner: &Arc<MultiProgress>,
        message: &str,
    ) {
        let elem = &progress_bars[index];
        let pb = multi_inner.insert(index, elem.progress_bar.clone());

        pb.enable_steady_tick(500);

        pb.set_message(&format!("{}: {}...", &elem.key, message));
    }
}