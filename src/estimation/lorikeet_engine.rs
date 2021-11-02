use abundance::abundance_calculator_engine::AbundanceCalculatorEngine;
use assembly::assembly_region_walker::AssemblyRegionWalker;
use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::CoverageEstimator;
use coverm::FlagFilter;
use estimation::bams::index_bams::*;
use haplotype::haplotype_clustering_engine::HaplotypeClusteringEngine;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use model::variant_context::VariantContext;
use model::variant_context_utils::VariantContextUtils;
use rayon::prelude::*;
use reference::reference_reader::ReferenceReader;
use reference::reference_reader_utils::ReferenceReaderUtils;
use reference::reference_writer::ReferenceWriter;
use scoped_threadpool::Pool;
use std::cmp::min;
use std::collections::HashMap;
use std::fs::create_dir_all;
use std::path::Path;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::NamedTempFile;
use utils::utils::get_cleaned_sample_names;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// The main lorikeet engine, takes any number of reference genomes and reads/bam files and performs
/// read mapping, variant calling, consensus genome calling, and strain genotyping
///
/// @author Rhys Newell <rhys.newell@.hdr.qut.edu.au>
pub struct LorikeetEngine<'a> {
    args: &'a clap::ArgMatches<'a>,
    short_read_bam_count: usize,
    long_read_bam_count: usize,
    coverage_estimators: Vec<CoverageEstimator>,
    flag_filters: FlagFilter,
    genomes_and_contigs: GenomesAndContigs,
    concatenated_genomes: Option<NamedTempFile>,
    tmp_bam_file_cache: Option<TempDir>,
    reference_map: HashMap<usize, String>,
    references: Vec<&'a str>,
    multi: Arc<MultiProgress>,
    multi_inner: Arc<MultiProgress>,
    tree: Arc<Mutex<Vec<&'a Elem>>>,
    progress_bars: &'a Vec<Elem>,
    threads: usize,
    mode: &'a str,
    run_in_parallel: bool,
}

impl<'a> LorikeetEngine<'a> {
    pub fn apply_per_reference(&self) {
        let parallel_genomes = self
            .args
            .value_of("parallel-genomes")
            .unwrap()
            .parse()
            .unwrap();
        let mut pool = Pool::new(parallel_genomes);
        let n_threads = std::cmp::max(
            self.threads / min(parallel_genomes as usize, self.references.len()),
            2,
        );
        let output_prefix = match self.args.is_present("output-directory") {
            true => {
                match std::fs::create_dir_all(
                    self.args.value_of("output-directory").unwrap().to_string(),
                ) {
                    Ok(_) => {}
                    Err(err) => panic!("Unable to create output directory {:?}", err),
                };
                self.args.value_of("output-directory").unwrap()
            }
            false => "./",
        };

        pool.scoped(|scope| {
            Self::begin_tick(0, &self.progress_bars, &self.multi_inner, "");
            Self::begin_tick(1, &self.progress_bars, &self.multi_inner, "");

            for (ref_idx, reference_stem) in self.reference_map.clone().into_iter() {
                let mode = self.mode;
                let multi_inner = &self.multi_inner;
                let tree = &self.tree;
                let progress_bars = &self.progress_bars;
                let flag_filters = &self.flag_filters;
                let reference_map = &self.reference_map;
                let references = &self.references;
                let tmp_bam_file_cache = match self.tmp_bam_file_cache.as_ref() {
                    Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
                    None => None,
                };
                let concatenated_genomes = match self.concatenated_genomes.as_ref() {
                    Some(file) => Some(file.path().to_str().unwrap().to_string()),
                    None => None,
                };
                let mut coverage_estimators = self.coverage_estimators.clone();
                let genomes_and_contigs = self.genomes_and_contigs.clone();

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
                    let cache = glob::glob(&format!(
                        "{}/*{}",
                        &output_prefix,
                        if mode == "call" {
                            ".vcf"
                        } else if mode == "genotype" {
                            "strain_coverages.tsv"
                        } else if mode == "consensus" {
                            "consensus_*.fna"
                        } else {
                            ".vcf"
                        }
                    ))
                    .expect("failed to interpret glob")
                    .map(|p| {
                        p.expect("Failed to read cached vcf path")
                            .to_str()
                            .unwrap()
                            .to_string()
                    });
                    if cache.count() > 0 {
                        {
                            let elem = &progress_bars[ref_idx + 2];
                            let pb = multi_inner.insert(ref_idx + 2, elem.progress_bar.clone());
                        }
                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];

                            pb.progress_bar.set_message(format!(
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
                                    .finish_with_message(format!("All genomes analyzed {}", "✔",));
                            }
                        }
                        {
                            let pb = &tree.lock().unwrap()[0];
                            pb.progress_bar.inc(
                                ((self.short_read_bam_count + self.long_read_bam_count) as u64) + 1,
                            );
                            pb.progress_bar.reset_eta();
                            let pos = pb.progress_bar.position();
                            let len = pb.progress_bar.length();
                            if pos >= len {
                                pb.progress_bar
                                    .finish_with_message(format!("All steps completed {}", "✔",));
                            }
                        }
                        continue;
                    }
                }

                scope.execute(move || {
                    let reference = &genomes_and_contigs.genomes[ref_idx];
                    Self::begin_tick(
                        ref_idx + 2,
                        &progress_bars,
                        &multi_inner,
                        "Preparing variants",
                    );

                    // Read BAMs back in as indexed
                    let mut indexed_bam_readers = recover_bams(
                        &self.args,
                        &concatenated_genomes,
                        self.short_read_bam_count,
                        self.long_read_bam_count,
                        &genomes_and_contigs,
                        n_threads as u32,
                        &tmp_bam_file_cache,
                        self.run_in_parallel,
                        ref_idx,
                    );

                    debug!("Indexed bam readers {:?}", &indexed_bam_readers);

                    // let mut reference_reader = ReferenceReader::new(
                    //     &Some(concatenated_genomes.as_ref().unwrap().to_string()),
                    //     genomes_and_contigs.clone(),
                    //     genomes_and_contigs.contig_to_genome.len(),
                    // );

                    let mut reference_reader = ReferenceReader::new(
                        &Some(references[ref_idx].to_string()),
                        genomes_and_contigs.clone(),
                        genomes_and_contigs.contig_to_genome.len(),
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
                        n_threads,
                    );

                    let mut shards = assembly_engine.collect_shards(
                        self.args,
                        &indexed_bam_readers,
                        &genomes_and_contigs,
                        &concatenated_genomes,
                        flag_filters,
                        n_threads,
                        tree,
                        &mut reference_reader,
                    );

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar.set_message(format!(
                            "{}: Performing variant calling on active regions...",
                            pb.key
                        ));
                    }
                    let contexts = assembly_engine.traverse(
                        shards,
                        flag_filters,
                        self.args,
                        &indexed_bam_readers,
                        &mut reference_reader,
                    );
                    debug!("example variant {:?}", &contexts.first());

                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(format!("{}: Variant calling complete...", pb.key));
                    }

                    let cleaned_sample_names = get_cleaned_sample_names(&indexed_bam_readers);

                    // ensure output path exists
                    create_dir_all(&output_prefix).expect("Unable to create output directory");
                    if mode == "call" {
                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];
                            pb.progress_bar
                                .set_message(format!("{}: Generating VCF file...", &reference,));
                        }
                        assembly_engine.evaluator.write_vcf(
                            &output_prefix,
                            &contexts,
                            &cleaned_sample_names,
                            &reference_reader,
                            false,
                        )
                    } else if mode == "genotype" {
                        // Calculate the geometric mean values and CLR for each variant, reference specific
                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];
                            pb.progress_bar.set_message(format!(
                                "{}: Starting clustering engine...",
                                &reference,
                            ));
                        }

                        // If a variant context contains more than one allele, we need to split
                        // this context into n different contexts, where n is number of variant
                        // alleles
                        let split_contexts = VariantContextUtils::split_contexts(
                            contexts,
                            self.args
                                .value_of("qual-by-depth-filter")
                                .unwrap()
                                .parse()
                                .unwrap(),
                        );

                        if split_contexts.len() >= 1 {
                            // Perform UMAP and HDBSCAN clustering followed by variant group
                            // read linkage clustering.
                            let mut clustering_engine = HaplotypeClusteringEngine::new(
                                output_prefix.as_str(),
                                split_contexts,
                                &reference_reader,
                                ref_idx,
                                indexed_bam_readers.len(),
                                n_threads,
                            );
                            let (n_strains, split_contexts) = clustering_engine.perform_clustering(
                                &indexed_bam_readers,
                                flag_filters,
                                n_threads,
                                tree,
                            );
                            debug!(
                                "example variant after clustering {:?}",
                                &split_contexts.first()
                            );

                            // Get strain abundances
                            {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar.set_message(format!(
                                    "{}: Calculating genotype abundances...",
                                    &reference,
                                ));
                            }
                            let mut abundance_calculator_engine = AbundanceCalculatorEngine::new(
                                split_contexts,
                                &reference_reader.genomes_and_contigs.genomes[ref_idx],
                                &output_prefix,
                                &cleaned_sample_names,
                            );

                            let (mut strain_ids_present, split_contexts) =
                                abundance_calculator_engine.run_abundance_calculator(
                                    n_strains,
                                    cleaned_sample_names.len(),
                                );
                            // let strain_ids_present = (0..n_strains).into_iter().collect::<Vec<usize>>();
                            {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar
                                    .set_message(
                                        format!("{}: Generating VCF file...", &reference,),
                                    );
                            }
                            assembly_engine.evaluator.write_vcf(
                                &output_prefix,
                                &split_contexts,
                                &cleaned_sample_names,
                                &reference_reader,
                                true,
                            );

                            // Write genotypes to disk, reference specific
                            {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar
                                    .set_message(format!("{}: Writing strains...", &reference,));
                            }
                            let mut reference_writer =
                                ReferenceWriter::new(reference_reader, &output_prefix);
                            reference_writer.generate_strains(
                                split_contexts,
                                ref_idx,
                                if strain_ids_present.len() > 0 {
                                    strain_ids_present
                                } else {
                                    vec![0]
                                },
                            );
                        } else {
                            // Write genotypes to disk, reference specific
                            {
                                let pb = &tree.lock().unwrap()[ref_idx + 2];
                                pb.progress_bar.set_message(format!(
                                    "{}: Writing reference strain...",
                                    &reference,
                                ));
                            }
                            let mut reference_writer =
                                ReferenceWriter::new(reference_reader, &output_prefix);
                            reference_writer.generate_strains(split_contexts, ref_idx, vec![0]);
                        }
                    } else if mode == "consensus" {
                        // Get sample distances
                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];
                            pb.progress_bar
                                .set_message(format!("{}: Generating VCF file...", &reference,));
                        }
                        assembly_engine.evaluator.write_vcf(
                            &output_prefix,
                            &contexts,
                            &cleaned_sample_names,
                            &reference_reader,
                            false,
                        );

                        {
                            let pb = &tree.lock().unwrap()[ref_idx + 2];
                            pb.progress_bar.set_message(format!(
                                "{}: Generating consensus genomes...",
                                &reference,
                            ));
                        }
                        // variant_matrix.generate_distances();
                        let mut reference_writer =
                            ReferenceWriter::new(reference_reader, &output_prefix);
                        reference_writer.generate_consensus(
                            contexts,
                            ref_idx,
                            &cleaned_sample_names,
                        );
                    };

                    // if self.args.is_present("bam-file-cache-directory") && self.run_in_parallel {
                    //     // clean up split bam files
                    //     if self.short_read_bam_count > 0 {
                    //         let cmd_string = format!(
                    //             "rm -rf {}/short/{}/",
                    //             self.args.value_of("bam-file-cache-directory").unwrap(),
                    //             &genomes_and_contigs.genomes[ref_idx]
                    //         );
                    //     }
                    // }
                    {
                        let pb = &tree.lock().unwrap()[ref_idx + 2];
                        pb.progress_bar
                            .set_message(format!("{}: All steps completed {}", &reference, "✔",));
                        pb.progress_bar.finish_and_clear();
                    }
                    {
                        let pb = &tree.lock().unwrap()[1];
                        pb.progress_bar.inc(1);
                        let pos = pb.progress_bar.position();
                        let len = pb.progress_bar.length();
                        if pos >= len {
                            pb.progress_bar
                                .finish_with_message(format!("All genomes analyzed {}", "✔",));
                        }
                    }
                    {
                        let pb = &tree.lock().unwrap()[0];
                        pb.progress_bar.inc(1);
                        let pos = pb.progress_bar.position();
                        let len = pb.progress_bar.length();
                        if pos >= len {
                            pb.progress_bar
                                .finish_with_message(format!("All steps completed {}", "✔",));
                        }
                    }
                });
            }

            // pb_main.finish_with_message("All genomes staged...");
            self.multi.join().unwrap();
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
                progress_bar: ProgressBar::new((short_sample_count + long_sample_count) as u64),
            };
            debug!("Reference {}", reference,);
            reference_map
                .entry(ref_idx)
                .or_insert_with(|| reference.to_string());
        }

        progress_bars[0] = Elem {
            key: "Operations remaining".to_string(),
            index: 0,
            progress_bar: ProgressBar::new(
                ((references.len() * (short_sample_count + long_sample_count)) + references.len())
                    as u64,
            ),
        };

        let sty_eta = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");

        let sty_aux = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");
        progress_bars
            .par_iter()
            .for_each(|pb| pb.progress_bar.set_style(sty_aux.clone()));
        progress_bars[0].progress_bar.set_style(sty_eta);

        return progress_bars;
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

        pb.set_message(format!("{}: {}...", &elem.key, message));
    }
}

pub fn start_lorikeet_engine<
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
    let threads = m.value_of("threads").unwrap().parse().unwrap();
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

    // All different counts of samples I need. Changes depends on when using concatenated genomes or not
    let short_read_bam_count = bam_readers.len();
    let mut long_read_bam_count = 0;
    let mut reference_count = references.len();

    let longreads = match longreads {
        Some(vec) => {
            long_read_bam_count += vec.len();
            vec
        }
        None => vec![],
    };

    let parallel_genomes: usize = m.value_of("parallel-genomes").unwrap().parse().unwrap();
    let mut run_in_parallel = false;
    if parallel_genomes > 1 && reference_count > 1 {
        run_in_parallel = true;
    };
    // Finish each BAM source
    if m.is_present("longreads") || m.is_present("longread-bam-files") {
        info!("Processing long reads...");
        finish_bams(
            longreads,
            threads,
            &genomes_and_contigs,
            run_in_parallel,
            !m.is_present("longread-bam-files"),
        );
    }

    if m.is_present("coupled")
        || m.is_present("interleaved")
        || m.is_present("read1")
        || m.is_present("read2")
        || m.is_present("single")
        || m.is_present("bam-files")
    {
        info!("Processing short reads...");
        finish_bams(
            bam_readers,
            threads,
            &genomes_and_contigs,
            run_in_parallel,
            !m.is_present("bam-files"),
        );
    }

    let mut reference_map = HashMap::new();

    // Set up multi progress bars
    let multi = Arc::new(MultiProgress::new());

    let multi_inner = Arc::clone(&multi);
    let mut progress_bars = LorikeetEngine::setup_progress_bars(
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

    {
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
            tree,
            progress_bars: &progress_bars,
            threads,
            mode,
            run_in_parallel,
        };

        lorikeet_engine.apply_per_reference();
    }
}
