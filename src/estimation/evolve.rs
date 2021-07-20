// use std;
// use std::collections::HashMap;
//
// use bird_tool_utils::command;
// use coverm::bam_generator::*;
// use estimation::bams::{index_bams::*, process_bam::*};
// use estimation::codon_structs::*;
// use estimation::variant_matrix::*;
// use estimation::vcfs::process_vcf::*;
// use external_command_checker;
// use utils::*;
// use haplotype::active_regions;
//
// use crate::*;
// use bio::io::gff;
// use coverm::genomes_and_contigs::GenomesAndContigs;
// use coverm::mosdepth_genome_coverage_estimators::*;
// use coverm::FlagFilter;
// use glob::glob;
// use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
// use rayon::prelude::*;
// use scoped_threadpool::Pool;
// use std::path::Path;
// use std::process::Stdio;
// use std::str;
// use std::sync::{Arc, Mutex};
// use tempdir::TempDir;
// use tempfile::NamedTempFile;
//
// #[allow(unused)]
// pub fn pileup_variants<
//     'a,
//     R: NamedBamReader,
//     S: NamedBamReaderGenerator<R>,
//     T: NamedBamReader,
//     U: NamedBamReaderGenerator<T>,
// >(
//     m: &clap::ArgMatches,
//     bam_readers: Vec<S>,
//     longreads: Option<Vec<U>>,
//     mode: &str,
//     coverage_estimators: &mut Vec<CoverageEstimator>,
//     flag_filters: FlagFilter,
//     mapq_threshold: u8,
//     min_var_depth: usize,
//     min: f32,
//     max: f32,
//     contig_end_exclusion: u64,
//     output_prefix: &str,
//     threads: usize,
//     method: &str,
//     coverage_fold: f32,
//     include_indels: bool,
//     include_soft_clipping: bool,
//     is_long_read: bool,
//     genomes_and_contigs: GenomesAndContigs,
//     tmp_bam_file_cache: Option<TempDir>,
//     concatenated_genomes: Option<NamedTempFile>,
// ) {
//     // TODO: We need to split up analyses per reference to help contain memory issues
//     let alpha: f64 = m.value_of("fdr-threshold").unwrap().parse().unwrap();
//     let parallel_genomes = m.value_of("parallel-genomes").unwrap().parse().unwrap();
//     let mut pool = Pool::new(parallel_genomes);
//     let n_threads = std::cmp::max(threads / parallel_genomes as usize, 2);
//
//     let references = parse_references(&m);
//     let references = references.iter().map(|p| &**p).collect::<Vec<&str>>();
//
//     retrieve_reference(&Some(
//         concatenated_genomes
//             .as_ref()
//             .unwrap()
//             .path()
//             .to_str()
//             .unwrap()
//             .to_string(),
//     ));
//
//     // All different counts of samples I need. Changes depends on when using concatenated genomes or not
//     let mut short_sample_count = bam_readers.len();
//     let mut long_sample_count = 0;
//     let mut reference_count = references.len();
//
//     let longreads = match longreads {
//         Some(vec) => {
//             long_sample_count += vec.len();
//             vec
//         }
//         None => vec![],
//     };
//
//
//
//     // Finish each BAM source
//     if m.is_present("longreads") || m.is_present("longread-bam-files") {
//         info!("Processing long reads...");
//         finish_bams(longreads, threads);
//     }
//
//     info!("Processing short reads...");
//     finish_bams(bam_readers, threads);
//
//
//     let mut reference_map = HashMap::new();
//
//     let mut progress_bars = utils::setup_progress_bars(
//         &references,
//         &mut reference_map,
//         &genomes_and_contigs,
//         short_samples_count,
//         long_sample_count,
//     );
//
//     debug!(
//         "{} Longread BAM files, {} Shortread BAM files {} Total BAMs over {} genome(s)",
//         long_sample_count,
//         short_sample_count,
//         (short_sample_count + long_sample_count),
//         reference_count
//     );
//
//
//
//     // Set up multi progress bars
//     let multi = Arc::new(MultiProgress::new());
//     let sty_eta = ProgressStyle::default_bar()
//         .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");
//
//     let sty_aux = ProgressStyle::default_bar()
//         .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");
//     progress_bars
//         .par_iter()
//         .for_each(|pb| pb.progress_bar.set_style(sty_aux.clone()));
//     progress_bars[0].progress_bar.set_style(sty_eta.clone());
//
//     // let pb_main = multi.add(ProgressBar::new(reference_map.keys().len() as u64));
//     // pb_main.set_style(sty_eta.clone());
//
//     let tree: Arc<Mutex<Vec<&Elem>>> =
//         Arc::new(Mutex::new(Vec::with_capacity(progress_bars.len())));
//     {
//         let mut tree = tree.lock().unwrap();
//         for pb in progress_bars.iter() {
//             tree.push(pb)
//         }
//     }
//     // let tree2 = Arc::clone(&tree);
//
//     let multi_inner = Arc::clone(&multi);
//
//     pool.scoped(|scope| {
//         {
//             // Total steps eta progress bar
//             let elem = &progress_bars[0];
//             let pb = multi_inner.insert(0, elem.progress_bar.clone());
//
//             pb.enable_steady_tick(500);
//
//             pb.set_message(&format!("{}...", &elem.key,));
//         }
//         {
//             // completed genomes progress bar
//             let elem = &progress_bars[1];
//             let pb = multi_inner.insert(1, elem.progress_bar.clone());
//
//             pb.enable_steady_tick(500);
//
//             pb.set_message(&format!("{}...", &elem.key,));
//         }
//         for (ref_idx, reference_stem) in reference_map.clone().into_iter() {
//             // let _ = std::thread::spawn(move || {
//             //     multi.join().unwrap();
//             // });
//
//             // let ref_idx = Arc::new(ref_idx);
//
//             // let ref_idx = Arc::new(ref_idx);
//             let multi_inner = &multi_inner;
//             let tree = &tree;
//             let progress_bars = &progress_bars;
//             let flag_filters = &flag_filters;
//             let reference_map = &reference_map;
//             let references = references.clone();
//             let tmp_bam_file_cache = match tmp_bam_file_cache.as_ref() {
//                 Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
//                 None => None,
//             };
//             let concatenated_genomes = match concatenated_genomes.as_ref() {
//                 Some(file) => Some(file.path().to_str().unwrap().to_string()),
//                 None => None,
//             };
//             let mut coverage_estimators = coverage_estimators.clone();
//             let genomes_and_contigs = &genomes_and_contigs;
//
//             let output_prefix = format!(
//                 "{}/{}",
//                 &output_prefix,
//                 Path::new(&reference_stem)
//                     .file_stem()
//                     .unwrap()
//                     .to_str()
//                     .unwrap(),
//             );
//
//             if Path::new(&output_prefix).exists() && !m.is_present("force") {
//                 let cache = glob(&format!("{}/*.vcf", &output_prefix))
//                     .expect("failed to interpret glob")
//                     .map(|p| {
//                         p.expect("Failed to read cached bam path")
//                             .to_str()
//                             .unwrap()
//                             .to_string()
//                     })
//                     .collect::<Vec<String>>();
//                 if cache.len() > 0 {
//                     {
//                         let elem = &progress_bars[ref_idx + 2];
//                         let pb = multi_inner.insert(ref_idx + 2, elem.progress_bar.clone());
//                     }
//                     {
//                         let pb = &tree.lock().unwrap()[ref_idx + 2];
//
//                         pb.progress_bar.set_message(&format!(
//                             "{}: Output already present. Run with --force to overwrite",
//                             &genomes_and_contigs.genomes[ref_idx]
//                         ));
//                         pb.progress_bar.finish_and_clear();
//                     }
//                     {
//                         let pb = &tree.lock().unwrap()[1];
//                         pb.progress_bar.inc(1);
//                         pb.progress_bar.reset_eta();
//                         let pos = pb.progress_bar.position();
//                         let len = pb.progress_bar.length();
//                         if pos >= len {
//                             pb.progress_bar
//                                 .finish_with_message(&format!("All genomes analyzed {}", "✔",));
//                         }
//                     }
//                     {
//                         let pb = &tree.lock().unwrap()[0];
//                         pb.progress_bar.inc(
//                             ((short_sample_count + long_sample_count)
//                                 as u64)
//                                 * 2
//                                 + 1,
//                         );
//                         pb.progress_bar.reset_eta();
//                         let pos = pb.progress_bar.position();
//                         let len = pb.progress_bar.length();
//                         if pos >= len {
//                             pb.progress_bar
//                                 .finish_with_message(&format!("All steps completed {}", "✔",));
//                         }
//                     }
//                     continue;
//                 }
//             }
//
//             scope.execute(move || {
//                 let reference = &genomes_and_contigs.genomes[ref_idx];
//
//                 {
//                     let elem = &progress_bars[ref_idx + 2];
//                     let pb = multi_inner.insert(ref_idx + 2, elem.progress_bar.clone());
//
//                     pb.enable_steady_tick(500);
//
//                     pb.set_message(&format!("{}: Preparing variants...", &elem.key, ));
//                     // multi.join().unwrap();
//
//                     // tree.lock().unwrap().insert(elem.index, &elem);
//                 }
//                 let mut codon_table = CodonTable::setup();
//
//                 // Read BAMs back in as indexed
//                 let mut indexed_bam_readers = recover_bams(
//                     m,
//                     &concatenated_genomes,
//                     short_sample_count,
//                     long_sample_count,
//                     &genomes_and_contigs,
//                     n_threads as u32,
//                     &tmp_bam_file_cache,
//                 );
//                 let mut per_reference_samples = 0;
//                 let mut per_reference_short_samples = 0;
//                 let mut variant_matrix = match concatenated_genomes {
//                     Some(ref concat) => {
//                         per_reference_samples =
//                             short_sample_count + long_sample_count;
//                         per_reference_short_samples = short_sample_count;
//                         debug!(
//                             "Per reference samples concatenated {}",
//                             &per_reference_samples
//                         );
//
//                         VariantMatrix::new_matrix(per_reference_samples)
//                     }
//                     None => {
//                         per_reference_samples =
//                             (short_sample_count + long_sample_count)
//                                 / references.len();
//                         per_reference_short_samples = short_sample_count / references.len();
//                         debug!(
//                             "Per reference samples not concatenated {}",
//                             &per_reference_samples
//                         );
//                         VariantMatrix::new_matrix(per_reference_samples)
//                     }
//                 };
//
//                 // Gff map lists coding regions
//                 let mut gff_map = HashMap::new();
//                 match mode {
//                     "evolve" => {
//                         codon_table.get_codon_table(11);
//                         // ani = 0.;
//
//                         let mut gff_reader;
//                         if m.is_present("gff") || !m.is_present("gff") {
//                             external_command_checker::check_for_prodigal();
//
//                             // create new fifo and give read, write and execute rights to the owner.
//                             let gff_dir = TempDir::new("lorikeet-prodigal")
//                                 .expect("unable to create prodigal directory");
//                             for reference in references.iter() {
//                                 let cmd_string = format!(
//                                     "set -e -o pipefail; \
//                                     prodigal -o {}/{}.gff -i {} -f gff {}",
//                                     // prodigal
//                                     gff_dir
//                                         .path()
//                                         .to_str()
//                                         .expect("Failed to convert tempfile path to str"),
//                                     Path::new(&reference).file_stem().unwrap().to_str().unwrap(),
//                                     &reference,
//                                     m.value_of("prodigal-params").unwrap(),
//                                 );
//                                 debug!("Queuing cmd_string: {}", cmd_string);
//                                 command::finish_command_safely(
//                                     std::process::Command::new("bash")
//                                         .arg("-c")
//                                         .arg(&cmd_string)
//                                         .stdout(Stdio::piped())
//                                         .stderr(Stdio::piped())
//                                         .spawn()
//                                         .expect("Unable to execute bash"),
//                                     "prodigal",
//                                 );
//
//                                 // Read in newly created gff
//                                 gff_reader = gff::Reader::from_file(
//                                     format!(
//                                         "{}/{}.gff",
//                                         gff_dir
//                                             .path()
//                                             .to_str()
//                                             .expect("Failed to convert tempfile path to str"),
//                                         Path::new(&reference)
//                                             .file_stem()
//                                             .unwrap()
//                                             .to_str()
//                                             .unwrap()
//                                     ),
//                                     bio::io::gff::GffType::GFF3,
//                                 )
//                                     .expect("Failed to read GFF file");
//
//                                 // Map to reference id
//                                 let ref_stem = Path::new(reference)
//                                     .file_stem()
//                                     .expect("problem determining file stem")
//                                     .to_str()
//                                     .unwrap()
//                                     .to_string();
//
//                                 let ref_idx = genomes_and_contigs.genome_index(&ref_stem).unwrap();
//                                 gff_reader.records().into_iter().for_each(|record| {
//                                     match record {
//                                         Ok(rec) => {
//                                             let gff_ref =
//                                                 gff_map.entry(ref_idx).or_insert(HashMap::new());
//                                             let contig_genes = gff_ref
//                                                 .entry(format!("{}~{}", &ref_stem, rec.seqname()))
//                                                 .or_insert(Vec::new());
//                                             contig_genes.push(rec);
//                                         }
//                                         _ => {}
//                                     };
//                                 });
//                             }
//
//                             gff_dir.close().expect("Failed to close temp directory");
//                         }
//                     }
//                     _ => {}
//                 }
//             }
//         }
//     }
// }