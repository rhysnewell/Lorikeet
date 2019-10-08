extern crate lorikeet;
use lorikeet::bam_generator::*;
use lorikeet::filter;
use lorikeet::external_command_checker;
use lorikeet::mapping_parameters::*;
use lorikeet::shard_bam_reader::*;
use lorikeet::FlagFilter;
use lorikeet::genome_exclusion::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate bio;
use bio::alignment::sparse::*;

use std::env;
use std::str;
use std::process;
use std::collections::BTreeMap;
use std::io::Write;
use std::fs::File;
use std::path::Path;

extern crate itertools;
use itertools::Itertools;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate tempfile;
use tempfile::NamedTempFile;
#[macro_use]
extern crate lazy_static;

const CONCATENATED_REFERENCE_CACHE_STEM: &str = "lorikeet-genome";

fn filter_full_help() -> &'static str {
    "lorikeet filter: Remove alignments with insufficient identity.

Only primary, non-supplementary alignments are considered, and output files
are grouped by reference, but not sorted by position.

Files (both required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM file(s)
   -o, --output-bam-files <PATH> ..      Path to corresponding output file(s)

Thresholds:
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only. [default: 0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Implies --proper-pairs-only. [default 0.0]
   --proper-pairs-only                   Require reads to be mapped as proper pairs

Other:
   -t, --threads <INT>                   Number of threads for output compression
                                         [default 1]
   --inverse                             Only keep reads which are unmapped or
                                         align below thresholds. Note that output
                                         records may still be marked as mapped
                                         if they do not meet the thresholds.
                                         [default false]
   --verbose                             Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Example usage:

  lorikeet filter -b in.bam -o out.bam --min-read-aligned-length 75

Rhys J.P. Newell <r.newell near uq.edu.au>"
}

fn polymorph_full_help() -> &'static str {
    "lorikeet contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly.
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs
   -t, --threads <INT>                   Number of threads for mapping / sorting
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Sharding i.e. multiple reference sets (optional):
   --sharded                             If -b/--bam-files was used:
                                           Input BAM files are read-sorted alignments
                                           of a set of reads mapped to multiple
                                           reference contig sets. Choose the best
                                           hit for each read pair.

                                         Otherwise if mapping was carried out:
                                           Map reads to each reference, choosing the
                                           best hit for each pair.

Alignment filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --allow-improper-pairs                Allows reads to be mapped as improper pairs

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   -k, --kmer-size <INT>                 K-mer size used to generate k-mer frequency
                                         table. [default: 4]
   -d, --depth-threshold <INT>           Minimum number of reads needed to verify
                                         a variant. [default: 50]
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 40]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --variant-fraction-threshold      Minimum percentage fraction of the depth
                                         threshold value a variant must occur at
                                         for it to be verified. I.e. If depth threshold
                                         is 50 and variant fraction is 0.1, then
                                         at least 5 reads need to contain a variant
                                         for it to be whitelisted [default: 0.1]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   -t, --threads                         Number of threads used. [default: 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>"
}

fn summarize_full_help() -> &'static str {
    "lorikeet contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly.
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs
   -t, --threads <INT>                   Number of threads for mapping / sorting
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Sharding i.e. multiple reference sets (optional):
   --sharded                             If -b/--bam-files was used:
                                           Input BAM files are read-sorted alignments
                                           of a set of reads mapped to multiple
                                           reference contig sets. Choose the best
                                           hit for each read pair.

                                         Otherwise if mapping was carried out:
                                           Map reads to each reference, choosing the
                                           best hit for each pair.

Alignment filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Conflicts --allow-improper-pairs. [default 0.0]
   --allow-improper-pairs                Allows reads to be mapped as improper pairs

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   -k, --kmer-size <INT>                 K-mer size used to generate k-mer frequency
                                         table. [default: 4]
   -d, --depth-threshold <INT>           Minimum number of reads needed to verify
                                         a variant. [default: 50]
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 40]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --variant-fraction-threshold      Minimum percentage fraction of the depth
                                         threshold value a variant must occur at
                                         for it to be verified. I.e. If depth threshold
                                         is 50 and variant fraction is 0.1, then
                                         at least 5 reads need to contain a variant
                                         for it to be whitelisted [default: 0.1]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   -t, --threads                         Number of threads used. [default: 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>"}

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {

        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            if m.is_present("full-help") {
                println!("{}", filter_full_help());
                process::exit(1);
            }
            set_log_level(m, true);

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let output_bam_files: Vec<&str> = m.values_of("output-bam-files").unwrap().collect();
            if bam_files.len() != output_bam_files.len() {
                panic!("The number of input BAM files must be the same as the number output")
            }

            let filter_params = FilterParameters::generate_from_clap(m);

            let num_threads = value_t!(m.value_of("threads"), u16).unwrap();

            for (bam, output) in bam_files.iter().zip(output_bam_files.iter()) {
                let reader = bam::Reader::from_path(bam).expect(
                    &format!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer = bam::Writer::from_path(
                    output,
                    &header
                ).expect(&format!("Failed to write BAM file {}", output));
                writer.set_threads(num_threads as usize).expect("Failed to set num threads in writer");
                let mut filtered = filter::ReferenceSortedBamFilter::new(
                    reader,
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    !m.is_present("inverse"));

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        },
        Some("polymorph") => {
            let m = matches.subcommand_matches("polymorph").unwrap();
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    let bam_readers = lorikeet::bam_generator::generate_filtered_bam_readers_from_bam_files(
                        bam_files,
                        filter_params.flag_filters.clone(),
                        filter_params.min_aligned_length_single,
                        filter_params.min_percent_identity_single,
                        filter_params.min_aligned_percent_single,
                        filter_params.min_aligned_length_pair,
                        filter_params.min_percent_identity_pair,
                        filter_params.min_aligned_percent_pair);
                    run_pileup(m,
                               bam_readers,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               bam_readers,
                               filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               bam_readers,
                               filter_params.flag_filters);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup(m,
                               all_generators,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m, &None, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               generator_sets,
                               filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                               all_generators,
                               filter_params.flag_filters.clone());
                }
            }
        },
        Some("summarize") => {
            let m = matches.subcommand_matches("summarize").unwrap();
            if m.is_present("full-help") {
                println!("{}", summarize_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    let bam_readers = lorikeet::bam_generator::generate_filtered_bam_readers_from_bam_files(
                        bam_files,
                        filter_params.flag_filters.clone(),
                        filter_params.min_aligned_length_single,
                        filter_params.min_percent_identity_single,
                        filter_params.min_aligned_percent_single,
                        filter_params.min_aligned_length_pair,
                        filter_params.min_percent_identity_pair,
                        filter_params.min_aligned_percent_pair);
                    run_pileup_contigs(m,
                                       bam_readers,
                                       filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup_contigs(m,
                                       bam_readers,
                                       filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup_contigs(m,
                                       bam_readers,
                                       filter_params.flag_filters);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_pileup_contigs(m,
                                       all_generators,
                                       filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m, &None, &NoExclusionGenomeFilter{});
                    run_pileup_contigs(m,
                                       generator_sets,
                                       filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup_contigs(m,
                                       all_generators,
                                       filter_params.flag_filters.clone());
                }
            }
        },
        Some("kmer") => {
            let m = matches.subcommand_matches("kmer").unwrap();
            if m.is_present("full-help") {
//                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let reference_path = Path::new(m.value_of("reference").unwrap());
            let fasta_reader = match bio::io::fasta::Reader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(e) => {
                    eprintln!("Missing or corrupt fasta file {}", e);
                    process::exit(1);
                },
            };
            let contigs = fasta_reader.records().collect_vec();
            // Initialize bound contig variable
            let mut tet_freq = BTreeMap::new();
            let contig_count = contigs.len();
            let mut contig_idx = 0 as usize;
            let mut contig_names = vec![String::new(); contig_count];
            for contig in contigs{
                let contig = contig.unwrap();
                contig_names[contig_idx] = contig.id().to_string();
                debug!("Parsing contig: {}", contig.id());
                let kmers = hash_kmers(contig.seq(), 4);
                // Get kmer counts in a contig
                for (kmer, pos) in kmers {
                    let k = tet_freq.entry(kmer.to_vec()).or_insert(vec![0; contig_count]);
                    k[contig_idx] = pos.len();
                }
                contig_idx += 1;
            }
            for (idx, contig) in contig_names.iter().enumerate() {
                print!("{}\t", contig);
                for (_kmer, counts) in &tet_freq{
                    print!("{}\t", counts[idx])
                }
                print!("\n");
            }

        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}


fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("methods") { return false; }
            let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
            if methods.contains(&"metabat") {
                if methods.len() > 1 {
                    panic!("Cannot specify the metabat method with any other coverage methods")
                } else {
                    return true
                }
            }
            return false
        },
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false
        }
    }
}

#[derive(Debug)]
struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        let mut f = FilterParameters {
            flag_filters: FlagFilter {
                include_improper_pairs: m.is_present("allow-improper-pairs"),
                include_secondary: false,
                include_supplementary: false,
            },
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => value_t!(m.value_of("min-read-aligned-length"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_single: match m.is_present("min-read-percent-identity") {
                true => value_t!(m.value_of("min-read-percent-identity"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_percent_single: match m.is_present("min-read-aligned-percent") {
                true => value_t!(m.value_of("min-read-aligned-percent"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => value_t!(m.value_of("min-read-aligned-length-pair"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_pair: match m.is_present("min-read-percent-identity-pair") {
                true => value_t!(m.value_of("min-read-percent-identity-pair"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_percent_pair: match m.is_present("min-read-aligned-percent-pair") {
                true => value_t!(m.value_of("min-read-aligned-percent-pair"), f32).unwrap(),
                false => 0.0
            }
        };
        if doing_metabat(&m) {
            debug!("Setting single read percent identity threshold at 0.97 for \
                    MetaBAT adjusted coverage.");
            // we use >= where metabat uses >. Gah.
            f.min_percent_identity_single = 0.97001;
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
            f.flag_filters.include_secondary = true;
        }
        debug!("Filter parameters set as {:?}", f);
        return f
    }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0 || self.min_percent_identity_pair > 0.0 ||
            self.min_aligned_percent_single > 0.0 || self.min_aligned_percent_pair > 0.0 ||
            self.min_aligned_length_single > 0 || self.min_aligned_length_pair > 0
    }
}

fn get_sharded_bam_readers<'a, 'b, T>(
    m: &'a clap::ArgMatches,
    reference_tempfile: &'a Option<NamedTempFile>,
    genome_exclusion: &'b T)
    -> Vec<ShardedBamReaderGenerator<'b, T>>
where T: GenomeExclusion {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");
    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
    let mut bam_readers = vec![];

    for reference_wise_params in params {
        let index = lorikeet::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference
                        },
                        naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::shard_bam_reader::generate_named_sharded_bam_readers_from_reads(
                        index.index_path(),
                        p.read1,
                        p.read2,
                        p.read_format.clone(),
                        p.threads,
                        bam_file_cache(p.read1).as_ref().map(String::as_ref),
                        discard_unmapped,
                        p.bwa_options,
                        reference_tempfile.is_none()));
        }

        debug!("Finished BAM setup");
    };
    let gen = ShardedBamReaderGenerator {
        stoit_name: "stoita".to_string(),
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec!(gen);
}

fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    reference_tempfile: &'a Option<NamedTempFile>)
    -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
    let mut generator_set = vec!();
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = lorikeet::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference
                        },
                        naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::bam_generator::generate_named_bam_readers_from_reads(
                    index.index_path(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.bwa_options,
                    reference_tempfile.is_none()));
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index
        };
        generator_set.push(to_return);
    };
    return generator_set;
}

fn generate_cached_bam_file_name(directory: &str, reference: &str, read1_path: &str) -> String {
    debug!("Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
           directory, reference, read1_path);
    std::path::Path::new(directory).to_str()
        .expect("Unable to covert bam-file-cache-directory name into str").to_string()+"/"+
        &std::path::Path::new(reference).file_name()
        .expect("Unable to convert reference to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+"."+
        &std::path::Path::new(read1_path).file_name()
        .expect("Unable to convert read1 name to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+".bam"
}

fn setup_bam_cache_directory(cache_directory: &str) {
    let path = std::path::Path::new(cache_directory);
    if path.is_dir() {
        if path.metadata().expect("Unable to read metadata for cache directory")
            .permissions().readonly() {
                panic!("Cache directory {} does not appear to be writeable, not continuing",
                       cache_directory);
            } else {
                info!("Writing BAM files to already existing directory {}", cache_directory)
            }
    } else {
        match path.parent() {
            Some(parent) => {
                let parent2 = match parent == std::path::Path::new("") {
                    true => std::path::Path::new("."),
                    false => parent
                };
                if parent2.canonicalize().expect(
                    &format!("Unable to canonicalize parent of cache directory {}", cache_directory)).is_dir() {
                    if parent2.metadata().expect(
                        &format!("Unable to get metadata for parent of cache directory {}",
                                 cache_directory))
                        .permissions().readonly() {
                            panic!(
                                "The parent directory of the (currently non-existent) \
                                 cache directory {} is not writeable, not continuing",
                                cache_directory);
                        } else {
                            info!("Creating cache directory {}", cache_directory);
                            std::fs::create_dir(path).expect("Unable to create cache directory");
                        }
                } else {
                    panic!("The parent directory of the cache directory {} does not \
                            yet exist, so not creating that cache directory, and not continuing.",
                           cache_directory)
                }
            },
            None => {
                panic!("Cannot create root directory {}", cache_directory)
            }
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        panic!("Failed to create test file in bam cache directory: {}", tf_result.err().unwrap())
    }
}

fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    reference_tempfile: &Option<NamedTempFile>)
    -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
    let mut generator_set = vec!();
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let filter_params = FilterParameters::generate_from_clap(m);
        let index = lorikeet::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference
                        },
                        naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::bam_generator::generate_filtered_named_bam_readers_from_reads(
                    index.index_path(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    p.bwa_options,
                    discard_unmapped,
                    reference_tempfile.is_none()));
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

fn run_pileup<'a,
    R: lorikeet::bam_generator::NamedBamReader,
    T: lorikeet::bam_generator::NamedBamReaderGenerator<R>>(
    m: &clap::ArgMatches,
    bam_readers: Vec<T>,
    flag_filters: FlagFilter) {

    let mut variant_consensus_file = "uninit.fna".to_string();
    let print_zeros = !m.is_present("no-zeros");
    let var_fraction = m.value_of("variant-fraction-threshold").unwrap().parse().unwrap();
    let print_consensus = m.is_present("variant-consensus-fasta");
    if print_consensus {
        variant_consensus_file = m.value_of("variant-consensus-fasta").unwrap().to_string();
    }
    let depth_threshold = m.value_of("depth-threshold").unwrap().parse().unwrap();
    let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();

    let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
    let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Missing or corrupt fasta file {}", e);
            process::exit(1);
        },
    };
    let threads = m.value_of("threads").unwrap().parse().unwrap();

    let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();

    if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
        eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
        process::exit(1)
    }
    let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
    let min = value_t!(m.value_of("trim-min"), f32).unwrap();
    let max = value_t!(m.value_of("trim-max"), f32).unwrap();
    if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
        panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
    }

    lorikeet::pileups::pileup_variants(
        bam_readers,
        fasta_reader,
        print_zeros,
        flag_filters,
        depth_threshold,
        mapq_threshold,
        var_fraction,
        min,
        max,
        min_fraction_covered,
        contig_end_exclusion,
        variant_consensus_file,
        print_consensus,
        threads);

}

fn run_pileup_contigs<'a,
    R: lorikeet::bam_generator::NamedBamReader,
    T: lorikeet::bam_generator::NamedBamReaderGenerator<R>>(
    m: &clap::ArgMatches,
    bam_readers: Vec<T>,
    flag_filters: FlagFilter) {

    let print_zeros = !m.is_present("no-zeros");
    let var_fraction = m.value_of("variant-fraction-threshold").unwrap().parse().unwrap();
    let depth_threshold = m.value_of("depth-threshold").unwrap().parse().unwrap();
    let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
    let kmer_size = m.value_of("kmer-size").unwrap().parse().unwrap();
    let reference_path = Path::new(m.value_of("reference").unwrap());
    let fasta_reader = match bio::io::fasta::Reader::from_file(&reference_path){
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Missing or corrupt fasta file {}", e);
            process::exit(1);
        },
    };

    let output_prefix = m.value_of("output-prefix").unwrap();
    let threads = m.value_of("threads").unwrap().parse().unwrap();

    let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();

    if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
        eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
        process::exit(1)
    }
    let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
    let min = value_t!(m.value_of("trim-min"), f32).unwrap();
    let max = value_t!(m.value_of("trim-max"), f32).unwrap();
    if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
        panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
    }

    let contigs = fasta_reader.records().collect_vec();
    // Initialize bound contig variable
    let mut tet_freq = BTreeMap::new();
    let contig_count = contigs.len();
    let mut contig_idx = 0 as usize;
    let mut contig_names = vec![String::new(); contig_count];
    for contig in contigs{
        let contig = contig.unwrap();
        contig_names[contig_idx] = contig.id().to_string();
        debug!("Parsing contig: {}", contig.id());
        let kmers = hash_kmers(contig.seq(), kmer_size);
        // Get kmer counts in a contig
        for (kmer, pos) in kmers {
            let k = tet_freq.entry(kmer.to_vec()).or_insert(vec![0; contig_count]);
            k[contig_idx] = pos.len();
        }
        contig_idx += 1;
    }

    let file_name = output_prefix.to_string() + &"_".to_owned()
        + &kmer_size.clone().to_string() + &"mer_counts".to_owned()
        + &".tsv".to_owned();
    let file_path = Path::new(&file_name);
    let mut file_open = match File::create(file_path) {
        Ok(fasta) => fasta,
        Err(e) => {
            println!("Cannot create file {:?}", e);
            std::process::exit(1)
        },
    };
    for (tid, name) in contig_names.iter().enumerate() {
        write!(file_open, "{}\t",
               name).unwrap();
        for (_kmer, counts) in tet_freq.iter(){
            write!(file_open, "{}\t", counts[tid]).unwrap();
        }
        write!(file_open, "\n").unwrap();
    }

    let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Missing or corrupt fasta file or no index file {}", e);
            process::exit(1);
        },
    };

    lorikeet::pileups::pileup_contigs(
        bam_readers,
        fasta_reader,
        print_zeros,
        flag_filters,
        depth_threshold,
        mapq_threshold,
        var_fraction,
        min,
        max,
        min_fraction_covered,
        contig_end_exclusion,
        output_prefix,
        threads);

}


fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("lorikeet version {}", crate_version!());
    }
}

fn build_cli() -> App<'static, 'static> {
    // specify static lazily because need to define it at runtime.
    lazy_static! {
        static ref POLYMORPH_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet polymorph --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet polymorph --method metabat --bam-files my.bam --reference assembly.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet polymorph --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet polymorph"),
            ansi_term::Colour::Green.paint(
                "Generate variant positions along contigs with a sample"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate variant positions from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate variant positions using MetaBAT adjusted coverage from a sorted BAM file, saving
the unfiltered BAM files in the saved_bam_files folder:")
        ).to_string();

        static ref SUMMARIZE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet summarize --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet summarize --method metabat --bam-files my.bam --reference assembly.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet summarize --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet summarize"),
            ansi_term::Colour::Green.paint(
                "Summarizes contigs stats including mean coverage, mean genotypes \
                and kmer frequencies"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate contig stats across samples"),
            ansi_term::Colour::Purple.paint(
                "Example: Summarizes contigs defined in reference from a sorted BAM file:"),
        ).to_string();

        static ref FILTER_HELP: String = format!(
            "
                            {}
                     {}

{}

  lorikeet filter --bam-files input.bam --output-bam filtered.bam
    --min-read-aligned-length 50

{}

  lorikeet filter -b input.bam -o inverse_filtered.bam --inverse
    --min-read-percent-identity 0.95 --threads 16

See lorikeet filter --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet filter"),
            ansi_term::Colour::Green.paint(
                "Filter BAM file alignments"),
            ansi_term::Colour::Purple.paint(
                "Example: Filter a BAM file by removing alignments shorter than 50bp:"),
            ansi_term::Colour::Purple.paint(
                "Example: Filter inverse: Keep alignments that have <95% alignment identity\n\
                 and those which do map at all. Note that the output BAM file will likely\n\
                 records that are still mapped, but align with < 95% identity. Use 16\n\
                 threads for output compression:"),
        ).to_string();
    }

    return App::new("lorikeet")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .help("
Mapping coverage analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
\tpolymorph\tCalculate variants along contig positions
\tsummarize\tSummarizes contig stats from multiple samples

Less used utility subcommands:
\tkmer\tCalculate kmer frequencies within contigs
\tfilter\tRemove (or only keep) alignments with insufficient identity

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <r.newell near uq.edu.au>
")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("polymorph")
                .about("Perform variant calling analysis")
                .help(POLYMORPH_HELP.as_str())

                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true))
                .arg(Arg::with_name("sharded")
                    .long("sharded")
                    .required(false))
                .arg(Arg::with_name("read1")
                    .short("-1")
                    .multiple(true)
                    .takes_value(true)
                    .requires("read2")
                    .required_unless_one(
                        &["bam-files","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                    .short("-2")
                    .multiple(true)
                    .takes_value(true)
                    .requires("read1")
                    .required_unless_one(
                        &["bam-files","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                    .short("-c")
                    .long("coupled")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                    .long("interleaved")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","coupled","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                    .long("single")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","coupled","interleaved","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                    .short("-r")
                    .long("reference")
                    .takes_value(true)
                    .required(true))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(Arg::with_name("bwa-params")
                    .long("bwa-params")
                    .long("bwa-parameters")
                    .takes_value(true)
                    .allow_hyphen_values(true)
                    .requires("reference"))
                .arg(Arg::with_name("discard-unmapped")
                    .long("discard-unmapped")
                    .requires("bam-file-cache-directory"))

                .arg(Arg::with_name("min-read-aligned-length")
                    .long("min-read-aligned-length")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                    .long("min-read-percent-identity")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                    .long("min-read-aligned-percent")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("variant-consensus-fasta")
                    .long("variant-consensus-fasta")
                    .takes_value(true)
                    .required(false))
                .arg(Arg::with_name("methods")
                    .short("m")
                    .long("method")
                    .long("methods")
                    .takes_value(true)
                    .multiple(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("variant-fraction-threshold")
                    .long("variant-fraction-threshold")
                    .short("f")
                    .default_value("0.1"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("40"))
                .arg(Arg::with_name("contig-end-exclusion")
                    .long("contig-end-exclusion")
                    .default_value("75"))
                .arg(Arg::with_name("trim-min")
                    .long("trim-min")
                    .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                    .long("trim-max")
                    .default_value("0.95"))
                .arg(Arg::with_name("no-zeros")
                    .long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs")
                    .long("allow-improper-pairs"))
                .arg(Arg::with_name("verbose")
                    .short("v")
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .long("quiet")))
        .subcommand(
            SubCommand::with_name("summarize")
                .about("Perform variant calling analysis and then binning")
                .help(SUMMARIZE_HELP.as_str())
                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true))
                .arg(Arg::with_name("sharded")
                    .long("sharded")
                    .required(false))
                .arg(Arg::with_name("read1")
                    .short("-1")
                    .multiple(true)
                    .takes_value(true)
                    .requires("read2")
                    .required_unless_one(
                        &["bam-files","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                    .short("-2")
                    .multiple(true)
                    .takes_value(true)
                    .requires("read1")
                    .required_unless_one(
                        &["bam-files","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                    .short("-c")
                    .long("coupled")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","interleaved","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                    .long("interleaved")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","coupled","single","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                    .long("single")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-files","read1","coupled","interleaved","full-help"])
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                    .short("-r")
                    .long("reference")
                    .takes_value(true)
                    .required(true))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(Arg::with_name("bwa-params")
                    .long("bwa-params")
                    .long("bwa-parameters")
                    .takes_value(true)
                    .allow_hyphen_values(true)
                    .requires("reference"))
                .arg(Arg::with_name("discard-unmapped")
                    .long("discard-unmapped")
                    .requires("bam-file-cache-directory"))

                .arg(Arg::with_name("min-read-aligned-length")
                    .long("min-read-aligned-length")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                    .long("min-read-percent-identity")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                    .long("min-read-aligned-percent")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("methods")
                    .short("m")
                    .long("method")
                    .long("methods")
                    .takes_value(true)
                    .multiple(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("variant-fraction-threshold")
                    .long("variant-fraction-threshold")
                    .short("f")
                    .default_value("0.1"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("40"))
                .arg(Arg::with_name("kmer-size")
                    .long("kmer-size")
                    .short("k")
                    .default_value("4"))
                .arg(Arg::with_name("contig-end-exclusion")
                    .long("contig-end-exclusion")
                    .default_value("75"))
                .arg(Arg::with_name("trim-min")
                    .long("trim-min")
                    .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                    .long("trim-max")
                    .default_value("0.95"))
                .arg(Arg::with_name("no-zeros")
                    .long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs")
                    .long("allow-improper-pairs"))
                .arg(Arg::with_name("output-prefix")
                    .long("output-prefix")
                    .short("o")
                    .default_value("output"))
                .arg(Arg::with_name("verbose")
                    .short("v")
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .long("quiet")))
        .subcommand(
        SubCommand::with_name("kmer")
            .about("Generate kmer count matrix for contigs")
//                .help(CONTIG_HELP.as_str())
            .arg(Arg::with_name("full-help")
                .long("full-help"))
            .arg(Arg::with_name("reference")
                .short("-r")
                .long("reference")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("verbose")
                .short("v")
                .long("verbose")))
        .subcommand(
            SubCommand::with_name("filter")
                .about("Remove alignments with insufficient identity")
                .help(FILTER_HELP.as_str())

                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("output-bam-files")
                    .short("o")
                    .long("output-bam-files")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("inverse")
                    .long("inverse"))

                .arg(Arg::with_name("min-read-aligned-length")
                    .long("min-read-aligned-length")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                    .long("min-read-percent-identity")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                    .long("min-read-aligned-percent")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("allow-improper-pairs"))

                .arg(Arg::with_name("allow-improper-pairs")
                    .long("allow-improper-pairs"))
                .arg(Arg::with_name("threads")
                    .long("threads")
                    .short("t")
                    .default_value("1"))

                .arg(Arg::with_name("verbose")
                    // .short("v") // Do not use since could be confused with
                    // inverse (a la grep -v)
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .short("q")
                    .long("quiet")));
}
