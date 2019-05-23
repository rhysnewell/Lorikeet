extern crate strainm;
use strainm::bam_generator::*;
use strainm::bam_generator::*;
use strainm::filter;
use strainm::external_command_checker;
use strainm::mapping_parameters::*;
use strainm::FlagFilter;
use strainm::genome::*;
use strainm::CONCATENATED_FASTA_FILE_SEPARATOR;
use strainm::genomes_and_contigs::GenomesAndContigs;


extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::env;
use std::str;
use std::process;

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


fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    let mut print_stream = &mut std::io::stdout();
    set_log_level(&matches, false);

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            if m.is_present("full-help") {
                println!("{}", genome_full_help());
                process::exit(1);
            }
            set_log_level(m, true);

            let filter_params = FilterParameters::generate_from_clap(m);
            let separator = parse_separator(m);

            let single_genome = m.is_present("single-genome");
            let genomes_and_contigs_option = match separator.is_some() || single_genome {
                true => None,
                false => {
                    let mut genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m);
                    info!("Reading contig names for {} genomes ..", genome_fasta_files.len());
                    Some(strainm::read_genome_fasta_files(
                        &genome_fasta_files.iter().map(|s| s.as_str()).collect()))
                }
            };

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    run_genome(
                        strainm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.flag_filters,
                            filter_params.min_aligned_length_single,
                            filter_params.min_percent_identity_single,
                            filter_params.min_aligned_percent_single,
                            filter_params.min_aligned_length_pair,
                            filter_params.min_percent_identity_pair,
                            filter_params.min_aligned_percent_pair),
                        m,
                        separator,
                        genomes_and_contigs_option);

                } else if m.is_present("read-sorted-shard-bam-files") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let mut bam_readers = strainm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads);
                    run_genome(
                        bam_readers,
                        m,
                        separator,
                        genomes_and_contigs_option);
                } else {
                    run_genome(
                        strainm::bam_generator::generate_named_bam_readers_from_bam_files(
                            bam_files),
                        m,
                        separator,
                        genomes_and_contigs_option);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();

                // Generate a temporary file of concatenated genomes if needed.
                let mut concatenated_genomes: Option<NamedTempFile> = None;
                if !m.is_present("reference") && !m.is_present("bam-files") {
                    info!("Generating reference FASTA file of concatenated genomes ..");
                    let list_of_genome_fasta_files = parse_list_of_genome_fasta_files(m);
                    concatenated_genomes = Some(
                        strainm::bwa_index_maintenance::generate_concatenated_fasta_file(
                            &list_of_genome_fasta_files));
                }

                if filter_params.doing_filtering() {
                    debug!("Mapping and filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m, &concatenated_genomes);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        separator,
                        genomes_and_contigs_option);
                } else if m.is_present("read-sorted-shard-bam-files") {
                    let generator_sets = get_sharded_bam_readers(m, &concatenated_genomes);
                    run_genome(
                        generator_sets,
                        m,
                        separator,
                        genomes_and_contigs_option);
                } else {
                    let generator_sets = get_streamed_bam_readers(m, &concatenated_genomes);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        separator,
                        genomes_and_contigs_option);
                };
            }
        },
        _ => {}
    }
}

fn parse_separator(m: &clap::ArgMatches) -> Option<u8> {
    let single_genome = m.is_present("single-genome");
    if single_genome {
        Some("0".as_bytes()[0])
    } else if m.is_present("separator") {
        let separator_str = m.value_of("separator").unwrap().as_bytes();
        if separator_str.len() != 1 {
            eprintln!(
                "error: Separator can only be a single character, found {} ({}).",
                separator_str.len(),
                str::from_utf8(separator_str).unwrap());
            process::exit(1);
        }
        Some(separator_str[0])
    } else if m.is_present("bam-files") || m.is_present("reference") {
        // Argument parsing enforces that genomes have been specified as FASTA
        // files.
        None
    } else {
        // Separator is set by strainm and written into the generated reference
        // fasta file.
        Some(CONCATENATED_FASTA_FILE_SEPARATOR.as_bytes()[0])
    }
}


fn parse_list_of_genome_fasta_files(m: &clap::ArgMatches) -> Vec<String> {
    match m.is_present("genome-fasta-files") {
        true => {
            m.values_of("genome-fasta-files").unwrap().map(|s| s.to_string()).collect()
        },
        false => {
            if m.is_present("genome-fasta-directory") {
                let dir = m.value_of("genome-fasta-directory").unwrap();
                let paths = std::fs::read_dir(dir).unwrap();
                let mut genome_fasta_files: Vec<String> = vec!();
                let extension = m.value_of("genome-fasta-extension").unwrap();
                for path in paths {
                    let file = path.unwrap().path();
                    match file.extension() {
                        Some(ext) => {
                            if ext == extension {
                                let s = String::from(file.to_string_lossy());
                                genome_fasta_files.push(s);
                            } else {
                                info!(
                                    "Not using directory entry '{}' as a genome FASTA file, as \
                                     it does not end with the extension '{}'",
                                    file.to_str().expect("UTF8 error in filename"),
                                    extension);
                            }
                        },
                        None => {
                            info!("Not using directory entry '{}' as a genome FASTA file",
                                  file.to_str().expect("UTF8 error in filename"));
                        }
                    }
                }
                if genome_fasta_files.len() == 0 {
                    panic!("Found 0 genomes from the genome-fasta-directory, cannot continue.");
                }
                genome_fasta_files // return
            } else {
                panic!("Either a separator (-s) or path(s) to genome FASTA files \
                        (with -d or -f) must be given");
            }
        }
    }
}


fn run_genome<'a,
    R: strainm::bam_generator::NamedBamReader,
    T: strainm::bam_generator::NamedBamReaderGenerator<R>>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    separator: Option<u8>,
    genomes_and_contigs_option: Option<GenomesAndContigs>) {

//    let print_zeros = !m.is_present("no-zeros");
    let proper_pairs_only = m.is_present("proper-pairs-only");
    let single_genome = m.is_present("single-genome");
    let reads_mapped = match separator.is_some() || single_genome {
        true => {
            strainm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator.unwrap(),
                proper_pairs_only,
                single_genome)
        },

        false => {
            strainm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                &genomes_and_contigs_option.unwrap(),
                proper_pairs_only)
        }
    };
}


fn build_cli() -> App<'static, 'static> {
    // specify static lazily because need to define it at runtime.
    lazy_static! {
        static ref CONTIG_HELP: String = format!(
            "
                            {}
              {}

{}

  strainm contig --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna

{}

  strainm contig --method metabat --bam-files my.bam
    --bam-file-cache-directory saved_bam_files

See strainm contig --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "strainm contig"),
            ansi_term::Colour::Green.paint(
                "Calculate coverage of individual contigs"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate mean coverage from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate MetaBAT adjusted coverage from a sorted BAM file, saving
the unfiltered BAM files in the saved_bam_files folder:")
        ).to_string();

        static ref GENOME_HELP: String = format!(
            "
                            {}
               {}

{}

  strainm genome --coupled read1.fastq.gz read2.fastq.gz
    --reference assembly.fna --separator '~'

{}

  strainm genome --bam-files my.bam --genome-fasta-directory genomes_directory/

See strainm genome --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "strainm genome"),
            ansi_term::Colour::Green.paint(
                "Calculate coverage of individual genomes"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference where the FASTA header separates
the genome name from the contig name with '~' e.g. >genome10~contig15"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate coverage of genomes defined as .fna files in
genomes_directory/ from a sorted BAM file:"),
        ).to_string();

        static ref FILTER_HELP: String = format!(
            "
                            {}
                     {}

{}

  strainm filter --bam-files input.bam --output-bam filtered.bam
    --min-read-aligned-length 50

{}

  strainm filter -b input.bam -o inverse_filtered.bam --inverse
    --min-read-percent-identity 0.95 --threads 16

See strainm filter --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "strainm filter"),
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

    let make_help: &'static str =
        "strainm make: Generate BAM files through mapping.

Output (required):
   -o, --output-directory <DIR>          Where generated BAM files will go

Mapping parameters:
   -r, --reference <PATH>                FASTA file(s) of contig(s) or BWA index stem.
                                         If multiple reference FASTA files are provided,
                                         reads will be mapped to each reference separately
                                         and create sharded BAM files.
   -t, --threads <INT>                   Number of threads to use for mapping
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
   --discard-unmapped                    Exclude unmapped reads from generated BAM files.

Example usage:

  strainm make -r combined_genomes.fna -1 read1.fq -2 read2.fq

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    return App::new("strainm")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .help("
Mapping coverage analysis for metagenomics

Usage: strainm <subcommand> ...

Main subcommands:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Less used utility subcommands:
\tmake\tGenerate BAM files through alignment
\tfilter\tRemove (or only keep) alignments with insufficient identity

Other options:
\t-V, --version\tPrint version information

Ben J. Woodcroft <benjwoodcroft near gmail.com>
")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")
                .help(GENOME_HELP.as_str())

                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true))
                .arg(Arg::with_name("read-sorted-shard-bam-files")
                    .long("read-sorted-shard-bam-files")
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
                    .multiple(true)
                    .conflicts_with("bam-files"))
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
                    .requires("reference")) // TODO: Relax this for autoconcatenation
                .arg(Arg::with_name("discard-unmapped")
                    .long("discard-unmapped")
                    .requires("bam-file-cache-directory"))

                .arg(Arg::with_name("separator")
                    .short("s")
                    .long("separator")
                    .conflicts_with("genome-fasta-files")
                    .conflicts_with("genome-fasta-directory")
                    .conflicts_with("single-genome")
                    .required_unless_one(
                        &["genome-fasta-files","genome-fasta-directory","single-genome","full-help"])
                    .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                    .short("f")
                    .long("genome-fasta-files")
                    .multiple(true)
                    .conflicts_with("separator")
                    .conflicts_with("genome-fasta-directory")
                    .conflicts_with("single-genome")
                    .required_unless_one(
                        &["separator","genome-fasta-directory","single-genome","full-help"])
                    .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                    .short("d")
                    .long("genome-fasta-directory")
                    .conflicts_with("separator")
                    .conflicts_with("genome-fasta-files")
                    .conflicts_with("single-genome")
                    .required_unless_one(
                        &["genome-fasta-files","separator","single-genome","full-help"])
                    .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                    .short("x")
                    .long("genome-fasta-extension")
                    // Unsure why, but uncommenting causes test failure - clap
                    // bug?
                    //.requires("genome-fasta-directory")
                    .default_value("fna")
                    .takes_value(true))
                .arg(Arg::with_name("single-genome")
                    .long("single-genome")
                    .conflicts_with("separator")
                    .conflicts_with("genome-fasta-files")
                    .conflicts_with("genome-fasta-directory"))

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
                    .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .requires("proper-pairs-only"))

                .arg(Arg::with_name("methods")
                    .short("m")
                    .long("method")
                    .long("methods")
                    .takes_value(true)
                    .multiple(true)
                    .possible_values(&[
                        "relative_abundance",
                        "mean",
                        "trimmed_mean",
                        "coverage_histogram",
                        "covered_fraction",
                        "covered_bases",
                        "variance",
                        "length",
                        "count",
                        "reads_per_base"])
                    .default_value("relative_abundance"))
                .arg(Arg::with_name("trim-min")
                    .long("trim-min")
                    .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                    .long("trim-max")
                    .default_value("0.95"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.10"))
                .arg(Arg::with_name("contig-end-exclusion")
                    .long("contig-end-exclusion")
                    .default_value("75"))
                .arg(Arg::with_name("no-zeros")
                    .long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
                .arg(Arg::with_name("output-format")
                    .long("output-format")
                    .possible_values(&["sparse","dense"])
                    .default_value("dense"))

                .arg(Arg::with_name("verbose")
                    .short("v")
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .short("q")
                    .long("quiet")))
}
