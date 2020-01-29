extern crate lorikeet;
use lorikeet::bam_generator::*;
use lorikeet::filter;
use lorikeet::external_command_checker;
use lorikeet::mapping_parameters::*;
use lorikeet::shard_bam_reader::*;
use lorikeet::FlagFilter;
use lorikeet::genome_exclusion::*;
use lorikeet::mosdepth_genome_coverage_estimators::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate seq_io;
use seq_io::fasta;

extern crate bio;
use bio::alignment::sparse::*;

use std::env;
use std::str;
use std::process;
use std::collections::BTreeMap;
use std::io::Write;
use std::process::Stdio;
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

const MAPPING_SOFTWARE_LIST: &[&str] = &["bwa-mem", "minimap2-sr", "minimap2-ont", "minimap2-pb","minimap2-no-preset"];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";
const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

const MAPPER_HELP: &'static str =
    "   -p, --mapper <NAME>                   Underlying mapping software used
                                         (\"minimap2-sr\", \"bwa-mem\", \"minimap2-ont\",
                                         \"minimap2-pb\", or \"minimap2-no-preset\").
                                         minimap2 -sr, -ont, -pb, -no-preset specify
                                         '-x' preset of minimap2 to be used
                                         (with map-ont, map-pb for -ont, -pb).
                                         [default: \"minimap2-sr\"]";

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
                                              bases must be aligned. [default 0.97]
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
   --proper-pairs-only                        Require reads to be mapped as proper pairs

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
    lazy_static! {
        static ref POLYMORPH_HELP: String = format!(
    "lorikeet polymorph: Print variant sites along a contig

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
{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not
                                         as a FASTA file.
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
                                              bases must be aligned. [default 0.97]
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
   --allow-improper-pairs                     Allows reads to be mapped as improper pairs
   --include-supplementary                    Includes read alignments flagged as supplementary
   --include-secondary                        Includes read alignments flagged as secondary

Other arguments (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         method is available at
                                         https://github.com/rhysnewell/lorikeet
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth               Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --include-indels                      Flag indicating whether to attempt to calculate INDEL sites
                                         Not recommended if using nanopore long read data.
   --coverage-fold                       Percentage value of coverage to look above and below
                                         when calculating variant locations. e.g. if coverage-fold
                                         is equal to 0.1, only areas of coverage * (1.0 - 0.1) and
                                         coverage * (1.0 + 0.1) will be considered.
                                         [default: 0.5]
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

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &POLYMORPH_HELP
}

fn evolve_full_help() -> &'static str {
    lazy_static! {
        static ref EVOLVE_HELP: String = format!(
    "lorikeet evolve: Calculate dN/dS values in coding regions based on variants found in read mappings

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
   --gff, -g <PATH>                      GFF3 file containing gene locations
                                         present in reference genome.
   --prodigal-params                     Paramaters passed onto prodigal to call
                                         gene locations on reference genome.
                                         -i and -f are already set. Only used if
                                         not GFF file is not premade.
{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not
                                         as a FASTA file.
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
                                              bases must be aligned. [default 0.97]
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
   --allow-improper-pairs                     Allows reads to be mapped as improper pairs
   --include-supplementary                    Includes read alignments flagged as supplementary
   --include-secondary                        Includes read alignments flagged as secondary

Other arguments (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         method is available at
                                         https://github.com/rhysnewell/lorikeet
   --include-indels                      Flag indicating whether to attempt to calculate INDEL sites
                                         Not recommended if using nanopore long read data.
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth               Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --coverage-fold                       Percentage value of coverage to look above and below
                                         when calculating variant locations. e.g. if coverage-fold
                                         is equal to 0.1, only areas of coverage * (1.0 - 0.1) and
                                         coverage * (1.0 + 0.1) will be considered.
                                         [default: 0.5]
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

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &EVOLVE_HELP
}


fn summarize_full_help() -> &'static str {
    lazy_static! {
        static ref SUMMARIZE_HELP: String = format!(
    "lorikeet summarize: Provides per contig variant statistics for a metagenome

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
{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not
                                         as a FASTA file.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]


Alignment filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.97]
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
   --include-supplementary               Includes read alignments flagged as supplementary
   --include-secondary                   Includes read alignments flagged as secondary

Other arguments (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth      Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --include-indels                      Flag indicating whether to attempt to calculate INDEL sites
                                         Not recommended if using nanopore long read data.
   --coverage-fold                       Percentage value of coverage to look above and below
                                         when calculating variant locations. e.g. if coverage-fold
                                         is equal to 0.1, only areas of coverage * (1.0 - 0.1) and
                                         coverage * (1.0 + 0.1) will be considered.
                                         [default: 0.5]
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

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &SUMMARIZE_HELP
}

fn genotype_full_help() -> &'static str {
    lazy_static! {
        static ref GENOTYPE_HELP: String = format!(
    "lorikeet genotype: Resolves strain-level genotypes and abundance from metagenomes

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
{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not
                                         as a FASTA file.
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
                                         bases must be aligned. [default 0.97]
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
   --include-supplementary                    Includes read alignments flagged as supplementary
   --include-secondary                        Includes read alignments flagged as secondary

Other arguments (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   -k, --kmer-size <INT>                 K-mer size used to generate k-mer frequency
                                         table. [default: 4]
   --include-indels                      Flag indicating whether to attempt to calculate INDEL sites
                                         Not recommended if using nanopore long read data.
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth      Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --coverage-fold                       Percentage value of coverage to look above and below
                                         when calculating variant locations. e.g. if coverage-fold
                                         is equal to 0.1, only areas of coverage * (1.0 - 0.1) and
                                         coverage * (1.0 + 0.1) will be considered.
                                         [default: 0.5]
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

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &GENOTYPE_HELP
}

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
                    &header,
                    rust_htslib::bam::Format::BAM)
                    .expect(&format!("Failed to write BAM file {}", output));
                writer
                    .set_threads(num_threads as usize)
                    .expect("Failed to set num threads in writer");
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
            let mode = "polymorph";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
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
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
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
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters.clone());
                }
            }
        },
        Some("summarize") => {
            let m = matches.subcommand_matches("summarize").unwrap();
            let mode = "summarize";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
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
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters);
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
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
                                       mode,
                                       &mut estimators,
                                       all_generators,
                                       filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       generator_sets,
                                       filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       all_generators,
                                       filter_params.flag_filters.clone());
                }
            }
        },
        Some("genotype") => {
            let m = matches.subcommand_matches("genotype").unwrap();
            let mode = "genotype";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", genotype_full_help());
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
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
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
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters.clone());
                }
            }
        },
        Some("evolve") => {
            let m = matches.subcommand_matches("evolve").unwrap();
            let mode = "evolve";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
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
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers = lorikeet::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                        bam_files, sort_threads, &NoExclusionGenomeFilter{});
                    run_pileup(m, mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       bam_readers,
                                       filter_params.flag_filters);
                }
            } else if m.is_present("read1") |
                m.is_present("interleaved") |
                m.is_present("coupled") |
                m.is_present("single"){
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
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
                                       mode,
                                       &mut estimators,
                                       all_generators,
                                       filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       generator_sets,
                                       filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                                       mode,
                                       &mut estimators,
                                       all_generators,
                                       filter_params.flag_filters.clone());
                }
            }
        },
        Some("polish") => {
            let m = matches.subcommand_matches("polish").unwrap();
            let mode = "polish";
            let mut estimators = EstimatorsAndTaker::generate_from_clap(m);
            if m.is_present("full-help") {
                println!("{}", polymorph_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

            if m.is_present("bam-file") {
                let bam_files: Vec<&str> = m.values_of("bam-file").unwrap().collect();
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
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                } else {
                    let bam_readers = lorikeet::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               bam_readers,
                               filter_params.flag_filters);
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
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
                               mode,
                               &mut estimators,
                               all_generators,
                               filter_params.flag_filters);
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_pileup(m,
                               mode,
                               &mut estimators,
                               generator_sets,
                               filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_pileup(m,
                               mode,
                               &mut estimators,
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

fn setup_mapping_index(
    reference_wise_params: &SingleReferenceMappingParameters,
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
) -> Option<Box<dyn lorikeet::mapping_index_maintenance::MappingIndex>> {
    match mapping_program {
        MappingProgram::BWA_MEM => Some(lorikeet::mapping_index_maintenance::generate_bwa_index(
            reference_wise_params.reference,
            None
        )),
        MappingProgram::MINIMAP2_SR |
        MappingProgram::MINIMAP2_ONT |
        MappingProgram::MINIMAP2_PB |
        MappingProgram::MINIMAP2_NO_PRESET => {
            if m.is_present("minimap2-reference-is-index") || reference_wise_params.len() == 1 {
                info!("Not pre-generating minimap2 index");
                if m.is_present("minimap2-reference-is-index") {
                    warn!("Minimap2 uses mapping parameters defined when the index was created, \
                    not parameters defined when mapping. Proceeding on the assumption that you \
                    passed the correct parameters when creating the minimap2 index.");
                }
                None
            } else {
                Some(lorikeet::mapping_index_maintenance::generate_minimap2_index(
                    reference_wise_params.reference,
                    Some(m.value_of("minimap2-params").unwrap_or("")),
                    mapping_program,
                ))
            }
        }
    }
}

fn parse_mapping_program(m: &clap::ArgMatches) -> MappingProgram {
    let mapping_program = match m.value_of("mapper") {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!(
            "Unexpected definition for --mapper: {:?}",
            m.value_of("mapper")
        ),
    };
    match mapping_program {
        MappingProgram::BWA_MEM => {
            external_command_checker::check_for_bwa();
        }
        MappingProgram::MINIMAP2_SR |
        MappingProgram::MINIMAP2_ONT |
        MappingProgram::MINIMAP2_PB |
        MappingProgram::MINIMAP2_NO_PRESET => {
            external_command_checker::check_for_minimap2();
        }
    }
    return mapping_program;
}

fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                process::exit(1);
            }
            info!("Using {} {}%", parameter, percentage * 100.0);
            percentage
        }
        false => 0.0,
    }
}


fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("method") { return false; }
            let methods: &str = m.value_of("method").unwrap();
            if methods.contains(&"metabat") {
                return true
            }
            return false
        },
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false
        }
    }
}

struct EstimatorsAndTaker {
    estimators: Vec<CoverageEstimator>,
    columns_to_normalise: Vec<usize>,
    rpkm_column: Option<usize>,
}

impl EstimatorsAndTaker {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> EstimatorsAndTaker {
        let mut estimators = vec![];
        let min_fraction_covered = parse_percentage(&m, "min-covered-fraction");
        let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();

        let methods: Vec<&str> = m.values_of("method").unwrap().collect();
        let mut columns_to_normalise: Vec<usize> = vec![];

//        let output_format = m.value_of("output-format").unwrap();
        let mut rpkm_column = None;

        if doing_metabat(&m) {
            estimators.push(CoverageEstimator::new_estimator_length());
            estimators.push(CoverageEstimator::new_estimator_mean(
                min_fraction_covered,
                contig_end_exclusion,
                false,
            ));
            estimators.push(CoverageEstimator::new_estimator_variance(
                min_fraction_covered,
                contig_end_exclusion,
            ));

            debug!("Cached regular coverage taker for metabat mode being used");

        } else {
            for (_i, method) in methods.iter().enumerate() {
                match method {
                    &"mean" => {
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        )); // TODO: Parameterise exclude_mismatches

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"trimmed_mean" => {
                        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                            error!(
                                "error: Trim bounds must be between 0 and 1, and \
                                 min must be less than max, found {} and {}",
                                min, max
                            );
                            process::exit(1);
                        }
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min,
                            max,
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    _ => unreachable!(),
                };
            }
        }

        // Check that min-covered-fraction is being used as expected
        if min_fraction_covered != 0.0 {
            let die = |estimator_name| {
                error!(
                    "The '{}' coverage estimator cannot be used when \
                     --min-covered-fraction is > 0 as it does not calculate \
                     the covered fraction. You may wish to set the \
                     --min-covered-fraction to 0 and/or run this estimator \
                     separately.",
                    estimator_name
                );
                process::exit(1)
            };
            for e in &estimators {
                match e {
                    CoverageEstimator::ReadCountCalculator { .. } => die("counts"),
                    CoverageEstimator::ReferenceLengthCalculator { .. } => die("length"),
                    CoverageEstimator::ReadsPerBaseCalculator { .. } => die("reads_per_base"),
                    _ => {}
                }
            }
        }

        return EstimatorsAndTaker {
            estimators: estimators,
            columns_to_normalise: columns_to_normalise,
            rpkm_column: rpkm_column,
        };
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
                include_secondary: m.is_present("include-secondary"),
                include_supplementary: m.is_present("include-supplementary"),
            },
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => value_t!(m.value_of("min-read-aligned-length"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => value_t!(m.value_of("min-read-aligned-length-pair"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_pair: parse_percentage(&m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(&m, "min-read-aligned-percent-pair"),

        };
        if m.is_present("nanopore") {
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
        }
        if doing_metabat(&m) {
            debug!(
                "Setting single read percent identity threshold at 0.97 for \
                 MetaBAT adjusted coverage."
            );
            // we use >= where metabat uses >. Gah.
            f.min_percent_identity_single = 0.97001;
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
            f.flag_filters.include_secondary = true;
        }
        debug!("Filter parameters set as {:?}", f);
        return f;
    }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0
            || self.min_percent_identity_pair > 0.0
            || self.min_aligned_percent_single > 0.0
            || self.min_aligned_percent_pair > 0.0
            || self.min_aligned_length_single > 0
            || self.min_aligned_length_pair > 0;
    }
}

fn get_sharded_bam_readers<'a, 'b, T>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    genome_exclusion: &'b T,
) -> Vec<ShardedBamReaderGenerator<'b, T>>
    where
        T: GenomeExclusion,
{
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");
    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut bam_readers = vec![];
    let mut concatenated_reference_name: Option<String> = None;
    let mut concatenated_read_names: Option<String> = None;

    for reference_wise_params in params {
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let reference_name = std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string();
        concatenated_reference_name = match concatenated_reference_name {
            Some(prev) => Some(format!("{}|{}", prev, reference_name)),
            None => Some(reference_name),
        };
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::shard_bam_reader::generate_named_sharded_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                ),
            );
            let name = &std::path::Path::new(p.read1)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string();
            concatenated_read_names = match concatenated_read_names {
                Some(prev) => Some(format!("{}|{}", prev, name)),
                None => Some(name.to_string()),
            };
        }

        debug!("Finished BAM setup");
    }
    let gen = ShardedBamReaderGenerator {
        stoit_name: format!(
            "{}/{}",
            concatenated_reference_name.unwrap(),
            concatenated_read_names.unwrap()
        ),
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec![gen];
}

fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

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
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::bam_generator::generate_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
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
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

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
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                lorikeet::bam_generator::generate_filtered_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
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
                    p.mapping_options,
                    discard_unmapped,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

fn generate_faidx(m: &clap::ArgMatches) -> bio::io::fasta::IndexedReader<File> {
    external_command_checker::check_for_samtools();
    info!("Generating reference index");
    let cmd_string = format!(
        "set -e -o pipefail; \
                     samtools faidx {}",
        m.value_of("reference").unwrap());
    debug!("Queuing cmd_string: {}", cmd_string);

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd_string)
        .stdout(Stdio::piped())
        .output()
        .expect("Unable to execute bash");

    let reference_path = Path::new(m.value_of("reference").unwrap());
    return bio::io::fasta::IndexedReader::from_file(&reference_path).expect(
        "Unable to generate index")
}

fn run_pileup<'a,
    R: lorikeet::bam_generator::NamedBamReader,
    T: lorikeet::bam_generator::NamedBamReaderGenerator<R>>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
    bam_readers: Vec<T>,
    flag_filters: FlagFilter) {
    match mode {
        "polymorph" => {
            let mut variant_consensus_file = "uninit.fna".to_string();
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let output_prefix = m.value_of("output-prefix").unwrap().to_string();
            let include_indels = m.is_present("include-indels");

            // Make sure we are dealing with a fresh file
            let file_name = output_prefix.to_string()
                + &".fna".to_owned();

            let file_path = Path::new(&file_name);

            let mut file_open = File::create(file_path)
                .expect("No Read or Write Permission in current directory");

            let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let threads = m.value_of("threads").unwrap().parse().unwrap();


            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }


            info!("Beginning polymorph with {} bam readers and {} threads", bam_readers.len(), threads);
            println!("sample\ttid\tpos\tvariant\treference\tvariant_depth\tdepth\tgenotypes\tvaf_cluster\tgenotype");
            lorikeet::pileups::pileup_variants(
                m,
                bam_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
                print_zeros,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                &output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                false);
        },
        "genotype" => {
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let kmer_size = m.value_of("kmer-size").unwrap().parse().unwrap();
            let reference_path = Path::new(m.value_of("reference").unwrap());
            let fasta_reader = match fasta::Reader::from_path(reference_path){
                Ok(reader) => reader,
                Err(e) => {
                    eprintln!("Missing or corrupt fasta file {}", e);
                    process::exit(1);
                },
            };

            let output_prefix = m.value_of("output-prefix").unwrap();
            let include_indels = m.is_present("include-indels");
            let include_soft_clipping = m.is_present("include-soft-clipping");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            let contigs = fasta_reader.into_records().collect_vec();
            // Initialize bound contig variable
            let mut tet_freq = BTreeMap::new();
            let contig_count = contigs.len();
            let mut contig_idx = 0 as usize;
            let mut contig_names = vec![String::new(); contig_count];
            info!("Calculating K-mer frequencies");
            for contig in contigs{
                let contig = contig.unwrap();
                contig_names[contig_idx] = String::from_utf8(contig.head).unwrap();
                let kmers = hash_kmers(&contig.seq, kmer_size);
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
                Err(_e) => generate_faidx(m),
            };

            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            lorikeet::pileups::pileup_variants(
                m,
                bam_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
                print_zeros,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                include_soft_clipping);
        },
        "summarize" => {
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let reference_path = Path::new(m.value_of("reference").unwrap());

            let output_prefix = m.value_of("output-prefix").unwrap();
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };

            info!("Beginning summarize with {} bam readers and {} threads", bam_readers.len(), threads);
            lorikeet::pileups::pileup_variants(
                m,
                bam_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
                print_zeros,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                coverage_fold,
                include_indels,
                false);
        },
        "evolve" => {
            let mut variant_consensus_file = "uninit.fna".to_string();
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path){
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let include_indels = m.is_present("include-indels");

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }

            info!("Beginning evolve with {} bam readers and {} threads", bam_readers.len(), threads);
            println!("gene\tstart\tend\tframe\tstrand\tdnds\tposition\tvariant\treference\tabundance\tdepth\tinfo");
            lorikeet::pileups::pileup_variants(
                m,
                bam_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
                print_zeros,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                "",
                threads,
                method,
                coverage_fold,
                include_indels,
                false);
        },
        "polish" => {
            let mut variant_consensus_file = "uninit.fna".to_string();
            let print_zeros = !m.is_present("no-zeros");
            let var_fraction = m.value_of("min-variant-depth").unwrap().parse().unwrap();
            let output_prefix = m.value_of("output-prefix").unwrap().to_string();

            // Make sure we are dealing with a fresh file
            let file_name = output_prefix.to_string()
                + &".fna".to_owned();

            let file_path = Path::new(&file_name);

            let mut file_open = File::create(file_path)
                .expect("No Read or Write Permission in current directory");

            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();
            let coverage_fold = m.value_of("coverage-fold").unwrap().parse().unwrap();
            let method = m.value_of("method").unwrap();

            let reference_path = Path::new(m.value_of("reference").unwrap());
//            let index_path = reference_path.clone().to_owned() + ".fai";
            let fasta_reader = match bio::io::fasta::IndexedReader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(_e) => generate_faidx(m),
            };
            let threads = m.value_of("threads").unwrap().parse().unwrap();


            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
            }


            info!("Beginning polishing with {} bam readers and {} threads", bam_readers.len(), threads);
            lorikeet::pileups::pileup_variants(
                m,
                bam_readers,
                mode,
                &mut estimators.estimators,
                fasta_reader,
                print_zeros,
                flag_filters,
                mapq_threshold,
                var_fraction,
                min,
                max,
                contig_end_exclusion,
                &output_prefix,
                threads,
                method,
                coverage_fold,
                true,
                false);
        },
        _ => panic!("Unknown lorikeet mode"),
    }

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

        static ref EVOLVE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet evolve --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet evolve --method metabat --bam-files my.bam --reference assembly.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet evolve --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet evolve"),
            ansi_term::Colour::Green.paint(
                "Calculate dN/dS values of genes based on read mappings"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate gene dN/dS values from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate gene dN/dS values using MetaBAT adjusted coverage from a sorted BAM file, saving
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
                "Summarizes contigs stats including mean variant abundance, total variants, \
                and standard deviations"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate contig stats across samples"),
            ansi_term::Colour::Purple.paint(
                "Example: Summarizes contigs defined in reference from a sorted BAM file:"),
        ).to_string();

        static ref GENOTYPE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet genotype --method metabat --bam-files my.bam --reference assembly.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet genotype --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet genotype"),
            ansi_term::Colour::Green.paint(
                "*EXPERIMENTAL* Report strain-level genotypes and abundances based on variant read mappings"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate genotypes"),
            ansi_term::Colour::Purple.paint(
                "Example: Generate strain-level genotypes from read mappings compared to reference from a sorted BAM file:"),
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
        .author("Rhys J.P. Newell <r.newell near uq.edu.au>")
        .about("Variant analysis of metagenomic datasets")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .help("
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
\tgenotype \tReport strain-level genotypes and abundances from metagenomes (*experimental*)
\tpolymorph\tReport variant sites along contigs
\tsummarize\tSummarizes contig stats from one or multiple samples
\tevolve   \tCalculate dN/dS values for genes from read mappings

Less used utility subcommands:
\tkmer     \tCalculate kmer frequencies within contigs
\tfilter   \tRemove (or only keep) alignments with insufficient identity

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
                    .takes_value(true)
                    .required_unless_one(
                        &["read1","read2","coupled","interleaved","single","full-help"]))
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
                    .multiple(true)
                    .required_unless_one(
                        &["full-help"]))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
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
                    .takes_value(true)
                    .default_value("0.97"))
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
                .arg(Arg::with_name("output-prefix")
                    .long("output-prefix")
                    .short("o")
                    .default_value("output"))
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .multiple(false)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("epsilon")
                    .long("epsilon")
                    .short("e")
                    .default_value("0.05"))
                .arg(Arg::with_name("min-cluster-size")
                    .long("min-cluster-size")
                    .short("s")
                    .default_value("10"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("0.5"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("10"))
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
                .arg(Arg::with_name("nanopore")
                    .long("nanopore"))
                .arg(Arg::with_name("include-secondary")
                    .long("include-secondary"))
                .arg(Arg::with_name("include-supplementary")
                    .long("include-supplementary"))
                .arg(Arg::with_name("include-indels")
                    .long("include-indels"))
                .arg(Arg::with_name("verbose")
                    .short("v")
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .long("quiet")))
        .subcommand(
            SubCommand::with_name("evolve")
                .about("Pinpoint sites of synonymous and non-synonymous mutations")
                .help(EVOLVE_HELP.as_str())
                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["read1","read2","coupled","interleaved","single","full-help"]))
                .arg(Arg::with_name("gff")
                    .short("g")
                    .long("gff")
                    .takes_value(true))
                .arg(Arg::with_name("prodigal-params")
                    .long("prodigal-params")
                    .takes_value(true)
                    .conflicts_with("gff"))
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
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
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
                    .default_value("0.97")
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
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("0.5"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("10"))
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
                .arg(Arg::with_name("nanopore")
                    .long("nanopore"))
                .arg(Arg::with_name("include-secondary")
                    .long("include-secondary"))
                .arg(Arg::with_name("include-supplementary")
                    .long("include-supplementary"))
                .arg(Arg::with_name("include-indels")
                    .long("include-indels"))
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
            SubCommand::with_name("summarize")
                .about("Perform variant calling analysis and then binning")
                .help(SUMMARIZE_HELP.as_str())
                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["read1","read2","coupled","interleaved","single","full-help"]))
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
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
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
                    .default_value("0.97")
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
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("epsilon")
                    .long("epsilon")
                    .short("e")
                    .default_value("0.05"))
                .arg(Arg::with_name("min-cluster-size")
                    .long("min-cluster-size")
                    .short("s")
                    .default_value("10"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("0.5"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("10"))
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
                .arg(Arg::with_name("nanopore")
                    .long("nanopore"))
                .arg(Arg::with_name("include-secondary")
                    .long("include-secondary"))
                .arg(Arg::with_name("include-supplementary")
                    .long("include-supplementary"))
                .arg(Arg::with_name("include-indels")
                    .long("include-indels"))
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
            SubCommand::with_name("genotype")
                .about("Perform variant calling analysis and then binning")
                .help(GENOTYPE_HELP.as_str())
                .arg(Arg::with_name("full-help")
                    .long("full-help"))

                .arg(Arg::with_name("bam-files")
                    .short("b")
                    .long("bam-files")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["read1","read2","coupled","interleaved","single","full-help"]))
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
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
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
                    .default_value("0.97")
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
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("epsilon")
                    .long("epsilon")
                    .short("e")
                    .default_value("0.05"))
                .arg(Arg::with_name("min-cluster-size")
                    .long("min-cluster-size")
                    .short("s")
                    .default_value("10"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("0.5"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("depth-threshold")
                    .long("depth-threshold")
                    .short("d")
                    .default_value("50"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("10"))
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
                .arg(Arg::with_name("nanopore")
                    .long("nanopore"))
                .arg(Arg::with_name("include-secondary")
                    .long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping")
                    .long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary")
                    .long("include-supplementary"))
                .arg(Arg::with_name("include-indels")
                    .long("include-indels"))
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
                    .default_value("0.97")
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
                    .long("quiet")))
        .subcommand(
            SubCommand::with_name("polish")
                .about("Polish an assembly using highly abundant variant calls")
                .help(POLYMORPH_HELP.as_str())
                .arg(Arg::with_name("full-help")
                    .long("full-help"))
                .arg(Arg::with_name("bam-file")
                    .short("b")
                    .long("bam-file")
                    .multiple(false)
                    .takes_value(true)
                    .required_unless_one(
                        &["read1","read2","coupled","interleaved","single","full-help"]))
                .arg(Arg::with_name("sharded")
                    .long("sharded")
                    .required(false))
                .arg(Arg::with_name("read1")
                    .short("-1")
                    .multiple(false)
                    .takes_value(true)
                    .requires("read2")
                    .required_unless_one(
                        &["bam-file","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("read2")
                    .short("-2")
                    .multiple(false)
                    .takes_value(true)
                    .requires("read1")
                    .required_unless_one(
                        &["bam-file","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("coupled")
                    .short("-c")
                    .long("coupled")
                    .multiple(false)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-file","read1","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("interleaved")
                    .long("interleaved")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-file","read1","coupled","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("single")
                    .long("single")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-file","read1","coupled","interleaved","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("reference")
                    .short("-r")
                    .long("reference")
                    .takes_value(true)
                    .required_unless_one(
                        &["full-help"]))
                .arg(Arg::with_name("bam-file-cache-directory")
                    .long("bam-file-cache-directory")
                    .takes_value(true)
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
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
                    .takes_value(true)
                    .default_value("0.97"))
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
                .arg(Arg::with_name("output-prefix")
                    .long("output-prefix")
                    .short("o")
                    .takes_value(true)
                    .required(true))
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .multiple(false)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("0.5"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("10"))
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
                    .long("quiet")));
}
