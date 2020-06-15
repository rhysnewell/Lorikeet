use clap::*;


const MAPPING_SOFTWARE_LIST: &[&str] = &["bwa-mem", "minimap2-sr", "minimap2-ont", "minimap2-pb","minimap2-no-preset"];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";


const MAPPER_HELP: &'static str =
    "   -p, --mapper <NAME>                   Underlying mapping software used
                                         (\"minimap2-sr\", \"bwa-mem\", \"minimap2-ont\",
                                         \"minimap2-pb\", or \"minimap2-no-preset\").
                                         minimap2 -sr, -ont, -pb, -no-preset specify
                                         '-x' preset of minimap2 to be used
                                         (with map-ont, map-pb for -ont, -pb).
                                         [default: \"minimap2-sr\"]";

pub fn filter_full_help() -> &'static str {
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

pub fn polymorph_full_help() -> &'static str {
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
   -d, --outdir                          Output directory
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
                                              Conflicts --proper-pairs-only. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                              identity e.g. 0.95 for 95%.
                                              Conflicts --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                              bases e.g. 0.95 means 95% of the read's
                                              bases must be aligned.
                                              Conflicts --proper-pairs-only. [default 0.0]
   --proper-pairs-only                     Allows reads to be mapped as improper pairs
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
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &POLYMORPH_HELP
}

pub fn evolve_full_help() -> &'static str {
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
   -d, --outdir                          Output directory.
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
                                              Conflicts --proper-pairs-only. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                              identity e.g. 0.95 for 95%.
                                              Conflicts --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                              bases e.g. 0.95 means 95% of the read's
                                              bases must be aligned.
                                              Conflicts --proper-pairs-only. [default 0.0]
   --proper-pairs-only                     Allows reads to be mapped as improper pairs
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
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &EVOLVE_HELP
}


pub fn summarize_full_help() -> &'static str {
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
   -d, --outdir                          Output directory
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
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --proper-pairs-only                Allows reads to be mapped as improper pairs
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
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &SUMMARIZE_HELP
}

pub fn genotype_full_help() -> &'static str {
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
   -d, --outdir                          Output directory
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
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --proper-pairs-only                Allows reads to be mapped as improper pairs
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
   -k, --kmer-size <INT>                 K-mer size used to generate k-mer frequency
                                         table. [default: 4]
   --include-indels                      Flag indicating whether to attempt to calculate INDEL sites
                                         Not recommended if using nanopore long read data.
   -q, mapq-threshold <INT>              Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth      Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --e-min                               Minimum epsilon value used in fuzzyDBSCAN algorithm.
                                         The minimum distance between two points required for clustering.
   --e-max                               Maximum epsilon value used in fuzzyDBSCAN algorithm.
                                         The maximum distance between two points for border clustering.
   --pts-min                             Minimum points as percentage in fuzzyDBSCAN algorithm.
                                         The fraction of points needed to be within e-max
                                         to begin core clustering.
   --pts-max                             Maximum points as percentage in fuzzyDBSCAN algorithm.
                                         The fraction of points needed to be within e-max
                                         to begin border clustering.
   --include-longread-svs                Include structural variants produced by SVIM in genotyping
                                         analysis. Can often overestimate number of variants present.
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

   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", MAPPER_HELP);
    }
    &GENOTYPE_HELP
}

pub fn build_cli() -> App<'static, 'static> {
    // specify _2 lazily because need to define it at runtime.
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
                .arg(Arg::with_name("outdir")
                    .long("bam-file-cache-directory")
                    .short("d")
                    .takes_value(true))
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
                    .default_value("0.0"))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
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
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("0"))
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
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
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
                .arg(Arg::with_name("outdir")
                    .long("bam-file-cache-directory")
                    .short("d")
                    .takes_value(true))
                .arg(Arg::with_name("threads")
                    .short("t")
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
                    .default_value("0.0")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
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
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("0"))
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
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
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
                .arg(Arg::with_name("longreads")
                    .long("longreads")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longread-bam-files", "sharded"]))
                .arg(Arg::with_name("longread-bam-files")
                    .short("l")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longreads", "sharded"]))
                .arg(Arg::with_name("reference")
                    .short("-r")
                    .long("reference")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("outdir")
                    .long("bam-file-cache-directory")
                    .short("d")
                    .takes_value(true))
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
                    .default_value("0.0")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
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
                .arg(Arg::with_name("strain-ani")
                    .long("strain-ani")
                    .short("a"))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .short("q")
                    .default_value("0"))
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
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
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
                .arg(Arg::with_name("longreads")
                    .long("longreads")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longread-bam-files", "sharded"]))
                .arg(Arg::with_name("longread-bam-files")
                    .short("l")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longreads", "sharded"]))
                .arg(Arg::with_name("reference")
                    .short("-r")
                    .long("reference")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("outdir")
                    .long("bam-file-cache-directory")
                    .short("d")
                    .takes_value(true))
                .arg(Arg::with_name("vcfs")
                    .long("vcfs")
                    .multiple(true)
                    .required(false))
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
                    .default_value("0.0")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("method")
                    .short("m")
                    .long("method")
                    .takes_value(true)
                    .possible_values(&[
                        "trimmed_mean",
                        "mean",
                        "metabat"])
                    .default_value("trimmed_mean"))
                .arg(Arg::with_name("e-min")
                    .long("e-min")
                    .default_value("0.05"))
                .arg(Arg::with_name("e-max")
                    .long("e-max")
                    .default_value("0.15"))
                .arg(Arg::with_name("pts-min")
                    .long("pts-min")
                    .default_value("10"))
                .arg(Arg::with_name("pts-max")
                    .long("pts-max")
                    .default_value("25"))
                .arg(Arg::with_name("phi")
                    .long("phi")
                    .default_value("0.0"))
                .arg(Arg::with_name("min-covered-fraction")
                    .long("min-covered-fraction")
                    .default_value("0.0"))
                .arg(Arg::with_name("coverage-fold")
                    .long("coverage-fold")
                    .default_value("1.0"))
                .arg(Arg::with_name("min-variant-depth")
                    .long("min-variant-depth")
                    .short("f")
                    .default_value("10"))
                .arg(Arg::with_name("strain-ani")
                    .long("strain-ani")
                    .short("a")
                    .takes_value(true))
                .arg(Arg::with_name("mapq-threshold")
                    .long("mapq-threshold")
                    .default_value("10"))
                .arg(Arg::with_name("base-quality-threshold")
                    .long("base-quality-threshold")
                    .short("q")
                    .default_value("13"))
                .arg(Arg::with_name("min-repeat-entropy")
                    .long("min-repeat-entropy")
                    .default_value("1.5"))
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
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
                .arg(Arg::with_name("nanopore")
                    .long("nanopore"))
                .arg(Arg::with_name("include-longread-svs")
                    .long("include-longread-svs"))
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
                .arg(Arg::with_name("ploidy")
                    .long("ploidy")
                    .default_value("2")
                    .required(false))

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
                    .default_value("0.0")
                    .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
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
                        &["bam-file","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("read2")
                    .short("-2")
                    .multiple(true)
                    .takes_value(true)
                    .requires("read1")
                    .required_unless_one(
                        &["bam-file","coupled","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("coupled")
                    .short("-c")
                    .long("coupled")
                    .multiple(true)
                    .takes_value(true)
                    .required_unless_one(
                        &["bam-file","read1","interleaved","single","full-help"])
                    .conflicts_with("bam-file"))
                .arg(Arg::with_name("longreads")
                    .long("longreads")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longread-bam-files", "sharded"]))
                .arg(Arg::with_name("longread-bam-files")
                    .short("l")
                    .multiple(true)
                    .takes_value(true)
                    .required(false)
                    .conflicts_with_all(&["longreads", "sharded"]))
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
                .arg(Arg::with_name("outdir")
                    .long("bam-file-cache-directory")
                    .short("d")
                    .takes_value(true))
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
                    .default_value("0.0"))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                    .long("min-read-aligned-length-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                    .long("min-read-percent-identity-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                    .long("min-read-aligned-percent-pair")
                    .takes_value(true)
                    .conflicts_with("proper-pairs-only"))
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
                    .default_value("0"))
                .arg(Arg::with_name("include-indels")
                    .long("include-indels"))
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
                .arg(Arg::with_name("proper-pairs-only")
                    .long("proper-pairs-only"))
                .arg(Arg::with_name("verbose")
                    .short("v")
                    .long("verbose"))
                .arg(Arg::with_name("quiet")
                    .long("quiet")));
}
