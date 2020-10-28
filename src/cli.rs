use clap::*;

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-no-preset",
    "ngmlr",
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

const LONGREAD_MAPPING_SOFTWARE_LIST: &[&str] =
    &["minimap2-ont", "minimap2-pb", "ngmlr-ont", "ngmlr-pb"];
const DEFAULT_LONGREAD_MAPPING_SOFTWARE: &str = "minimap2-ont";

const MAPPER_HELP: &'static str =
    "    -p, --mapper <NAME>             Underlying mapping software used
                                         (\"minimap2-sr\", \"bwa-mem\",
                                         \"ngmlr-ont\", \"ngmlr-pb\", \"minimap2-ont\",
                                         \"minimap2-pb\", or \"minimap2-no-preset\").
                                         minimap2 -sr, -ont, -pb, -no-preset specify
                                         '-x' preset of minimap2 to be used
                                         (with map-ont, map-pb for -ont, -pb).
                                         [default: \"minimap2-sr\"] \n
         --minimap2-params PARAMS        Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"] \n
         --minimap2-reference-is-index   Treat reference as a minimap2 database, not
                                         as a FASTA file.\n
         --bwa-params PARAMS             Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]\n
         --ngmlr-params PARAMS           Extra parameters to provide to NGMLR.
                                         --bam-fix, -x ont, -t are already set. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.\n";

const VARIANT_CALLING_HELP: &'static str =
    "        --mapq-threshold <INT>                Mapping quality threshold used to verify
                                              a variant. [default: 10]\n
        -q, --base-quality-threshold <INT>    The minimum PHRED score for base in a read for it to be
                                              considered in the variant calling process.\n
        --fdr-threshold <FLOAT>               False discovery rate threshold for filtering variants
                                              based on the quality scores and accounting for the
                                              presence in all available samples.\n
        --ploidy <INT>                        Sets the default ploidy for the analysis to N.  [default: 1]\n
        --min-repeat-entropy <FLOAT>          To detect interrupted repeats, build across sequence until it has
                                              entropy > N bits per bp. Set to 0 to turn off. [default: 1.3]\n
        -o, --output-prefix <STRING>          Output prefix for files. [default: output]\n
        -f, --min-variant-depth <INT>         Minimum depth threshold value a variant must occur at
                                              for it to be considered. [default: 10]\n
        --min-variant-quality <INT>           Minimum QUAL value required for a variant to be included in
                                              analysis. [default: 10]\n
        --include-longread-svs                Include structural variants produced by SVIM in genotyping
                                              analysis. Can often overestimate number of variants present.\n
        --ulimit                              Sets the ulimit stack size to help prevent segmentation faults
                                              in freebayes recursive calls. Lower this on smaller systems.
                                              Increase this if your bam files contain thousands of contigs and
                                              freebayes is segfaulting. [default: 81920]\n
        --force                               Forcefully overwrite previous runs.\n";

const ALIGNMENT_OPTIONS: &'static str = "Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).
  -l, --longread-bam-files <PATH> ..     Path to BAM files(s) generated from longreads.
                                         Must be reference sorted.

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly.
                                         If multiple reference FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to reference
                                         separately as sharded BAMs
   -d, --genome-fasta-directory <PATH>   Directory containing FASTA files to be analyzed
   -x, --genome-fasta-extension <STR>    FASTA file extension in --genome-fasta-directory
                                         [default \"fna\"]
   -t, --threads <INT>                   Number of threads for mapping / sorting
                                         [default 1]
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
                                         [default 1]
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --longreads <PATH> ..                 pacbio or oxford nanopore long reads FASTA/Q files(s).
   --bam-file-cache-directory            Directory to store cached BAM files. BAM files are stored
                                         in /tmp by default.
   -d, --output-directory                Output directory";

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
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
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

pub fn polish_full_help() -> &'static str {
    lazy_static! {
        static ref POLISH_HELP: String = format!(
            "lorikeet polish: Produce consensus genomes for input genomes per sample

{}
{}
{}

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
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --plot                                Produce SNP density plots
   -w, --window-size <FLOAT>             Window size in kilobase pairs at which to calculate SNP and
                                         SV density.
   -t, --threads                         Number of threads used. [default: 1]
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
                                         [default 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage

   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>",
            ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP
        );
    }
    &POLISH_HELP
}

pub fn evolve_full_help() -> &'static str {
    lazy_static! {
        static ref EVOLVE_HELP: String = format!(
    "lorikeet evolve: Calculate dN/dS values in coding regions based on variants found in read mappings

{}
{}
{}

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
   mapq-threshold <INT>                  Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -q, --base-quality-threshold <INT>    The minimum PHRED score for base in a read for it to be
                                         considered in the variant calling process.
   --fdr-threshold <FLOAT>               False discovery rate threshold for filtering variants
                                         based on the quality scores and accounting for the
                                         presence in all available samples.
   --heterozygosity <FLOAT>              The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.01 implies on average a SNP should be detected
                                         every 1 in 100 bp [default: 0.01]
   --indel-heterozygosity <FLOAT>        The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.001 implies on average an INDEL should be detected
                                         every 1 in 1000 bp [default: 0.001]
   --ploidy <INT>                        Sets the default ploidy for the analysis to N.  (default: 1)
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth <INT>         Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --min-variant-quality <INT>           Minimum QUAL value required for a variant to be included in
                                         analysis. [default: 10]
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
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
                                         [default 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP);
    }
    &EVOLVE_HELP
}

pub fn summarize_full_help() -> &'static str {
    lazy_static! {
        static ref SUMMARIZE_HELP: String = format!(
    "lorikeet summarize: Provides per contig variant statistics for a metagenome

{}
{}
{}

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
   -w, --window-size <FLOAT>             Window size in kilobase pairs at which to calculate SNP and
                                         SV density.
   --include-longread-svs                Include structural variants produced by SVIM in genotyping
                                         analysis. Can often overestimate number of variants present.
   --mapq-threshold <INT>                Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -q, --base-quality-threshold <INT>    The minimum PHRED score for base in a read for it to be
                                         considered in the variant calling process.
   --fdr-threshold <FLOAT>               False discovery rate threshold for filtering variants
                                         based on the quality scores and accounting for the
                                         presence in all available samples.
   --heterozygosity <FLOAT>              The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.01 implies on average a SNP should be detected
                                         every 1 in 100 bp [default: 0.01]
   --indel-heterozygosity <FLOAT>        The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.001 implies on average an INDEL should be detected
                                         every 1 in 1000 bp [default: 0.001]
   --ploidy <INT>                        Sets the default ploidy for the analysis to N.  (default: 1)
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth <INT>         Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --min-variant-quality <INT>           Minimum QUAL value required for a variant to be included in
                                         analysis. [default: 10]
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --include-longread-svs                Flag indicating whether to use SVIM to calculate structural
                                         variants.
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   -t, --threads                         Number of threads used. [default: 1]
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
                                         [default 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP);
    }
    &SUMMARIZE_HELP
}

pub fn genotype_full_help() -> &'static str {
    lazy_static! {
        static ref GENOTYPE_HELP: String = format!(
    "lorikeet genotype: Resolves strain-level genotypes and abundance from metagenomes

{}
{}
{}

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
   --mapq-threshold <INT>                  Mapping quality threshold used to verify
                                         a variant. [default: 10]
   -q, --base-quality-threshold <INT>    The minimum PHRED score for base in a read for it to be
                                         considered in the variant calling process.
   --fdr-threshold <FLOAT>               False discovery rate threshold for filtering variants
                                         based on the quality scores and accounting for the
                                         presence in all available samples.
   --heterozygosity <FLOAT>              The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.01 implies on average a SNP should be detected
                                         every 1 in 100 bp [default: 0.01]
   --indel-heterozygosity <FLOAT>        The expected heterozygosity value used to compute prior
                                         probability that a locus is non-reference. A value of
                                         0.001 implies on average an INDEL should be detected
                                         every 1 in 1000 bp [default: 0.001]
   --ploidy <INT>                        Sets the default ploidy for the analysis to N.  (default: 1)
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]
   -f, --min-variant-depth <INT>         Minimum depth threshold value a variant must occur at
                                         for it to be considered. [default: 10]
   --min-variant-quality <INT>           Minimum QUAL value required for a variant to be included in
                                         analysis. [default: 10]
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
   -n, --minimum-seed-size <INT>         Minimum seed size for initiating fuzzy DBSCAN with. Corresponds
                                         to the number of variants within the seed. [default: 20]
   -s, --maximum-seed-similarity <FLOAT> The maximum Jaccard's similarity value allowed between two
                                         seeds before they get concatenated into one seed. Prevents
                                         overclustering. [default: 0.97]
   --minimum-reads-in-link <INT>         Minimum amount of reads required to be shared between two
                                         variants before they are counted as 'linked'. [default: 5]
   --include-longread-svs                Include structural variants produced by SVIM in genotyping
                                         analysis. Can often overestimate number of variants present.
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --plot                                Produce SNP density plots
   -w, --window-size <FLOAT>             Window size in kilobase pairs at which to calculate SNP and
                                         SV density.
   -t, --threads                         Number of threads used. [default: 1]
   --parallel-genomes                    Number of genomes to run in parallel.
                                         Increases memory usage linearly.
                                         [default 1]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage

   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Rhys J. P. Newell <r.newell near uq.edu.au>", ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP);
    }
    &GENOTYPE_HELP
}

pub fn build_cli() -> App<'static, 'static> {
    // specify _2 lazily because need to define it at runtime.
    lazy_static! {
        static ref POLISH_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet polish --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet polish --bam-files my.bam --reference genome.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet polish --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet polish"),
            ansi_term::Colour::Green.paint(
                "Generate consensus genomes from multiple samples"),
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

  lorikeet evolve --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10

See lorikeet evolve --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet evolve"),
            ansi_term::Colour::Green.paint(
                "Calculate dN/dS values of genes based on read mappings"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate gene dN/dS values from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate gene dN/dS values from reads against many genomes in a directory and cache bams and save output to directory:")
        ).to_string();


        static ref SUMMARIZE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet summarize --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --window-size 1

{}

  lorikeet summarize --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10

See lorikeet summarize --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet summarize"),
            ansi_term::Colour::Green.paint(
                "Summarizes contigs stats across given window size"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate contig stats across samples using a 1kb window size"),
            ansi_term::Colour::Purple.paint(
                "Example: Summarizes genomic variation across contigs in genomes in directory using long and short reads:"),
        ).to_string();

        static ref GENOTYPE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10 --plot

See lorikeet genotype --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet genotype"),
            ansi_term::Colour::Green.paint(
                "*EXPERIMENTAL* Report strain-level genotypes and abundances based on variant read mappings"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate genotypes"),
            ansi_term::Colour::Purple.paint(
                "Example: Generate strain-level genotypes from read mappings compared to reference from a sorted BAM file and plots the results:"),
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
        .args_from_usage(
            "-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'",
        )
        .help(
            "
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
\tgenotype \tReport strain-level genotypes and abundances from metagenomes (*experimental*)
\tsummarize\tSummarizes contig stats from one or multiple samples
\tevolve   \tCalculate dN/dS values for genes from read mappings

Less used utility subcommands:
\tkmer     \tCalculate kmer frequencies within contigs
\tfilter   \tRemove (or only keep) alignments with insufficient identity

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <r.newell near uq.edu.au>
",
        )
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("evolve")
                .about("Pinpoint sites of synonymous and non-synonymous mutations")
                .help(EVOLVE_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(Arg::with_name("gff").long("gff").takes_value(true))
                .arg(
                    Arg::with_name("prokka-params")
                        .long("prokka-params")
                        .takes_value(true)
                        .conflicts_with("gff"),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short("r")
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short("d")
                        .takes_value(true)
                        .required_unless_one(&["reference", "genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short("x")
                        .takes_value(true)
                        .default_value("fna"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .long("parallel-genomes")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
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
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("coverage-fold")
                        .long("coverage-fold")
                        .default_value("0.5"),
                )
                .arg(
                    Arg::with_name("min-variant-depth")
                        .long("min-variant-depth")
                        .short("f")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-variant-quality")
                        .long("min-variant-quality")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("fdr-threshold")
                        .long("fdr-threshold")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("heterozygosity")
                        .long("heterozygosity")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("kmer-size")
                        .long("kmer-size")
                        .short("k")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(
                    Arg::with_name("window-size")
                        .long("window-size")
                        .short("w")
                        .default_value("1"),
                )
                .arg(Arg::with_name("plot").long("plot"))
                .arg(Arg::with_name("nanopore").long("nanopore"))
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping").long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(Arg::with_name("include-indels").long("include-indels"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("min-repeat-entropy")
                        .long("min-repeat-entropy")
                        .default_value("1.3")
                        .required(false),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet"))
                .arg(
                    Arg::with_name("ulimit")
                        .long("ulimit")
                        .default_value("81920"),
                ),
        )
        .subcommand(
            SubCommand::with_name("summarize")
                .about("Perform variant calling analysis and then binning")
                .help(SUMMARIZE_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-file")
                        .short("r")
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short("d")
                        .takes_value(true)
                        .required_unless_one(&["reference", "genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short("x")
                        .takes_value(true)
                        .default_value("fna"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .long("parallel-genomes")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
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
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(
                    Arg::with_name("window-size")
                        .long("window-size")
                        .short("w")
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("epsilon")
                        .long("epsilon")
                        .short("e")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("min-cluster-size")
                        .long("min-cluster-size")
                        .short("s")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("coverage-fold")
                        .long("coverage-fold")
                        .default_value("0.5"),
                )
                .arg(
                    Arg::with_name("min-variant-depth")
                        .long("min-variant-depth")
                        .short("f")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-variant-quality")
                        .long("min-variant-quality")
                        .default_value("10"),
                )
                .arg(Arg::with_name("strain-ani").long("strain-ani").short("a"))
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("fdr-threshold")
                        .long("fdr-threshold")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("heterozygosity")
                        .long("heterozygosity")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("kmer-size")
                        .long("kmer-size")
                        .short("k")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("min-repeat-entropy")
                        .long("min-repeat-entropy")
                        .default_value("1.3")
                        .required(false),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("nanopore").long("nanopore"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet"))
                .arg(
                    Arg::with_name("ulimit")
                        .long("ulimit")
                        .default_value("81920"),
                ),
        )
        .subcommand(
            SubCommand::with_name("genotype")
                .about("Perform variant calling analysis and then binning")
                .help(GENOTYPE_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short("r")
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short("d")
                        .takes_value(true)
                        .required_unless_one(&["reference", "genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short("x")
                        .takes_value(true)
                        .default_value("fna"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("vcfs")
                        .long("vcfs")
                        .multiple(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .long("parallel-genomes")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
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
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(Arg::with_name("e-min").long("e-min").default_value("0.1"))
                .arg(Arg::with_name("e-max").long("e-max").default_value("0.25"))
                .arg(
                    Arg::with_name("pts-min")
                        .long("pts-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("pts-max")
                        .long("pts-max")
                        .default_value("0.1"),
                )
                .arg(Arg::with_name("phi").long("phi").default_value("0.0"))
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("coverage-fold")
                        .long("coverage-fold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("min-variant-depth")
                        .long("min-variant-depth")
                        .short("f")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-variant-quality")
                        .long("min-variant-quality")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("strain-ani")
                        .long("strain-ani")
                        .short("a")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("minimum-seed-size")
                        .long("minimum-seed-size")
                        .short("n")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("minimum-reads-in-link")
                        .long("minimum-reads-in-link")
                        .default_value("5"),
                )
                .arg(
                    Arg::with_name("maximum-seed-similarity")
                        .long("maximum-seed-similarity")
                        .short("s")
                        .default_value("0.97"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("fdr-threshold")
                        .long("fdr-threshold")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("heterozygosity")
                        .long("heterozygosity")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("kmer-size")
                        .long("kmer-size")
                        .short("k")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(
                    Arg::with_name("window-size")
                        .long("window-size")
                        .short("w")
                        .default_value("1"),
                )
                .arg(Arg::with_name("plot").long("plot"))
                .arg(Arg::with_name("nanopore").long("nanopore"))
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping").long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("min-repeat-entropy")
                        .long("min-repeat-entropy")
                        .default_value("1.3")
                        .required(false),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet"))
                .arg(
                    Arg::with_name("ulimit")
                        .long("ulimit")
                        .default_value("81920"),
                ),
        )
        .subcommand(
            SubCommand::with_name("kmer")
                .about("Generate kmer count matrix for contigs")
                //                .help(CONTIG_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .takes_value(true)
                        .required(true),
                )
                .arg(Arg::with_name("verbose").short("v").long("verbose")),
        )
        .subcommand(
            SubCommand::with_name("filter")
                .about("Remove alignments with insufficient identity")
                .help(FILTER_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help"]),
                )
                .arg(
                    Arg::with_name("output-bam-files")
                        .short("o")
                        .long("output-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help"]),
                )
                .arg(Arg::with_name("inverse").long("inverse"))
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(
                    Arg::with_name("threads")
                        .long("threads")
                        .short("t")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .long("parallel-genomes")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(
                    Arg::with_name("verbose")
                        // .short("v") // Do not use since could be confused with
                        // inverse (a la grep -v)
                        .long("verbose"),
                )
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("polish")
                .about("Polish an assembly using highly abundant variant calls")
                .help(POLISH_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short("r")
                        .long("reference")
                        .alias("reference")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short("d")
                        .takes_value(true)
                        .required_unless_one(&["reference", "genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short("x")
                        .takes_value(true)
                        .default_value("fna"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .long("parallel-genomes")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
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
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .takes_value(true)
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .multiple(false)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("coverage-fold")
                        .long("coverage-fold")
                        .default_value("0.5"),
                )
                .arg(
                    Arg::with_name("min-variant-depth")
                        .long("min-variant-depth")
                        .short("f")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-variant-quality")
                        .long("min-variant-quality")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("fdr-threshold")
                        .long("fdr-threshold")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("heterozygosity")
                        .long("heterozygosity")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("plot").long("plot"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("min-repeat-entropy")
                        .long("min-repeat-entropy")
                        .default_value("1.3")
                        .required(false),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(
                    Arg::with_name("window-size")
                        .long("window-size")
                        .short("w")
                        .default_value("1"),
                )
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet"))
                .arg(
                    Arg::with_name("ulimit")
                        .long("ulimit")
                        .default_value("81920"),
                ),
        );
}
