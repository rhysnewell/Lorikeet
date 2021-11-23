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

const MAPPER_HELP: &str = "
Read mapping options:
  --mapper <NAME>                 Underlying mapping software used for short reads
                                  (\"minimap2-sr\", \"bwa-mem\",
                                  \"ngmlr-ont\", \"ngmlr-pb\", \"minimap2-ont\",
                                  \"minimap2-pb\", or \"minimap2-no-preset\").
                                  minimap2 -sr, -ont, -pb, -no-preset specify
                                  '-x' preset of minimap2 to be used
                                  (with map-ont, map-pb for -ont, -pb).
                                  [default: \"minimap2-sr\"] \n
  --longread-mapper <NAME>        Underlying mapping software used for long reads
                                  (\"minimap2-sr\", \"bwa-mem\",
                                  \"ngmlr-ont\", \"ngmlr-pb\", \"minimap2-ont\",
                                  \"minimap2-pb\", or \"minimap2-no-preset\").
                                  minimap2 -sr, -ont, -pb, -no-preset specify
                                  '-x' preset of minimap2 to be used
                                  (with map-ont, map-pb for -ont, -pb).
                                  [default: \"minimap2-ont\"] \n
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

const VARIANT_CALLING_HELP: &str = "
Variant calling options (Basic):
  -k, --kmer-sizes <INT>                        K-mer sizes used to generate DeBruijn Graphs.
                                                Multiple values at once are accepted and encouraged
                                                e.g. 10 25 [default: 25] \n
  --ploidy <INT>                                Sets the default ploidy for the analysis to N.
                                                [default: 1]\n
  -f, --features-vcf                            The set of alleles to force-call regardless
                                                of evidence. Note: The sight containing these alleles
                                                has to be called as 'active' in order for them to appear
                                                in the final VCF. Addtionally, Provided file must be
                                                compressed using bgzip and indexed using bcftools index. If no index
                                                is present, and index will be attempted to be created.
                                                If the file is not properly compressed, Lorikeet will
                                                segfault with no error message.
  -q, --min-base-quality                        Minimum base quality required to consider a
                                                base for calling. [default: 10]
  --base-quality-score-threshold                Base qualities below this threshold will
                                                be reduced to the minimum (6). [default: 18]
  --max-input-depth                             The maximum number of reads included within an
                                                assembly region across all samples. Larger numbers
                                                increase run time. If the depth of an assembly region
                                                exceeds this value, then the reads will be filtered
                                                by mean base quality. [default: 200000]
  --min-contig-size                             The minimum contig size to call variants on. Smaller
                                                contigs can often contain highly variable regions that
                                                mostly represent noise. Call variants on them can often
                                                be slow and not produce anything fruitful. If you
                                                wish to call variants on all available contigs,
                                                then set this to 0. [default: 2500]

Variant calling options (Advanced):
  --phred-scaled-global-read-mismapping-rate    The global assumed mismapping rate for reads. [default: 45]
  --pair-hmm-gap-continuation-penalty           Flat gap continuation penalty for use in the Pair HMM. [default: 10]
  --pcr-indel-model                             The PCR indel model to use. [default: conservative]
  --heterozygosity                              Heterozygosity value used to compute prior
                                                likelihoods for any locus. [default: 0.001]
  --heterozygosity-stdev                        Standard deviation of heterozygosity for SNP and
                                                indel calling. [default: 0.01]
  --indel-heterozygosity                        Heterozygosity for indel calling. [default: 0.000125]
  -C, --standard-min-confidence-threshold-for-calling
                                                The minimum phred-scaled confidence threshold at
                                                which variants should be called. [default: 30.0]
  --use-posteriors-to-calculate-qual            if available, use the genotype posterior
                                                probabilities to calculate the site QUAL.
  --annotate-with-num-discovered-alleles        If provided, we will annotate records with the
                                                number of alternate alleles that were discovered
                                                (but not necessarily genotyped) at a given site.
  --active-probability-threshold                Minimum probability for a locus to be
                                                considered active. [default: 0.002]
  --min-assembly-region-size                    Minimum size of an assembly region. [default: 50]
  --max-assembly-region-size                    Maximum size of an assembly region. [default: 300]
  --assembly-region-padding                     Number of additional bases of context to
                                                include around each assembly region. [default: 100]
  --dont-increase-kmer-sizes-for-cycles         Disable iterating over kmer sizes when
                                                graph cycles are detected.
  --allow-non-unique-kmers-in-ref               Allow graphs that have non-unique kmers in the reference.
  --do-not-run-physical-phasing                 Disable physical phasing.
  --recover-all-dangling-branches               Recover all dangling branches.
  --min-dangling-branch-length                  Minimum length of a dangling branch to
                                                attempt recovery. [default: 4]
  --graph-output                                Write debug assembly graph information to this file.
  --num-pruning-samples                         Number of samples that must pass the
                                                min_pruning threshold [default: 1]
  --dont-use-soft-clipped-bases                 Do not analyze soft clipped bases in the reads.
  --initial-error-rate-for-pruning              Initial base error rate estimate for adaptive
                                                pruning. [default: 0.001]
  --pruning-log-odds-threshold                  Likelihood ratio threshold for adaptive
                                                pruning algorithm. This value will be converted to
                                                log odds value. [default: 1.0]
  --max-unpruned-variants                       Maximum number of variants in graph the
                                                adaptive pruner will allow. [default: 100]
  --max-prob-propagation-distance               Upper limit on how many bases away probability mass
                                                can be moved around when calculating the boundaries
                                                between active and inactive assembly regions. [default: 50]
  --max-mnp-distance                            Two or more phased substitutions separated by
                                                this distance or less are merged into MNPs. [default: 0]
  --disable-optimizations                       Don't skip calculations in ActiveRegions with no variants
  --disable-avx                                 Disable the use of the GKL-rs AVX acceleration components
                                                for PairHMM and Smith-Waterman calculations.
  --limiting-interval                           Mainly used for debugging purposes. Only call variants
                                                within this given span on all contigs. E.g. providing
                                                '1000-2000' would only call variants between the 1000
                                                and 2000 bp span on each provided contig.
  --force                                       Forcefully overwrite previous runs.\n";

const ALIGNMENT_OPTIONS: &str = "
Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).
   -l, --longread-bam-files <PATH> ..    Path to BAM files(s) generated from longreads.
                                         Must be reference sorted.
   --query-assembly-bam-files <PATH>     The results of mapping one or more metagenome assemblies
                                         back on to your MAGs.

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
   -d, --output-directory                Output directory
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]\n


Sharding i.e. multiple reference sets (optional):
  --sharded                              If -b/--bam-files was used:
                                           Input BAM files are read-sorted alignments
                                           of a set of reads mapped to multiple
                                           reference contig sets. Choose the best
                                           hit for each read pair.

                                         Otherwise if mapping was carried out:
                                           Map reads to each reference, choosing the
                                             best hit for each pair.

Alignment filtering (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                           metabat (\"MetaBAT adjusted coverage\")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
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
   --allow-improper-pairs                Allows reads to be mapped as improper pairs
   --include-supplementary               Includes read alignments flagged as supplementary
   --include-secondary                   Includes read alignments flagged as secondary
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
";

const GENERAL_HELP: &str = "
General options:
  -t, --threads                         Maximum number of threads used. [default: 8]
  -p, --parallel-genomes                Number of genomes to run in parallel.
                                        Increases memory usage linearly.
                                        Thread usage qill not exceed the value
                                        provided by --threads [default 4]
  -v, --verbose                         Print extra debugging information
  -q, --quiet                           Unless there is an error, do not print
                                        log messages

Author: Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>
";

pub fn genotype_full_help() -> String {
    format!(
        "lorikeet genotype: Resolves strain-level genotypes and abundance from metagenomes

{}
{}
{}

Genotyping arguments (optional):

  --qual-by-depth-filter                The minimum QD value for a variant to have for it to be
                                        included in the genotyping analysis. [default: 20]
  --min-variant-depth-for-genotyping    The minimum total depth of a variant - across all samples -
                                        for it to be included in the strain genotyping process.
                                        Lower values tend to confuse and break the UMAP embedding
                                        and strain abundance calculation. [default: 5]
{}
        ",
        ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP, GENERAL_HELP
    )
}

pub fn call_full_help() -> String {
    format!(
        "lorikeet call: Call variants using local reassembly across multiple genomes and samples

        {}
        {}
        {}
        {}",
        ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP, GENERAL_HELP
    )
}

pub fn consensus_full_help() -> String {
    format!(
        "lorikeet consensus: Generate consensus genomes for each provided sample and genome

        {}
        {}
        {}
        {}",
        ALIGNMENT_OPTIONS, MAPPER_HELP, VARIANT_CALLING_HELP, GENERAL_HELP
    )
}

pub fn build_cli() -> App<'static, 'static> {
    // specify _2 lazily because need to define it at runtime.
    lazy_static! {
        static ref CONSENSUS_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet consensus --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  lorikeet consensus --bam-files my.bam --reference genome.fna
    --bam-file-cache-directory saved_bam_files --threads 10

See lorikeet consensus --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet consensus"),
            ansi_term::Colour::Green.paint(
                "Generate consensus genomes from multiple samples"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate variant positions from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate variant positions using MetaBAT adjusted coverage from a sorted BAM file, saving
    the unfiltered BAM files in the saved_bam_files folder:")
        );

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
        );


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
        );

        static ref GENOTYPE_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25

{}

  lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10

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
        );

        static ref CALL_HELP: String = format!(
            "
                            {}
              {}

{}

  lorikeet call --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25

{}

  lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10 --kmer-sizes 10 25

See lorikeet genotype --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet call"),
            ansi_term::Colour::Green.paint(
                "Perform read mapping and variant calling using local reassembly of active regions"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference and generate genotypes"),
            ansi_term::Colour::Purple.paint(
                "Example: Perform read read mapping and variant calling on an entire directory of genomes and save the bam files:"),
        );
    }

    return App::new("lorikeet")
        .version(crate_version!())
        .author("Rhys J.P. Newell <rhys.newell near hdr.qut.edu.au>")
        .about("Variant analysis of metagenomic datasets")
        .args_from_usage(
            "-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'",
        )
        .help(
            "
Variant calling and strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
\tgenotype  \tReport strain-level genotypes and abundances from metagenomes (*experimental*)
\tconsensus \tCreates consensus genomes for each input reference and for each sample
\tcall      \tPerforms variant calling on the provides genomes
\tevolve    \tCalculate dN/dS values for genes from read mappings

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>
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
                .arg(
                    Arg::with_name("assembly-bam-files")
                        .long("query-assembly-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["assembly"]),
                )
                .arg(Arg::with_name("gff").long("gff").takes_value(true))
                .arg(
                    Arg::with_name("prodigal-params")
                        .long("prodigal-params")
                        .takes_value(true)
                        .default_value("-p meta")
                        .conflicts_with("gff"),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("1")
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
                        .short("2")
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
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short("p")
                        .long("parallel-genomes")
                        .default_value("4")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
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
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
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
                    Arg::with_name("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .default_value("45"),
                )
                .arg(
                    Arg::with_name("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .possible_values(&[
                            "none",
                            "None",
                            "NONE",
                            "hostile",
                            "Hostile",
                            "HOSTILE",
                            "aggresive",
                            "Agressive",
                            "AGGRESSIVE",
                            "conservative",
                            "Conservative",
                            "CONSERVATIVE",
                        ]),
                )
                .arg(
                    Arg::with_name("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::with_name("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short("C")
                        .default_value("30.0"),
                )
                .arg(
                    Arg::with_name("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .possible_values(&[
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual"),
                )
                .arg(
                    Arg::with_name("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles"),
                )
                .arg(
                    Arg::with_name("active-probability-threshold")
                        .long("active-probability-threshold")
                        .default_value("0.002"),
                )
                .arg(
                    Arg::with_name("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .default_value("300"),
                )
                .arg(
                    Arg::with_name("kmer-sizes")
                        .long("kmer-sizes")
                        .short("k")
                        .multiple(true)
                        .default_value("25"), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("256")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles"),
                )
                .arg(
                    Arg::with_name("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref"),
                )
                .arg(
                    Arg::with_name("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing"),
                )
                .arg(
                    Arg::with_name("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches"),
                )
                .arg(
                    Arg::with_name("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::with_name("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("num-pruning-samples")
                        .long("num-pruning-samples")
                        .multiple(true)
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("min-prune-factor")
                        .long("min-prune-factor")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-adaptive-pruning")
                        .long("dont-use-adaptive-pruning")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases"),
                )
                .arg(
                    Arg::with_name("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .default_value("4.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("max-input-depth")
                        .long("max-input-depth")
                        .short("i")
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("error-correct-reads")
                        .long("error-correct-reads")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-mnp-distance")
                        .long("max-mnp-distance")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .default_value("-1")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("assembly-region-padding")
                        .long("assembly-region-padding")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .hidden(true),
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
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short("q")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .default_value("18"),
                )
                .arg(
                    Arg::with_name("enable-dynamic-read-disqualification-for-genotyping")
                        .long("enable-dynamic-read-disqualification-for-genotyping")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .default_value("1.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .default_value("0.02")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .hidden(true),
                )
                .arg(Arg::with_name("disable-optimizations").long("disable-optimizations"))
                .arg(Arg::with_name("disable-avx").long("disable-avx"))
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs").long("allow-improper-pairs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
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
                        .short("1")
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
                        .short("2")
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
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short("f")
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short("p")
                        .long("parallel-genomes")
                        .default_value("4")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
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
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
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
                    Arg::with_name("pts-min")
                        .long("pts-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("pts-max")
                        .long("pts-max")
                        .default_value("0.1"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("min-contig-size")
                        .long("min-contig-size")
                        .default_value("2500"),
                )
                .arg(
                    Arg::with_name("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .default_value("45"),
                )
                .arg(
                    Arg::with_name("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .possible_values(&[
                            "none",
                            "None",
                            "NONE",
                            "hostile",
                            "Hostile",
                            "HOSTILE",
                            "aggresive",
                            "Agressive",
                            "AGGRESSIVE",
                            "conservative",
                            "Conservative",
                            "CONSERVATIVE",
                        ]),
                )
                .arg(
                    Arg::with_name("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::with_name("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short("C")
                        .default_value("30.0"),
                )
                .arg(
                    Arg::with_name("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .possible_values(&[
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual"),
                )
                .arg(
                    Arg::with_name("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles"),
                )
                .arg(
                    Arg::with_name("active-probability-threshold")
                        .long("active-probability-threshold")
                        .default_value("0.002"),
                )
                .arg(
                    Arg::with_name("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .default_value("300"),
                )
                .arg(
                    Arg::with_name("kmer-sizes")
                        .long("kmer-sizes")
                        .short("k")
                        .multiple(true)
                        .default_value("25"), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("256")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles"),
                )
                .arg(
                    Arg::with_name("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref"),
                )
                .arg(
                    Arg::with_name("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing"),
                )
                .arg(
                    Arg::with_name("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches"),
                )
                .arg(
                    Arg::with_name("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::with_name("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("num-pruning-samples")
                        .long("num-pruning-samples")
                        .multiple(true)
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("min-prune-factor")
                        .long("min-prune-factor")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-adaptive-pruning")
                        .long("dont-use-adaptive-pruning")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases"),
                )
                .arg(
                    Arg::with_name("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .default_value("4.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("max-input-depth")
                        .long("max-input-depth")
                        .short("i")
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("error-correct-reads")
                        .long("error-correct-reads")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-mnp-distance")
                        .long("max-mnp-distance")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .default_value("-1")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("assembly-region-padding")
                        .long("assembly-region-padding")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .hidden(true),
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
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short("q")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .default_value("18"),
                )
                .arg(
                    Arg::with_name("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-variant-depth-for-genotyping")
                        .long("min-variant-depth-for-genotyping")
                        .default_value("5")
                )
                .arg(
                    Arg::with_name("enable-dynamic-read-disqualification-for-genotyping")
                        .long("enable-dynamic-read-disqualification-for-genotyping")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .default_value("1.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .default_value("0.02")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .hidden(true),
                )
                .arg(Arg::with_name("disable-optimizations").long("disable-optimizations"))
                .arg(Arg::with_name("disable-avx").long("disable-avx"))
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs").long("allow-improper-pairs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("call")
                .about("Perform variant calling across the given genomes and samples")
                .help(CALL_HELP.as_str())
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
                        .short("1")
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
                        .short("2")
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
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short("f")
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short("p")
                        .long("parallel-genomes")
                        .default_value("4")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
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
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("min-contig-size")
                        .long("min-contig-size")
                        .default_value("2500"),
                )
                .arg(
                    Arg::with_name("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .default_value("45"),
                )
                .arg(
                    Arg::with_name("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .possible_values(&[
                            "none",
                            "None",
                            "NONE",
                            "hostile",
                            "Hostile",
                            "HOSTILE",
                            "aggresive",
                            "Agressive",
                            "AGGRESSIVE",
                            "conservative",
                            "Conservative",
                            "CONSERVATIVE",
                        ]),
                )
                .arg(
                    Arg::with_name("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::with_name("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short("C")
                        .default_value("30.0"),
                )
                .arg(
                    Arg::with_name("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .possible_values(&[
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual"),
                )
                .arg(
                    Arg::with_name("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles"),
                )
                .arg(
                    Arg::with_name("active-probability-threshold")
                        .long("active-probability-threshold")
                        .default_value("0.002"),
                )
                .arg(
                    Arg::with_name("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .default_value("300"),
                )
                .arg(
                    Arg::with_name("kmer-sizes")
                        .long("kmer-sizes")
                        .short("k")
                        .multiple(true)
                        .default_value("25"), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("256")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles"),
                )
                .arg(
                    Arg::with_name("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref"),
                )
                .arg(
                    Arg::with_name("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing"),
                )
                .arg(
                    Arg::with_name("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches"),
                )
                .arg(
                    Arg::with_name("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::with_name("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("num-pruning-samples")
                        .long("num-pruning-samples")
                        .multiple(true)
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("min-prune-factor")
                        .long("min-prune-factor")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-adaptive-pruning")
                        .long("dont-use-adaptive-pruning")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases"),
                )
                .arg(
                    Arg::with_name("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .default_value("4.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("max-input-depth")
                        .long("max-input-depth")
                        .short("i")
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("error-correct-reads")
                        .long("error-correct-reads")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-observations-for-kmers-to-be-solid")
                        .long("min-observations-for-kmers-to-be-solid")
                        .default_value("20")
                        .hidden(true)
                )
                .arg(
                    Arg::with_name("max-mnp-distance")
                        .long("max-mnp-distance")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .default_value("-1")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("assembly-region-padding")
                        .long("assembly-region-padding")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .hidden(true),
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
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short("q")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .default_value("18"),
                )
                .arg(
                    Arg::with_name("enable-dynamic-read-disqualification-for-genotyping")
                        .long("enable-dynamic-read-disqualification-for-genotyping")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .default_value("1.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .default_value("0.02")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .hidden(true),
                )
                .arg(Arg::with_name("disable-optimizations").long("disable-optimizations"))
                .arg(Arg::with_name("disable-avx").long("disable-avx"))
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs").long("allow-improper-pairs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("consensus")
                .about("Generate consensus genomes across all provided samples")
                .help(CONSENSUS_HELP.as_str())
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
                        .short("1")
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
                        .short("2")
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
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short("f")
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short("p")
                        .long("parallel-genomes")
                        .default_value("4")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
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
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("min-contig-size")
                        .long("min-contig-size")
                        .default_value("2500"),
                )
                .arg(
                    Arg::with_name("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .default_value("45"),
                )
                .arg(
                    Arg::with_name("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .possible_values(&[
                            "none",
                            "None",
                            "NONE",
                            "hostile",
                            "Hostile",
                            "HOSTILE",
                            "aggresive",
                            "Agressive",
                            "AGGRESSIVE",
                            "conservative",
                            "Conservative",
                            "CONSERVATIVE",
                        ]),
                )
                .arg(
                    Arg::with_name("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::with_name("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short("C")
                        .default_value("30.0"),
                )
                .arg(
                    Arg::with_name("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .possible_values(&[
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual"),
                )
                .arg(
                    Arg::with_name("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles"),
                )
                .arg(
                    Arg::with_name("active-probability-threshold")
                        .long("active-probability-threshold")
                        .default_value("0.002"),
                )
                .arg(
                    Arg::with_name("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .default_value("300"),
                )
                .arg(
                    Arg::with_name("kmer-sizes")
                        .long("kmer-sizes")
                        .short("k")
                        .multiple(true)
                        .default_value("25"), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("256")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles"),
                )
                .arg(
                    Arg::with_name("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref"),
                )
                .arg(
                    Arg::with_name("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing"),
                )
                .arg(
                    Arg::with_name("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches"),
                )
                .arg(
                    Arg::with_name("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::with_name("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("num-pruning-samples")
                        .long("num-pruning-samples")
                        .multiple(true)
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("min-prune-factor")
                        .long("min-prune-factor")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-adaptive-pruning")
                        .long("dont-use-adaptive-pruning")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases"),
                )
                .arg(
                    Arg::with_name("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .default_value("0.001"),
                )
                .arg(
                    Arg::with_name("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .default_value("4.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("max-input-depth")
                        .long("max-input-depth")
                        .short("i")
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .default_value("50"),
                )
                .arg(
                    Arg::with_name("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("error-correct-reads")
                        .long("error-correct-reads")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-mnp-distance")
                        .long("max-mnp-distance")
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .default_value("-1")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("assembly-region-padding")
                        .long("assembly-region-padding")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .default_value("75")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .default_value("25")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .hidden(true),
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
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short("q")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .default_value("18"),
                )
                .arg(
                    Arg::with_name("enable-dynamic-read-disqualification-for-genotyping")
                        .long("enable-dynamic-read-disqualification-for-genotyping")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .default_value("1.0")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .default_value("0.02")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .default_value("2")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .hidden(true),
                )
                .arg(Arg::with_name("disable-optimizations").long("disable-optimizations"))
                .arg(Arg::with_name("disable-avx").long("disable-avx"))
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("allow-improper-pairs").long("allow-improper-pairs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        );
}
