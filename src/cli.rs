use bird_tool_utils::clap_utils::{add_clap_verbosity_flags, default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Example, Flag, Manual, Opt, Section};
use clap::*;
use clap_complete::*;
use galah::cluster_argument_parsing::GalahClustererCommandDefinition;
use roff::bold as roff_bold;
use roff::Roff;

// See https://github.com/rust-cli/roff-rs/issues/19
fn bold(s: &str) -> String {
    Roff::new().text([roff_bold(s)]).to_roff()
}

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "bwa-mem2",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-hifi",
    "minimap2-no-preset",
    "ngmlr-ont",
    "ngmlr-pb"
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

const LONGREAD_MAPPING_SOFTWARE_LIST: &[&str] =
    &["minimap2-ont", "minimap2-pb", "minimap2-hifi", "ngmlr-ont", "ngmlr-pb"];
const DEFAULT_LONGREAD_MAPPING_SOFTWARE: &str = "minimap2-ont";

fn add_mapping_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Mapping algorithm options")
            .option(Opt::new("NAME").long("--mapper").help(&format!(
                "Underlying mapping software used for short reads {}. One of: {}",
                default_roff("minimap2-sr"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("minimap2-sr"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x sr"))
                    ],
                    &[
                        &monospace_roff("bwa-mem"),
                        &format!("bwa mem using default parameters")
                    ],
                    &[
                        &monospace_roff("bwa-mem2"),
                        &format!("bwa-mem2 using default parameters")
                    ],
                    &[
                        &monospace_roff("minimap2-ont"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("minimap2-pb"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("minimap2-hifi"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
                    ],
                    &[
                        &monospace_roff("ngmlr-ont"),
                        &format!("ngmlr with '{}' option", &monospace_roff("-x ont"))
                    ],
                    &[
                        &monospace_roff("ngmlr-pb"),
                        &format!("ngmlr with '{}' option", &monospace_roff("-x pb"))
                    ],
                ])
            )))
            .option(Opt::new("NAME").long("--longread-mapper").help(&format!(
                "Underlying mapping software used for long reads {}. One of: {}",
                default_roff("minimap2-sr"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("minimap2-ont"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("minimap2-pb"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("minimap2-hifi"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
                    ],
                    &[
                        &monospace_roff("ngmlr-ont"),
                        &format!("ngmlr with '{}' option", &monospace_roff("-x ont"))
                    ],
                    &[
                        &monospace_roff("ngmlr-pb"),
                        &format!("ngmlr with '{}' option", &monospace_roff("-x pb"))
                    ],
                ])
            )))
            .option(Opt::new("PARAMS").long("--minimap2-params").help(&format!(
                "Extra parameters to provide to minimap2, \
        both indexing command (if used) and for \
        mapping. Note that usage of this parameter \
        has security implications if untrusted input \
        is specified. '{}' is always specified to minimap2. \
        [default: none]",
                &monospace_roff("-a")
            )))
            .flag(Flag::new().long("--minimap2-reference-is-index").help(
                "Treat reference as a minimap2 database, not as a FASTA file. [default: not set]",
            ))
            .option(Opt::new("PARAMS").long("--bwa-params").help(
                "Extra parameters to provide to BWA or BWA-MEM2. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--ngmlr-params").help(
                "Extra parameters to provide to NGMLR. \
        --bam-fix, -x ont, -t are already set. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
            )),
    )
}

fn add_thresholding_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Alignment thresholding")
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length")
                    .help(&format!(
                        "Exclude reads with smaller numbers of \
        aligned bases. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity")
                    .help(&format!(
                        "Exclude reads by overall percent \
        identity e.g. 95 for 95%. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent")
                    .help(&format!(
                        "Exclude reads by percent aligned \
        bases e.g. 95 means 95% of the read's \
        bases must be aligned. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length-pair")
                    .help(&format!(
                        "Exclude pairs with smaller numbers of \
        aligned bases. \
        Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity-pair")
                    .help(&format!(
                        "Exclude pairs by overall percent \
                identity e.g. 95 for 95%. \
                Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent-pair")
                    .help(&format!(
                        "Exclude reads by percent aligned \
                bases e.g. 95 means 95% of the read's \
                bases must be aligned. \
                Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .flag(
                Flag::new()
                    .long("--proper-pairs-only")
                    .help("Require reads to be mapped as proper pairs. [default: not set]"),
            )
            .flag(
                Flag::new()
                    .long("--exclude-supplementary")
                    .help("Exclude supplementary alignments. [default: not set]"),
            )
            .flag(
                Flag::new()
                    .long("--include-secondary")
                    .help("Include secondary alignments. [default: not set]"),
            )
            .option(
                Opt::new("INT")
                    .long("--contig-end-exclusion")
                    .help(
                        "Exclude bases at the ends of reference \n
                         sequences from calculation [default: 0]"
                    ),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--trim-min")
                    .help(
                        "Remove this smallest fraction of positions \n
                         when calculating trimmed_mean [default: 0.00]"
                    ),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--trim-max")
                    .help(
                        "Maximum fraction for trimmed_mean \n
                         calculations [default: 1.00]"
                    ),
            )
            .flag(
                Flag::new()
                    .long("--split-bams")
                    .help(
                        "Split the mapped read files up per reference.
                         Useful if you think run time is being hampered
                         by I/O. Most of the time this will not improve
                         performance and instead just increase disk usage."
                    ),
            )
    )
}

fn read_mapping_params_section() -> Section {
    Section::new("Read mapping parameters")
        .option(
            Opt::new("PATH ..")
                .short("-1")
                .help("Forward FASTA/Q file(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH ..")
                .short("-2")
                .help("Reverse FASTA/Q file(s) for mapping. These may be gzipped or not."),
        )
        .option(Opt::new("PATH ..").short("-c").long("--coupled").help(
            "One or more pairs of forward and reverse \
        possibly gzipped FASTA/Q files for mapping in order \
        <sample1_R1.fq.gz> <sample1_R2.fq.gz> \
        <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..",
        ))
        .option(
            Opt::new("PATH ..")
                .long("--interleaved")
                .help("Interleaved FASTA/Q files(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH ..")
                .long("--single")
                .help("Unpaired FASTA/Q files(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH ..")
                .long("--longreads")
                .help("Longread FASTA/Q files(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH")
                .short("-b")
                .long("--bam-files")
                .help(&format!(
                    "Path to BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken.",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        )
        .option(
            Opt::new("PATH")
                .short("-l")
                .long("--longread-bam-files")
                .help(&format!(
                    "Path to longread BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken.",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        )
}

fn reference_options() -> Section {
    Section::new("Input reference options")
        .option(
            Opt::new("PATH")
                .short("-r")
                .long("--reference")
                .help(
                    &format!(
                        "FASTA files of contigs e.g. concatenated \
                    genomes or metagenome assembly
                    [required unless {} is specified]",
                        monospace_roff("-d/--genome-fasta-directory")
                    )
                )
        )
        .option(
            Opt::new("PATH")
                .short("-d")
                .long("--genome-fasta-directory")
                .help(
                    &format!(
                        "Directory containing FASTA files of contigs e.g. \
                    genomes or metagenome assembly
                    [required unless {} is specified]",
                        monospace_roff("-r/--reference")
                    )
                )

        )
        .option(
            Opt::new("STR")
                .short("-x")
                .long("--genome-fasta-extension")
                .help(
                    &format!(
                        "FASTA file extension in --genome-fasta-directory \
                        [default \"fna\"]"
                    )
                )
        )
}

fn threads_options() -> Section {
    Section::new("Threading options")
        .option(
            Opt::new("INT")
                .long("--threads")
                .short("-t")
                .help(
                    "Maximum number of threads used. [default: 8]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--parallel-genomes")
                .short("-p")
                .help(
                    "Number of genomes to run in parallel. \
                     Increases memory usage linearly. \
                     Thread usage qill not exceed the value \
                     provided by --threads [default 4]"
                )
        )
}

fn add_help_options(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-h")
                .long("--help")
                .help("Output a short usage message. [default: not set]"),
        )
        .flag(
            Flag::new()
                .long("--full-help")
                .help("Output a full help message and display in 'man'. [default: not set]"),
        )
        .flag(Flag::new().long("--full-help-roff").help(
            "Output a full help message in raw ROFF format for \
        conversion to other formats. [default: not set]",
        ))
}

fn add_help_options_to_section(section: Section) -> Section {
    section
        .flag(
            Flag::new()
                .short("-h")
                .long("--help")
                .help("Output a short usage message. [default: not set]"),
        )
        .flag(
            Flag::new()
                .long("--full-help")
                .help("Output a full help message and display in 'man'. [default: not set]"),
        )
        .flag(Flag::new().long("--full-help-roff").help(
            "Output a full help message in raw ROFF format for \
        conversion to other formats. [default: not set]",
        ))
}

fn sharding_section() -> Section {
    Section::new("Sharding").flag(Flag::new().long("--sharded").help(&format!(
        "If {} was used: \
        Input BAM files are read-sorted alignments \
        of a set of reads mapped to multiple \
        reference contig sets. Choose the best \
        hit for each read pair. Otherwise if mapping was carried out: \
        Map reads to each reference, choosing the \
        best hit for each pair. [default: not set]",
        monospace_roff("-b/--bam-files")
    )))
}

fn faq_section() -> Section {
    Section::new("Frequently asked questions (FAQ)").paragraph(&format!(
        "{} Lorikeet makes use of \
        the system temporary directory (often {}) to store intermediate files. This can cause \
        problems if the amount of storage available there is small or used by many programs. \
        To fix, set the {} environment variable e.g. to set it to use the current directory: {}\n",
        bold("Can the temporary directory used be changed?"),
        monospace_roff("/tmp"),
        monospace_roff("TMPDIR"),
        monospace_roff("TMPDIR=. lorikeet call <etc>"),
    ))
}

fn add_verbosity_flags(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information. [default: not set]"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

fn variant_calling_section_basic() -> Section {
    Section::new("Variant calling options (Basic)")
        .option(
            Opt::new("INT ..")
                .long("--kmer-sizes")
                .short("-k")
                .help(
                    "K-mer sizes used to generate DeBruijn Graphs. \
                     Multiple values at once are accepted and encouraged \
                     e.g. 10 25 [default: 10 25] \n"
                )
        )
        .option(
            Opt::new("INT")
                .long("--ploidy")
                .help(
                    "Sets the default ploidy for the analysis to N. \
                    [default: 1]"
                )
        )
        .flag(
            Flag::new()
                .long("--calculate-fst")
                .help("Calculate Fst values between samples and variants.")
        )
        .flag(
            Flag::new()
                .long("--calculate-dnds")
                .help(
                    "Calculate coding regions and perform dN/dS calculations \
                    along them using called variants. *Microbial only*."
                )
        )
        .option(
            Opt::new("PATH")
                .short("-f")
                .long("--features-vcf")
                .help(
                    "The set of alleles to force-call regardless \
                     of evidence. Note: The sight containing these alleles \
                     has to be called as 'active' in order for them to appear \
                     in the final VCF. Addtionally, Provided file must be \
                     compressed using bgzip and indexed using bcftools index. If no index \
                     is present, and index will be attempted to be created. \
                     If the file is not properly compressed, Lorikeet will \
                     unfortunately SEGFAULT with no error message."
                )
        )
        .option(
            Opt::new("INT")
                .long("--qual-by-depth-filter")
                .help(
                    "The minimum QD value for a variant to have for it to be \
                     included in the genotyping or ANI analyses. [default: 25]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--qual-threshold")
                .help(
                    "The PHRED-scaled quality score threshold for use \
                     with ANI calculations. [default: 150]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--depth-per-sample-filter")
                .help(
                    "Minimum depth of a variant in a sample for that \
                     sample to be included in ANI & Fst calculations for that \
                     variant. [default: 5]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-long-read-size")
                .help(
                    "The minimum size for long reads to be used for analysis \
                    [default: 1500]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-long-read-average-base-qual")
                .help(
                    "The minimum average base quality of a long read \
                     for it to be used for analysis [default: 20]"
                )
        )
        .option(
            Opt::new("INT")
                .short("-q")
                .long("--min-base-quality")
                .help(
                    "Minimum base quality required to consider a \
                     base for calling. [default: 10]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-mapq")
                .help(
                    "Minimum MAPQ score for reads to be considered \
                     during variant calling. [default: 20]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--base-quality-score-threshold")
                .help(
                    "Base qualities below this threshold will \
                     be reduced to the minimum (6). [default: 18]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--max-input-depth")
                .help(
                    "The maximum number of reads included within an \
                     assembly region across all samples. Larger numbers \
                     increase run time. If the depth of an assembly region \
                     exceeds this value, then the reads will be filtered \
                     by mean base quality. [default: 200000]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-contig-size")
                .help(
                    "The minimum contig size to call variants on. Smaller \
                    contigs can often contain highly variable regions that \
                    mostly represent noise. Call variants on them can often \
                    be slow and not produce anything fruitful. If you \
                    wish to call variants on all available contigs, \
                    then set this to 0. [default: 2500]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-sv-qual")
                .help(
                    "Minimum structural variants quality returned by svim \
                     and used by lorikeet. Not PHRED-scaled quality, value \
                     determined by number of supporting reads. Consult \
                     svim documentation for details. [default: 3]"
                )
        )
        .flag(
            Flag::new()
                .long("--do-not-call-svs")
                .help(
                    "Opts not to use svim to call structural variants \
                     using provided longreads. If no longreads are provided \
                     this has no effect."
                )
        )
}

fn variant_calling_options_advanced() -> Section {
    Section::new("Variant calling options (Advanced)")
        .option(
            Opt::new("INT")
                .long("--phred-scaled-global-read-mismapping-rate")
                .help(
                    "The global assumed mismapping rate for reads. [default: 45]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--pair-hmm-gap-continuation-penalty ")
                .help(
                    "Flat gap continuation penalty for use in the Pair HMM. \
                    [default: 10]"
                )
        )
        .option(
            Opt::new("STR")
                .long("--pcr-indel-model")
                .help(
                    "The PCR indel model to use. [default: conservative]"
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--heterozygosity")
                .help(
                    "Heterozygosity value used to compute prior \
                     likelihoods for any locus. [default: 0.001]"
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--heterozygosity-stdev")
                .help(
                    "Standard deviation of heterozygosity for SNP and \
                     indel calling. [default: 0.01]"
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--indel-heterozygosity")
                .help(
                    "Heterozygosity for indel calling. [default: 0.000125]"
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--standard-min-confidence-threshold-for-calling")
                .short("-C")
                .help(
                    "The minimum phred-scaled confidence threshold at \
                     which variants should be called. [default: 30.0]"
                )
        )
        .flag(
            Flag::new()
                .long("--use-posteriors-to-calculate-qual")
                .help(
                    "if available, use the genotype posterior \
                     probabilities to calculate the site QUAL."
                )
        )
        .flag(
            Flag::new()
                .long("--annotate-with-num-discovered-alleles")
                .help(
                    "If provided, we will annotate records with the \
                     number of alternate alleles that were discovered \
                     (but not necessarily genotyped) at a given site."
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--active-probability-threshold")
                .help(
                    "Minimum probability for a locus to be \
                     considered active. [default: 0.002]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-assembly-region-size")
                .help(
                    "Minimum size of an assembly region. [default: 50]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--max-assembly-region-size")
                .help(
                    "Maximum size of an assembly region. [default: 300]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--assembly-region-padding")
                .help(
                    "Number of additional bases of context to \
                     include around each assembly region. [default: 100]"
                )
        )
        .flag(
            Flag::new()
                .long("--dont-increase-kmer-sizes-for-cycles")
                .help(
                    "Disable iterating over kmer sizes when \
                     graph cycles are detected."
                )
        )
        .flag(
            Flag::new()
                .long("--allow-non-unique-kmers-in-ref")
                .help(
                    "Allow graphs that have non-unique kmers in the reference."
                )
        )
        .flag(
            Flag::new()
                .long("--do-not-run-physical-phasing")
                .help(
                    "Disable physical phasing."
                )
        )
        .flag(
            Flag::new()
                .long("--recover-all-dangling-branches")
                .help(
                    "Recover all dangling branches."
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-dangling-branch-length")
                .help(
                    "Minimum length of a dangling branch to \
                     attempt recovery. [default: 4]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--min-prune-factor")
                .help(
                    "Minimum support to not prune paths in the graph. \
                     [default: 2]"
                )
        )
        .flag(
            Flag::new()
                .long("--use-adaptive-pruning")
                .help(
                    "Use more advanced pruning algorithm to prune paths in
                     graph. Better suited when performing variant calling
                     on when depth along a genome is variable e.g. RNA
                     and exome data."
                )
        )
        .option(
            Opt::new("INT")
                .long("--num-pruning-samples")
                .help(
                    "Number of samples that must pass the \
                     min_pruning threshold [default: 1]"
                )
        )
        .option(
            Opt::new("PATH")
                .long("--graph-output")
                .help(
                    "Write debug assembly graph information to this file."
                )
        )
        .flag(
            Flag::new()
                .long("--dont-use-soft-clipped-bases")
                .help(
                    "Do not analyse soft clipped bases in the reads."
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--initial-error-rate-for-pruning")
                .help(
                    "Initial base error rate estimate for adaptive \
                     pruning. [default: 0.001]"
                )
        )
        .option(
            Opt::new("FLOAT")
                .long("--pruning-log-odds-threshold")
                .help(
                    "Likelihood ratio threshold for adaptive \
                     pruning algorithm. This value will be converted to \
                     log odds value. [default: 1.0]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--max-unpruned-variants")
                .help(
                    "Maximum number of variants in graph the \
                     adaptive pruner will allow. [default: 100]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--max-prob-propagation-distance")
                .help(
                    "Upper limit on how many bases away probability mass \
                     can be moved around when calculating the boundaries \
                     between active and inactive assembly regions. [default: 50]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--max-mnp-distance")
                .help(
                    "Two or more phased substitutions separated by \
                     this distance or less are merged into MNPs. [default: 0]"
                )
        )
        .flag(
            Flag::new()
                .long("--disable-optimizations")
                .help(
                    "Don't skip calculations in ActiveRegions with no variants"
                )
        )
        .flag(
            Flag::new()
                .long("--disable-avx")
                .help(
                    "Disable the use of the GKL-rs AVX acceleration components \
                     for PairHMM and Smith-Waterman calculations."
                )
        )
        .option(
            Opt::new("STR")
                .long("--limiting-interval")
                .help(
                    "Mainly used for debugging purposes. Only call variants \
                     within this given span on all contigs. E.g. providing \
                     '1000-2000' would only call variants between the 1000 \
                     and 2000 bp span on each provided contig."
                )
        )
        .flag(
            Flag::new()
                .long("--force")
                .help(
                    "Forcefully overwrite previous runs."
                )
        )
}

fn add_verbosity_flags_to_section(section: Section) -> Section {
    section
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information. [default: not set]"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

pub fn genotype_full_help() -> Manual {
    let mut manual = Manual::new("lorikeet genotype")
        .about(
            &format!(
                "Call variants and cluster into potential strain haplotypes (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email("rhys.newell94 near gmail.com"))
        .description(
            "
            ======= EXPERIMENTAL =======\n
            lorikeet genotype discovers variants within a given set of reads and genomes and \n
            clusters the variants into candidate strain haplotypes. Lorikeet uses UMAP and HDBSCAN \n
            to cluster variants and an Expectation-Maximization algorithm to determine strain \n
            haplotype abudnances within each samples. Additionally, calculate strain \n
            diversity metrics like conANI, popANI, subpopANI, dN/dS, and the highly robust Hudson's Fst\n
            \n
            This process can be undertaken in several ways, for instance by specifying BAM files \n
            or raw reads as input, using different mapping programs, thresholding read alignments \n
            ============================\
            "
        );


    manual = manual.custom(
        threads_options()
    );
    manual = manual.custom(
        reference_options()
    );
    manual = manual.custom(
        read_mapping_params_section()
    );
    manual = manual.custom(
        sharding_section()
    );
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(
        variant_calling_section_basic()
    );
    manual = manual.custom(
        variant_calling_options_advanced()
    );
    manual = manual.custom(
        Section::new("Output options")
            .option(Opt::new("DIRECTORY").short("-o").long("--output-directory").help(
                "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used]",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--discard-unmapped")
                    .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
            )
    );

    manual = manual.example(
        Example::new()
            .text("Map paired reads to a reference and generate genotypes")
            .command(
                "lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25 51",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Generate strain-level genotypes from read mappings compared to reference from a sorted BAM file and plots the results",
            )
            .command(
                "lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam \
                --genome-fasta-directory genomes/ -x fna --bam-file-cache-directory saved_bam_files \
                --output-directory lorikeet_out/ --threads 10",
            ),
    );

    manual = add_verbosity_flags(manual);

    manual = manual.custom(faq_section());

    return manual;
}

pub fn call_full_help() -> Manual {
    let mut manual = Manual::new("lorikeet call")
        .about(
            &format!(
                "Call variants within a given set of genomes and samples using local reassembly (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email("rhys.newell94 near gmail.com"))
        .description(
            "
            ===========================\n
            lorikeet call discovers variants within a given set of reads and genomes using a local \n
            reassembly algorithm based on the GATK HaplotypeCaller. Additionally, calculate strain \n
            diversity metrics like conANI, popANI, subpopANI, dN/dS, and the highly robust Hudson's Fst \n\
            \n
            This process can be undertaken in several ways, for instance by specifying BAM files \n
            or raw reads as input, using different mapping programs, thresholding read alignments \n
            ============================\
            "
        );


    manual = manual.custom(
        threads_options()
    );
    manual = manual.custom(
        reference_options()
    );
    manual = manual.custom(
        read_mapping_params_section()
    );
    manual = manual.custom(
        sharding_section()
    );
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(
        variant_calling_section_basic()
    );
    manual = manual.custom(
        variant_calling_options_advanced()
    );
    manual = manual.custom(
        Section::new("Output options")
            .option(Opt::new("DIRECTORY").short("-o").long("--output-directory").help(
                "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used]",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--discard-unmapped")
                    .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
            )
    );

    manual = manual.example(
        Example::new()
            .text("Map paired reads to a reference call variants and calculate Fst")
            .command(
                "lorikeet call --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25 51 --calculate-fst",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Call variants from read mappings compared to reference from a sorted BAM file and plots the results",
            )
            .command(
                "lorikeet call --bam-files my.bam --longread-bam-files my-longread.bam \
                --genome-fasta-directory genomes/ -x fna --bam-file-cache-directory saved_bam_files \
                --output-directory lorikeet_out/ --threads 10",
            ),
    );

    manual = add_verbosity_flags(manual);

    manual = manual.custom(faq_section());

    return manual;
}

pub fn consensus_full_help() -> Manual {
    let mut manual = Manual::new("lorikeet consensus")
        .about(
            &format!(
                "Call variants within a given set of genomes and samples using local reassembly (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email("rhys.newell94 near gmail.com"))
        .description(
            "
            ===========================\n
            lorikeet consensus discovers variants within a given set of reads and genomes using a local \n
            reassembly algorithm based on the GATK HaplotypeCaller. Additionally, calculate strain \n
            diversity metrics like conANI, popANI, subpopANI, dN/dS, and the highly robust Hudson's Fst \n\
            \n\
            Lorikeet consensus also generates the consensus strain haplotypes for each sample and prints \
            them as a FASTA file for each input genome in the output directory. \n\
            \n\
            This process can be undertaken in several ways, for instance by specifying BAM files \n
            or raw reads as input, using different mapping programs, thresholding read alignments \n
            ============================\
            "
        );


    manual = manual.custom(
        threads_options()
    );
    manual = manual.custom(
        reference_options()
    );
    manual = manual.custom(
        read_mapping_params_section()
    );
    manual = manual.custom(
        sharding_section()
    );
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(
        variant_calling_section_basic()
    );
    manual = manual.custom(
        variant_calling_options_advanced()
    );
    manual = manual.custom(
        Section::new("Output options")
            .option(Opt::new("DIRECTORY").short("-o").long("--output-directory").help(
                "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used]",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--discard-unmapped")
                    .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
            )
    );

    manual = manual.example(
        Example::new()
            .text("Map paired reads to a reference call variants and calculate Fst")
            .command(
                "lorikeet consensus --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25 51 --calculate-fst",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Call variants from read mappings compared to reference from a sorted BAM file and plots the results",
            )
            .command(
                "lorikeet consensus --bam-files my.bam --longread-bam-files my-longread.bam \
                --genome-fasta-directory genomes/ -x fna --bam-file-cache-directory saved_bam_files \
                --output-directory lorikeet_out/ --threads 10",
            ),
    );

    manual = add_verbosity_flags(manual);

    manual = manual.custom(faq_section());

    return manual;
}

pub fn summarise_full_help() -> Manual {
    let mut manual = Manual::new("lorikeet summarise")
        .about(
            &format!(
                "Calculate ANI and Fst metrics on a given set of VCF files (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email("rhys.newell94 near gmail.com"))
        .description(
            "
            ===========================\n
            lorikeet summarise uses a set of VCF files as input and calculates conANI, popANI, \n\
            subpopANI, and Fst metrics for the variants in each file.

            NOTE: ANI metrics require coverage information to determine the number of shared bases
            in each sample. VCF files do not provide this information, so the shared base size is just
            the total size of the genome. In our experience, this doesn't really matter that much as
            the ANI metrics are quite insensitive when provided low coverage samples anyway.
            Fst tends to perform better for low and high coverage samples and does not require whole
            genome coverage information.
            ============================\
            "
        );

    manual = manual.option(
        Opt::new("PATH ..")
            .short("-i")
            .long("--vcfs")
            .help(
                "Paths to input VCF files. Can provide one or more."
            )
        )
        .option(
            Opt::new("DIRECTORY")
                .short("-i")
                .long("--vcfs")
                .help(
                    "Paths to input VCF files. Can provide one or more."
                )
        )
        .option(Opt::new("DIRECTORY").short("-o").long("--output").help(
            "Output directory. Folder will contain subfolders for each input VCF \n
             [default: ./]",
        ))
        .option(
            Opt::new("INT")
                .long("--threads")
                .short("-t")
                .help(
                    "Maximum number of threads used. [default: 8]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--qual-by-depth-filter")
                .help(
                    "The minimum QD value for a variant to have for it to be \
                     included in the genotyping or ANI analyses. [default: 25]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--qual-threshold")
                .help(
                    "The PHRED-scaled quality score threshold for use \
                     with ANI calculations. [default: 150]"
                )
        )
        .option(
            Opt::new("INT")
                .long("--depth-per-sample-filter")
                .help(
                    "Minimum depth of a variant in a sample for that \
                     sample to be included in ANI & Fst calculations for that \
                     variant. [default: 5]"
                )
        );

    manual = add_verbosity_flags(manual);
    return manual;
}

pub fn build_cli() -> Command<'static> {
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

  lorikeet summarise --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --window-size 1

{}

  lorikeet summarise --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10

See lorikeet summarise --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "lorikeet summarise"),
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
        .author(crate::AUTHOR_AND_EMAIL)
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
\tcall      \tPerforms variant calling on the provides genomes
\tconsensus \tCreates consensus genomes for each input reference and for each sample

Utility subcommands:
\tsummarise \tCalculate microdiversity statistics for a given set of VCF files
\tshell-completion  \tGenerate shell completion scripts

Experimental subcommands:
\tgenotype  \tReport strain-level genotypes and abundances from metagenomes

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>
",
        )
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("genotype")
                .about("Perform variant calling analysis and then binning")
                .help(GENOTYPE_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(Arg::with_name("full-help-roff").long("full-help-roff"))
                .arg(
                    Arg::with_name("bam-files")
                        .short('b')
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short('1')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short('2')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short('c')
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longread-bam-files"
                        ])
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short('l')
                        .long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                        ])
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short('r')
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .takes_value(true)
                        .required_unless_one(&["genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
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
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short('t')
                        .long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short('p')
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
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("high-memory")
                        .long("high-memory")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::with_name("split-bams").long("split-bams"))
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
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short('m')
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
                        .short('C')
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
                        .short('k')
                        .takes_value(true)
                        .multiple(true)
                        .default_values(&["10", "25"]), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("128")
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
                        .default_value("2"),
                )
                .arg(Arg::with_name("use-adaptive-pruning").long("use-adaptive-pruning"))
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
                        .short('i')
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("0"),
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
                        .default_value("0.00"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("1.00"),
                )
                .arg(
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-sv-qual")
                        .long("min-sv-qual")
                        .default_value("3"),
                )
                .arg(Arg::with_name("do-not-call-svs").long("do-not-call-svs"))
                .arg(
                    Arg::with_name("min-mapq")
                        .long("min-mapq")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-long-read-size")
                        .long("min-long-read-size")
                        .default_value("1500"),
                )
                .arg(
                    Arg::with_name("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
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
                        .default_value("25"),
                )
                .arg(
                    Arg::with_name("qual-threshold")
                        .long("qual-threshold")
                        .default_value("150"),
                )
                .arg(
                    Arg::with_name("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .default_value("5"),
                )
                .arg(
                    Arg::with_name("min-variant-depth-for-genotyping")
                        .long("min-variant-depth-for-genotyping")
                        .default_value("5"),
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
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("calculate-dnds")
                        .long("calculate-dnds")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("calculate-fst")
                        .long("calculate-fst")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("prodigal-params")
                        .long("prodigal-params")
                        .takes_value(true)
                        .default_value("-p meta"),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short('v').long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("call")
                .about("Perform variant calling across the given genomes and samples")
                .help(CALL_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short('b').long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short('1')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short('2')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short('c').long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longread-bam-files"
                        ])
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short('l').long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                        ])
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short('r').long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .takes_value(true)
                        .required_unless_one(&["genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
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
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short('t').long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short('p').long("parallel-genomes")
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
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("high-memory")
                        .long("high-memory")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::with_name("split-bams").long("split-bams"))
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
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short('m').long("method")
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
                        .short('C')
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
                        .short('k')
                        .takes_value(true)
                        .multiple(true)
                        .default_values(&["10", "25"]), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("128")
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
                        .default_value("2"),
                )
                .arg(Arg::with_name("use-adaptive-pruning").long("use-adaptive-pruning"))
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
                        .short('i')
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("0"),
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
                        .default_value("0.00"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("1.00"),
                )
                .arg(
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-sv-qual")
                        .long("min-sv-qual")
                        .default_value("3"),
                )
                .arg(Arg::with_name("do-not-call-svs").long("do-not-call-svs"))
                .arg(
                    Arg::with_name("min-mapq")
                        .long("min-mapq")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-long-read-size")
                        .long("min-long-read-size")
                        .default_value("1500"),
                )
                .arg(
                    Arg::with_name("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
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
                        .default_value("25"),
                )
                .arg(
                    Arg::with_name("qual-threshold")
                        .long("qual-threshold")
                        .default_value("150"),
                )
                .arg(
                    Arg::with_name("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .default_value("5"),
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
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("calculate-dnds")
                        .long("calculate-dnds")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("calculate-fst")
                        .long("calculate-fst")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("prodigal-params")
                        .long("prodigal-params")
                        .default_value("-p meta")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short('v').long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("consensus")
                .about("Generate consensus genomes across all provided samples")
                .help(CONSENSUS_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("bam-files")
                        .short('b')
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ]),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short('1')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("read2")
                        .short('2')
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short('c')
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
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
                            "full-help", "full-help-roff",
                            "longreads",
                            "longread-bam-files"
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longread-bam-files"
                        ])
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short('l')
                        .long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help", "full-help-roff",
                            "longreads",
                        ])
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short('r')
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["genome-fasta-directory", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .takes_value(true)
                        .required_unless_one(&["genome-fasta-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
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
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .takes_value(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short('t')
                        .long("threads")
                        .default_value("8")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("parallel-genomes")
                        .short('p')
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
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("high-memory")
                        .long("high-memory")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::with_name("split-bams").long("split-bams"))
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
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("method")
                        .short('m')
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
                        .short('C')
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
                        .short('k')
                        .takes_value(true)
                        .multiple(true)
                        .default_values(&["10", "25"]), //TODO: Wait for clap v3 and change this to default_values
                )
                .arg(
                    Arg::with_name("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .default_value("128")
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
                        .default_value("2"),
                )
                .arg(Arg::with_name("use-adaptive-pruning").long("use-adaptive-pruning"))
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
                        .short('i')
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("0"),
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
                        .default_value("0.00"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("1.00"),
                )
                .arg(
                    Arg::with_name("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .default_value("20")
                        .hidden(true),
                )
                .arg(
                    Arg::with_name("min-sv-qual")
                        .long("min-sv-qual")
                        .default_value("3"),
                )
                .arg(Arg::with_name("do-not-call-svs").long("do-not-call-svs"))
                .arg(
                    Arg::with_name("min-mapq")
                        .long("min-mapq")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-long-read-size")
                        .long("min-long-read-size")
                        .default_value("1500"),
                )
                .arg(
                    Arg::with_name("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .default_value("20"),
                )
                .arg(
                    Arg::with_name("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
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
                        .default_value("25"),
                )
                .arg(
                    Arg::with_name("qual-threshold")
                        .long("qual-threshold")
                        .default_value("150"),
                )
                .arg(
                    Arg::with_name("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .default_value("5"),
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
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("calculate-dnds")
                        .long("calculate-dnds")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("calculate-fst")
                        .long("calculate-fst")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("prodigal-params")
                        .long("prodigal-params")
                        .takes_value(true)
                        .default_value("-p meta"),
                )
                .arg(
                    Arg::with_name("limiting-interval")
                        .long("limiting-interval")
                        .required(false)
                        .takes_value(true),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short('v').long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("summarise")
                .about("Summarizes ANI values of a given set of VCF files")
                .arg(
                    Arg::with_name("vcfs")
                        .long("vcfs")
                        .short('i')
                        .takes_value(true)
                        .multiple(true)
                        .required(true),
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::with_name("threads")
                        .long("threads")
                        .short('t')
                        .default_value("8"),
                )
                .arg(
                    Arg::with_name("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .default_value("25"),
                )
                .arg(
                    Arg::with_name("qual-threshold")
                        .long("qual-threshold")
                        .default_value("150"),
                )
                .arg(
                    Arg::with_name("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .default_value("5"),
                )
                .arg(Arg::with_name("verbose").short('v').long("verbose")),
        )
        .subcommand(
            add_clap_verbosity_flags(Command::new("shell-completion"))
                .about("Generate a shell completion script for lorikeet")
                .arg(
                    Arg::new("output-file")
                        .short('o')
                        .long("output-file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("shell")
                        .long("shell")
                        .takes_value(true)
                        .required(true)
                        .allow_invalid_utf8(true)
                        .value_parser(value_parser!(Shell)),
                ),
        )
}
