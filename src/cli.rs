use bird_tool_utils::clap_utils::{add_clap_verbosity_flags, default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Example, Flag, Manual, Opt, Section};
use clap::*;
use clap_complete::*;
// use galah::cluster_argument_parsing::GalahClustererCommandDefinition;
use roff::bold as roff_bold;
use roff::Roff;
use crate::utils::utils::table_roff;

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
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

const LONGREAD_MAPPING_SOFTWARE_LIST: &[&str] = &[
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-hifi",
];
const DEFAULT_LONGREAD_MAPPING_SOFTWARE: &str = "minimap2-ont";

fn add_mapping_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Mapping algorithm options")
            .option(Opt::new("NAME").long("--mapper").help(&format!(
                "Underlying mapping software used for short reads {}. One of: {}",
                default_roff("minimap2-sr"),
                table_roff(&[
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
                ])
            )))
            .option(Opt::new("NAME").long("--longread-mapper").help(&format!(
                "Underlying mapping software used for long reads {}. One of: {}",
                default_roff("minimap2-ont"),
                table_roff(&[
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
                ])
            )))
            .option(Opt::new("PARAMS").long("--minimap2-params").help(&format!(
                "Extra parameters to provide to minimap2, \
        both indexing command (if used) and for \
        mapping. Note that usage of this parameter \
        has security implications if untrusted input \
        is specified. '{}' is always specified to minimap2. \
        [default: none] \n",
                &monospace_roff("-a")
            )))
            .flag(Flag::new().long("--minimap2-reference-is-index").help(
                "Treat reference as a minimap2 database, not as a FASTA file. [default: not set]",
            ))
            .option(Opt::new("PARAMS").long("--bwa-params").help(
                "Extra parameters to provide to BWA or BWA-MEM2. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none] \n",
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
        aligned bases. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity")
                    .help(&format!(
                        "Exclude reads by overall percent \
        identity e.g. 95 for 95%. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent")
                    .help(&format!(
                        "Exclude reads by percent aligned \
        bases e.g. 95 means 95% of the read's \
        bases must be aligned. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length-pair")
                    .help(&format!(
                        "Exclude pairs with smaller numbers of \
        aligned bases. \
        Conflicts with --allow-improper-pairs. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity-pair")
                    .help(&format!(
                        "Exclude pairs by overall percent \
                identity e.g. 95 for 95%. \
                Conflicts with --allow-improper-pairs. {} \n",
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
                Conflicts with --allow-improper-pairs. {} \n",
                        default_roff("0")
                    )),
            )
            .flag(
                Flag::new()
                    .long("--allow-improper-pairs")
                    .help("Require reads to be mapped as proper pairs. [default: not set] \n"),
            )
            .flag(
                Flag::new()
                    .long("--exclude-supplementary")
                    .help("Exclude supplementary alignments. [default: not set] \n"),
            )
            .flag(
                Flag::new()
                    .long("--include-secondary")
                    .help("Include secondary alignments. [default: not set] \n"),
            )
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                "Exclude bases at the ends of reference \n
                         sequences from calculation [default: 0]",
            ))
            .option(Opt::new("FLOAT").long("--trim-min").help(
                "Remove this smallest fraction of positions \n
                         when calculating trimmed_mean [default: 0.00]",
            ))
            .option(Opt::new("FLOAT").long("--trim-max").help(
                "Maximum fraction for trimmed_mean \n
                         calculations [default: 1.00]",
            ))
            .flag(Flag::new().long("--split-bams").help(
                "Split the mapped read files up per reference.
                         Useful if you think run time is being hampered
                         by I/O. Most of the time this will not improve
                         performance and instead just increase disk usage. \n",
            )),
    )
}

fn read_mapping_params_section() -> Section {
    Section::new("Read mapping parameters")
        .option(
            Opt::new("PATH ..")
                .short("-1")
                .help("Forward FASTA/Q file(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .short("-2")
                .help("Reverse FASTA/Q file(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(Opt::new("PATH ..").short("-c").long("--coupled").help(
            "One or more pairs of forward and reverse \
        possibly gzipped FASTA/Q files for mapping in order \
        <sample1_R1.fq.gz> <sample1_R2.fq.gz> \
        <sample2_R1.fq.gz> <sample2_R2.fq.gz> .. \n",
        ))
        .option(
            Opt::new("PATH ..")
                .long("--interleaved")
                .help("Interleaved FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .long("--single")
                .help("Unpaired FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .long("--longreads")
                .help("Longread FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
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
                with {}). When specified, no read mapping algorithm is undertaken. \n",
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
                with {}). When specified, no read mapping algorithm is undertaken. \n",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        )
}

fn reference_options() -> Section {
    Section::new("Input reference options")
        .option(
            Opt::new("PATH")
                .short("-r,-f")
                .long("--reference,--genome-fasta-files")
                .help(&format!(
                    "FASTA files of contigs e.g. concatenated \
                    genomes or metagenome assembly
                    [required unless {} is specified] \n",
                    monospace_roff("-d/--genome-fasta-directory")
                )),
        )
        .option(
            Opt::new("PATH")
                .short("-d")
                .long("--genome-fasta-directory")
                .help(&format!(
                    "Directory containing FASTA files of contigs e.g. \
                    genomes or metagenome assembly
                    [required unless {} is specified] \n",
                    monospace_roff("-r/--reference/-f/--genome-fasta-files")
                )),
        )
        .option(
            Opt::new("STR")
                .short("-x")
                .long("--genome-fasta-extension")
                .help(&format!(
                    "FASTA file extension in --genome-fasta-directory \
                        [default \"fna\"] \n"
                )),
        )
}

fn threads_options() -> Section {
    Section::new("Threading options")
        .option(
            Opt::new("INT")
                .long("--threads")
                .short("-t")
                .help("Maximum number of threads used. [default: 10] \n"),
        )
        .option(Opt::new("INT").long("--parallel-genomes").short("-p").help(
            "Number of genomes to run in parallel. \
                     Increases memory usage linearly. \
                     Thread usage qill not exceed the value \
                     provided by --threads [default 10] \n",
        ))
}

// fn add_help_options(manual: Manual) -> Manual {
//     manual
//         .flag(
//             Flag::new()
//                 .short("-h")
//                 .long("--help")
//                 .help("Output a short usage message. [default: not set] \n"),
//         )
//         .flag(
//             Flag::new()
//                 .long("--full-help")
//                 .help("Output a full help message and display in 'man'. [default: not set] \n"),
//         )
//         .flag(Flag::new().long("--full-help-roff").help(
//             "Output a full help message in raw ROFF format for \
//         conversion to other formats. [default: not set] \n",
//         ))
// }

// fn add_help_options_to_section(section: Section) -> Section {
//     section
//         .flag(
//             Flag::new()
//                 .short("-h")
//                 .long("--help")
//                 .help("Output a short usage message. [default: not set] \n"),
//         )
//         .flag(
//             Flag::new()
//                 .long("--full-help")
//                 .help("Output a full help message and display in 'man'. [default: not set] \n"),
//         )
//         .flag(Flag::new().long("--full-help-roff").help(
//             "Output a full help message in raw ROFF format for \
//         conversion to other formats. [default: not set] \n",
//         ))
// }

fn sharding_section() -> Section {
    Section::new("Sharding").flag(Flag::new().long("--sharded").help(&format!(
        "If {} was used: \
        Input BAM files are read-sorted alignments \
        of a set of reads mapped to multiple \
        reference contig sets. Choose the best \
        hit for each read pair. Otherwise if mapping was carried out: \
        Map reads to each reference, choosing the \
        best hit for each pair. [default: not set] \n",
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
                .help("Print extra debugging information. [default: not set] \n"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

fn variant_calling_section_basic() -> Section {
    Section::new("Variant calling options (Basic)")
        .option(Opt::new("STR").long("--profile").help(
            "Assembly profile to use for variant calling. \
            Overrides --kmer-sizes, --min-prune-factor, --allow-non-unique-kmers-in-ref, and --recover-all-dangling-branches. \
            Possible options include: 'fast', 'very-fast', 'sensitive', 'precise', 'super-sensitive' \
                     [default: not_set] \n",
        ))
        .option(Opt::new("INT ..").long("--kmer-sizes").short("-k").help(
            "K-mer sizes used to generate DeBruijn Graphs. \
            Multiple values at once are accepted and encouraged at a cost to increased runtime \
                     e.g. 17 21 25 [default: 21 33] \n",
        ))
        .option(Opt::new("INT").long("--ploidy").help(
            "Sets the default ploidy for the analysis to N. \
                    [default: 2] \n",
        ))
        .flag(
            Flag::new()
                .long("--calculate-fst")
                .help("Calculate Fst values between samples and variants. \n"),
        )
        .flag(Flag::new().long("--calculate-dnds").help(
            "Calculate coding regions and perform dN/dS calculations \
                    along them using called variants. *Microbial only*. \n",
        ))
        .option(Opt::new("PATH").short("-f").long("--features-vcf").help(
            "The set of alleles to force-call regardless \
                     of evidence. Note: The sight containing these alleles \
                     has to be called as 'active' in order for them to appear \
                     in the final VCF. Addtionally, Provided file must be \
                     compressed using bgzip and indexed using bcftools index. If no index \
                     is present, and index will be attempted to be created. \
                     If the file is not properly compressed, Lorikeet will \
                     unfortunately SEGFAULT with no error message. \n",
        ))
        .option(Opt::new("INT").long("--qual-by-depth-filter").help(
            "The minimum QD value for a variant to have for it to be \
                     included in the genotyping or ANI analyses. [default: 25] \n",
        ))
        .option(Opt::new("INT").long("--qual-threshold").help(
            "The PHRED-scaled quality score threshold for use \
                     with ANI calculations. [default: 150] \n",
        ))
        .option(Opt::new("INT").long("--depth-per-sample-filter").help(
            "Minimum depth of a variant in a sample for that \
                     sample to be included in ANI & Fst calculations for that \
                     variant. [default: 5] \n",
        ))
        .option(Opt::new("INT").long("--min-long-read-size").help(
            "The minimum size for long reads to be used for analysis \
                    [default: 1500] \n",
        ))
        .option(
            Opt::new("INT")
                .long("--min-long-read-average-base-qual")
                .help(
                    "The minimum average base quality of a long read \
                     for it to be used for analysis [default: 20] \n",
                ),
        )
        .option(Opt::new("INT").short("-q").long("--min-base-quality").help(
            "Minimum base quality required to consider a \
                     base for calling. [default: 10] \n",
        ))
        .option(Opt::new("INT").long("--min-mapq").help(
            "Minimum MAPQ score for reads to be considered \
                     during variant calling. [default: 20] \n",
        ))
        .option(Opt::new("INT").long("--base-quality-score-threshold").help(
            "Base qualities below this threshold will \
                     be reduced to the minimum (6). [default: 18] \n",
        ))
        .option(Opt::new("INT").long("--max-input-depth").help(
            "The maximum number of reads included within an \
                     assembly region across all samples. Larger numbers \
                     increase run time. If the depth of an assembly region \
                     exceeds this value, then the reads will be filtered \
                     by mean base quality. [default: 200000] \n",
        ))
        .option(Opt::new("INT").long("--min-contig-size").help(
            "The minimum contig size to call variants on. Smaller \
                    contigs can often contain highly variable regions that \
                    mostly represent noise. Call variants on them can often \
                    be slow and not produce anything fruitful. [default: 0] \n",
        ))
        .option(Opt::new("INT").long("--min-sv-qual").help(
            "Minimum structural variants quality returned by svim \
                     and used by lorikeet. Not PHRED-scaled quality, value \
                     determined by number of supporting reads. Consult \
                     svim documentation for details. [default: 3] \n",
        ))
        .flag(Flag::new().long("--do-not-call-svs").help(
            "Opts not to use svim to call structural variants \
                     using provided longreads. If no longreads are provided \
                     this has no effect. \n",
        ))
}

fn variant_calling_options_advanced() -> Section {
    Section::new("Variant calling options (Advanced)")
        .option(
            Opt::new("INT")
                .long("--phred-scaled-global-read-mismapping-rate")
                .help("The global assumed mismapping rate for reads. [default: 45] \n"),
        )
        .option(
            Opt::new("INT")
                .long("--pair-hmm-gap-continuation-penalty")
                .help(
                    "Flat gap continuation penalty for use in the Pair HMM. \
                    [default: 10] \n",
                ),
        )
        .option(
            Opt::new("STR")
                .long("--pcr-indel-model")
                .help("The PCR indel model to use. [default: conservative] \n"),
        )
        .option(Opt::new("FLOAT").long("--heterozygosity").help(
            "Heterozygosity value used to compute prior \
                     likelihoods for any locus. [default: 0.001] \n",
        ))
        .option(Opt::new("FLOAT").long("--heterozygosity-stdev").help(
            "Standard deviation of heterozygosity for SNP and \
                     indel calling. [default: 0.01] \n",
        ))
        .option(
            Opt::new("FLOAT")
                .long("--indel-heterozygosity")
                .help("Heterozygosity for indel calling. [default: 0.000125] \n"),
        )
        .option(
            Opt::new("FLOAT")
                .long("--standard-min-confidence-threshold-for-calling")
                .short("-C")
                .help(
                    "The minimum phred-scaled confidence threshold at \
                     which variants should be called. [default: 25.0] \n",
                ),
        )
        .flag(Flag::new().long("--use-posteriors-to-calculate-qual").help(
            "if available, use the genotype posterior \
                     probabilities to calculate the site QUAL. \n",
        ))
        .flag(
            Flag::new()
                .long("--annotate-with-num-discovered-alleles")
                .help(
                    "If provided, we will annotate records with the \
                     number of alternate alleles that were discovered \
                     (but not necessarily genotyped) at a given site. \n",
                ),
        )
        .option(
            Opt::new("FLOAT")
                .long("--active-probability-threshold")
                .help(
                    "Minimum probability for a locus to be \
                     considered active. [default: 0.002] \n",
                ),
        )
        .option(
            Opt::new("INT")
                .long("--min-assembly-region-size")
                .help("Minimum size of an assembly region. [default: 50] \n"),
        )
        .option(
            Opt::new("INT")
                .long("--max-assembly-region-size")
                .help("Maximum size of an assembly region. [default: 300] \n"),
        )
        .option(Opt::new("INT").long("--assembly-region-padding").help(
            "Number of additional bases of context to \
                     include around each assembly region. [default: 100] \n",
        ))
        .flag(
            Flag::new()
                .long("--dont-increase-kmer-sizes-for-cycles")
                .help(
                    "Disable iterating over kmer sizes when \
                     graph cycles are detected. \n",
                ),
        )
        .flag(
            Flag::new()
                .long("--allow-non-unique-kmers-in-ref")
                .help(
                    "Allow graphs that have non-unique kmers in the reference. \
                    Can increase runtime and sensitivity in repeat regions. \n"),
        )
        .flag(
            Flag::new()
                .long("--do-not-run-physical-phasing")
                .help("Disable physical phasing. \n"),
        )
        .flag(
            Flag::new()
                .long("--recover-all-dangling-branches")
                .help("Recover all dangling branches. \n"),
        )
        .option(Opt::new("INT").long("--min-dangling-branch-length").help(
            "Minimum length of a dangling branch to \
                     attempt recovery. [default: 1] \n",
        ))
        .option(Opt::new("INT").long("--min-prune-factor").help(
            "Minimum read support to not prune paths in the graph. \
            Lower values like 1 or 0 will drastically increase runtime, but recover more variants. \
            For faster performance, set this to 2 or 3 at a cost to recall. \
                     [default: 1] \n",
        ))
        .flag(Flag::new().long("--use-adaptive-pruning").help(
            "Use more advanced pruning algorithm to prune paths in
                     graph. Better suited when performing variant calling
                     on when depth along a genome is variable e.g. RNA
                     and exome data. \n",
        ))
        .option(Opt::new("INT").long("--num-pruning-samples").help(
            "Number of samples that must pass the \
                     min_pruning threshold [default: 1] \n",
        ))
        .option(
            Opt::new("PATH")
                .long("--graph-output")
                .help("Write debug assembly graph information to this file. \n"),
        )
        .flag(
            Flag::new()
                .long("--dont-use-soft-clipped-bases")
                .help("Do not analyse soft clipped bases in the reads. \n"),
        )
        .option(
            Opt::new("FLOAT")
                .long("--initial-error-rate-for-pruning")
                .help(
                    "Initial base error rate estimate for adaptive \
                     pruning. [default: 0.001] \n",
                ),
        )
        .option(Opt::new("FLOAT").long("--pruning-log-odds-threshold").help(
            "Likelihood ratio threshold for adaptive \
                     pruning algorithm. This value will be converted to \
                     log odds value. [default: 1.0] \n",
        ))
        .option(Opt::new("INT").long("--max-unpruned-variants").help(
            "Maximum number of variants in graph the \
                     adaptive pruner will allow. [default: 100] \n",
        ))
        .option(
            Opt::new("INT")
                .long("--max-prob-propagation-distance")
                .help(
                    "Upper limit on how many bases away probability mass \
                     can be moved around when calculating the boundaries \
                     between active and inactive assembly regions. [default: 50] \n",
                ),
        )
        .option(Opt::new("INT").long("--max-mnp-distance").help(
            "Two or more phased substitutions separated by \
                     this distance or less are merged into MNPs. [default: 0] \n",
        ))
        .flag(
            Flag::new()
                .long("--disable-optimizations")
                .help("Don't skip calculations in ActiveRegions with no variants \n"),
        )
        .flag(Flag::new().long("--disable-avx").help(
            "Disable the use of the GKL-rs AVX acceleration components \
                     for PairHMM and Smith-Waterman calculations. \n",
        ))
        .option(Opt::new("STR").long("--limiting-interval").help(
            "Mainly used for debugging purposes. Only call variants \
                     within this given span on all contigs. E.g. providing \
                     '1000-2000' would only call variants between the 1000 \
                     and 2000 bp span on each provided contig. \n",
        ))
        .flag(
            Flag::new()
                .long("--force")
                .help("Forcefully overwrite previous runs. \n"),
        )
}

// fn add_verbosity_flags_to_section(section: Section) -> Section {
//     section
//         .flag(
//             Flag::new()
//                 .short("-v")
//                 .long("--verbose")
//                 .help("Print extra debugging information. [default: not set]"),
//         )
//         .flag(Flag::new().short("-q").long("--quiet").help(
//             "Unless there is an error, do not print \
//     log messages. [default: not set] \n",
//         ))
// }

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

    manual = manual.custom(threads_options());
    manual = manual.custom(reference_options());
    manual = manual.custom(read_mapping_params_section());
    manual = manual.custom(sharding_section());
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(variant_calling_section_basic());
    manual = manual.custom(variant_calling_options_advanced());
    manual = manual.custom(
        Section::new("Output options")
            .option(
                Opt::new("DIRECTORY")
                    .short("-o")
                    .long("--output-directory")
                    .help(
                        "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./]",
                    ),
            )
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used] \n",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--keep-unmapped")
                    .help("Include unmapped reads from cached BAM files. [default: not set] \n"),
            ),
    );

    manual = manual.example(
        Example::new()
            .text("Map paired reads to a reference and generate genotypes")
            .command(
                "lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 17 25 51",
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
            ============================\n
            "
        );

    manual = manual.custom(threads_options());
    manual = manual.custom(reference_options());
    manual = manual.custom(read_mapping_params_section());
    manual = manual.custom(sharding_section());
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(variant_calling_section_basic());
    manual = manual.custom(variant_calling_options_advanced());
    manual = manual.custom(
        Section::new("Output options")
            .option(
                Opt::new("DIRECTORY")
                    .short("-o")
                    .long("--output-directory")
                    .help(
                        "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./]",
                    ),
            )
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used] \n",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--keep-unmapped")
                    .help("Include unmapped reads from cached BAM files. [default: not set]"),
            ),
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

    manual = manual.custom(threads_options());
    manual = manual.custom(reference_options());
    manual = manual.custom(read_mapping_params_section());
    manual = manual.custom(sharding_section());
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = manual.custom(variant_calling_section_basic());
    manual = manual.custom(variant_calling_options_advanced());
    manual = manual.custom(
        Section::new("Output options")
            .option(
                Opt::new("DIRECTORY")
                    .short("-o")
                    .long("--output-directory")
                    .help(
                        "Output directory. Folder will contain subfolders for each input genome \n
                [default: ./] \n",
                    ),
            )
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used] \n",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--Keep-unmapped")
                    .help("Include unmapped reads from cached BAM files. [default: not set]"),
            ),
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
            lorikeet summarise uses a set of VCF files as input and calculates conANI, popANI,
            subpopANI, and Fst metrics for the variants in each file. \n
            \n
            ANI metrics require coverage information to determine the number of shared bases
            in each sample. VCF files do not provide this information, so the shared base size is just
            the total size of the genome. In our experience, this doesn't really matter that much as
            the ANI metrics are quite insensitive when provided low coverage samples anyway.
            Fst tends to perform better for low and high coverage samples and does not require whole
            genome coverage information.
            ============================\n
            "
        );

    manual = manual
        .option(
            Opt::new("PATH ..")
                .short("-i")
                .long("--vcfs")
                .help("Paths to input VCF files. Can provide one or more. \n"),
        )
        .option(
            Opt::new("DIRECTORY")
                .short("-i")
                .long("--vcfs")
                .help("Paths to input VCF files. Can provide one or more. \n"),
        )
        .option(Opt::new("DIRECTORY").short("-o").long("--output").help(
            "Output directory. Folder will contain subfolders for each input VCF \n
             [default: ./] \n",
        ))
        .option(
            Opt::new("INT")
                .long("--threads")
                .short("-t")
                .help("Maximum number of threads used. [default: 8] \n"),
        )
        .option(Opt::new("INT").long("--qual-by-depth-filter").help(
            "The minimum QD value for a variant to have for it to be \
                     included in the genotyping or ANI analyses. [default: 25] \n",
        ))
        .option(Opt::new("INT").long("--qual-threshold").help(
            "The PHRED-scaled quality score threshold for use \
                     with ANI calculations. [default: 150] \n",
        ))
        .option(Opt::new("INT").long("--depth-per-sample-filter").help(
            "Minimum depth of a variant in a sample for that \
                     sample to be included in ANI & Fst calculations for that \
                     variant. [default: 5] \n",
        ));

    manual = add_verbosity_flags(manual);
    return manual;
}

pub fn build_cli() -> Command {
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

  lorikeet call --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 17 25

{}

  lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10 --kmer-sizes 17 25

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

    return Command::new("lorikeet")
        .version(crate_version!())
        .author(crate::AUTHOR_AND_EMAIL)
        .about("Variant analysis of metagenomic datasets")
        .args(&[
            arg!(-v --verbose "Print extra debug logging information"),
            arg!(-q --quiet "Unless there is an error, do not print logging information"),
        ])
        .override_help(
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
        .arg_required_else_help(true)
        .subcommand(
            Command::new("genotype")
                .about("Perform variant calling analysis and then binning")
                .override_help(GENOTYPE_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("sharded")
                        .long("sharded")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-genomes-from-deshard")
                        .long("exclude-genomes-from-deshard")
                        .requires("sharded"),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .long("read1")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .long("read2")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("longreads")
                        .long("longreads")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("longread-bam-files")
                        .short('l').long("longread-bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("genome-fasta-files")
                        .short('f')
                        .short_alias('r')
                        .alias("reference")
                        .long("genome-fasta-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&["genome-fasta-directory", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .required_unless_present_any(&["genome-fasta-files", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
                        .default_value("fna"),
                )
                .arg(
                    Arg::new("bam-file-cache-directory")
                        .long("bam-file-cache-directory"),
                )
                .arg(
                    Arg::new("output-directory")
                        .long("output-directory")
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::new("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .required(false),
                )
                .arg(
                    Arg::new("threads")
                        .short('t').long("threads")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("parallel-genomes")
                        .short('p').long("parallel-genomes")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("longread-mapper")
                        .long("longread-mapper")
                        .value_parser(LONGREAD_MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("high-memory")
                        .long("high-memory")
                        .hide(true),
                )
                .arg(
                    Arg::new("keep-unmapped")
                        .long("keep-unmapped")
                        .action(clap::ArgAction::SetTrue)
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::new("split-bams").long("split-bams").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-covered-fraction")
                        .long("min-covered-fraction")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-contig-size")
                        .long("min-contig-size")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("45"),
                )
                .arg(
                    Arg::new("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .value_parser(vec![
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
                    Arg::new("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.01"),
                )
                .arg(
                    Arg::new("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::new("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short('C')
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .value_parser(vec![
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hide(true),
                )
                .arg(
                    Arg::new("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("active-probability-threshold")
                        .long("active-probability-threshold")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.002"),
                )
                .arg(
                    Arg::new("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("300"),
                )
                .arg(
                    Arg::new("kmer-sizes")
                        .long("kmer-sizes")
                        .short('k')
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .value_parser(clap::value_parser!(usize))
                        .default_values(&["21", "33"]),
                )
                .arg(
                    Arg::new("disable-automatic-kmer-adjustment")
                        .long("disable-automatic-kmer-adjustment")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("128")
                        .hide(true),
                )
                .arg(
                    Arg::new("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::new("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hide(true),
                )
                .arg(
                    Arg::new("num-pruning-samples")
                        .long("num-pruning-samples")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("min-prune-factor")
                        .long("min-prune-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("disable-prune-factor-correction")
                        .long("disable-prune-factor-correction")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("use-adaptive-pruning").long("use-adaptive-pruning").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("4.0")
                        .hide(true),
                )
                .arg(
                    Arg::new("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0"),
                )
                .arg(
                    Arg::new("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("max-input-depth")
                        .long("max-input-depth")
                        .short('i')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("200000"),
                )
                .arg(
                    Arg::new("min-variant-depth-for-genotyping")
                        .long("min-variant-depth-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("error-correct-reads")
                        .long("error-correct-reads")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("min-observations-for-kmers-to-be-solid")
                        .long("min-observations-for-kmers-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-mnp-distance")
                        .long("max-mnp-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("-1")
                        .hide(true),
                )
                .arg(
                    Arg::new("assembly-region-padding")
                        .long("assembly-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("trim-min")
                        .long("trim-min")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.00"),
                )
                .arg(
                    Arg::new("trim-max")
                        .long("trim-max")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.00"),
                )
                .arg(
                    Arg::new("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-sv-qual")
                        .long("min-sv-qual")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("3"),
                )
                .arg(Arg::new("do-not-call-svs").long("do-not-call-svs").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-long-read-size")
                        .long("min-long-read-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1500"),
                )
                .arg(
                    Arg::new("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("18"),
                )
                .arg(
                    Arg::new("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("qual-threshold")
                        .long("qual-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("150.0"),
                )
                .arg(
                    Arg::new("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .value_parser(clap::value_parser!(i64))
                        .default_value("5"),
                )
                .arg(
                    Arg::new("disable-dynamic-read-disqualification-for-genotyping")
                        .long("disable-dynamic-read-disqualification-for-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(false),
                )
                .arg(
                    Arg::new("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0")
                        .hide(false),
                )
                .arg(
                    Arg::new("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.02")
                        .hide(false),
                )
                .arg(
                    Arg::new("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2")
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(Arg::new("disable-optimizations").long("disable-optimizations").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("disable-avx").long("disable-avx").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("no-zeros").long("no-zeros").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("allow-improper-pairs").long("allow-improper-pairs").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("include-secondary").long("include-secondary").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("exclude-supplementary").long("exclude-supplementary").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("ploidy")
                        .long("ploidy")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2"),
                )
                .arg(
                    Arg::new("calculate-dnds")
                        .long("calculate-dnds")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("calculate-fst")
                        .long("calculate-fst")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("prodigal-params")
                        .long("prodigal-params")
                        .default_value("-p meta"),
                )
                .arg(
                    Arg::new("limiting-interval")
                        .long("limiting-interval")
                        .required(false),
                )
                .arg(
                    Arg::new("profile")
                        .long("profile")
                        .value_parser(["fast", "very-fast", "sensitive", "precise", "super-sensitive"])
                        .required(false)
                )
                .arg(Arg::new("force").long("force").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("verbose").short('v').long("verbose").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("quiet").long("quiet").action(clap::ArgAction::SetTrue)),
        )
        .subcommand(
            Command::new("call")
                .about("Perform variant calling across the given genomes and samples")
                .override_help(CALL_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("sharded")
                        .long("sharded")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-genomes-from-deshard")
                        .long("exclude-genomes-from-deshard")
                        .requires("sharded"),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .long("read1")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .long("read2")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("longreads")
                        .long("longreads")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("longread-bam-files")
                        .short('l').long("longread-bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("genome-fasta-files")
                        .short('f')
                        .short_alias('r')
                        .alias("reference")
                        .long("genome-fasta-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&["genome-fasta-directory", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .required_unless_present_any(&["genome-fasta-files", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
                        .default_value("fna"),
                )
                .arg(
                    Arg::new("bam-file-cache-directory")
                        .long("bam-file-cache-directory"),
                )
                .arg(
                    Arg::new("output-directory")
                        .long("output-directory")
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::new("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .required(false),
                )
                .arg(
                    Arg::new("threads")
                        .short('t').long("threads")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("parallel-genomes")
                        .short('p').long("parallel-genomes")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("longread-mapper")
                        .long("longread-mapper")
                        .value_parser(LONGREAD_MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("high-memory")
                        .long("high-memory")
                        .hide(true),
                )
                .arg(
                    Arg::new("keep-unmapped")
                        .long("keep-unmapped")
                        .action(clap::ArgAction::SetTrue)
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::new("split-bams").long("split-bams").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-covered-fraction")
                        .long("min-covered-fraction")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-contig-size")
                        .long("min-contig-size")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("45"),
                )
                .arg(
                    Arg::new("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .value_parser(vec![
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
                    Arg::new("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.01"),
                )
                .arg(
                    Arg::new("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::new("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short('C')
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .value_parser(vec![
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hide(true),
                )
                .arg(
                    Arg::new("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("active-probability-threshold")
                        .long("active-probability-threshold")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.002"),
                )
                .arg(
                    Arg::new("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("300"),
                )
                .arg(
                    Arg::new("kmer-sizes")
                        .long("kmer-sizes")
                        .short('k')
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .value_parser(clap::value_parser!(usize))
                        .default_values(&["21", "33"]),
                )
                .arg(
                    Arg::new("disable-automatic-kmer-adjustment")
                        .long("disable-automatic-kmer-adjustment")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("128")
                        .hide(true),
                )
                .arg(
                    Arg::new("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::new("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hide(true),
                )
                .arg(
                    Arg::new("num-pruning-samples")
                        .long("num-pruning-samples")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("min-prune-factor")
                        .long("min-prune-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("disable-prune-factor-correction")
                        .long("disable-prune-factor-correction")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("use-adaptive-pruning").long("use-adaptive-pruning").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("4.0")
                        .hide(true),
                )
                .arg(
                    Arg::new("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0"),
                )
                .arg(
                    Arg::new("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("max-input-depth")
                        .long("max-input-depth")
                        .short('i')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("200000"),
                )
                .arg(
                    Arg::new("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("error-correct-reads")
                        .long("error-correct-reads")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("min-observations-for-kmers-to-be-solid")
                        .long("min-observations-for-kmers-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-mnp-distance")
                        .long("max-mnp-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("-1")
                        .hide(true),
                )
                .arg(
                    Arg::new("assembly-region-padding")
                        .long("assembly-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("trim-min")
                        .long("trim-min")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.00"),
                )
                .arg(
                    Arg::new("trim-max")
                        .long("trim-max")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.00"),
                )
                .arg(
                    Arg::new("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-sv-qual")
                        .long("min-sv-qual")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("3"),
                )
                .arg(Arg::new("do-not-call-svs").long("do-not-call-svs").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-long-read-size")
                        .long("min-long-read-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1500"),
                )
                .arg(
                    Arg::new("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("18"),
                )
                .arg(
                    Arg::new("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("qual-threshold")
                        .long("qual-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("150.0"),
                )
                .arg(
                    Arg::new("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .value_parser(clap::value_parser!(i64))
                        .default_value("5"),
                )
                .arg(
                    Arg::new("disable-dynamic-read-disqualification-for-genotyping")
                        .long("disable-dynamic-read-disqualification-for-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(false),
                )
                .arg(
                    Arg::new("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0")
                        .hide(false),
                )
                .arg(
                    Arg::new("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.02")
                        .hide(false),
                )
                .arg(
                    Arg::new("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2")
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(Arg::new("disable-optimizations").long("disable-optimizations").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("disable-avx").long("disable-avx").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("no-zeros").long("no-zeros").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("allow-improper-pairs").long("allow-improper-pairs").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("include-secondary").long("include-secondary").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("exclude-supplementary").long("exclude-supplementary").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("ploidy")
                        .long("ploidy")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2"),
                )
                .arg(
                    Arg::new("calculate-dnds")
                        .long("calculate-dnds")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("calculate-fst")
                        .long("calculate-fst")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("prodigal-params")
                        .long("prodigal-params")
                        .default_value("-p meta"),
                )
                .arg(
                    Arg::new("limiting-interval")
                        .long("limiting-interval")
                        .required(false),
                )
                .arg(
                    Arg::new("profile")
                        .long("profile")
                        .value_parser(["fast", "very-fast", "sensitive", "precise", "super-sensitive"])
                        .required(false)
                )
                .arg(Arg::new("force").long("force").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("verbose").short('v').long("verbose").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("quiet").long("quiet").action(clap::ArgAction::SetTrue)),
        )
        .subcommand(
            Command::new("consensus")
                .about("Generate consensus genomes across all provided samples")
                .override_help(CONSENSUS_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("sharded")
                        .long("sharded")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-genomes-from-deshard")
                        .long("exclude-genomes-from-deshard")
                        .requires("sharded"),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .long("read1")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .long("read2")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any(&[
                            "bam-files",
                            "single",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "single",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "longreads",
                            "longread-bam-files",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("longreads")
                        .long("longreads")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("longread-bam-files")
                        .short('l').long("longread-bam-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&[
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
                    Arg::new("genome-fasta-files")
                        .short('f')
                        .short_alias('r')
                        .alias("reference")
                        .long("genome-fasta-files")
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&["genome-fasta-directory", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .short('d')
                        .required_unless_present_any(&["genome-fasta-files", "full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short('x')
                        .default_value("fna"),
                )
                .arg(
                    Arg::new("bam-file-cache-directory")
                        .long("bam-file-cache-directory"),
                )
                .arg(
                    Arg::new("output-directory")
                        .long("output-directory")
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::new("features-vcf")
                        .long("features-vcf")
                        .short('f')
                        .required(false),
                )
                .arg(
                    Arg::new("threads")
                        .short('t').long("threads")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("parallel-genomes")
                        .short('p').long("parallel-genomes")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("longread-mapper")
                        .long("longread-mapper")
                        .value_parser(LONGREAD_MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index"),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("high-memory")
                        .long("high-memory")
                        .hide(true),
                )
                .arg(
                    Arg::new("keep-unmapped")
                        .long("keep-unmapped")
                        .action(clap::ArgAction::SetTrue)
                        .requires("bam-file-cache-directory"),
                )
                .arg(Arg::new("split-bams").long("split-bams").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .conflicts_with("allow-improper-pairs"),
                )
                .arg(
                    Arg::new("min-covered-fraction")
                        .long("min-covered-fraction")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.0"),
                )
                .arg(
                    Arg::new("min-contig-size")
                        .long("min-contig-size")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("phred-scaled-global-read-mismapping-rate")
                        .long("phred-scaled-global-read-mismapping-rate")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("45"),
                )
                .arg(
                    Arg::new("pair-hmm-gap-continuation-penalty")
                        .long("pair-hmm-gap-continuation-penalty")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("pcr-indel-model")
                        .long("pcr-indel-model")
                        .default_value("conservative")
                        .value_parser(vec![
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
                    Arg::new("heterozygosity-stdev")
                        .long("heterozygosity-stdev")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.01"),
                )
                .arg(
                    Arg::new("snp-heterozygosity")
                        .long("snp-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.000125"),
                )
                .arg(
                    Arg::new("standard-min-confidence-threshold-for-calling")
                        .long("standard-min-confidence-threshold-for-calling")
                        .short('C')
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("genotype-assignment-method")
                        .long("genotype-assignment-method")
                        .default_value("UsePLsToAssign")
                        .value_parser(vec![
                            "UsePLsToAssign",
                            "UsePosteriorProbabilities",
                            "BestMatchToOriginal",
                            "DoNotAssignGenotypes",
                        ])
                        .hide(true),
                )
                .arg(
                    Arg::new("use-posteriors-to-calculate-qual")
                        .long("use-posteriors-to-calculate-qual")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("annotate-with-num-discovered-alleles")
                        .long("annotate-with-num-discovered-alleles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("active-probability-threshold")
                        .long("active-probability-threshold")
                        .value_parser(clap::value_parser!(f32))
                        .default_value("0.002"),
                )
                .arg(
                    Arg::new("min-assembly-region-size")
                        .long("min-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("max-assembly-region-size")
                        .long("max-assembly-region-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("300"),
                )
                .arg(
                    Arg::new("kmer-sizes")
                        .long("kmer-sizes")
                        .short('k')
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .value_parser(clap::value_parser!(usize))
                        .default_values(&["21", "33"]),
                )
                .arg(
                    Arg::new("disable-automatic-kmer-adjustment")
                        .long("disable-automatic-kmer-adjustment")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("max-allowed-path-for-read-threading-assembler")
                        .long("max-allowed-path-for-read-threading-assembler")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("128")
                        .hide(true),
                )
                .arg(
                    Arg::new("dont-increase-kmer-sizes-for-cycles")
                        .long("dont-increase-kmer-sizes-for-cycles")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("allow-non-unique-kmers-in-ref")
                        .long("allow-non-unique-kmers-in-ref")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("debug-graph-transformations")
                        .long("debug-graph-transformations")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-recover-dangling-branches")
                        .long("do-not-recover-dangling-branches")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("do-not-run-physical-phasing")
                        .long("do-not-run-physical-phasing")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("recover-all-dangling-branches")
                        .long("recover-all-dangling-branches")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-dangling-branch-length")
                        .long("min-dangling-branch-length")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("graph-output")
                        .long("graph-output")
                        .default_value("lorikeet_haplotype_caller"),
                )
                .arg(
                    Arg::new("debug-graph-output")
                        .long("debug-graph-output")
                        .default_value("lorikeet_haplotype_caller_debug")
                        .hide(true),
                )
                .arg(
                    Arg::new("num-pruning-samples")
                        .long("num-pruning-samples")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("min-prune-factor")
                        .long("min-prune-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("disable-prune-factor-correction")
                        .long("disable-prune-factor-correction")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("use-adaptive-pruning").long("use-adaptive-pruning").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("dont-use-soft-clipped-bases")
                        .long("dont-use-soft-clipped-bases")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("initial-error-rate-for-pruning")
                        .long("initial-error-rate-for-pruning")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.001"),
                )
                .arg(
                    Arg::new("pruning-seeding-log-odds-threshold")
                        .long("pruning-seeding-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("4.0")
                        .hide(true),
                )
                .arg(
                    Arg::new("pruning-log-odds-threshold")
                        .long("pruning-log-odds-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0"),
                )
                .arg(
                    Arg::new("max-unpruned-variants")
                        .long("max-unpruned-variants")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("max-input-depth")
                        .long("max-input-depth")
                        .short('i')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("200000"),
                )
                .arg(
                    Arg::new("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("max-prob-propagation-distance")
                        .long("max-prob-propagation-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("50"),
                )
                .arg(
                    Arg::new("use-linked-debruijn-graph")
                        .long("use-linked-debruijn-graph")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("error-correct-reads")
                        .long("error-correct-reads")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("kmer-length-for-read-error-correction")
                        .long("kmer-length-for-read-error-correction")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("min-observations-for-kmers-to-be-solid")
                        .long("min-observations-for-kmers-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-mnp-distance")
                        .long("max-mnp-distance")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0"),
                )
                .arg(
                    Arg::new("min-observation-for-kmer-to-be-solid")
                        .long("min-observation-for-kmer-to-be-solid")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("enable-legacy-graph-cycle-detection")
                        .long("enable-legacy-graph-cycle-detection")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("min-matching-bases-to-dangling-end-recovery")
                        .long("min-matching-bases-to-dangling-end-recovery")
                        .value_parser(clap::value_parser!(i32))
                        .default_value("-1")
                        .hide(true),
                )
                .arg(
                    Arg::new("assembly-region-padding")
                        .long("assembly-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("indel-padding-for-genotyping")
                        .long("indel-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("str-padding-for-genotyping")
                        .long("str-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("75")
                        .hide(true),
                )
                .arg(
                    Arg::new("snp-padding-for-genotyping")
                        .long("snp-padding-for-genotyping")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20")
                        .hide(true),
                )
                .arg(
                    Arg::new("max-extension-into-region-padding")
                        .long("max-extension-into-region-padding")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("25")
                        .hide(true),
                )
                .arg(
                    Arg::new("soft-clip-low-quality-ends")
                        .long("soft-clip-low-quality-ends")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("trim-min")
                        .long("trim-min")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.00"),
                )
                .arg(
                    Arg::new("trim-max")
                        .long("trim-max")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.00"),
                )
                .arg(
                    Arg::new("mapping-quality-threshold-for-genotyping")
                        .long("mapping-quality-threshold-for-genotyping")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-sv-qual")
                        .long("min-sv-qual")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("3"),
                )
                .arg(Arg::new("do-not-call-svs").long("do-not-call-svs").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-long-read-size")
                        .long("min-long-read-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("1500"),
                )
                .arg(
                    Arg::new("min-long-read-average-base-qual")
                        .long("min-long-read-average-base-qual")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("20"),
                )
                .arg(
                    Arg::new("min-base-quality")
                        .long("min-base-quality")
                        .short('q')
                        .value_parser(clap::value_parser!(u8))
                        .default_value("10"),
                )
                .arg(
                    Arg::new("base-quality-score-threshold")
                        .long("base-quality-score-threshold")
                        .value_parser(clap::value_parser!(u8))
                        .default_value("18"),
                )
                .arg(
                    Arg::new("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("qual-threshold")
                        .long("qual-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("150.0"),
                )
                .arg(
                    Arg::new("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .value_parser(clap::value_parser!(i64))
                        .default_value("5"),
                )
                .arg(
                    Arg::new("disable-dynamic-read-disqualification-for-genotyping")
                        .long("disable-dynamic-read-disqualification-for-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(false),
                )
                .arg(
                    Arg::new("dynamic-read-disqualification-threshold")
                        .long("dynamic-read-disqualification-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("1.0")
                        .hide(false),
                )
                .arg(
                    Arg::new("expected-mismatch-rate-for-read-disqualification")
                        .long("expected-mismatch-rate-for-read-disqualification")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("0.02")
                        .hide(false),
                )
                .arg(
                    Arg::new("allele-informative-reads-overlap-margin")
                        .long("allele-informative-reads-overlap-margin")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2")
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-symmetric-hmm-normalizing")
                        .long("disable-symmetric-hmm-normalizing")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-cap-base-qualities-to-map-quality")
                        .long("disable-cap-base-qualities-to-map-quality")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(
                    Arg::new("disable-spanning-event-genotyping")
                        .long("disable-spanning-event-genotyping")
                        .action(clap::ArgAction::SetTrue)
                        .hide(true),
                )
                .arg(Arg::new("disable-optimizations").long("disable-optimizations").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("disable-avx").long("disable-avx").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("no-zeros").long("no-zeros").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("allow-improper-pairs").long("allow-improper-pairs").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("include-secondary").long("include-secondary").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("exclude-supplementary").long("exclude-supplementary").action(clap::ArgAction::SetTrue))
                .arg(
                    Arg::new("ploidy")
                        .long("ploidy")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("2"),
                )
                .arg(
                    Arg::new("calculate-dnds")
                        .long("calculate-dnds")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("calculate-fst")
                        .long("calculate-fst")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("prodigal-params")
                        .long("prodigal-params")
                        .default_value("-p meta"),
                )
                .arg(
                    Arg::new("limiting-interval")
                        .long("limiting-interval")
                        .required(false),
                )
                .arg(
                    Arg::new("profile")
                        .long("profile")
                        .value_parser(["fast", "very-fast", "sensitive", "precise", "super-sensitive"])
                        .required(false)
                )
                .arg(Arg::new("force").long("force").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("verbose").short('v').long("verbose").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("quiet").long("quiet").action(clap::ArgAction::SetTrue)),
        )
        .subcommand(
            Command::new("summarise")
                .about("Summarizes ANI values of a given set of VCF files")
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("vcfs")
                        .long("vcfs")
                        .short('i')
                        .action(ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(&["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .default_value("./"),
                )
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .short('t')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("8"),
                )
                .arg(
                    Arg::new("qual-by-depth-filter")
                        .long("qual-by-depth-filter")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("25.0"),
                )
                .arg(
                    Arg::new("qual-threshold")
                        .long("qual-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .default_value("150.0"),
                )
                .arg(
                    Arg::new("depth-per-sample-filter")
                        .long("depth-per-sample-filter")
                        .value_parser(clap::value_parser!(i64))
                        .default_value("5"),
                )
                .arg(Arg::new("verbose").short('v').long("verbose").action(ArgAction::SetTrue)),
        )
        .subcommand(
            add_clap_verbosity_flags(Command::new("shell-completion"))
                .about("Generate a shell completion script for lorikeet")
                .arg(
                    Arg::new("output-file")
                        .short('o')
                        .long("output-file")
                        .required(true),
                )
                .arg(
                    Arg::new("shell")
                        .long("shell")
                        .required(true)
                        .value_parser(value_parser!(Shell)),
                ),
        );
}
