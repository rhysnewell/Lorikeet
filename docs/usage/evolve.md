---
title: Evolve
---

Evolve
========

The `evolve` subcommand of Lorikeet is used to calculate dN/dS ratios for coding regions within a given set of genomes
and samples. This method either reuses variants found via `call` or calls new variants directly.

## Help

```
lorikeet evolve --full-help
lorikeet call: Call variants using local reassembly across multiple genomes and samples


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
                                         [default "fna"]
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
   -o, --output-prefix <STRING>          Output directory prefix [default: output]



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
                                           metabat ("MetaBAT adjusted coverage")
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Require --discard-improper-pairs. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Require --discard-improper-pairs. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Require --discard-improper-pairs. [default 0.0]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 0]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.00]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 1.00]
   --discard-improper-pairs              Discard improperly mapped read pairs for variant calling.
   --discard-supplementary               Discard read alignments flagged as supplementary
   --include-secondary                   Includes read alignments flagged as secondary
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   --high-memory                         Run in high memory mode. Can be slightly faster sometimes
                                         but consumes much more RAM than standard mode.
   --split-bams                          Split the mapped read files up per reference.
                                         Useful if you think run time is being hampered
                                         by I/O. Most of the time this will not improve
                                         performance and instead just increase disk usage.


Read mapping options:
  --mapper <NAME>                 Underlying mapping software used for short reads
                                  ("minimap2-sr", "bwa-mem",
                                  "ngmlr-ont", "ngmlr-pb", "minimap2-ont",
                                  "minimap2-pb", or "minimap2-no-preset").
                                  minimap2 -sr, -ont, -pb, -no-preset specify
                                  '-x' preset of minimap2 to be used
                                  (with map-ont, map-pb for -ont, -pb).
                                  [default: "minimap2-sr"]

  --longread-mapper <NAME>        Underlying mapping software used for long reads
                                  ("minimap2-sr", "bwa-mem",
                                  "ngmlr-ont", "ngmlr-pb", "minimap2-ont",
                                  "minimap2-pb", or "minimap2-no-preset").
                                  minimap2 -sr, -ont, -pb, -no-preset specify
                                  '-x' preset of minimap2 to be used
                                  (with map-ont, map-pb for -ont, -pb).
                                  [default: "minimap2-ont"]

  --minimap2-params PARAMS        Extra parameters to provide to minimap2,
                                  both indexing command (if used) and for
                                  mapping. Note that usage of this parameter
                                  has security implications if untrusted input
                                  is specified. '-a' is always specified.
                                  [default ""]

  --minimap2-reference-is-index   Treat reference as a minimap2 database, not
                                  as a FASTA file.

  --bwa-params PARAMS             Extra parameters to provide to BWA. Note
                                  that usage of this parameter has security
                                  implications if untrusted input is specified.
                                  [default ""]

  --ngmlr-params PARAMS           Extra parameters to provide to NGMLR.
                                  --bam-fix, -x ont, -t are already set. Note
                                  that usage of this parameter has security
                                  implications if untrusted input is specified.


Variant calling options (Basic):
  -k, --kmer-sizes <INT>                        K-mer sizes used to generate DeBruijn Graphs.
                                                Multiple values at once are accepted and encouraged
                                                e.g. 10 25 [default: 25]

  --ploidy <INT>                                Sets the default ploidy for the analysis to N.
                                                [default: 1]

  -f, --features-vcf                            The set of alleles to force-call regardless
                                                of evidence. Note: The sight containing these alleles
                                                has to be called as 'active' in order for them to appear
                                                in the final VCF. Addtionally, Provided file must be
                                                compressed using bgzip and indexed using bcftools index. If no index
                                                is present, and index will be attempted to be created.
                                                If the file is not properly compressed, Lorikeet will
                                                unfortunately SEGFAULT with no error message.
  --qual-by-depth-filter                        The minimum QD value for a variant to have for it to be
                                                included in the genotyping or ANI analyses. [default: 20]
  --qual-threshold                              The PHRED-scaled quality score threshold for use
                                                with ANI calculations. [default: 150]
  --depth-per-sample-filter                     Minimum depth of a variant in a sample for that
                                                sample to be included in ANI calculations for that
                                                variant. [default: 5]
  -q, --min-base-quality                        Minimum base quality required to consider a
                                                base for calling. [default: 10]
  --min-mapq                                    Minimum MAPQ score for longreads to be considered
                                                during variant calling. [default: 60]
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
  --do-not-call-svs                             Opts not to use svim to call structural variants
                                                using provided longreads. If no longreads are provided
                                                this has no effect.
  --min-sv-qual                                 Minimum structural variants quality returned by svim
                                                and used by lorikeet. Not PHRED-scaled quality, value
                                                determined by number of supporting reads. Consult
                                                svim documentation for details. [default: 3]

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
  --min-prune-factor                            Minimum support to not prune paths in the graph.
                                                [default: 2]
  --use-adaptive-pruning                        Use more advanced pruning algorithm to prune paths in
                                                graph. Better suited when performing variant calling
                                                on when depth along a genome is variable e.g. RNA
                                                and exome data.
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
  --force                                       Forcefully overwrite previous runs.


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

```