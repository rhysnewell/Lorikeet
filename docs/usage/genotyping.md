---
title: Genotype
---

Genotype
========

The `genotype` subcommand uses the lorikeet strain genome recovery algorithm. It attempts to phase the haplotypes 
produced by the variant calling algorithm into candidate strain genomes. Additionally, Lorikeet attempts to provide relative
abundance values for the recovered strain haplotypes across the samples.

## 1) Variant embedding and clustering 
Once all the variants for a given genome have been discovered then the read depths for both the reference and variant 
alleles are collected for each sample. These depths are written into a v × 2n numpy array where v is the number of 
alternate alleles and n is the number of samples. The array is passed into the python component of Lorikeet which aims 
to cluster variants into variant groupings. To achieve this a UMAP embeddings 
is calculated for the variant depth array. The UMAP algorithm is an alternative dimensionality reduction technique 
that attempts to reduce dimensionality of a given dataset whilst maintaining global and local topology. Parameters for 
UMAP are automatically decided by Lorikeet based on the size of the variant depth array. Prior to being passed to UMAP, 
the data is first transformed using the Centred Log-Ratio transformation as the input depth is inherently a proportional 
measurement. Two UMAP embedding are initially calculated using two different distance metrics. The first makes use of the 
proportionality metric, ρ, and the second makes use of the Aitchinson distance. Since the underlying data structure 
of UMAP is a graph-based layout, we can then intersect the two UMAP embeddings and generate a finalize embedding that 
captures only the relevant connections between points seen in both original embeddings. The intersected UMAP embeddings 
are then passed to the HDBSCAN algorithm which clusters the variants in UMAP space into fine-grained variant groups. 
These variant groups as well as an array containing the pairwise density separation between each group is returned into 
the main Rust component of Lorikeet.

## 2) Variant group linking
Variant groups by themselves most likely do not constitute an entire strain haplotype. This is because strains can share 
specific variants which can alter the coverage distribution of that specific variant across samples. Therefore, variant 
groups are more likely to represent small clusters of variants that are shared between one or more strains. Thus, to reconstruct 
the strain haplotypes we need to be able to thread these variant groups back together. 

## 3) Strain abundance calculation
The final part of the Lorikeet strain recovery algorithm aims to assign a relative abundance value to each of the 
candidate strains recovered in the previous two steps. Additionally, this stage will filter out any candidate strains 
that might not actually represent an extant haplotype. Since multiple strains can share the same variant, or groups of 
variants, reads that map to those variants will be assigned to multiple strains. As such, we need a way to appropriately 
split those reads across multiple strains to get an estimate of the relative abundance. To this end, Lorikeet makes use 
of an Expectation-Maximization algorithm similar to those developed for the Centrifuge, Cufflinks, and Sailfish algorithms. 


## Help

```
                                    lorikeet genotype
              Report strain-level genotypes and abundances based on variant read mappings

Example: Map paired reads to a reference and generate genotypes

  lorikeet genotype --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10 --kmer-sizes 10 25

Example: Generate strain-level genotypes from read mappings compared to reference from a sorted BAM file and plots the results:

  lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
    --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10

See lorikeet genotype --full-help for further options and further detail.
```

## Full Help
```
lorikeet genotype: Resolves strain-level genotypes and abundance from metagenomes


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
   -o, --output-prefix <STRING>          Output prefix for files. [default: output]



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


Read mapping options:
  -p, --mapper <NAME>             Underlying mapping software used
                                  ("minimap2-sr", "bwa-mem",
                                  "ngmlr-ont", "ngmlr-pb", "minimap2-ont",
                                  "minimap2-pb", or "minimap2-no-preset").
                                  minimap2 -sr, -ont, -pb, -no-preset specify
                                  '-x' preset of minimap2 to be used
                                  (with map-ont, map-pb for -ont, -pb).
                                  [default: "minimap2-sr"]

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
                                                segfault with no error message.
  -q, --min-base-quality                        Minimum base quality required to consider a
                                                base for calling. [default: 10]
  --base-quality-score-threshold                Base qualities below this threshold will
                                                be reduced to the minimum (6). [default: 18]
  --max-input-depth                             The maximum number of reads included within an
                                                assembly region across all samples. Larger numbers
                                                increase run time. If the depth of an assembly region
                                                exceeds this value, then the reads will be filtered
                                                by mean base quality. [default: 1000]
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
  --limiting-interval                           Mainly used for debugging purposes. Only call variants
                                                within this given span on all contigs. E.g. providing
                                                '1000-2000' would only call variants between the 1000
                                                and 2000 bp span on each provided contig.
  --force                                       Forcefully overwrite previous runs.


Genotyping arguments (optional):

  --qual-by-depth-filter                The minimum QD value for a variant to have for it to be
                                        included in the genotyping analysis. [default: 20]
  --min-variant-depth-for-genotyping    The minimum total depth of a variant - across all samples -
                                        for it to be included in the strain genotyping process.
                                        Lower values tend to confuse and break the UMAP embedding
                                        and strain abundance calculation. [default: 5]

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
