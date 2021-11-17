---
title: Call
---

Call
========

The `call` subcommand of Lorikeet is the primary method for users who wish to just call variants across one or more
reference genomes. The variant calling process uses an algorithm based on [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) which 
calls Single Nucleotide Polymorphisms (SNPs) and Insertions and Deletions (INDELs) via local re-assembly of haplotypes.

## Method

Lorikeet first scans each input genome for regions showing signs of variation, it then discards the 
existing mapping information and completely reassembles the reads in that region. This allows the Lorikeet to be 
more accurate when calling regions that are traditionally difficult to call, for example when they contain different 
types of variants close to each other. 


### 1)	Define active regions
Variant calling focusses on regions across the input genome where there appears to be significant evidence that the 
aligned reads do not match the reference genome. Such evidence can take the form of base mismatches in the read, the 
presence of insertion of deletion markers, or soft clipped portions of the mapped reads with high base quality scores. 
Each base is assigned an active probability based on the genotype likelihoods of the reference vs. any other possible 
allele at this position (See figure 1.1). A heterozygosity prior is used for the genotype likelihood calculation, 0.001 
for SNPs and 0.0001 for INDELs. The per base probabilities are smoothed across a region via process of convolution using 
a gaussian kernel. These smoothed values are used to determine whether a contiguous span along the genome is active by 
checking to see if the per base probability exceeds a given threshold (0.002 by default).

### 2)	Perform local reassembly of active regions
If a genomic region is deemed to be active, then all the reads that mapped to that region are collected. The reads are 
then split in k-mers and reassembled using a “de Bruijn-like” graph structure using the reference sequence in this region 
as a template (See figure 1.2). The graph is initialized by building a standard de Bruijn graph from the reference sequence. 
This graph is ensured to only contain a single path by choosing a k-mer size such that there are no non-unique vertices present. 
The paths in these graphs are then search for the most likely candidate haplotypes. The reads associated with the active 
region are then threaded into the de Bruijn graph, creating new vertices as new k-mers are introduced. Edges in the graph 
are assigned weights depending on how many reads were seen connecting those two vertices. Edges that have limited k-mer 
support are systematically pruned from the graph using an adaptive chain pruning algorithm. Paths through the graph from 
the starting vertex to the terminal vertex are treated as candidate haplotypes. Paths with fewer than two supporting reads are discarded.

### 3)	Reassign reads to haplotypes and determine likelihoods
The previous step can produce a high number of candidate haplotypes (See figure 1.3), so it is necessary to filter out unlikely haplotypes. 
A pair-wise Hidden Markov Model (PairHMM) is used in order to achieve this. In the PairHMM algorithm each read is aligned 
against each haplotype and produces a likelihood score of that read being assigned to that haplotype (See figure 1.4). 
These per read scores are then used to filter out the least likely haplotypes and return the reference haplotype and any associated variants at this location.

A complete breakdown of the HaplotypeCaller algorithm can be found in their paper18. However, there are a few steps in 
the algorithm where Lorikeet differs. For instance, Lorikeet does not make use of the Reference Confidence Model devised 
by the GATK team to speed up their pipeline when thousands of samples are included. Instead, Lorikeet relies on the speed 
of the Rust programming language and takes a highly parallel approach to calculating variant calls. 


![](/figures/variant_calling_diagram.png)
Figure 1) Steps involveded in the variant calling procedure

## Help

```
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