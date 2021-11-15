---
title: Usage
---

Getting started
========

Input can either be reads and reference genome, or MAG. Or a BAM file and associated genome.

Using `lorikeet -h` will provide the following message:

```
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
    genotype    Resolve strain-level genotypes of MAGs from microbial communities
    consensus   Creates consensus genomes for each input reference and for each sample
    call        Performs variant calling with no downstream analysis
    evolve      *unavailable* Calculate dN/dS values for genes from read mappings

Other options:
    -V, --version   Print version information

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>
```

# Quick Start

As a simple example, imagine we have a single sample where the reads have previously mapped to our metagenome or set
of references using [CoverM](https://github.com/wwood/coverm) or Lorikeet:

```
lorikeet call --bam-files my.bam --genome-fasta-directory genomes/ -x fna --output-directory lorikeet_out/ --threads 10
```

One of the parts of what makes Lorikeet faster than other available metagenomic variant calling tools is that it is
capable of handling multiple reference and samples at a time. If you provide Lorikeet with multiple references and mutliple samples
it will handle the mapping of all those samples on to all of the reference for you. This is generally the slowest part of the 
algorithm as read mapping is an expensive task, as such it is recommended that you save any bams that are produced by using the
`--bam-file-cache-directory` option. That way you can reuse the BAM files if you should want to rerun the analysis. Once read mapping
is completed Lorikeet then parallelizes the entire variant calling process across all references drastically increasing performance. 
Additionally you can provide both long and short read samples to Lorikeet with ease using the associate longread flags.:

```
lorikeet call -r input_genomes/*.fna -1 forward_reads/*_1.fastq -2 reverse_reads/*_2.fastq -l longreads/*.bam --parallel-genome 8 --threads 24
```

# Output

Lorikeet will create an output for each input reference genome within the supplied output folder:
```
lorikeet_output --
                 | - Genome1
                 | - Genome2
                ...
                 | - GenomeN --
                              |
                              | - VCF
                              | - Consensus, Population, Subpopulation ANI
                              | - Optional Reconstructed Strain Genomes and Strain Relative Abundances    
```