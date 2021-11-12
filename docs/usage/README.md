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

As a simple example, imagine we have a single sample where the reads have previously mapped to our metagenome or set
of references using [CoverM](https://github.com/wwood/coverm) or Lorikeet:

```
lorikeet call --bam-files my.bam --genome-fasta-directory genomes/ -x fna --output-directory lorikeet_out/ --threads 10
```

Call variants from short reads and longread bam files:

`lorikeet call -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam`