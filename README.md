[![Build Status](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)](https://travis-ci.com/rhysnewell/Lorikeet)
[![codecov](https://codecov.io/gh/rhysnewell/Lorikeet/branch/master/graph/badge.svg)](https://codecov.io/gh/rhysnewell/Lorikeet)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


![](docs/static/images/lorikeet_logo_crop.png)

A strain resolver for metagenomics.

## Installation
```
cargo install lorikeet
```


## Usage

Input can either be reads and reference genome, or MAG. Or a BAM file and associated genome.

```
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
    genotype    *Experimental* Resolve strain-level genotypes of MAGs from microbial communities
    polymorph   Calculate variants along contig positions
    summarize   Summarizes contig stats from multiple samples
    evolve  Calculate dN/dS values for genes from read mappings

Less used utility subcommands:
    kmer    Calculate kmer frequencies within contigs
    filter    Remove (or only keep) alignments with insufficient identity

Other options:
    -V, --version   Print version information

Rhys J. P. Newell <r.newell near uq.edu.au>
```


`lorikeet polymorph -b input.bam -r input_genome.fna`

OR

`lorikeet polymorph -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq`

# Output 
The output of polymorph is a tab-delimited file with 9 fields

`tid  pos variant reference abundance depth genotypes type  connected_bases`
