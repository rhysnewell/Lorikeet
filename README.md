[![Build Status](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)](https://travis-ci.com/rhysnewell/Lorikeet)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


![](docs/static/images/lorikeet_logo_crop.png)

A strain resolver for metagenomics.

## Installation

#### Option 1: Conda - Recommended
```
conda install lorikeet-genome
```

#### Option 2: Cargo
```
conda create -n lorikeet -c conda-forge -c bioconda -c defaults -y snippy=4.4.5 svim rust clangdev pkg-config zlib gsl \ 
openblas zlib gsl starcode openblas bwa minimap2 fastani dashing && \ 
cargo install lorikeet-genome
```

#### Option 3: Build manually
You may need to manually set the paths for `C_INCLUDE_PATH`, `LIBRARY_PATH`, `LIBCLANG_PATH`, and `OPENSSL_DIR` to their corresponding
paths in the your conda environment if they can't properly be found on your system.
```
conda create -n lorikeet -c conda-forge -c bioconda -c defaults -y snippy=4.4.5 svim rust clangdev pkg-config zlib gsl \ 
openblas zlib gsl starcode openblas bwa minimap2 fastani dashing && \ 
git clone https://github.com/rhysnewell/Lorikeet/git && \ 
cd Lorikeet && cargo build --release
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

Genotype from bam:

`lorikeet genotype -b input.bam -r input_genome.fna --e-min 0.1 --e-max 0.5 --pts-min 0.1 --pts-max 0.5`

Genotype from short reads and longread bam:

`lorikeet genotype -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam`

## Output

#### Genotype 
Genotype will produce multiple .fna files representative of the expected strain level genotypes

#### Polymorph
Polymorph produces a tab delimited file containing possible variants and their positions within the reference

#### Evolve
Evolve will produce dN/dS values within coding regions based on the possible variants found along the reference.
These dN/dS values only take single nucleotide polymorphisms into account but INDELs can still be reported.