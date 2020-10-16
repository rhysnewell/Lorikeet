![](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/license.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/version.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/platforms.svg)


![](docs/images/lorikeet_logo_crop_left.png)

## Introduction

Lorikeet is a within-species variant analysis pipeline for metagenomic communities that utilizes both long and short read datasets.
Lorikeet combines short read variant calling with Freebayes and long read structural variant calling with SVIM to generate a 
complete variant landscape present across samples within a microbial community. SNV are also linked together using read 
information to help provide likely genotypes based on observed physical linkages.

## Installation

#### Option 1: Conda - Recommended

Install into current conda environment:
```
conda install lorikeet-genome
```

Create fresh conda environment and install lorikeet there:
```
conda create -n lorikeet -c bioconda lorikeet-genome && \
conda activate lorikeet
```

#### Option 2: Cargo
```
conda create -n lorikeet -y -c conda-forge -c bioconda -c defaults -y python=3.7 parallel pysam=0.16 svim \ 
freebayes prodigal samtools=1.9 vt rust clangdev pkg-config zlib gsl starcode openblas bwa minimap2 \ 
fastani dashing r-base && \ 
conda activate lorikeet && \ 
cargo install lorikeet-genome
```

#### Option 3: Install manually
You may need to manually set the paths for `C_INCLUDE_PATH`, `LIBRARY_PATH`, `LIBCLANG_PATH`, and `OPENSSL_DIR` to their corresponding
paths in the your conda environment if they can't properly be found on your system.
```
conda create -n lorikeet -y -c conda-forge -c bioconda -c defaults -y python=3.7 parallel pysam=0.16 svim \ 
freebayes prodigal samtools=1.9 vt rust clangdev pkg-config zlib gsl starcode openblas bwa minimap2 \ 
fastani dashing r-base && \ 
conda activate lorikeet && \ 
git clone https://github.com/rhysnewell/Lorikeet.git && \ 
cd Lorikeet && \ 
bash install.sh # or e.g. `cargo run -- genotype`
```


## Usage

Input can either be reads and reference genome, or MAG. Or a BAM file and associated genome.

```
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
    genotype    *Experimental* Resolve strain-level genotypes of MAGs from microbial communities
    polish      Creates consensus genomes for each input reference and for each sample
    summarize   Summarizes contig stats from multiple samples
    evolve      Calculate dN/dS values for genes from read mappings

Less used utility subcommands:
    kmer      Calculate kmer frequencies within contigs
    filter    Remove (or only keep) alignments with insufficient identity

Other options:
    -V, --version   Print version information

Rhys J. P. Newell <r.newell near uq.edu.au>
```

Genotype from bam:

`lorikeet genotype --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
     --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10 --plot`

Genotype from short reads and longread bam:

`lorikeet genotype -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam`

## Workflow

![](docs/images/Lorikeet-workflow.png)


## Output

#### Genotype 
Genotype will produce:
- Variants that clustered into representations of the expected strain level genotypes
- Sample adjacency matrix displaying the number of shared variants seen between samples.
- VCF file detailing all observed variants
- Expected coverage values for each of the produced genotypes in each sample
- Per reference and per sample summary statistics displaying the mean SNPs and structural
  variations per provided base pair window
*optional*
- SNP density plot. Can take a long time to generate
#### Polish
Polish produces:
- Consensus genomes for each input reference across each sample.
- Sample adjacency matrix displaying the number of shared variants seen between samples.
- VCF file detailing all observed variants
- Per reference and per sample summary statistics displaying the mean SNPs and structural
  variations per provided base pair window
*optional*
- SNP density plot. Can take a long time to generate

#### Evolve
Evolve will produce:
- GFF file for each input reference with dN/dS values within coding regions based on the possible variants 
  found along the reference
    - These dN/dS values only take single nucleotide polymorphisms into account but INDELs can still be reported.
- VCF file detailing all observed variants
- Per reference and per sample summary statistics displaying the mean SNPs and structural
  variations per provided base pair window
*optional*
- SNP density plot. Can take a long time to generate

#### Summarize
- Per reference and per sample summary statistics displaying the mean SNPs and structural
  variations per provided base pair window
*optional*
- SNP density plot. Can take a long time to generate