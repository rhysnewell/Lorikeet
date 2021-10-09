![](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/license.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/version.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/platforms.svg)


![](docs/images/lorikeet_logo.png)

## Lorikeet

Lorikeet is a within-species variant analysis pipeline for metagenomic communities that utilizes both long and short reads.
Lorikeet utilizes a re-implementaion of the GATK HaplotypeCaller algorithm, performing local re-assembly of potentially active
regions within candidate genomes. Called variants can be clustered into likely strains using a combination of UMAP and HDBSCAN.

## Index
1. [Installation](#installation)
    - [Conda installation](#option-1-conda)
    - [Manual installation](#option-2-install-manually)
2. [Requirements](#requirements)
3. [Usage](#usage)
4. [Workflow](#workflow)
5. [Output](#output)
    - [Genotype](#genotype)
    - [Polish](#polish)
    - [Evolve](#evolve)
    - [Summarize](#summarize)
6. [FAQs](#faqs)
7. [Resource](#resources)
8. [Citation](#citation)

## Installation

#### Option 1: Conda

*NOTE:* The conda version is often a few commits and/or versions behind the development version. If you want the most
up to date version, follow the instruction in option 2. 

Install into current conda environment:
```
conda install lorikeet-genome
```

Create fresh conda environment and install lorikeet there:
```
conda create -n lorikeet -c bioconda lorikeet-genome && \
conda activate lorikeet
```

#### Option 2: Install manually
You may need to manually set the paths for `C_INCLUDE_PATH`, `LIBRARY_PATH`, `LIBCLANG_PATH`, and `OPENSSL_DIR` to their corresponding
paths in the your conda environment if they can't properly be found on your system.
```
git clone --recursive https://github.com/rhysnewell/Lorikeet.git \ 
cd Lorikeet \
conda env create -n lorikeet -f lorikeet.yml \ 
conda activate lorikeet \ 
bash install.sh # or run without installing e.g. `cargo run --release -- genotype -h` \
lorikeet genotype -h
```

Depending on your local network configuration, you may have problems obtaining Lorikeet via git.
If you see something like this you may be behind a proxy that blocks access to standard git:// port (9418).

```
$ git clone --recursive git://github.com/rhysnewell/Lorikeet.git
Cloning into 'Lorikeet'...
fatal: Unable to look up github.com (port 9418) (Name or service not known)
```

Luckily, thanks to this handy tip from the developer of [Freebayes](https://github.com/ekg/freebayes) we can work around it.
If you have access to https:// on port 443, then you can use this 'magic' command as a workaround to enable download of the submodules:

```
git config --global url.https://github.com/.insteadOf git://github.com/
```

## Usage

Input can either be reads and reference genome, or MAG. Or a BAM file and associated genome.

```
Strain genotyping analysis for metagenomics

Usage: lorikeet <subcommand> ...

Main subcommands:
    genotype    *Experimental* Resolve strain-level genotypes of MAGs from microbial communities
    consensus   Creates consensus genomes for each input reference and for each sample
    call        Performs variant calling with no downstream analysis
    evolve      Calculate dN/dS values for genes from read mappings

Other options:
    -V, --version   Print version information

Rhys J. P. Newell <r.newell near uq.edu.au>
```

Call variants from bam:

`lorikeet call --bam-files my.bam --longread-bam-files my-longread.bam --genome-fasta-directory genomes/ -x fna
     --bam-file-cache-directory saved_bam_files --output-directory lorikeet_out/ --threads 10 --plot`

Call variants from short reads and longread bam files:

`lorikeet call -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam`

## Workflow

An updated workflow image is on its way.

## Output

#### Genotype 
Genotype will produce:
- Variants that clustered into representations of the expected strain level genotypes
- Sample adjacency matrix displaying the number of shared variants seen between samples.
- VCF file detailing all observed variants
- Expected coverage values for each of the produced genotypes in each sample

#### Polish
Polish will produce:
- Consensus genomes for each input reference across each sample.
- Sample adjacency matrix displaying the number of shared variants seen between samples.
- VCF file detailing all observed variants

#### Evolve
Evolve will produce:
- GFF file for each input reference with dN/dS values within coding regions based on the possible variants 
  found along the reference
    - These dN/dS values only take single nucleotide polymorphisms into account but INDELs can still be reported.
- VCF file detailing all observed variants

#### Call
Call only produces the called variants in a VCF file without any downstream analysis.

## FAQs

Your feedback and input helps me as a developer create better tools for you. So please, if you have any questions
raise them as an issue on this GitHub. If I feel like it is an issue others might frequently have then I'll place my 
response here.

## Resources

The variant calling algorithm is basically a one-to-one re-implementation of the algorithm used in GATK HaplotypeCaller.
As such, many of the FAQs and documentation for HaplotypeCaller can be useful in understanding how Lorikeet actually
finds variants. An overview of the HaplotypeCaller pipeline can be found here: [HaplotypeCaller Docs](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).

Lorikeet makes use of a couple new and daunting algorithms. UMAP in particular is an amazing algorithm but might be cause 
for concern since it is difficult to understand how it works and what it is doing. So please look over this amazing article 
by Andy Coenen and Adam Pearce: [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

## Citation

*Watch this space* A paper is on its way. If you use rosella and like the results before the paper, then please cite this GitHub

## License

Code is [GPL-3.0](LICENSE)