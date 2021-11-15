---
title: Call
---

Call
========

The `call` subcommand of Lorikeet is the primary method for users who wish to just call variants across one or more
reference genomes. The variant calling process uses an algorithm based on [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) which 
calls Single Nucleotide Polymorphisms (SNPs) and Insertions and Deletions (INDELs) via local re-assembly of haplotypes.

## Quick start


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