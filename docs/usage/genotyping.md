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
alternate alleles and n is the number of samples. The array is passed into the 
python component of Lorikeet which aims to cluster variants into variant groupings. To achieve this a UMAP embeddings 
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
