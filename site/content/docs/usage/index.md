---
title: 'Usage'
date: 2019-02-11T19:30:08+10:00
draft: false
weight: 4
summary: Syntax highlighting and menus can be configured via `config.toml`.
---

## Quick Start

### Genotyping

If provided a complete or draft microbial reference genome, and one or more sets of reads or BAM files Lorikeet can find variant locations within the reference genome
based on the read mappings. Ideally, multiple samples will be provided so log-ratio variances can be calculated between
each pairwise combination of variants. If only a single sample is provided, a euclidean distance matrix is calculated
based on the abundance of each variant location. Non-negative matrix factorization is then performed on the resulting pairiwse matrix
reducing the variants into a lower latent space where each variant is assigned to its resulting dominant factor.


```commandline
lorikeet genotype -r draft_genome.fasta -1 sample_1.1.fq.gz sample_2.1.fq.gz -2 sample_1.2.fq.gz sample_2.2.fq.gz -f 10 -t 10 --eps-min 0.05 --eps-max 0.15
```

### Variant Calling


```commandline
lorikeet polymorph -r draft_genome.fasta -1 sample_1.1.fq.gz sample_2.1.fq.gz -2 sample_1.2.fq.gz sample_2.2.fq.gz -f 10 -t 10
```

