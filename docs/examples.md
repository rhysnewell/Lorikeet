---
title: Examples
---

Examples
========

## Forcing Variant Calls

Through use of the `-f, --features-vcf` argument it is possible to provide lorikeet with a list of known variants or variants of concern that will be forcibly called by
the algorithm if there is any activity occurring at that genomic location in the provided samples. It is possible that lorikeet
will also call any potential variation events surrounding a given location as well. This can be useful if you are more 
interested in the activity occurring around some given locations rather than specific variation events.

The list of variants must be provided in VCF format with some caveats on how the variant locations are written. E.g.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=lorikeet-v0.6.0
##contig=<ID=random10000~random_sequence_length_10000_1,length=10000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
random10000~random_sequence_length_10000_1	223	.	G	C	1445.63	.
random10000~random_sequence_length_10000_1	435	.	G	C	1941.63	.
random10000~random_sequence_length_10000_1	949	.	T	A	383.629	.
```

Note that you need to specific both the genome name and contig name in the `CHROM` column separated by the `~` character
 e.g. `random10000~random_sequence_length_10000_1`. Additionally, the standard `FORMAT` and `INFO` columns are optional when
 provding a VCF file.
 
Finally, the provided VCF file must be compressed using `bgzip` and indexed using `bcftools index` e.g.
```
bgzip -c random10000.vcf > random10000.vcf;
bcftools index random10000.vcf.gz;
lorikeet call -r random10000.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam -f random10000.vcf.gz
```