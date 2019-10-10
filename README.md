# Lorikeet
*Note:* StrainM will be changing its name to lorikeet

A strain resolver for metagenomics. Currently under active development, StrainM currently acts as a variant caller for reads mapped to a metagenome assembled genome.

Input can either be reads and reference genome, or MAG. Or a BAM file and associated genome.

# Usage

`lorikeet polymorph -b input.bam -r input_genome.fna`

OR

`lorikeet polymorph -r input_genome.fna -1 forward_reads.fastq -2 reverse_reads.fastq`

# Output 
The output of polymorph is a tab-delimited file with 9 fields

`tid  pos variant reference abundance depth genotypes type  connected_bases`
