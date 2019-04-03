# Tumor Mutational Burden Calculator

## Requirements:
* samtools depth in path
* VCF called by sentieon and annotated with VEP
* Tumor concentration


## Manual:

`./filter.pl -h`

`OPT	Description`

`-v	input VCF REQUIRED`

`-b	input BAM`

`-c	coverage cutoff REQUIRED if -b in use`

`-t	tumour concentration REQUIRED`

`-s	include Synonymous mutations`

`-e	ExAC MAF cutoff. Default 0.01`

`-h	this help message`

`-o	output filtered vcf`

