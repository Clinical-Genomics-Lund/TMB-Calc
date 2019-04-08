# Tumor Mutational Burden Calculator

## Requirements:
* VCF called by sentieon and annotated with VEP
* VCF called by Freebayes and annotated with VEP
* Tumor concentration
* [SGZ-installation](https://github.com/jsunfmi/SGZ)


## **tmb_calc_wrapper.pl**
Runs filter.pl for both Sentieon and Freebayes VCF, takes ID, VCF_SENT, VCF_FREEBAYES, TUMOR CONC and coding size in mb as input.
Takes the intersect of the two filtered VCFs and filter by SGZ. TMB is printed to STDOUT for each sample

### Input file
**for each sample:**

ID:path/to/sentieon/vcf:path/to/freebayes/vcf:tumorconc

### Options:

`./tmb_calc_wrapper.pl *coding-size* *out/path/*`

## **filter.pl**

### Options:

`./filter.pl -h`

**-flag** | **description**
--- | ---
**-v** | input VCF ***REQUIRED***
**-d** | coverage cutoff. Integer
**-t** | tumour concentration. FLoat. ***REQUIRED***
**-s** | exclude synonymous mutations
**-e** | ExAC MAF cutoff. Float. ***Default 0.0001***
**-h** | this help message
**-o** | output filtered BED

