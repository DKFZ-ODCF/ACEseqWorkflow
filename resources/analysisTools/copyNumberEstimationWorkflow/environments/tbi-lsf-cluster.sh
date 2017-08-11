#!/usr/bin/env bash

modules load htslib/0.2.5
modules load perl/5.20.2
modules load python/2.7.9
modules load R/3.3.1
modules load bedtools/2.16.2
modules load samtools/0.1.19
modules load vcftools/0.1.10

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export BCFTOOLS_BINARY=bcftools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
