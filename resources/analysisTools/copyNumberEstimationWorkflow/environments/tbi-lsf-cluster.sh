#!/usr/bin/env bash

modules load "htslib/$HTSLIB_VERSION"
modules load "perl/$PERL_VERSION"
modules load "python/$PYTHON_VERSION"
modules load "R/$RSCRIPT_VERSION"
modules load "bedtools/$BEDTOOLS_VERSION"
modules load "samtools/$SAMTOOLS_VERSION"
modules load "vcftools/$VCFTOOLS_VERSION"

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export BCFTOOLS_BINARY=bcftools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
