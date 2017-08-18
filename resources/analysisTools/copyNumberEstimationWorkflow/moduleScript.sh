#!/usr/bin/bash

source ${CONFIG_FILE}

#Rscript
module load R/$RSCRIPT_VERSION

#python
module load python/$PYTHON_VERSION

#perl
module load perl/$PERL_VERSION

#bgzip
#module load bgzip/$BGZIP_VERSION
#
#tabix
#module load tabix/$TABIX_VERSION
module load htslib/$TABIX_VERSION

#vcf sort
#module load vcf-sort/$VCFTOOLS_SORT_VERSION

#samtools (loads bcftools too)1
module load samtools/$SAMTOOLS_VERSION
#module load bcftools/$BCFTOOLS_VERSION

#intersectbed
module load bedtools/$INTERSECTBED_VERSION

BGZIP_BINARY=bgzip
TABIX_BINARY=tabix
PERL_BINARY=perl
PYTHON_BINARY=python
RSCRIPT_BINARY=Rscript
INTERSECTBED_BINARY=intersectBed
BCFTOOLS_BINARY=bcftools
VCFTOOLS_SORT_BINARY=vcf-sort
SAMTOOLS_BINARY=samtools
