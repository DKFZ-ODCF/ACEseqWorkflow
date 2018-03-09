#!/usr/bin/env bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

module load "htslib/$HTSLIB_VERSION"
module load "perl/$PERL_VERSION"
module load "python/$PYTHON_VERSION"
module load "R/$RSCRIPT_VERSION"
module load "bedtools/$BEDTOOLS_VERSION"
module load "samtools/$SAMTOOLS_VERSION"
module load "vcftools/$VCFTOOLS_VERSION"

source /ibios/tbi_cluster/virtualenvs/warsow/python_2.7.9_SNVCalling_1.2.166-1/bin/activate

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export BCFTOOLS_BINARY=bcftools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
