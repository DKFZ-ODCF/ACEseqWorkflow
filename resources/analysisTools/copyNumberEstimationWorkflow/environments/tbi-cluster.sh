#!/usr/bin/env bash
# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
# Load a Conda environment.

module load "htslib/$HTSLIB_VERSION"
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix

module load "perl/$PERL_VERSION"
export PERL_BINARY=perl

module load "python/$PYTHON_VERSION"
export PYTHON_BINARY=python

module load "R/$RSCRIPT_VERSION"
export RSCRIPT_BINARY=Rscript

module load "bedtools/$BEDTOOLS_VERSION"
export INTERSECTBED_BINARY=intersectBed

module load "vcftools/$VCFTOOLS_VERSION"
export VCFTOOLS_SORT_BINARY=vcf-sort

module load "samtools/$SAMTOOLS_VERSION"
export BCFTOOLS_BINARY=bcftools
export SAMTOOLS_BINARY=samtools
