#!/usr/bin/env bash
# Load a Conda environment.

createCleanCondaEnvironment () {
    unset LD_LIBRARY_PATH
    unset MODULE_PATH
    unset PKG_CONFIG_PATH
    unset PYTHONHOME
    unset PYTHONPATH
    unset PYTHONSTARTUP
    unset PYTHON_LIB
    unset PBS_SHARED_BIN
    unset PERL5LIB
    unset PERL_LOCAL_LIB_ROOT
    unset PERL_MB_OPT
    unset PERL_MM_OPT
    export PATH="$HOME/miniconda3/bin:/opt/torque/bin:/usr/lib64/mpi/gcc/openmpi/bin:/opt/maui/bin:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/X11R6/bin:/usr/lib/mit/bin"
}

createCleanCondaEnvironment

source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export PERL_BINARY=perl
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export INTERSECTBED_BINARY=intersectBed
export BCFTOOLS_BINARY=bcftools
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMTOOLS_BINARY=samtools
