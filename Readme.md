# ACEseq Workflow 

Author: Kortine Kleinheinz
k.kleinheinz@dkfz-heidelberg.de

## Description
ACEseq (Allele-specific copy number estimation with whole genome sequencing) is a tool to estimate allele-specific copy numbers from WGS data and comes along with a variety of features:
* GC/replication timin Bias correction
* quality check
* SV breakpoint inclusion
* automated estimation of ploidy and tumore cell content
* HRD/TAI/LST score estimation 
* with/without matched control processing 

The algorithm is further documented under aceseq.readthedocs.io

## Installation

### High-Throughput Cluster with Shared Storage

#### Conda

The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. You should set up the Conda environment at a centralized position available from all compute hosts. 

First install the BioConda channels:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then install the environment

```
conda env create -n ACEseqWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/copyNumberEstimationWorkflow/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable in `resources/configurationFiles/analysisCopyNumberEstimation.xml`.

#### Other

If you do not want to use Conda, you can get a complete list of all packages and package versions Conda would install from the `resources/configurationFiles/conda.yaml`.



## Prepackaged files

There are two different packaged versions of ACEseq and Roddy available. One is a prepackaged zip file with a full Roddy / Plugin setup and different scripts to install everything properly. The second one is a dockerized version of the workflow. You can download both formats here: (http://bfg-nfs3.ipmb.uni-heidelberg.de).

Further information about the usage of both versions can be found at our Read The Docs website: (http://aceseq.readthedocs.io).

