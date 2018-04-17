# ACEseq Workflow

Original Author: Kortine Kleinheinz
k.kleinheinz@dkfz-heidelberg.de

Current Author: Gregor Warsow
g.warsow@dkfz-heidelberg.de

## Description
ACEseq (Allele-specific copy number estimation with whole genome sequencing) is a tool to estimate allele-specific copy numbers from WGS data and comes along with a variety of features:
* GC/replication timin Bias correction
* quality check
* SV breakpoint inclusion
* automated estimation of ploidy and tumore cell content
* HRD/TAI/LST score estimation
* with/without matched control processing

## Installation

### Prepackaged files

Note that an older version of ACEseq is available as Docker version for downloaded at http://bfg-nfs3.ipmb.uni-heidelberg.de. The same ACEseq-version and a matching Roddy version for running on an HTC cluster are available there as well. The prepackaged zip file contains a full Roddy / Plugin setup and includes different scripts to install all necessary software and download the required reference files. Currently, we do not intent to update these prepackaged installation files or the Docker version. Please see the standard way to install recent workflow versions.

### Standard

To run the ACEseq-workflow multiple components are needed:

  * ACEseq workflow plugin
  * Roddy as the workflow execution framework
  * Software stack
  * Reference data

The standard way to install the workflow is to just download the zip-archive from [Github-Releases](https://github.com/eilslabs/ACEseqWorkflow/releases). The file `buildinfo.txt` in the ACEseq zip shows you the Roddy API version that you need for the chosen ACEseq workflow version. Please see the [Roddy repository](https://github.com/TheRoddyWMS/Roddy) for installation instructions for Roddy. Note that the release archive of ACEseq already contains a Jar-archive with the compiled Java/Groovy code (JAR-file) for the given Roddy API version. No compilation of the plugin is therefore required.

### Software Stack (Conda)

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
conda env create -n ACEseqWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable in `resources/configurationFiles/analysisCopyNumberEstimation.xml`. But you can also set the `condaEnvironmentName` variable in any of you other higher-priorized configuration files (see [Roddy documentation](http://roddy-documentation.readthedocs.io/)) or even directly in your Roddy call via `--cvalues="condaEnvironmentName:$value"`.

If you do not want to use Conda, you can get a complete list of all packages and package versions Conda would install from the `resources/configurationFiles/conda.yaml` in the ACEseq repository.

### Reference Data

The workflow uses various files as reference files, such as a reference genome or annotation files. Depending on the contents of these files also the outcome of your analysis may change. We provide installation scripts in the `installation/` directory (currently only in the `github` branch of the repository). To download and prepare the reference files please check out the ACEseq repository and do

```bash
bash $ACEseqPluginDirectory/installation/downloadReferences.sh $targetDirectory
```

with `$targetDirectory` being the directory into which you want to install the files.



