Requirements
=============


Hardware
^^^^^^^^
ACEseq is requires the execution of multiple jobs that are highly parallelized in the beginning but linearize towards the end of the workflow.
It requires a maximum of 50g RAM in few of the Jobs.
On a HPC cluster with multiple cores available it will usually finish within 24h (100-160 CPU h). The final output usually requires between 4 and 6g memory.

Software
^^^^^^^^
ACEseq can be run as stand alone pipeline with the Roddy pipeline framework (roddy-documentation.readthedocs.io, https://github.com/eilslabs/Roddy). A complete zipped version required to run ACEseq
is available under http://bfg-nfs3.ipmb.uni-heidelberg.de. Currently SGE, PBS and LSF schedulers are supported but a local execution is also possible. 

Additionally a Docker version of ACEseq can be downloaded from http://bfg-nfs3.ipmb.uni-heidelberg.de.

All software required to run ACEseq is stored in Bioconda and can be downloaded to set up a conda environment. Specifications about the packages are given in `$PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml`

The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. You should set up the Conda environment at a centralized position available from all compute hosts. 

First install the BioConda channels:

::

    conda config --add channels r
::

    conda config --add channels defaults

::

    conda config --add channels conda-forge

::

    conda config --add channels bioconda

Then install the environment

::

    conda env create -n ACEseqWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable in `resources/configurationFiles/analysisCopyNumberEstimation.xml`.


If you do not want to use Conda, you can get a complete list of all packages and package versions Conda would install from the `resources/configurationFiles/conda.yaml`.


Reference files
^^^^^^^^^^^^^^^^^
To obtain all necessary reference files run $PATH_TO_PLUGIN_DIRECTORY/installation/downloadRefrences. 
It is tuned to be run on hg19 but a liftover of all files should probably enable a run on GRch38.

Please convert the bigwig file in databases/UCSC to a BedGraph (https://genome.ucsc.edu/goldenpath/help/bigWig.html) and save it under wgEncodeCrgMapabilityAlign100mer_chr.bedGraph, 
compress it with bgzip and index with tabix.
The variable baseDirectoryReference in the project.xml  needs to be set to the path from which the downloader script was run.



