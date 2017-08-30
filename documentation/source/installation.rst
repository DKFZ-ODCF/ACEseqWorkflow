Installation & Run instructions
==============

Reference files
^^^^^^^^^^^^^^^^
To get all necessary reference files run the script $PATH_TO_PLUGIN_DIRECTORY/installation/downloadReferences.sh from the destination path for all files.
Please convert the bigwig file in databases/UCSC to a BedGraph (https://genome.ucsc.edu/goldenpath/help/bigWig.html) and save it under wgEncodeCrgMapabilityAlign100mer_chr.bedGraph, 
compress it with bgzip and index with tabix.
The variable baseDirectoryReference in the project.xml  needs to be set to the path from which the downloader script was run.

Software
^^^^^^^^^
All software required to run ACEseq is stored in Bioconda and can be downloaded to set up a conda environment. Specifications about the packages are given in `$PATH_TO_PLUGIN_DIRECTORY/resources/configurationFiles/analysisCopyNumberEstimation.xml`.

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

	conda env create -n ACEseqWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/copyNumberEstimationWorkflow/environments/conda.yml


Roddy-based version
^^^^^^^^^^^^^^^^^^^^^
The Roddy-based version can be downloaded from XYZ. The output directory can be specified within configurations/projectsACEseqTest.xml. Here the variables "baseDirectoryReference", "inputBaseDirectory", "outputBaseDirectory", "outputAnalysisBaseDirectory" need to be set. Additionally the SV=no and runWithSv=false needs to be set if SV breakpoints should not be included. Otherwise "svOutputDirectory" and the SV bedpe filename in the filenames section needs to be set.

Technical specifications need to be set within the file configurations/applicationProperties.ini. The path to the project.xml and the path to the plugins ($PATH_TO_PLUGIN_DIRECTORY/Roddy/dist/plugins/ needs to be set under configurationDirectories and pluginDirectories. Finally the job manager and execution host need to be set.

The execution can be Please have a look at the following default application properties ini
file:

.. code-block:: ini

    [COMMON]
    useRoddyVersion=current                     # Use the most current version for tests

    [DIRECTORIES]
    configurationDirectories=[FOLDER_WITH_CONFIGURATION_FILES]
    pluginDirectories=[FOLDER_WITH_PLUGINS]

    [COMMANDS]
    jobManagerClass=de.dkfz.roddy.execution.jobs.direct.synchronousexecution.DirectSynchronousExecutionJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.pbs.PBSJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.sge.SGEJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.slurm.SlurmJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.lsf.rest.LSFRestJobManager
    commandFactoryUpdateInterval=300
    commandLogTruncate=80                       # Truncate logged commands to this length. If <= 0, then no truncation.

    [COMMANDLINE]
    CLI.executionServiceUser=USERNAME
    CLI.executionServiceClass=de.dkfz.roddy.execution.io.LocalExecutionService
    #CLI.executionServiceClass=de.dkfz.roddy.execution.io.SSHExecutionService
    CLI.executionServiceHost=[YOURHOST]
    CLI.executionServiceAuth=keyfile
    #CLI.executionServiceAuth=password
    CLI.executionServicePasswd=
    CLI.executionServiceStorePassword=false
    CLI.executionServiceUseCompression=false
    CLI.fileSystemInfoProviderClass=de.dkfz.roddy.execution.io.fs.FileSystemInfoProvider


::

    sh $PATH_TO_PLUGIN_DIRECTORY/Roddy/roddy.sh rerun ACEseq@copyNumberEstimation $pid \
    --useconfig=$PATH_TO_PLUGIN_DIRECTORY/applicationProperties.ini \
    --cvalues="bamfile_list:$pathToControlBamFile;$pathToTumorBamFile,sample_list:control;tumor,possibleControlSampleNamePrefixes:control,possibleTumorSampleNamePrefixes:tumor"


More information on Roddy can be found `here <https://roddy-documentation.readthedocs.io/>`_.

Docker version
^^^^^^^^^^^^^^^
The Docker-based version can be downloaded form XYZ.
