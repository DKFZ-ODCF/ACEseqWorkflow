Installation & Run instructions
================================

The Roddy based version for direct execution on hpc clusters as well as the Docker version of ACEseq can be found under http://bfg-nfs3.ipmb.uni-heidelberg.de. For both version an initial download of Reference files is necessary.

New versions of the ACEseq plugin can be obtained from https://github.com/eilslabs/ACEseqWorkflow and can be used in the Roddy-based version.


Roddy-based version
^^^^^^^^^^^^^^^^^^^^^
To run the Roddy-based version of ACEseq please download the pre-packed zip file from http://bfg-nfs3.ipmb.uni-heidelberg.de. Three steps are reuired to ensure running of ACEseq.

1. Run the "prepareRoddyInstallation.sh" script.
2. Download all reference files as specified in the section below. 
3. Set up the conda environment or install the necessary software as specified in the section below.

Prior to running ACEseq a few parameters need to be adjusted in the configuration files. The output directory is specified within $PATH_TO_ACEseq_RODDY_VERSION/configurations/projectsACEseqTest.xml. Here the variables "baseDirectoryReference", "inputBaseDirectory", "outputBaseDirectory", "outputAnalysisBaseDirectory" need to be set. If no SVs should be included the following cvalues should be included:

.. code-block:: ini

    <cvalue name='runWithSv' value='true' type="boolean"/>
    <cvalue name='SV' value='yes' type="boolean"/>  

Otherwise "svOutputDirectory" and the SV bedpe filename in the filenames section needs to be set.


.. code-block:: ini

    <configurationvalues>
  
      <cvalue name='svOutputDirectory' value='${outputAnalysisBaseDirectory}/nameOfDirectoryWithSVResults' type="path"/>
    </configurationvalues>
  
    <filenames package='de.dkfz.b080.co.files' filestagesbase='de.dkfz.b080.co.files.COFileStage'>
       <filename class="TextFile" onMethod="de.dkfz.b080.co.aceseq.ACESeqMethods.mergeSv"
                selectiontag="svFileTag"
                pattern='${svOutputDirectory}/${pid}_svs.bedpe'/>
    </filenames>

Technical specifications need to be set within the file $PATH_TO_ACEseq_RODDY_VERSION/configurations/applicationProperties.ini. The path to the project.xml and the path to the plugins ($PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/) need to be set under configurationDirectories and pluginDirectories. Finally the job manager and execution host need to be set.

Please have a look at the following default applicationProperties.ini file:

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


To execute ACEseq run

::

    sh $PATH_TO_ACEseq_RODDY_VERSION//Roddy/roddy.sh rerun ACEseq@copyNumberEstimation $pid \
    --useconfig=$PATH_TO_ACEseq_RODDY_VERSION/configuration/applicationProperties.ini \
    --cvalues="bamfile_list:$pathToControlBamFile;$pathToTumorBamFile,sample_list:control;tumor,possibleControlSampleNamePrefixes:control,possibleTumorSampleNamePrefixes:tumor"


More information on Roddy can be found `here <https://roddy-documentation.readthedocs.io/>`_.

Docker version
^^^^^^^^^^^^^^^
1. Download all reference files as specified in the section below. 
2. Download the Base and ACEseq Docker images from the website: http://bfg-nfs3.ipmb.uni-heidelberg.de
3. Import both files with:

::

	docker load < BaseDockerContainer.tar.gz

::

	docker load < ACEseqDockerContainer.tar.gz

4. Download the control files archive and extract them. The directory contains the file "roddy.sh". Please call this script with: bash roddy.sh

Reference files
^^^^^^^^^^^^^^^^
To get all necessary reference files run the script $PATH_TO_PLUGIN_DIRECTORY/installation/downloadReferences.sh from the destination path for all files.
Please convert the bigwig file in databases/UCSC to a BedGraph (https://genome.ucsc.edu/goldenpath/help/bigWig.html) and save it under wgEncodeCrgMapabilityAlign100mer_chr.bedGraph, 
compress it with bgzip and index with tabix.
The variable baseDirectoryReference in the project.xml  needs to be set to the path from which the downloader script was run.

Software
^^^^^^^^^
All software required to run ACEseq is stored in Bioconda and can be downloaded to set up a conda environment. Specifications about the packages are given in `$PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml` (for the zipped Roddy version the $PATH_TO_PLUGIN_DIRECTORY is $PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/).

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




