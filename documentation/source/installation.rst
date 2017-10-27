Installation & Run instructions
================================

Two strategies for the deployment of Roddy and the ACEseq workflow are described here. The files for both, the Roddy-based deployment for direct
execution on HPC clusters as well as a Docker version of ACEseq can be found under http://bfg-nfs3.ipmb.uni-heidelberg.de. For both version an initial
download of reference files is necessary.

New versions of the ACEseq plugin can be obtained from https://github.com/eilslabs/ACEseqWorkflow and can be used in the Roddy-based version.

Roddy-based Deployment
^^^^^^^^^^^^^^^^^^^^^^^
The workflow management tool Roddy supports the batch processing systems Torque/PBS and IBM LSF. Other systems are planned but currently not supported.

To run the Roddy-based version of ACEseq please download the pre-packed zip file from http://bfg-nfs3.ipmb.uni-heidelberg.de. Three steps are required to ensure running of ACEseq.

1. Run the "prepareRoddyInstallation.sh" script.
2. Download all reference files as specified in the section "Reference files" (below).
3. Set up the Conda environment or install the necessary software as specified in the section "Software" (below).

Before running ACEseq a few parameters need to be adjusted in the configuration files. The output directory is specified in $PATH_TO_ACEseq_RODDY_VERSION/configurations/projectsACEseqTest.xml. Here the variables "baseDirectoryReference", "inputBaseDirectory", "outputBaseDirectory", "outputAnalysisBaseDirectory" need to be set. If no SVs should be included the following configuration values (cvalues) should be included:

.. code-block:: ini

    <cvalue name='runWithSv' value='true' type="boolean"/>
    <cvalue name='SV' value='yes' type="boolean"/>


Otherwise "svOutputDirectory" and the SV bedpe filename in the filenames section need to be set.

.. code-block:: ini

    <configurationvalues>
      <cvalue name='svOutputDirectory' value='${outputAnalysisBaseDirectory}/nameOfDirectoryWithSVResults' type="path"/>
    </configurationvalues>

    <filenames package='de.dkfz.b080.co.files' filestagesbase='de.dkfz.b080.co.files.COFileStage'>
       <filename class="TextFile" onMethod="de.dkfz.b080.co.aceseq.ACESeqMethods.mergeSv"
                selectiontag="svFileTag"
                pattern='${svOutputDirectory}/${pid}_svs.bedpe'/>
    </filenames>

Technical specifications are set in the file $PATH_TO_ACEseq_RODDY_VERSION/configurations/applicationProperties.ini. The path to the project.xml and the path to the plugins ($PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/) need to be set under configurationDirectories and pluginDirectories. Finally the job manager and execution host need to be set.

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
3. Import both files with (names might differ based on supplied version):

::

	docker load < BaseDockerContainer.tar.gz

::

	docker load < ACEseqDockerContainer.tar.gz

4. Download the control files archive and extract them. The directory contains the file "roddy.sh". Please call this script with: bash roddy.sh. You will see:

::

        #!/bin/bash
        # 1: Run mode, which might be "run" or "testrun"
        # 2: Configuration identifier, normally "ACEseq"
        # 3: Configuration directory
        # 4: Dataset identifier / PID
        # 5: Control bam file
        # 6: Tumor bam file
        # 7: Control bam sample name
        # 8: Tumor bam sample name
        # 9: Reference files path
        # 10: Output folder
        # 11: Optional: The SV file

An example call is:

::

        bash roddy.sh run ACEseq ./config/ stds /home/roddy/someproject/control_MB99_merged.mdup.bam /home/roddy/someproject/tumor_MB99_merged.mdup.bam control tumor /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes ./output

Here you tell roddy to run the ACEseq configuration using the config folder in the current directory with a control and tumor bam. Also you tell Roddy the samples for both files namely control and tumor. Finally, you supply the path to the reference files and the folder where you will store your output data.

Software
^^^^^^^^^
The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. You should set up the Conda environment at a centralized position available from all compute hosts. The full specification of the required packages is given in `$PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml` (for the zipped Roddy version the $PATH_TO_PLUGIN_DIRECTORY is $PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/).

Unless you have already done so, you should first install [Conda](https://conda.io/docs/). Then you need to set up the BioConda channel that contains many of the required software packages:

::

    conda config --add channels r

::

    conda config --add channels defaults

::

    conda config --add channels conda-forge

::

    conda config --add channels bioconda


Eventually, you can import the conda environment with

::

	conda env create -n ACEseqWorkflow -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml


Reference files
^^^^^^^^^^^^^^^^
To get all necessary reference files run the script $PATH_TO_PLUGIN_DIRECTORY/installation/downloadReferences.sh in the destination path for all files. The script will create directories directly beneath the path where it executed.

Please convert the bigwig file in databases/UCSC to BedGraph format (https://genome.ucsc.edu/goldenpath/help/bigWig.html) and save it under wgEncodeCrgMapabilityAlign100mer_chr.bedGraph, compress it with bgzip and index it with tabix. Please use the tabix from htslib 0.2.5. We suggest you simply use the previously installed Conda environment to do that. This is more or less the code to do the conversion:

::

    bigWigToBedGraph wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer_chr.bedGraph
    bgzip wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz
    tabix -0 wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz

Finally, set the variable baseDirectoryReference in the project.xml to the path from which the downloader script was run.
