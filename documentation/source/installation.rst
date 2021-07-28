.. _old-server: http://bfg-nfs3.ipmb.uni-heidelberg.de/
.. _ppcg-server: https://ppcg.dkfz.de/pipelines/


.. role:: bash(code)
   :language: bash
.. role:: xml(code)
   :language: xml

.. _installation:

Installation & Run instructions
===============================

To run the ACEseq-workflow multiple components are needed:

  * The `Roddy workflow management framework <https://github.com/TheRoddyWMS/Roddy>`_
  * `ACEseqWorkflow <https://github.com/DKFZ-ODCF/ACEseqWorkflow>`_
  * `COWorkflowsBasePlugin <https://github.com/TheRoddyWMS/COWorkflowsBasePlugin>`_
  * `PluginBase <https://github.com/TheRoddyWMS/Roddy-Base-Plugin>`_
  * `DefaultPlugin <https://github.com/TheRoddyWMS/Roddy-Default-Plugin>`_
  * Bioinformatic software stack
  * Reference data

The most flexible way to install the workflow is the :ref:`manual way`. Alternatively, there are Docker containers for some versions, which contain all components and that will run the complete workflow within the single container. See :ref:`old-containers` and :ref:`new-containers`. Furthermore, there is an installation specifically for ACEseqWorkflow 1.2.10 (:ref:`aceseq_1.2.10-installation`)

.. _manual way:

The Manual Way
----------------

Roddy Installation
^^^^^^^^^^^^^^^^^^

The first step is to have a Roddy installation. Please, first read a bit about how Roddy requires the plugins to be installed. See `here <https://roddy-documentation.readthedocs.io/>`_. In short, the Roddy core component and the "PluginBase" and the "DefaultPlugin" are installed with a specific directory layout. Unless you use a very old version of ACEseq (< 1.2) you should use the newest available Roddy, PluginBase and DefaultPlugin versions.

1. Download the most recent RoddyEnv ZIP from https://github.com/TheRoddyWMS/Roddy/releases.
2. Download matching Roddy Roddy ZIP from https://github.com/TheRoddyWMS/Roddy/releases and extract it into some directory. The `dist/plugins/` directory will be created in that directory.
3. Download the DefaultPlugin https://github.com/TheRoddyWMS/Roddy-Default-Plugin/releases and extract it into `dist/plugins/DefaultPlugin_x.y.z`
4. Download the PluginBase from https://github.com/TheRoddyWMS/Roddy-Base-Plugin/releases and extract it into `dist/plugins/PluginBease_x.y.z`

Note that the '_' that separates the plugin name from its version number can be replaced by a minus symbol '-'.

Plugin Installation
^^^^^^^^^^^^^^^^^^^

1. Download the ACEseqWorkflow https://github.com/DKFZ-ODCF/ACEseqWorkflow/tags. Version 1.2.8-4 is used in our production systems. GRCh38/hg38 support is in development and will be available only with version 6. Extract the archive into a directory called "ACEseqWorkflow_$version" with `$version` being the correct version number.
2. Look up the version of the COWorkflowsBasePlugin listed in the "buildinfo.txt" at the root of the ACEseqWorkflow plugin that you just downloaded. Download that version from https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/tags. Create a "COWorkflowsBasePlugin_$version" directory.

You can create the COWorkflowsBasePlugin and ACEseqWorkflow directories in the `dist/plugins` directory, which was extracted with the RoddyEnv ZIP.

Software Stack (Conda)
^^^^^^^^^^^^^^^^^^^^^^

The software stack is the software that the workflow jobs need to do the bioinformatic work. If you run the workflow on a cluster you need to ensure that the software stack is also available on all cluster nodes that execute jobs.

For reproducibility the software stack has been stored with `Conda Pack <https://conda.github.io/conda-pack/>`_ in an archive. To restore the environment you need to

1. Download the appropriate archive from https://ppcg.dkfz.de/pipelines/conda-pack
2. Unpack and set up the environment
    .. code-block:: bash

        mkdir $envBaseDir/ACEseqWorkflow    # The directory name is the environment name.
        cd $envBaseDir/ACEseqWorkflow
        source bin/activate
        conda unpack

Now, if you do `conda env list` the environment should be listed. If not make sure you installed the environment into `$CONDA_PREFIX/envs/` or another `environments directory <https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#specify-environment-directories-envs-dirs>`_ configured for your Conda installation.

Try out `conda activate ACEseqWorkflow` to activate the environment. You may also want to start Python or R to test whether the packages can be imported and libraries loaded.

To configure the workflow to use the environment you need to set two configuration variables (e.g. in the command-line or a project XML; see `Roddy documentation <http://roddy-documentation.readthedocs.io/>`_). For instance on the CLI

.. code-block:: bash

    roddy.sh ... \
        --cvalues="workflowEnvironmentScript:workflowEnvironment_conda;condaEnvironmentName:$value"

The Conda installation, the environment and the `.condarc` (if you have modified it to set the `envs_dirs`) need to be available on your cluster nodes.


.. _aceseq_1.2.10-installation:

Prepackaged files (ACEseq 1.2.10 only)
--------------------------------------

On old-server_ you can find archives for the 1.2.10 plugin version. The prepackaged zip files contains a full Roddy / Plugin setup and include different scripts to install all necessary software and download the required reference files. Currently, we do not intent to update these prepackaged installation files or the Docker version. Note that the Roddy version packaged not capable of submitting to LSF. Only Torque/PBS is supported.

.. _old-containers:

Stand-alone Roddy for Execution on HTC Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the Roddy-based version of ACEseq please download the pre-packed zip file from old-server_. Three steps are required to ensure running of ACEseq.

1. Run the "prepareRoddyInstallation.sh" script.
2. Download all reference files as specified in the section `Reference files <install-reference-files_>`_.
3. Set up the Conda environment or install the necessary software as specified in the section "Software".

Before running ACEseq a few parameters need to be adjusted in the configuration files. The output directory is specified in $PATH_TO_ACEseq_RODDY_VERSION/configurations/projectsACEseqTest.xml. Here the variables "baseDirectoryReference", "inputBaseDirectory", "outputBaseDirectory", "outputAnalysisBaseDirectory" need to be set. If no SVs should be included the following configuration values (cvalues) should be included:

.. code-block:: xml

    <cvalue name='runWithSv' value='false' type="boolean"/>
    <cvalue name='SV' value='false' type="boolean"/>

If you set these values to "true", then "svOutputDirectory" and the SV bedpe filename in the filenames section need to be set.

.. code-block:: xml

    <configurationvalues>
      <cvalue name='svOutputDirectory' value='${outputAnalysisBaseDirectory}/nameOfDirectoryWithSVResults' type="path"/>
    </configurationvalues>


Technical specifications are set in the file $PATH_TO_ACEseq_RODDY_VERSION/configurations/applicationProperties.ini. The path to the project.xml and the path to the plugins ($PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/) need to be set with "configurationDirectories" and "pluginDirectories". Finally the job manager and execution host need to be set. Please have a look at the following default applicationProperties.ini file:

.. code-block:: ini

    [COMMON]
    useRoddyVersion=current                     # Use the most current version for tests

    [DIRECTORIES]
    configurationDirectories=[FOLDER_WITH_CONFIGURATION_FILES]
    pluginDirectories=[FOLDER_WITH_PLUGINS]

    [COMMANDS]
    # Choose your job-manager. The first one executes the jobs locally.
    jobManagerClass=de.dkfz.roddy.execution.jobs.direct.synchronousexecution.DirectSynchronousExecutionJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.pbs.PBSJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.sge.SGEJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.slurm.SlurmJobManager
    #jobManagerClass=de.dkfz.roddy.execution.jobs.cluster.lsf.rest.LSFRestJobManager
    commandFactoryUpdateInterval=300
    commandLogTruncate=0                    # Truncate logged commands to this length. If <= 0, then no truncation.

    [COMMANDLINE]
    CLI.executionServiceUser=USERNAME
    # The execution service determines how commands are exectuted. Locally, or via SSH.
    # SSHExecution service is needed if the host on which you run Roddy is different from the
    # submission host that allows executing the bsub/qsub command.
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

.. code-block:: bash

    sh $PATH_TO_ACEseq_RODDY_VERSION/Roddy/roddy.sh rerun ACEseq@copyNumberEstimation $pid \
        --useconfig=$PATH_TO_ACEseq_RODDY_VERSION/configuration/applicationProperties.ini \
        --cvalues="bamfile_list:$pathToControlBamFile;$pathToTumorBamFile,sample_list:control;tumor,possibleControlSampleNamePrefixes:control,possibleTumorSampleNamePrefixes:tumor"


More information on Roddy can be found `here <https://roddy-documentation.readthedocs.io/>`_.

Docker version
^^^^^^^^^^^^^^

1. Download all reference files as specified in the section below.
2. Download the Base and ACEseq Docker images from the website: old-server_
3. Import both files with (names might differ based on supplied version):

.. code-block:: bash

	docker load < BaseDockerContainer.tar.gz

.. code-block:: bash

	docker load < ACEseqDockerContainer.tar.gz

4. Download the control files archive and extract them. The directory contains the file "roddy.sh". Please call this script with: :bash:`bash roddy.sh`. You will see:

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

.. code-block:: bash

        bash roddy.sh run ACEseq \
            ./config/ \
            stds \
            /home/roddy/someproject/control_MB99_merged.mdup.bam \
            /home/roddy/someproject/tumor_MB99_merged.mdup.bam \
            control \
            tumor \
            /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes \
            ./output

Here you tell roddy to run the ACEseq configuration using the config folder in the current directory with a control and tumor bam. Also you tell Roddy the samples for both files namely control and tumor. Finally, you supply the path to the reference files and the folder where you will store your output data.


.. _new-containers:

Docker version
--------------

Different versions of the ACE-seq workflow have been packaged with other workflows and the reference data. These containers have the required software-stacks installed but run all workflow jobs within the same container (so they are not optimized for throughput). To download these pipelines and reference data see ppcg-server_. Instructions for the installation are given in the archives.

.. _Github-Releases: https://github.com/eilslabs/ACEseqWorkflow/releases

.. _install-reference-files:

Reference files
---------------

The workflow uses various files as reference files, such as a reference genome or annotation files. We provide installation scripts in the `installation/` directory. To download and prepare the reference files please check out the ACEseq repository and do

.. code-block:: bash

   bash $PATH_TO_PLUGIN_DIRECTORY/installation/downloadReferences $targetDirectory

with `$targetDirectory` being the directory into which you want to install the files. The variable `baseDirectoryReference` in your configurations needs to be set to the `$targetDirectory` path.

Note that the current plugin version is tuned to be run on the hg19 human assembly, but a liftover of all files should probably enable a run on GRch38.

Alternative reference files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reference data can also be downloaded from *TODO SERVER NEEDS TO BE SET*.

.. WARNING:: The "database" file https://ppcg.dkfz.de/pipelines/ is huge (> 20 GB). Please only download it, if you need it.

