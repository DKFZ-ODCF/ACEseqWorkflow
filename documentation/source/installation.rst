.. _old-server: https://hub.dkfz.de/s/XjxXEgCTjjfyMJD
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

Section :ref:`manual way` describes how these components are manually installed, which gives you the most flexibility. This is probably mostly relevant if you have to run multiple different Roddy-based workflows and in a development setup.

There is an installation with pre-packaged files of ACEseqWorkflow 1.2.10 (see :ref:`here <aceseq_1.2.10-installation>`). This installation is mostly interesting for historical reasons and considered "legacy". We can probably not help you much with this installation, if you encounter problems.

Finally, it is possible to run ACEseq in a container. Currently, all available containers wrap the complete workflow including the workflow manager Roddy and even a batch processing system (Sun GridEngine) in a single container. This is done because the workflow manager Roddy does not support submitting containerized cluster jobs. The disadvantage is that the one-size-fits all container for all the cluster jobs may have the wrong size for much of the ACEseq analysis.

There are different variants of these containers:

  * There is a set of containers using versions 1.2.8 and 1.2.10, as described in the ACEseq article. See section :ref:`old-containers`.
  * An updated set of containers based on the Conda environment (see :ref:`here <new-containers>`).
  * ACEseq v1.0.189 was used in the `PanCancer Analysis of Whole Genomes (PCAWG) <https://doi.org/10.1038/s41586-020-1969-6>`_ and is part of the `PCAWG container <https://dockstore.org/containers/quay.io/pancancer/pcawg-dkfz-workflow:2.2.0>`_.


.. _manual way:

The Manual Way
----------------

Roddy Installation
^^^^^^^^^^^^^^^^^^

The file ``buildinfo.txt`` in the ACEseq repository shows you the Roddy API version that you need. We suggest you use the Roddy version with a matching major version number and the highest released minor number. For instance, if the plugin requires Roddy 3.0, use e.g. Roddy 3.6.0. To install that Roddy version please follow the instructions at the `Roddy repository <https://github.com/TheRoddyWMS/Roddy>`_ or the documentation. Note that the installation also requires you to install the "PluginBase" plugin and the "DefaulPlugin" in the versions stated in the ``buildinfo.txt`` files of the plugins in the ``dist/plugins/`` directory of Roddy installation.

With a Roddy installation you can install the ACEseq workflow and its dependencies. To manually resolve the plugin's versions you again need to look at the ``buildinfo.txt``. You start from the ACEseqWorkflow plugin and recurse to its dependencies. Usually, you will probably want to install the released versions of these plugins from the zip-archives that can be found in the Github repositories (e.g. for `ACEseq <https://github.com/DKFZ-ODCF/ACEseqWorkflow/releases>`_) of the respective plugin or component, but you can also use the cloned repositories with the correct release tag (or just commit) checked out.

Note that the ACEseqWorkflow and the COWorkflowBasePlugin can be installed in any directory, as long as all subdirectories there match the pattern for plugin directories of Roddy. So ideally this directory should only contain installations of plugins.

Plugin Installation
^^^^^^^^^^^^^^^^^^^

1. Download the ACEseqWorkflow https://github.com/DKFZ-ODCF/ACEseqWorkflow/tags. Version 1.2.8-4 is used in our production systems. GRCh38/hg38 support is in development and will be available only with version 6. Extract the archive into a directory called "ACEseqWorkflow_$version" with `$version` being the correct version number.
2. Look up the version of the COWorkflowsBasePlugin listed in the "buildinfo.txt" at the root of the ACEseqWorkflow plugin that you just downloaded. Download that version from https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/tags. Create a "COWorkflowsBasePlugin_$version" directory.

You can create the COWorkflowsBasePlugin and ACEseqWorkflow directories in the `dist/plugins` directory, which was extracted with the RoddyEnv ZIP.

Software Stack (Conda)
^^^^^^^^^^^^^^^^^^^^^^

The software stack is the software that the workflow jobs need to do the bioinformatic work. If you run the workflow on a cluster you need to ensure that the software stack is also available on all cluster nodes that execute jobs. Unfortunately, we do not have ACEseq in a workflow manager that supports job execution in containers.

You can try to rebuild the Conda environment based on the ``resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml`` file. Note, however, that some packages may not be available in the referenced Conda channels anymore. Consequently, the following setup based on the YAML file is really more for the developer:


1. Set up the required channels:

    .. code-block:: bash

        conda config --add channels r
        conda config --add channels defaults
        conda config --add channels conda-forge
        conda config --add channels bioconda

2. Rebuild the environment

   .. code-block:: bash

        cd $PATH_TO_PLUGIN_DIRECTORY
        conda env create \
            -n ACEseqWorkflow_1.2.8-4 \
            -f resources/analysisTools/copyNumberEstimationWorkflow/environments/conda.yml

The name of the Conda environment is arbitrary but needs to be consistent with the ``condaEnvironmentName`` variable. You can set the ``condaEnvironmentName`` variable in any of the loaded configuration files (see `Roddy documentation <http://roddy-documentation.readthedocs.io/>`_) or even directly in your Roddy call via ``--cvalues="condaEnvironmentName:$value"``.


Given the problems with old packages in some Conda channels, we offer a work-around. For reproducibility the software stack has been stored with `Conda Pack <https://github.com/conda/conda-pack>`_ in an `archive <old-server>`_.

1. Download the appropriate archive from old-server_ (e.g. ``ACEseqWorkflow_1.2.8-4_conda_4.10.3_x86.tgz``)
2. Unpack and set up the environment
    .. code-block:: bash

        mkdir $envBaseDir/ACEseqWorkflow_1.2.8-4    # The directory name is the environment name.
        cd $envBaseDir/ACEseqWorkflow_1.2.8-4
        source bin/activate
        conda unpack

Now, if you do ``conda env list`` the environment should be listed. If not make sure you installed the environment into ``$CONDA_PREFIX/envs/`` or another `environments directory <https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#specify-environment-directories-envs-dirs>`_ can be configured for your Conda installation.


.. _new-containers:

Docker version
--------------

Different versions of the ACE-seq workflow have been packaged with other workflows and the reference data. These containers have the required software-stacks installed but run all workflow jobs within the same container. To download these pipelines and reference data see ppcg-server_. Instructions for the installation are given in the archives.


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

The reference data can also be downloaded from the ppcg-server_.



Legacy Installations
--------------------

The following installation approaches are kept in the documentation for historical reasons.


.. _aceseq_1.2.10-installation:

Prepackaged files (ACEseq 1.2.10 only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On old-server_ you can find archives for the 1.2.10 plugin version. The prepackaged zip files contains a full Roddy / Plugin setup and include different scripts to install all necessary software and download the required reference files. Currently, we do not intent to update these prepackaged installation files or the Docker version. Note that the Roddy version packaged is not capable of submitting to LSF. Only Torque/PBS is supported.


Stand-alone Roddy for Execution on HTC Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the Roddy-based version of ACEseq please download the pre-packed zip file from the old-server_. Three steps are required to ensure running of ACEseq.

1. Run the ``prepareRoddyInstallation.sh`` script.
2. Download all reference files as specified in the section "Reference files" (below).
3. Set up the Conda environment or install the necessary software as specified in the section "Software" (below).

Before running ACEseq a few parameters need to be adjusted in the configuration files. The output directory is specified in ``$PATH_TO_ACEseq_RODDY_VERSION/configurations/projectsACEseqTest.xml``. Here the variables ``baseDirectoryReference``, ``inputBaseDirectory``, ``outputBaseDirectory``, ``outputAnalysisBaseDirectory`` need to be set. If no SVs should be included the following configuration values (``cvalues``) should be included:

.. code-block:: xml

    <cvalue name='runWithSv' value='false' type="boolean"/>
    <cvalue name='SV' value='false' type="boolean"/>

If you set these values to "true", then "svOutputDirectory" and the SV bedpe filename in the filenames section need to be set, but the SV calls are then used for improved CNV calling.


.. code-block:: xml

    <configurationvalues>
      <cvalue name='svOutputDirectory' value='${outputAnalysisBaseDirectory}/nameOfDirectoryWithSVResults' type="path"/>
    </configurationvalues>

    <filenames package='de.dkfz.b080.co.files' filestagesbase='de.dkfz.b080.co.files.COFileStage'>
       <filename class="TextFile" onMethod="de.dkfz.b080.co.aceseq.ACESeqMethods.mergeSv"
                selectiontag="svFileTag"
                pattern='${svOutputDirectory}/${pid}_svs.bedpe'/>
    </filenames>

Technical specifications are set in the file ``$PATH_TO_ACEseq_RODDY_VERSION/configurations/applicationProperties.ini``. The path to the project.xml and the path to the plugins (``$PATH_TO_ACEseq_RODDY_VERSION/Roddy/dist/plugins/``) need to be set under ``configurationDirectories`` and ``pluginDirectories``. Finally the job manager and execution host need to be set.

Please have a look at the following default ``applicationProperties.ini`` file:

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


.. _old-containers:

Legacy Docker versions
^^^^^^^^^^^^^^^^^^^^^^

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

