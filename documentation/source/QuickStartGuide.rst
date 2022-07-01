QuickStart
=============

To start ACEseq download package from `here <https://LinkToGitHub.html/>`_ and install the reference files and conda package as described under :ref:`installation`.

::

    sh $PATH_TO_PLUGIN_DIRECTORY/Roddy/roddy.sh \
        rerun \
        ACEseq@copyNumberEstimation \
        $pid \
        --useconfig=$PATH_TO_PLUGIN_DIRECTORY/applicationProperties.ini \
        --cvalues="bamfile_list:$pathToControlBamFile;$pathToTumorBamFile,sample_list:control;tumor,possibleControlSampleNamePrefixes:control,possibleTumorSampleNamePrefixes:tumor"

Following parameters should be changed in the project.xml:

- baseDirectoryReference
- outputBaseDirectory
- outputFileGroup (in case all outputfiles should have different group than primary group)

Alternative running modes:

- runWithoutControl/isNoControlWorkflow (in case it should be run without control)
- runWithFakeControl (in case the coverage should be taken from a different control)
