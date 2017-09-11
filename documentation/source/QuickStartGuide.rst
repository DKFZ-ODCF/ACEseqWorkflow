QuickStart
=============

To start ACEseq download package from `here <https://LinkToGitHub.html/>`_ an install the reference files and conda package as described under requirements.

::

    sh $PATH_TO_PLUGIN_DIRECTORY/Roddy/roddy.sh rerun ACEseq@copyNumberEstimation $pid \
    --useconfig=$PATH_TO_PLUGIN_DIRECTORY/applicationProperties.ini \
    --cvalues="bamfile_list:$pathToControlBamFile;$pathToTumorBamFile,sample_list:control;tumor,possibleControlSampleNamePrefixes:control,possibleTumorSampleNamePrefixes:tumor"

Following parameters need to be changed in the project.xml:

- baseDirectoryReference
- outputBaseDirectory
- outputFileGroup (in case all outputfiles should have different group than primary group)

Alternative running modes:

- runWithoutControl (in case it should be run without control)
- runwithFakeControl (in case the coverage should be taken from a different control)
