Alternative Running Modes
===========================

Multiple alternative running modes are enabled with ACEseq. 


Run With "Fake" Control
^^^^^^^^^^^^^^^^^^^^^^^^

We often observed extremely noisy coverage profiles in matched controls from projects
outside the ICGC MMML-Seq, possibly due to wrong handling of blood samples, preventing
accurate copy number calls based on tumor/control ratios. For such samples ACEseq offers
an option to replace the coverage signal from the matched control with an independent
control whilst still maintaining the BAFs of the matched control. This control replacement
option enables full analysis of these sample pairs including reliable discrimination between
runs of homozygosity (ROH) in the germline and somatic loss of heterozygosity (LOH).
Furthermore ACEseq can be run without matched control enlarging the spectrum of samples
that can be processed.

To run the workflow in this mode, the ``runWithFakeControl`` option should be set to "true".

.. code-block:: ini

     <cvalue name="runWithFakeControl" value="true" type="boolean/>
     <cvalue name='MALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path'
             description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for male patients" />
     <cvalue name='FEMALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path'
             description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for female patients" />
     <cvalue name='FAKE_CONTROL_POST' value=".cnv.anno.tab.gz" type='string'
             description="suffix for chromosome wise 1kb coverage files used for fake control workflow"/>

The fake control files should thus be located at ``${*_FAKE_CONTROL_PRE}${chromosome}${FAKE_CONTROL_POST}``.
Each file should be a gzip-compressed TSV with a commented (``#``) header:

::

    #chr    pos     end     normal  tumor   map

Of these columns the ``chr`` and ``pos`` columns are used to combine the analysis results of the tumor
with the "fake" control file. The ``normal`` value from the "fake" control is inserted into the
tumor results file (see ``resources/analysisTools/copyNumberEstimationWorkflow/replaceControlACEseq.R``).

If you are operating at the DKFZ you will find a path prefix to a suitable generic control in the
default configuration of the workflow.


Run Without Control
^^^^^^^^^^^^^^^^^^^^

If no control sample is available, but ACEseq was already used to process
other tumor sample pairs, one of their control coverage profiles can be
used for normalization. In this case, no BAFs can be used from a matching control sample
and also the patient's sex is not inferred.

For the configuration you need to specify the path and prefix to a control coverage profile
for a male and a female patient so it can be matched to the processed sample. To activate this
option the configuration value ``runWithoutControl`` (for versions < 3)
or ``isNoControlWorkflow`` (for versions >= 3) needs to be set to 'true',
either via the command line execution under cvalues or in the project.xml. Furthermore, the
patient's sex needs to be set explicitly with ``PATIENTSEX="male|female|klinefelter"``.


.. code-block:: ini

     <cvalue name="isNoControlWorkflow" value="true" type="boolean
             description="since version 3"/>
     <cvalue name="runWithoutControl" value="true" type="boolean"
             description="up to and including version 2; better use a more recent version!"/>
     <cvalue name="PATIENTSEX" value="male|female|klinefelter" />
     <cvalue name='MALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path'  
             description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for male patients" />
     <cvalue name='FEMALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path' 
             description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for female patients" />
     <cvalue name='FAKE_CONTROL_POST' value=".cnv.anno.tab.gz" type='string'
             description="suffix for chromosome wise 1kb coverage files used for fake control workflow"/>

Note that if run in no-control mode with SV input (you have ``svOutputDirectory`` set), then ACEseq does not expect the SV file to be named ``svs_${PID}_filtered_somatic_minEventScore3.tsv``, like for the tumor/control case, but ``svs_${PID}_filtered_minEventScore3.tsv``. If you use an output directory of the `Sophia workflow <https://github.com/DKFZ-ODCF/SophiaWorkflow>`_ in no-control mode, you can simply symlink the ``svs_${PID}_filtered_minEventScore3.tsv`` to create the "somatic".

Run quality check only
^^^^^^^^^^^^^^^^^^^^^^^

In case you do not want to run the full ACEseq pipeline immediately, 
but would rather access the sample's quality first you can start 
ACEseq with the option "runQualityCheckOnly" set to "true". 

Replace low quality control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a control sample is very noisy and masks CNAs it can be replaced with the coverage profile from a different control of the same sex.
For this run ACEseq with "runWithFakeControl" set to "true" and specify the values "FEMALE_FAKE_CONTROL_PRE" and "MALE_FAKE_CONTROL_PRE" as described in the section for analysis without matched control.


Run with/without SV breakpoint incorporation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To process samples with incorporation of SV breakpoints set the following in the project.xml:

.. code-block:: ini

      <configurationvalues>
        <cvalue name='svOutputDirectory' value='${outputAnalysisBaseDirectory}/nameOfDirectoryWithSVResults' type="path"/>
        <cvalue name='runWithSv' value='true' type="boolean"/>
      </configurationvalues>
    
      <filenames package='de.dkfz.b080.co.files' filestagesbase='de.dkfz.b080.co.files.COFileStage'>
            <filename class="TextFile" onMethod="de.dkfz.b080.co.aceseq.ACESeqMethods.mergeSv"
                      selectiontag="svFileTag"
                      pattern='${svOutputDirectory}/${pid}_svs.bedpe'/>
      </filenames>

If the bedpe file does not exist ACEseq will submit all
steps until the bedpe file is required. A rerun once 
the SV file is generated will start the pipeline up from
the point where SV breakpoints are incorporated.

To process a samples without SVs please set the following in the project.xml:

.. code-block:: ini

    <cvalue name='runWithSv' value='false' type="boolean"/>
    <cvalue name='SV' value='no' type="string"/>


