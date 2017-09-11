Alternative Running Modes
===========================

Multiple alternative running modes are enabled with ACEseq. 


Run Without Control
^^^^^^^^^^^^^^^^^^^^

If no control sample is available, but ACEseq was already used to process 
other tumor sample pairs one of their control coverage profile can be
used for normalization. In this case the patient's sex needs to be set
with PATIENTSEX="male|female|klinefelter".

Please specify the path and prefix to a control coverage profile for a male (MALE_FAKE_CONTROL_PRE)
and a female patient (FEMALE_FAKE_CONTROL_PRE) so it can be matched to the processed sample. To 
activate this option the value runWithout control needs to be set to 'true',
either via the command line execution under cvalues or in the project.xml.


.. code-block:: ini

     <cvalue name="runWithoutControl" value="true" type="boolean" />
     <cvalue name="PATIENTSEX" value="male|female|klinefelter" type="boolean" />
     <cvalue name='MALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path'  
                description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for male patients" />  
     <cvalue name='FEMALE_FAKE_CONTROL_PRE' value="pathToPID/${pid}/ACEseq/cnv_snp/${pid}.chr" type='path' 
                description="path and prefix to chromosome-wise 1kb coverage file used for fake control workflow for female patients" /> 
     

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


