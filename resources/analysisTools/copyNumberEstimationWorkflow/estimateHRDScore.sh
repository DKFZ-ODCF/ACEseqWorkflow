#!/usr/bin/sh

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


PIPELINE_DIR=`dirname ${TOOL_HRD_ESTIMATION}`
	if [[ "$?" != 0 ]]
	then
		echo "TOOL_HRD_ESTIMATION not found!" 
		exit 2
	fi

$PYTHON_BINARY ${TOOL_PARSE_JSON} -f ${FILENAME_PARAMETER_JSON} >${FILENAME_PARAMETER_JSON}.tmp
cat ${FILENAME_PARAMETER_JSON}.tmp | while read line
do
	solutions=$line
	for item in ${solutions[@]}; do
		eval $item
	done
	combProFile=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.txt


	mostImportantFile=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.txt
	##remove artifact regions
	combProFileNew=$(echo $combProFile | sed 's/.txt/.smoothed.txt/')
	combProFileNoArtifacts=$(echo $combProFile | sed 's/.txt/.noartifacts.txt/')

	COMBPROFILE_FIFO=${aceseqOutputDirectory}/combProFile_FIFO
	if [[ -p ${COMBPROFILE_FIFO} ]]; then rm ${COMBPROFILE_FIFO}; fi
	mkfifo ${COMBPROFILE_FIFO}

    if [[ ${legacyMode} == 'true' ]]; then
        # replace 'crest' column name by the new name 'SV.Type'
        cat ${COMBPROFILE_FIFO} | awk 'BEGIN{ FS="\t"; OFS="\t" } FNR==1{ for (col=1; col<=NF; ++col) if ($col == "crest") $col = "SV.Type"; } {print}' >$combProFile.tmp &
        pid_legacyMode=$!
    else
        cat ${COMBPROFILE_FIFO} >$combProFile.tmp &
        pid_legacyMode=$!
	fi

	${INTERSECTBED_BINARY} -header -v -f 0.7 \
				     -a $combProFile -b $PIPELINE_DIR/$blacklistFileName \
				     >${COMBPROFILE_FIFO}


	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code intersecting with bedfile" 
		exit 2
	fi
    wait ${pid_legacyMode}; [[ ! $? -eq 0 ]] && echo "Error in legacyMode transcoding process" && exit 10
    rm ${COMBPROFILE_FIFO}

	#smooth Data
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp
    [[ "$?" != 0 ]] && echo "There was a non-zero exit code while removing breakpoints (first time)" && exit 2
	mv $combProFile.tmp.tmp $combProFile.tmp
	${PYTHON_BINARY} ${TOOL_MERGE_ARTIFACTS} -f $combProFile.tmp -o $combProFile.tmp.tmp
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code while merging artifacts" && exit 2
	mv $combProFile.tmp.tmp $combProFile.tmp
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp
	[[ "$?" != 0 ]] && echo "There was a non-zero exit code while removing breakpoints (second time)" && exit 2
	mv $combProFile.tmp.tmp $combProFile.tmp


	(head -1 $combProFile.tmp ; tail -n +2 $combProFile.tmp | sort -k 1,1 -V -k 2,2n ) >$combProFileNoArtifacts

	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while sorting segment file" 
		exit 2
	fi

	#this file could be written out and sorted according to chromosomes
	${PYTHON_BINARY} ${TOOL_SMOOTH_DATA} -f $combProFile.tmp  -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp  && \
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp
	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while smoothing segments" 
		exit 2
	fi


    patientsex=`cat ${FILENAME_SEX:-iDoNotExist.txt}`

	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code getting patient sex" 
		exit 2
	fi

	HRDFile=${aceseqOutputDirectory}/${pid}_HRDscore_${ploidyFactor}_${tcc}.txt
	HRD_DETAILS_FILE=${aceseqOutputDirectory}/${pid}_HRDscore_contributingSegments_${ploidyFactor}_${tcc}.txt
	LST_DETAILS_FILE=${aceseqOutputDirectory}/${pid}_LSTscore_contributingSegments_${ploidyFactor}_${tcc}.CentromerReduced.txt
	MERGED_REDUCED_FILE=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.smoothed.CentromerReduced.txt

	${RSCRIPT_BINARY} ${TOOL_HRD_ESTIMATION} \
		 $combProFileNoArtifacts \
		 ${combProFile}.tmp \
		 $patientsex \
		 $ploidy \
		 $tcc \
		 $pid \
		 ${HRDFile}.tmp \
		 ${HRD_DETAILS_FILE}.tmp \
		 ${LST_DETAILS_FILE}.tmp \
		 ${MERGED_REDUCED_FILE}.tmp \
		 ${FILENAME_CENTROMERES} \
		 ${cytobandsFile} \
		 ${PIPELINE_DIR}


	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while estimating HRD score" 
		exit 2
	fi

	(head -1 $combProFile.tmp ; tail -n +2 $combProFile.tmp | sort -k 1,1 -V -k 2,2n ) \
		>$combProFileNew

	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while sorting segment file" 
		exit 2
	fi

	mv ${HRDFile}.tmp ${HRDFile}
	mv ${HRD_DETAILS_FILE}.tmp ${HRD_DETAILS_FILE}
	mv ${LST_DETAILS_FILE}.tmp ${LST_DETAILS_FILE}
	mv ${MERGED_REDUCED_FILE}.tmp ${MERGED_REDUCED_FILE}
	rm ${combProFile}.tmp
done
if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while processing solutions" 
	exit 2
fi
rm ${FILENAME_PARAMETER_JSON}.tmp 
touch ${FILENAME_CHECKPOINT_HRD}
