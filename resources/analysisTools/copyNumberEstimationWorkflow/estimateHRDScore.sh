#!/usr/bin/sh

source ${CONFIG_FILE}

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
	HRDFile=${aceseqOutputDirectory}/${pid}_HRDscore_${ploidyFactor}_${tcc}.txt
	echo $combProFile

	mostImportantFile=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.txt
	##remove artifact regions
	combProFileNew=$(echo $combProFile | sed 's/.txt/.smoothed.txt/')
	combProFileNoArtifacts=$(echo $combProFile | sed 's/.txt/.noartifacts.txt/')
	${INTERSECTBED_BINARY} -header -v -f 0.7 \
				     -a $combProFile -b $PIPELINE_DIR/$blacklistFileName \
				     >$combProFile.tmp
	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code intersecting with bedfile" 
		exit 2
	fi


	#smooth Data
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	${PYTHON_BINARY} ${TOOL_MERGE_ARTIFACTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while removing breakpoints" 
		exit 2
	fi


	(head -1 $combProFile.tmp ; tail -n +2 $combProFile.tmp | sort -k 1,1 -V -k 2,2n ) \
		>$combProFileNoArtifacts

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


	${RSCRIPT_BINARY} ${TOOL_HRD_ESTIMATION} \
		 $combProFileNoArtifacts \
		 $combProFile.tmp \
		 $gender \
		 $ploidy \
		 $tcc \
		 $pid \
		 $HRDFile.tmp \
		 $PIPELINE_DIR/${centromerFilename} \
		 $PIPELINE_DIR 

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

	mv $HRDFile.tmp $HRDFile
	rm $combProFile.tmp 
done
if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while processing solutions" 
	exit 2
fi
rm ${FILENAME_PARAMETER_JSON}.tmp 
touch $FILENAME_CHECKPOINT_HRD 
