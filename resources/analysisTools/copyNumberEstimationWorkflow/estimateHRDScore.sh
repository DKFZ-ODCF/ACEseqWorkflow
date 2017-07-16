#!/usr/bin/sh

ACEseq=${aceseqOutputDirectory}
blacklists=${ACEseqBlacklistRegions} #blacklistRegions


FILENAME_PARAMETER_JSON=$1
#$PYTHON_BINARY ${TOOL_PARSE_JSON} -f ${FILENAME_PARAMETER_JSON} | while read line
python parseJson.py -f ${FILENAME_PARAMETER_JSON} | while read line
do
	solutions=$(line)
	for item in ${solutions[@]}; do
		eval $item
	done
	combProFile=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.txt
	mostImportantFile=${aceseqOutputDirectory}/${pid}_comb_pro_extra${ploidyFactor}_${tcc}.txt
	##remove artifact regions
	combProFile.new=`echo $combProFile | sed 's/.txt/.smoothed.txt/'`
	combProFile.noArtifacts=`echo $combProFile | sed 's/.txt/.noartifacts.txt/'`
	${BEDTOOLS_BINARY} intersect -header -v -f 0.7 \
				     -a $combProFile -b $blacklists/artifact.homoDel1500.txt $file | \
	${BEDTOOLS_BINARY} intersect -header -v -f 0.7 \
				     -a - -b $blacklists/potentialArtifacts.txt \
				     >$combProFile.tmp

	#smooth Data
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	${PYTHON_BINARY} ${TOOL_MERGE_ARTIFACTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \

	cp $combProFile.tmp $combProFile.noArtifacts
	#this file could be written out and sorted according to chromosomes
	${PYTHON_BINARY} ${TOOL_SMOOTH_DATA} -f $combProFile.tmp  -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
	${PYTHON_BINARY} ${TOOL_REMOVE_BREAKPOINTS} -f $combProFile.tmp -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp

	patientsex=`cat $SEX_FILE`

	${RSCRIPT_BINARY} ${TOOL_ESTIMATE_HRD} \
		--patientsex $gender \
		--ploidy $ploidy \
		--tcc $tcc \
		--id $pid \
		--segmentfile $combProFile.noArtifacts \
		--mergedfile $combProFile.tmp 

	mv $combProFile.tmp $combProFileNew
done
