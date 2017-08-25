#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmpSegmentsPeaks=${FILENAME_SEGMENTS_W_PEAKS}_tmp

${RSCRIPT_BINARY}  --vanilla "${TOOL_ESTIMATE_PEAKS}" \
	      --file     "${FILENAME_ALL_SNP_UPDATE3}" \
	      --gender   "${FILENAME_SEX}" \
	      --segments "${FILENAME_CLUSTERED_SEGMENTS}" \
	      --segOut   "${tmpSegmentsPeaks}" \
	      --out	     "${aceseqOutputDirectory}"
 
if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while generating combi_level.txt, combi.txt and peaks.txt files;" 
	exit 2
fi

mv ${tmpSegmentsPeaks} ${FILENAME_SEGMENTS_W_PEAKS}
