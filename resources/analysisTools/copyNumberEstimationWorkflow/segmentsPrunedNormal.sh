#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmpUpdatedSnps=${FILENAME_ALL_SNP_UPDATE3}_tmp

${PYTHON_BINARY} "${TOOL_ADD_SEGMENTS_TO_SNP_DATA}" \
	--pscbs  "${FILENAME_ALL_SNP_UPDATE2}" \
	--input  "${FILENAME_CLUSTERED_SEGMENTS}" \
	--output "${tmpUpdatedSnps}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while creating all_seg_2.txt"
	exit 2
fi

mv ${tmpUpdatedSnps} ${FILENAME_ALL_SNP_UPDATE3}


${TABIX_BINARY} -f -s 1 -b 3 -e 4  "${FILENAME_ALL_SNP_UPDATE3}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code using bgzip or tabix on all_seg_2.txt" 
	exit 1 
fi
