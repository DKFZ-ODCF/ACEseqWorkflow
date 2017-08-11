#!/bin/bash


tmpSnpsUpdate=${FILENAME_ALL_SNP_UPDATE1}_tmp

${PYTHON_BINARY} "${TOOL_ADD_SEGMENTS_TO_SNP_DATA}" \
     --pscbs  "${FILE_PSCBS_DATA}" \
     --input  "${FILENAME_SEGMENTS_W_HOMDEL}" \
     --output "${tmpSnpsUpdate}"

mv ${tmpSnpsUpdate} ${FILENAME_ALL_SNP_UPDATE1}

${TABIX_BINARY} -f -s 1 -b 3 -e 4 ${FILENAME_ALL_SNP_UPDATE1}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while creating knownSegments.txt and pscbs_data.txt;" 
	exit 2
fi
