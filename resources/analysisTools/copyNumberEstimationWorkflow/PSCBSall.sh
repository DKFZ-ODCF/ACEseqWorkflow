#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmpSegments=${FILENAME_SEGMENTS}_tmp
nocontrol=${isNoControlWorkflow^^}

${RSCRIPT_BINARY} --vanilla "${TOOL_PSCBS_SEGMENTATION}" \
	--file_data         "${FILE_PSCBS_DATA}" \
	--file_breakpoints  "${FILENAME_BREAKPOINTS}" \
	--chrLengthFile     "${CHROMOSOME_LENGTH_FILE}" \
	--file_fit          "${tmpSegments}" \
	--minwidth          $min_seg_width \
	--undo.SD           $undo_SD \
	-h                  $pscbs_prune_height \
	--sv                $SV \
	--libloc            "${libloc_PSCBS}" \
	--nocontrol	    ${nocontrol}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while generating fit.txt file;" 
	exit 2
fi

mv ${tmpSegments} ${FILENAME_SEGMENTS}
