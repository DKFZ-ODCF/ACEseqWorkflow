#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
            --inputpath    "${SNP_VCF_CNV_PATH}" \
            --inputsuffix  ".${CNV_ANNO_SUFFIX}" \
            --output       "${FILENAME_COV_WINDOWS_WG}" \
	        --coverage     ${cnv_min_coverage} \
            --mappability  ${mapping_quality} \
            --NoOfWindows  ${min_windows} \
            > /dev/stderr

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

