#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


#do not set minimum coverage in case no control is used
if [[ ${isNoControlWorkflow} == true ]]
then
	snp_min_coverage=0
fi


tmp_snp_filename=${FILENAME_SNP_POSITIONS_WG}.tmp

${PYTHON_BINARY} "${TOOL_MERGE_FILTER_SNP}" \
            --inputpath    "${SNP_VCF_CNV_PATH}" \
            --inputsuffix  ".${SNP_SUFFIX}" \
            --output       "$tmp_snp_filename" \
            --coverage     ${snp_min_coverage}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the snp merge and filter process;" 
	exit 2
fi

${TABIX_BINARY} -s 1 -b 2 -e 2 $tmp_snp_filename

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the tabix process;" 
	exit 2
fi


mv $tmp_snp_filename ${FILENAME_SNP_POSITIONS_WG}
mv ${tmp_snp_filename}.tbi ${FILENAME_SNP_POSITIONS_WG}.tbi
