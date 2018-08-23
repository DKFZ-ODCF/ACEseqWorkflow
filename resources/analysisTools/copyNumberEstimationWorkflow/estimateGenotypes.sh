#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

#change header and index file

#get pseudo control counts
${RSCRIPT_BINARY} $TOOL_GET_HET_SNPS -s ${FILENAME_SNP_POSITIONS_WG} | bgzip >${FILENAME_SNP_POSITIONS_WG_FAKE}.tmp

if [[ $? != 0 ]]
then
	echo "genotype estimation failed"
	exit 2
fi

${TABIX_BINARY} -h -f -s 1 -b 2 -e 2 ${FILENAME_SNP_POSITIONS_WG_FAKE}.tmp

if [[ $? != 0 ]]
then
	echo "genotype estimation failed"
	exit 2
fi

mv ${FILENAME_SNP_POSITIONS_WG_FAKE}.tmp ${FILENAME_SNP_POSITIONS_WG_FAKE}
mv ${FILENAME_SNP_POSITIONS_WG_FAKE}.tmp.tbi ${FILENAME_SNP_POSITIONS_WG_FAKE}.tbi

