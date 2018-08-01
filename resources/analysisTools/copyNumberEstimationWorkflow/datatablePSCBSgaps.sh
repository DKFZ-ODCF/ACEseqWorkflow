#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmp_knownsegments=${FILENAME_KNOWNSEGMENTS}_tmp
tmp_pscbsData=${FILE_PSCBS_DATA}_tmp

${RSCRIPT_BINARY}  --vanilla "${TOOL_DEFINE_BREAKPOINTS}" \
	      --file_cnv           "${FILENAME_GC_CORRECTED_WINDOWS}" \
	      --file_snp           "${FILENAME_SNP_POSITIONS_HAPLO_WG}" \
	      --file_beta          "${FILE_DENSITYBETA}" \
	      --file_sex	       "${FILENAME_SEX}" \
	      --file_knownSegments "${tmp_knownsegments}" \
	      --file_data          "${tmp_pscbsData}" \
	      --libloc             "${libloc_PSCBS}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while creating knownSegments.txt and pscbs_data.txt;" 
	exit 2
fi

mv ${tmp_pscbsData} ${FILE_PSCBS_DATA}

${TABIX_BINARY} -f -s 2 -b 1 ${FILE_PSCBS_DATA}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code with tabix;"
	exit 2
fi

mv ${tmp_knownsegments} ${FILENAME_KNOWNSEGMENTS}
