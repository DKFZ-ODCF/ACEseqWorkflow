#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

FILENAME_PARAMETER_JSON_tmp=${FILENAME_PARAMETER_JSON}.tmp

${RSCRIPT_BINARY} 	--vanilla "${TOOL_GENERATE_PLOTS}" \
		--SNPfile "${FILENAME_ALL_SNP_UPDATE3}" \
		--svFile "${FILENAME_SV_POINTS}" \
		--segments "${FILENAME_SEGMENTS_W_PEAKS}" \
		--outfile "${PLOT_PRE}" \
		--chrLengthFile "${CHROMOSOME_LENGTH_FILE}" \
		--outDir "${aceseqOutputDirectory}" \
		--pp "${FILENAME_PURITY_PLOIDY}" \
		--file_sex "${FILENAME_SEX}" \
		--sv_YN $SV \
		--ID ${PID} \
		--pipelineDir `dirname ${TOOL_GENERATE_PLOTS}`

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while generating Plots" 
	exit 2
fi

#create json file for OTP
${PYTHON_BINARY}  ${TOOL_GET_FINAL_PURITY_PLOIDY} \
				--pid ${PID} \
				--path $aceseqOutputDirectory  \
				--out ${FILENAME_PARAMETER_JSON_tmp} \
				--solutionFile ${FILENAME_PURITY_PLOIDY}
if [[ "$?" != 0 ]]
then
	echo "Selection of optimal solution failed"
	exit 2
fi

mv ${FILENAME_PARAMETER_JSON_tmp} ${FILENAME_PARAMETER_JSON}

touch ${FILENAME_CHECKPOINT_PLOTS}
