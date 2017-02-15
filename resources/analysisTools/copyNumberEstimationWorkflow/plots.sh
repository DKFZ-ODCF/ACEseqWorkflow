#!/bin/bash

set -o pipefail
source ${CONFIG_FILE}
set -x 

${RSCRIPT_BINARY} 	--vanilla "${TOOL_GENERATE_PLOTS}" \
		--SNPfile "${FILENAME_ALL_SNP_UPDATE3}" \
		--crestFile "${FILENAME_SV_POINTS}" \
		--segments "${FILENAME_SEGMENTS_W_PEAKS}" \
		--outfile "${PLOT_PRE}" \
		--chrLengthFile "${CHROMOSOME_LENGTH_FILE}" \
		--outDir "${aceseqOutputDirectory}" \
		--pp "${FILENAME_PURITY_PLOIDY}" \
		--file_sex "${FILENAME_SEX}" \
		--crest_YN $CREST \
		--ID ${PID} \
		--pipelineDir `dirname ${TOOL_GENERATE_PLOTS}`

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while generating Plots" 
	exit 2
fi

#create json file for OTP
optimalSolution=`python ${TOOL_GET_FINAL_PURITY_PLOIDY} $PID $aceseqOutputDirectory`
${FILENAME_PURITY_PLOIDY} $optimalSolution

touch ${FILENAME_CHECKPOINT_PLOTS}
