#!/bin/bash

set -o pipefail
source ${CONFIG_FILE}
set -x


#testing=FILES_SNP
#array=${CHROMOSOME_INDICES[@]}
#type=file_array

#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1


#testing=SNP_OUT
#type=create

#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2


${PYTHON_BINARY} "${TOOL_MERGE_FILTER_SNP}" \
            --inputpath    "${SNP_VCF_CNV_PATH}" \
            --inputsuffix  ".${SNP_SUFFIX}" \
            --output       "${FILENAME_SNP_POSITIONS_WG}" \
            --coverage     ${snp_min_coverage}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the snp merge and filter process;" 
	exit 2
fi
