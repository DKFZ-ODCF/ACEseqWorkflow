#!/bin/bash
#
source ${CONFIG_FILE}
#
#PID=$PID
#
#testing=CNV_ARRAY
#array=${CHROMOSOME_INDICES[@]}
#type=file_array
#
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of bam files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#
#testing=CNV_OUT
#type=create
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1


${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
            --inputpath    "${SNP_VCF_CNV_PATH}" \
            --inputsuffix  ".${VCF_SUFFIX}" \
            --output       "${FILENAME_COV_WINDOWS_WG}" \
	        --coverage     ${cnv_min_coverage} \
            --mappability  ${mapping_quality} \
            --NoOfWindows  ${min_windows} 

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

