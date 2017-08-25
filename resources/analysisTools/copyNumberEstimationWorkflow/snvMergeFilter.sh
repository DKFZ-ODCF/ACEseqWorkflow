#!/bin/bash


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

#do not set minimum coverage in case no control is used
if [[ ${runWithoutControl} == true ]]
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
