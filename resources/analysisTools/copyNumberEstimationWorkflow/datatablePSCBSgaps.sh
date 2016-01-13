#!/bin/bash

set -o pipefail
set -x
#
#PID=$PID
#
#testing=CNV_OUT
#type=eval_job
#
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#
#testing=PSCBS_DATATABLE_OUT
#type=create
#
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2
#
#if [[ "${haplotype}"=='yes' ]]
#then
#	eval "SNP_FILE=${TARGET_DIR}/all.snp.haplo.tab.gz"
#else
#	eval "SNP_FILE=${TARGET_DIR}/all.snp.tab.gz"
#fi
#
#if [[ ! -f ${SNP_FILE}  || ! -s ${SNP_FILE} || ! -r ${SNP_FILE} ]]
# then
#		echo -e "${SNP_FILE} is not a regular, non-empty readable file, maybe haplotype is set to 'yes' but no imputation was run?"
#		exit 1
#fi

${RSCRIPT_BIN}  --vanilla "${TOOL_DATATABLE_AND_PSCBS_GAPS}" \
	      --file_cnv           "${FILENAME_COV_WINDOWS_WG}" \
	      --file_snp           "${FILENAME_SNP_POSITIONS_HAPLO_WG}" \
	      --file_beta          "${FILE_DENSITYBETA}" \
	      --file_sex	       "${FILENAME_SEX}" \
	      --file_knownSegments "${FILE_KNOWNSEGMENTS}" \
	      --file_data          "${FILE_PSCBS_DATA}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while creating knownSegments.txt and pscbs_data.txt;"
	exit 2
fi

${TABIX_BINARY} -f -s 2 -b 5 -e 5 ${FILE_PSCBS_DATA}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in with tabix;"
	exit 2
fi
#
#echo "starting to merge PSCBS breakpoints and crest points"
#
#testing=FILES_TO_EVALUATE
#type=eval_job
#
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of bam files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#
#testing=OUT_MERGE_PSCBS
#type=create
#
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2


if [[ "${CREST}" == 'yes' ]]
then

	${PYTHON_BINARY} "${TOOL_PSCBS_PLUS_CREST}" \
		    --crest_deldupinv "${FILE_CREST_DDI}" \
		    --crest_tx        "${FILE_CREST_TX}" \
		    --known_segments "${FILE_KNOWNSEGMENTS}" \
		    --output          "${FILE_BREAKPOINTS}" \
		    --crest_out       "${FILE_SV_POINTS}" \
		    --DDI_length      "$min_DDI_length"


else
	cp ${FILE_KNOWNSEGMENTS} ${FILE_BREAKPOINTS}
	sed -i '1s/^chr/#chr/' ${FILE_BREAKPOINTS}
	echo "" > "${FILE_SV_POINTS}"

fi
