#!/bin/bash

set -o pipefail
source ${CONFIG_FILE}
set -x

#tmp_breakpoints=${FILENAME_BREAKPOINTS}_tmp
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

#echo "starting to merge PSCBS breakpoints and crest points"
#
#if [[ "${CREST}" == 'yes' ]]
#then
#
#	${PYTHON_BINARY} "${TOOL_ADD_DELLY_TO_PSCBS_GAPS}" \
#		    --dele 	     "${FILE_DELLY_DEL}" \
#		    --tx         "${FILE_DELLY_TX}" \
#		    --dup 	     "${FILE_DELLY_DUP}" \
#		    --inv        "${FILE_DELLY_INV}" \
#		    --known_segments "${FILE_KNOWNSEGMENTS}" \
#		    --output         "${tmp_breakpoints}" \
#		    --sv_out         "${FILENAME_SV_POINTS}" \
#		    --DDI_length     "$min_DDI_length"
#
#else
#	cp ${FILE_KNOWNSEGMENTS} ${tmp_breakpoints}
#	sed -i '1s/^chr/#chr/' ${tmp_breakpoints}
#	echo "" > "${FILE_SV_POINTS}"
#
#fi
#
#if [[ "$?" != 0 ]]
#then
#	echo "There was a non-zero exit code while merging knownSegments.txt and SVs;"
#	exit 2
#fi
#
#mv ${tmp_breakpoints} ${FILENAME_BREAKPOINTS}
#
