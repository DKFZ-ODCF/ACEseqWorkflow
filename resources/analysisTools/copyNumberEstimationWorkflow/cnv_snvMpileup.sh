#!/bin/bash

set -o pipefail
source ${CONFIG_FILE}
set -x

[[ -z ${PARM_CHR_INDEX-} ]] && echo "Variable is missing" && exit -5

tmpFileSnpPos=${FILENAME_SNP_POS}_tmp
tmpFileCovWin=${FILENAME_COV_WIN}_tmp

source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${FILE_TUMOR_BAM} # Sets CHR_PREFIX and REFERENCE_GENOME

CHR_NR=${CHR_PREFIX}${CHR_NAME}

$SAMTOOLS_BINARY mpileup ${MPILEUP_OPTS} \
	-f "${REFERENCE_GENOME}" \
	-r ${CHR_NR} \
	"${FILE_CONTROL_BAM}" "${FILE_TUMOR_BAM}" \
	| ${PYTHON_BINARY} ${TOOL_SNP_POS_CNV_WIN_GENERATOR} \
	-Q $mpileup_qual \
	"${dbSNP_FILE}" \
	- \
	${tmpFileSnpPos} \
        ${tmpFileCovWin}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the snv/cnv mpileup; exiting..."
	exit 2
fi

mv $tmpFileSnpPos ${FILENAME_SNP_POS}
mv $tmpFileCovWin ${FILENAME_COV_WIN}