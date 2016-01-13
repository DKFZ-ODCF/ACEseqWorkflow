#!/bin/bash

set -o pipefail
source ${CONFIG_FILE}
set -x

#estimate sex of patient from X and Y coverage
Y_FILE="${SNP_VCF_CNV_PATH}Y.${CNV_SUFFIX}"
X_FILE="${SNP_VCF_CNV_PATH}X.${CNV_SUFFIX}"

tmp_sex_file=${FILENAME_SEX}_tmp

${RSCRIPT_BINARY} ${TOOL_ESTIMATE_SEX} \
	 --file_dataY ${Y_FILE} \
	 --file_dataX ${X_FILE} \
	 --file_size ${CHROMOSOME_LENGTH_FILE} \
	 --cnv_files "${SNP_VCF_CNV_PATH}*.${CNV_SUFFIX}" \
	 --min_Y_ratio ${min_Y_ratio} \
	 --min_X_ratio ${min_X_ratio} \
	 --file_out ${tmp_sex_file}


if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the gender estimation" 
	exit 2
fi

mv $tmp_sex_file ${FILENAME_SEX}

for CHR_NR in "${CHROMOSOME_INDICES[@]}"
 do

  A_FILE="${SNP_VCF_CNV_PATH}${CHR_NR}.${CNV_SUFFIX}"
  echo ${SNP_VCF_CNV_PATH}
  O_FILE="${SNP_VCF_CNV_PATH}${CHR_NR}.${VCF_SUFFIX}"
  tmp_out="${O_FILE}_tmp"

  ${PERL_BINARY} "${TOOL_ANNOTATE_CNV_VCF}" \
	-a ${A_FILE} \
	--aFileType=custom \
	--aChromColumn chr \
	--aPosColumn pos \
	--aEndColumn end \
	-b "${MAPPABILITY_FILE}" \
	--bFileType=bed \
	--reportBFeatCoord \
	--columnName map |
  ${PYTHON_BINARY} ${TOOL_ADD_MAPPABILITY} \
	-o "${tmp_out}"

  mv $tmp_out ${O_FILE}

 done

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while adding the mappability values."
	exit 2
fi

 

