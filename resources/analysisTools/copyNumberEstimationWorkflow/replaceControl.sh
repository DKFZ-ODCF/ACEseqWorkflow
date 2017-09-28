#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


sex=`cat ${FILENAME_SEX}`


declare -a CHROMOSOME_INDICES="$CHROMOSOME_INDICES"

for CHR_NR in "${CHROMOSOME_INDICES[@]}"
 do

	if [[ ${sex} == "female" ]]
	then
		FAKE_CONTROL="${FEMALE_FAKE_CONTROL_PRE}${CHR_NR}${FAKE_CONTROL_POST}"
	elif [[ ${sex} == "male" ]]
	then
		FAKE_CONTROL="${MALE_FAKE_CONTROL_PRE}${CHR_NR}${FAKE_CONTROL_POST}"
	else
		echo "No valid control file for klinefelter patients!"
		exit 2
	fi

	O_FILE="${SNP_VCF_CNV_PATH}${CHR_NR}.${CNV_ANNO_SUFFIX}"
	tmp_out="${O_FILE}_tmp"

	${RSCRIPT_BINARY} ${TOOL_FAKE_CONTROL} \
		-b ${O_FILE}  \
		-f ${FAKE_CONTROL} \
		-o ${tmp_out}

	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code while incorporating a fake control for chromosome $CHR_NR."
		exit 2
	fi

	mv $tmp_out ${O_FILE}
done

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while incorporating a fake control."
	exit 2
fi

touch ${FILENAME_CHECKPOINT_FAKECONTROL}
