#!/usr/bin/sh

set -o pipefail
source ${CONFIG_FILE}
set -x

cd ${aceseqOutputDirectory}

name=${FILE_MOST_IMPORTANT_INFO_SEG_PRE}
find -name "$name*" | while read line
do 
	purity_ploidy=`echo $line | sed "s/\/$name//" | sed 's/\.txt//' `
	echo $purity_ploidy
	tmp_file_segment_vcf=${FILE_SEGMENT_VCF_PRE}.tmp
	targetFile=${FILE_SEGMENT_VCF_PRE}${purity_ploidy}${FILE_SEGMENT_VCF_POST}

	$PYTHON_BINARY ${TOOL_CONVERT_TO_VCF} \
		--file ${line} \
		--out  ${tmp_file_segment_vcf}
#		--id   ${pid} \ #taken out for pancan

	if [[ $? != 0 ]]
	then
		echo "Something went while creating the gc file for ${purity_ploidy}...\n\n"
		exit 2
	fi

	${VCFTOOLS_SORT_BINARY} ${tmp_file_segment_vcf} >${tmp_file_segment_vcf}_sorted
	mv ${tmp_file_segment_vcf}_sorted ${targetFile}

	${BGZIP_BINARY} -f ${targetFile}
	$TABIX_BINARY -p vcf ${targetFile}.gz
done

if [[ $? != 0 ]]
then
	echo "Something went while creating the gc file...\n\n"
	exit 2
fi	

touch ${FILENAME_CHECKPOINT_VCF}
