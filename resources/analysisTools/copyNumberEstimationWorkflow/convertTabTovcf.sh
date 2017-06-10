#!/usr/bin/sh


source ${CONFIG_FILE}


cd ${aceseqOutputDirectory}

cat ${FILENAME_PURITY_PLOIDY} | grep -v ploidy | while read line
do 
	purity_ploidy=`echo $line | cut -d " " -f 2,3 | sed 's/ /_/' `
	echo $purity_ploidy
	tmp_file_segment_vcf=${FILE_SEGMENT_VCF_PRE}_${purity_ploidy}.tmp
	infoFile=${FILE_MOST_IMPORTANT_INFO_SEG_PRE}${purity_ploidy}${FILE_MOST_IMPORTANT_INFO_SEG_POST}
	targetFile=${FILE_SEGMENT_VCF_PRE}${purity_ploidy}${FILE_SEGMENT_VCF_POST}

#	$PYTHON_BINARY ${TOOL_CONVERT_TO_VCF} \
#		--file ${infoFile} \
#		--out  ${tmp_file_segment_vcf}
##		--id   ${pid} \ #taken out for pancan
#
#	if [[ $? != 0 ]]
#	then
#		echo "Something went while creating the vcf file for ${purity_ploidy}...\n\n"
#		exit 2
#	fi
#
#	${VCFTOOLS_SORT_BINARY} ${tmp_file_segment_vcf} >${tmp_file_segment_vcf}_sorted
#	mv ${tmp_file_segment_vcf}_sorted ${targetFile}
#
#	${BGZIP_BINARY} -f ${targetFile}
#	$TABIX_BINARY -p vcf ${targetFile}.gz
#	rm $tmp_file_segment_vcf
done

if [[ $? != 0 ]]
then
	echo "Something went while creating the gc file...\n\n"
	exit 2
fi	

touch ${FILENAME_CHECKPOINT_VCF}
