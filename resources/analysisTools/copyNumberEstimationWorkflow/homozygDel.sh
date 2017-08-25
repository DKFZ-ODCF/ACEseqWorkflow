#!/bin/bash


source ${CONFIG_FILE}


#annotate_vcf & addMappability: add mappability values to fit.txt segments
#homozygous_del: add sv information (DEL,CTX...) and mappability classification (mappable, unmappable or homozygDel (if no tcn defined) )

#add "#" to header of fit.txt
sed -i '1s/^chr/#chr/' ${FILENAME_SEGMENTS}

tmpSegsHomDel=${FILENAME_SEGMENTS_W_HOMDEL}_tmp

${PERL_BINARY} "${TOOL_ANNOTATE_CNV_VCF}" \
	-a "${FILENAME_SEGMENTS}" \
	--aFileType=custom \
	--aChromColumn chromosome \
	--aPosColumn "start" \
	--aEndColumn end \
	-b "${MAPPABILITY_FILE}" \
	--bFileType=bed \
	--reportBFeatCoord \
	--columnName map \
	--chromXtr 23:X \
	--chromYtr 24:Y \
        | \
${PYTHON_BINARY} "${TOOL_ADD_MAPPABILITY}" \
	-s start \
	-e end \
	-m map \
	| \
${PERL_BINARY} "${TOOL_ADD_HOMOZYGOUS_DELETION}" \
	--a "${FILENAME_SV_POINTS}" \
	--b $min_segment_map \
    | \
${BGZIP_BINARY} > "${tmpSegsHomDel}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while obtaining homozygous deletions;" 
	exit 2
fi

mv ${tmpSegsHomDel} ${FILENAME_SEGMENTS_W_HOMDEL}

$TABIX_BINARY -f -s 1 -b 2 -e 3 "${FILENAME_SEGMENTS_W_HOMDEL}"
 
#
#if [[ "$?" != 0 ]]
#then
#	echo "There was a non-zero exit code while obtaining homozygous deletions;"
#	exit 2
#fi
