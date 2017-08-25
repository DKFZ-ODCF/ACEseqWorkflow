#!/bin/bash


zcat $dbSNP_FILE | perl $TOOL_ANNOTATE_CNV_VCF  \
		-a - -b ${FILENAME_SNP_POSITIONS_WG_FAKE} \
		--bFileType vcflike \
		--chromXtr X:23 \
		--chromYtr Y:24 \
		--columnName genotype \
		--aColNameLineStart "#CHROM" | \
		grep -v "#"  | \
		$PERL_BINARY $TOOL_PARSE_VCF ${imputeOutputDirectory} $unphasedGenotypesFilePrefix $unphasedGenotypesFileSuffix 


if [[ $? != 0 ]]
then
	echo "creation of unphased files failed"
	exit 2
fi
