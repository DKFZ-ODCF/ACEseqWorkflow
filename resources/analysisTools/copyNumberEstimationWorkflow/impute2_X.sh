#!/bin/bash

source ${CONFIG_FILE}




source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${FILE_CONTROL_BAM} # Sets CHR_PREFIX and REFERENCE_GENOME


CHR_NAME=X
CHR_NR=${CHR_PREFIX}${CHR_NAME}

#DEFINE FILENAMES
PHASED_GENOTYPE=${FILE_PHASED_GENOTYPE}${CHR_NAME}.${FILE_TXT_SUF}
PHASED_HAPS=${PHASED_GENOTYPE}_haps
PHASED_HAPS_CONF=${PHASED_GENOTYPE}_haps_confidence
PHASED_INFO=${PHASED_GENOTYPE}_info
PHASED_INFO_SAMPLE=${PHASED_GENOTYPE}_info_by_sample
PHASED_SUMMARY=${PHASED_GENOTYPE}_summary
PHASED_WARNINGS=${PHASED_GENOTYPE}_warnings
tmpphased=${FILENAME_PHASED_GENOTYPES}_tmp #These two files should have 23 as chromosomes name rather than 'X'
tmphaploblocks=${FILENAME_HAPLOBLOCK_GROUPS}_tmp
UNPHASED="${FILE_UNPHASED_PRE}${CHR_NAME}.${FILE_VCF_SUF}"
FILENAME_UNPHASED_GENOTYPE=${FILE_UNPHASED_GENOTYPE}${CHR_NAME}.${FILE_TXT_SUF}


#check whether the patient is female or male
if grep -Pv 'female|klinefelter'  "${FILENAME_SEX}"
 then
   echo "Patient $PID is male."
   echo " " >"${FILENAME_HAPLOBLOCK_GROUPS}"
   echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_control_${PID}"  >"${FILENAME_PHASED_GENOTYPES}"
   exit 0
 fi

if [[ ${runWithoutControl} == false ]]
then


        ${SAMTOOLS_BINARY} mpileup ${MPILEUP_OPTS} -u \
        	    -f "${REFERENCE_GENOME}" \
        	    -r ${CHR_NR} \
        	    "${FILE_CONTROL_BAM}" \
        	    | \
        	    ${BCFTOOLS_BINARY} view ${BCFTOOLS_OPTS} - \
        	    > "${UNPHASED}"
        
        if [[ "$?" != 0 ]]
        then
        	echo "Non zero exit status for mpileup in impute2.sh"
        	exit 2
        fi
fi

${PYTHON_BINARY} "${TOOL_EXTRACT_GENOTYPE_VCF}" \
    --vcf_file "${UNPHASED}" \
    --outfile  "${FILENAME_UNPHASED_GENOTYPE}"

if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while extracting genotype in impute2.sh"
	exit 2
fi

	#create sample_g file
	echo "ID_1 ID_2 missing sex" > "${FILE_SAMPLE_G}"
	echo "0 0 0 D" >> "${FILE_SAMPLE_G}"
	echo "${PID} ${PID} 0 2" >> "${FILE_SAMPLE_G}"


	# This runs the impute2 haplotype imputation program
	# strand_g option necessary?
	echo -n> "${PHASED_HAPS}"
	echo -n> "${PHASED_HAPS_CONF}"
	echo -n> "${PHASED_INFO}"
	echo -n> "${PHASED_INFO_SAMPLE}"
	echo -n> "${PHASED_SUMMARY}"
	echo -n> "${PHASED_WARNINGS}"

	# For smaller chromosomes, impute2 will simply not print output for the
	# last few segments.
	for SEGMENT in $(seq 0 50)
	do
		PHASED_GENOTYPE_PART=${FILE_PHASED_GENOTYPE}${CHR_NAME}.part${SEGMENT}.${FILE_TXT_SUF}
		PHASED_HAPS_PART=${PHASED_GENOTYPE_PART}_haps
		PHASED_HAPS_CONF_PART=${PHASED_GENOTYPE_PART}_haps_confidence
		PHASED_INFO_PART=${PHASED_GENOTYPE_PART}_info
		PHASED_INFO_SAMPLE_PART=${PHASED_GENOTYPE_PART}_info_by_sample
		PHASED_SUMMARY_PART=${PHASED_GENOTYPE_PART}_summary
		PHASED_WARNINGS_PART=${PHASED_GENOTYPE_PART}_warnings

		echo -n> "${PHASED_GENOTYPE_PART}"
		echo -n> "${PHASED_HAPS_PART}"
		echo -n> "${PHASED_HAPS_CONF_PART}"
		echo -n> "${PHASED_INFO_PART}"
		echo -n> "${PHASED_INFO_SAMPLE_PART}"
		echo -n> "${PHASED_SUMMARY_PART}"
		echo -n> "${PHASED_WARNINGS_PART}"

		${TOOL_IMPUTE} \
		    -seed 25041988 \
		    -phase \
		    -m $(echo "${GENETIC_MAP_FILE_X}" | sed "s/\${CHR_NR}/${CHR_NAME}/g") \
		    -h $(echo "${KNOWN_HAPLOTYPES_FILE_X}" | sed "s/\${CHR_NR}/${CHR_NAME}/g") \
		    -l $(echo "${KNOWN_HAPLOTYPES_LEGEND_FILE_X}" | sed "s/\${CHR_NR}/${CHR_NAME}/g") \
		    -g "${FILENAME_UNPHASED_GENOTYPE}" \
		    -int $[5000000*${SEGMENT}] $[5000000*${SEGMENT} + 4999999] \
		    -Ne 20000 \
		    -o "${PHASED_GENOTYPE_PART}"

		if [[ "$?" != 0 ]]
			then
			echo "WARNING: Non zero exit status during segmentation of segment ${SEGMENT} on chr ${CHR_NAME} in impute2.sh"
			exit 2
		fi


		cat "${PHASED_HAPS_PART}" \
		    >> "${PHASED_HAPS}"
		cat "${PHASED_HAPS_CONF_PART}" \
		    >> "${PHASED_HAPS_CONF}"
		cat "${PHASED_INFO_PART}" \
		    >> "${PHASED_INFO}"
		cat "${PHASED_INFO_SAMPLE_PART}" \
		    >> "${PHASED_INFO_SAMPLE}"
		cat "${PHASED_SUMMARY_PART}" \
		    >> "${PHASED_SUMMARY}"
		cat "${PHASED_WARNINGS_PART}" \
		    >> "${PHASED_WARNINGS}"

		rm  "${PHASED_GENOTYPE_PART}" \
		    "${PHASED_HAPS_PART}" \
		    "${PHASED_HAPS_CONF_PART}" \
		    "${PHASED_INFO_PART}" \
		    "${PHASED_INFO_SAMPLE_PART}" \
		    "${PHASED_SUMMARY_PART}" \
		    "${PHASED_WARNINGS_PART}"
	done

${PYTHON_BINARY} "${TOOL_EMBED_HAPLOTYPES_VCF}" \
	    --hap_file "${PHASED_HAPS}" \
	    --vcf_file "${UNPHASED}" \
	    --outfile  "${tmpphased}"

if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while embedding haplotypes in impute2.sh"
	exit 2
fi

${PYTHON_BINARY} "${TOOL_GROUP_HAPLOTYPES}" \
	--infile "${tmpphased}" \
	--out "${tmphaploblocks}" \
	--minHT ${minHT}
	
if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while grouping haplotypes in impute2.sh"
	exit 2
fi

mv ${tmpphased} ${FILENAME_PHASED_GENOTYPES}
mv ${tmphaploblocks} ${FILENAME_HAPLOBLOCK_GROUPS}
