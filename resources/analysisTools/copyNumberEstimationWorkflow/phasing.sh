#!/bin/bash/

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

#DEFINE FILE NAMES
UNPHASED="${FILE_UNPHASED_PRE}${CHR_NAME}.${FILE_VCF_SUF}"
UNPHASED_TWOSAMPLES="${FILE_UNPHASED_PRE}${CHR_NAME}_2samples.${FILE_VCF_SUF}"
PHASED_TWOSAMPLES="${FILE_PHASED_GENOTYPE}${CHR_NAME}_2samples"
tmpPhased="${FILENAME_PHASED_GENOTYPES}_tmp"
tmpHaploblocks="${FILENAME_HAPLOBLOCK_GROUPS}_tmp"

if [[ ${isNoControlWorkflow} == false ]]
then

    source ${TOOL_ANALYZE_BAM_HEADER}
    getRefGenomeAndChrPrefixFromHeader ${FILE_CONTROL_BAM} # Sets CHR_PREFIX and REFERENCE_GENOME

    CHR_NR=${CHR_PREFIX}${CHR_NAME:?CHR_NAME is not set}

         ${SAMTOOLS_BINARY} mpileup ${CNV_MPILEUP_OPTS} -u \
         	    -f "${REFERENCE_GENOME}" \
         	    -r ${CHR_NR} \
         	    "${FILE_CONTROL_BAM}" \
         	    | \
         	    ${BCFTOOLS_BINARY} view ${BCFTOOLS_OPTS} - \
         	    > "${UNPHASED}"
         
         if [[ "$?" != 0 ]]
         then
         	echo "Non zero exit status for mpileup in phasing.sh" >> /dev/stderr
         	exit 2
         fi
fi

echo -n > "${UNPHASED_TWOSAMPLES}"
echo -n > "${tmpPhased}"
echo -n > "${tmpHaploblocks}"
echo -n > "${FILENAME_HAPLOBLOCK_GROUPS}"
echo -n > "${FILENAME_PHASED_GENOTYPES}"

${PYTHON_BINARY} "${TOOL_BEAGLE_CREATE_FAKE_SAMPLES}" \
    --in_file "${UNPHASED}" \
    --out_file "${UNPHASED_TWOSAMPLES}"

if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while creating 2nd sample in vcf-file in phasing.sh" >> /dev/stderr
	exit 2
fi

${JAVA_BINARY} \
    -jar ${TOOL_BEAGLE} \
    gt="${UNPHASED_TWOSAMPLES}" \
    ref="${BEAGLE_REFERENCE_FILE}" \
    out="${PHASED_TWOSAMPLES}" \
    map="${BEAGLE_GENETIC_MAP}" \
    impute=false \
    seed=25041988

if [[ "$?" != 0 ]]
then
    echo "Non zero exit status while phasing with Beagle in phasing.sh" >> /dev/stderr 
    exit 2
fi

${PYTHON_BINARY} "${TOOL_BEAGLE_EMBED_HAPLOTYPES_VCF}" \
    --hap_file "${PHASED_TWOSAMPLES}.vcf.gz" \
    --vcf_file "${UNPHASED}" \
    --out_file  "${tmpPhased}"

if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while embedding haplotypes in phasing.sh" >> /dev/stderr
	exit 2
fi


${PYTHON_BINARY} "${TOOL_GROUP_HAPLOTYPES}" \
	--infile "${tmpPhased}" \
	--out "${tmpHaploblocks}" \
	--minHT ${minHT}
	
if [[ "$?" != 0 ]]
then
	echo "Non zero exit status while grouping haplotypes in phasing.sh" >> /dev/stderr
	exit 2
fi

mv $tmpPhased ${FILENAME_PHASED_GENOTYPES}
mv $tmpHaploblocks ${FILENAME_HAPLOBLOCK_GROUPS}
