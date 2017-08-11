#!/bin/bash


HAPLO_PATH=`dirname ${PHASED_GENOTYPE_X_FILE}`/${phasedGenotypesFilePrefix}

tmpSnpsHaplo="${FILENAME_SNP_POSITIONS_HAPLO_WG}_tmp"

$PYTHON_BINARY "${TOOL_ADD_HAPLOTYPE}" \
            --inputpath    "${HAPLO_PATH}" \
            --inputsuffix  ".${phasedGenotypesFileSuffix}" \
            --snps         "${FILENAME_SNP_POSITIONS_WG}" \
            --out          "${tmpSnpsHaplo}"

if [[ "$?" != 0 ]]
then
        echo "There was a non-zero exit code in the snp merge and filter process;" 
        exit 2
fi

mv ${tmpSnpsHaplo} ${FILENAME_SNP_POSITIONS_HAPLO_WG}

${TABIX_BINARY} -f -s 1 -b 2 -e 2 ${FILENAME_SNP_POSITIONS_HAPLO_WG}

if [[ "$?" != 0 ]]
then
        echo "There was a non-zero exit code in with tabix;"
        exit 2
fi
