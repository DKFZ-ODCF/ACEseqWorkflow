#!/bin/bash


${RSCRIPT_BINARY} ${TOOL_CREATE_BAF_PLOTS} --file_snp  ${FILENAME_SNP_POSITIONS_HAPLO_WG} \
					   --file_sex  ${FILENAME_SEX} \
					   --chrLengthFile ${CHROMOSOME_LENGTH_FILE} \
					   --pid ${PID} \
					   --plot_Dir ${plotOutputDirectory}

if [[ "$?" != 0 ]]
then
        echo "There was an error during BAF plot creation;"
        exit 2
fi

touch ${FILENAME_BAF_PLOT_CHECKPOINT}
