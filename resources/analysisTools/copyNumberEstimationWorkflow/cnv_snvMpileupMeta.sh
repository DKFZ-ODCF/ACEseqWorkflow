#!/bin/bash

source $CONFIG_FILE


outpath=${aceseqOutputDirectory}

# TODO Assumption #C_IND % 2 == 0!
# Take the nth and n - 1
CNT=${#CHROMOSOME_INDICES[@]} # 24 by default
MAX=`expr ${CNT} / 2 - 1` # Arrays are indexed from 0 to (including) 23
for n in `seq 0 $MAX`
do
    export forwardIndex=$n
    export reverseIndex=`expr $CNT - 1 - $n` # Get n - 1

    export CHR_INDEX_FRONT=${CHROMOSOME_INDICES_OPTIMIZED_CNV_SNV_MPILEUP[$forwardIndex]}
    export CHR_INDEX_BACK=${CHROMOSOME_INDICES_OPTIMIZED_CNV_SNV_MPILEUP[$reverseIndex]}
    export FILE_FRONT_SNP_POS=${outpath}/${pid}.chr${CHR_INDEX_FRONT}.${SNP_SUFFIX}
    export FILE_FRONT_COV_WIN=${outpath}/${pid}.chr${CHR_INDEX_FRONT}.${CNV_SUFFIX}
    export FILE_BACK_SNP_POS=${outpath}/${pid}.chr${CHR_INDEX_BACK}.${SNP_SUFFIX}
    export FILE_BACK_COV_WIN=${outpath}/${pid}.chr${CHR_INDEX_BACK}.${CNV_SUFFIX}
    export FILE_FRONT_CP=${FILE_FRONT_SNP_POS}.cp
    export FILE_BACK_CP=${FILE_BACK_SNP_POS}.cp

    if [[ ! -f ${FILE_FRONT_SNP_POS} ]] && [[ ! -f ${FILE_BACK_SNP_POS} ]]
    then
    # load snvCalling script twice with a wait after another.
        (FILENAME_SNP_POS=${FILE_FRONT_SNP_POS} FILENAME_COV_WIN=${FILE_FRONT_COV_WIN} FILENAME_VCF_SNVS_CHECKPOINT=${FILE_FRONT_CP} PARM_CHR_INDEX=${CHR_INDEX_FRONT} bash ${TOOL_CNV_SNP_GENERATION} ; \
         FILENAME_SNP_POS=${FILE_BACK_SNP_POS}  FILENAME_COV_WIN=${FILE_BACK_COV_WIN}  FILENAME_VCF_SNVS_CHECKPOINT=${FILE_BACK_CP}  PARM_CHR_INDEX=${CHR_INDEX_BACK}  bash ${TOOL_CNV_SNP_GENERATION}) &
    else
        echo "Both files exist, skipping jobs for ${CHR_INDEX_FRONT}, ${CHR_INDEX_BACK}: (${FILE_FRONT_SNP_POS}, ${FILE_BACK_SNP_POS}, ${FILE_FRONT_COV_WIN}, ${FILE_BACK_COV_WIN})";
    fi
done

wait
