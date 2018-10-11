#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmpClusteredSeg=${FILENAME_CLUSTERED_SEGMENTS}_tmp
tmpSnpsOut=${FILENAME_ALL_SNP_UPDATE2}_tmp


${RSCRIPT_BINARY} --vanilla 		  "${TOOL_MANUAL_PRUNING_SEGMENTS}" \
            	  --file              "${FILENAME_ALL_SNP_UPDATE1}" \
            	  --segments          "${FILENAME_SEGMENTS_W_HOMDEL}" \
	              --functions		  "${TOOL_CLUSTER_FUNCTIONS}" \
            	  --out               "${aceseqOutputDirectory}" \
	              --segOut		      "${tmpClusteredSeg}" \
            	  --min_seg_length    ${min_seg_length_prune} \
            	  --clustering_YN     $clustering \
            	  --min_num_cluster   $min_cluster_number \
	              --min_num_SNPs	  $min_num_SNPs \
            	  --min_membership    $min_membership \
            	  --min_distance      $min_distance \
    	    	  --blockPre		  ${haplogroupFilePath}	\
	              --blockSuf		  ${haplogroupFileSuffix}	\
	              --adjustAlleles	  ${TOOL_ADJUST_ALLELE_ASSIGNMENT} \
	              --newFile		      ${tmpSnpsOut} \
	              --sex		          ${FILENAME_SEX} \
	              --gcCovWidthFile    ${FILENAME_GC_CORRECTED_QUALITY} \
	              --chrLengthFile     ${CHROMOSOME_LENGTH_FILE} \
	              --pid               ${pid} \
	              --libloc            "${libloc_flexclust}" \
	              --runInDebugMode    ${runInDebugMode}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while pruning manually;" 
	exit 2
fi

mv ${tmpClusteredSeg} ${FILENAME_CLUSTERED_SEGMENTS}
mv ${tmpSnpsOut} ${FILENAME_ALL_SNP_UPDATE2}

$TABIX_BINARY -f -s 1 -b 2 -e 2 ${FILENAME_ALL_SNP_UPDATE2}

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while pruning manually;" 
	exit 2
fi
