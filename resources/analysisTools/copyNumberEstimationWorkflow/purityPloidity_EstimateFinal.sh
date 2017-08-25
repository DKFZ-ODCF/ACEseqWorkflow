#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmpPurityPloidy=${FILENAME_PURITY_PLOIDY}_tmp

${RSCRIPT_BINARY} --vanilla "${TOOL_ESTIMATE_PURITY_PLOIDY}" \
	     --segments          "${FILENAME_SEGMENTS_W_PEAKS}" \
	     --file_sex		     "${FILENAME_SEX}" \
	     --purity_ploidy	 "${tmpPurityPloidy}" \
	     --out               "${aceseqOutputDirectory}" \
	     --min_length_purity  $min_length_purity \
	     --min_hetSNPs_purity $min_hetSNPs_purity \
	     --dh_Stop            $dh_stop \
	     --min_length_dh_stop $min_length_dh_stop \
	     --dh_zero            $dh_zero \
	     --purity_min         $purity_min \
	     --purity_max         $purity_max \
	     --ploidy_min         $ploidy_min \
	     --ploidy_max         $ploidy_max \
	     --pid		  $PID


if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while generating ploidy_purity_2D.txt;" 
	exit 2
fi

mv ${tmpPurityPloidy} ${FILENAME_PURITY_PLOIDY}
