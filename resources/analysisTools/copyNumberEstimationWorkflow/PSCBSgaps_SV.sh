#!/bin/bash

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


source ${CONFIG_FILE}


echo "starting to merge PSCBS breakpoints and sv points"

tmp_breakpoints=${FILENAME_BREAKPOINTS}_tmp
tmp_svPoints=${FILENAME_SV_POINTS}_tmp

if [[ "${SV}" == 'yes' ]]
then

	${PYTHON_BINARY} "${TOOL_ADD_SV_TO_PSCBS_GAPS}" \
		    --variants       "${FILENAME_SV}" \
		    --known_segments "${FILENAME_KNOWNSEGMENTS}" \
		    --output         "${tmp_breakpoints}" \
		    --sv_out         "${tmp_svPoints}" \
		    --DDI_length     "$min_DDI_length" \
		    --selectCol	     "${selSVColumn}"

else
	cp ${FILENAME_KNOWNSEGMENTS} ${tmp_breakpoints}
	sed -i '1s/^chr/#chr/' ${tmp_breakpoints}
	echo "" > "${tmp_svPoints}"

fi


if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while merging knownSegments.txt and SVs;"
	exit 2
fi

sleep 5s

mv ${tmp_breakpoints} ${FILENAME_BREAKPOINTS}
mv ${tmp_svPoints} ${FILENAME_SV_POINTS}
