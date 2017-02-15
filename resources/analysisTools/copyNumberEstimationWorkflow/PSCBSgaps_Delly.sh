#!/bin/bash


source ${CONFIG_FILE}


echo "starting to merge PSCBS breakpoints and crest points"

tmp_breakpoints=${FILENAME_BREAKPOINTS}_tmp
tmp_svPoints=${FILENAME_SV_POINTS}_tmp

if [[ "${CREST}" == 'yes' ]]
then

	${PYTHON_BINARY} "${TOOL_ADD_DELLY_TO_PSCBS_GAPS}" \
		    --variants       "${FILENAME_DELLY_SV}" \
		    --known_segments "${FILENAME_KNOWNSEGMENTS}" \
		    --output         "${tmp_breakpoints}" \
		    --sv_out         "${tmp_svPoints}" \
		    --DDI_length     "$min_DDI_length" \
		    --selectCol	     "${selSVColumn}"

else
	cp ${FILENAME_KNOWNSEGMENTS} ${tmp_breakpoints}
	sed -i '1s/^chr/#chr/' ${tmp_breakpoints}
	echo "" > "${svPoints_tmp}"

fi


if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while merging knownSegments.txt and SVs;"
	exit 2
fi

sleep 5s

mv ${tmp_breakpoints} ${FILENAME_BREAKPOINTS}
mv ${tmp_svPoints} ${FILENAME_SV_POINTS}
