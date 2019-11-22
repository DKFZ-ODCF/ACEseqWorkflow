#!/usr/bin/sh

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).


tmp_corrected_windowfile="$FILENAME_GC_CORRECTED_WINDOWS.tmp"

$RSCRIPT_BINARY "$TOOL_CORRECT_GC_BIAS_R \
	--windowFile "$cnvSnpOutputDirectory/$PID".all.cnv.chr*.tab \
	--timefile   "$REPLICATION_TIME_FILE" \
	--chrLengthFile $CHROMOSOME_LENGTH_FILE" \
	--pid	     "$PID" \
	--email      "$EMAIL" \
	--outfile    "$tmp_corrected_windowfile" \
	--corPlot    "$FILENAME_GC_CORRECT_PLOT" \
	--gcFile     "$GC_CONTENT_FILE" \
	--outDir     "$aceseqOutputDirectory" \
	--lowess_f   "$LOWESS_F" \
	--scaleFactor "$SCALE_FACTOR"



if [[ $? != 0 ]]
then
	echo "Something went wrong during GC correction. Program had non-zero exit status, exiting pipeline...\n\n"
	exit 2
fi	

mv "$tmp_corrected_windowfile" "$FILENAME_GC_CORRECTED_WINDOWS"
