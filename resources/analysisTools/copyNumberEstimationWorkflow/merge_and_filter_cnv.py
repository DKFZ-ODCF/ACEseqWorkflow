#!/usr/bin/env python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

# This script replaces all.cnv.pl.
#
# usage: merge_and_filter_cnv.py --inputpath [PATH] --inputsuffix [SUFFIX] --output [FILE] --coverage [INT]
#
# This script generates a file name ( inputpath + chromosome + inputsuffix ) for
# each chromosome. 
# For example "~/Data/patient.chr" and ".snp" result in
# "~/Data/patient.chr1.snp", ... "~/Data/patient.chr22.snp",
# "~/Data/patient.chrX.snp", "~/Data/patient.chrY.snp".
# It than merges these files into one output file while filtering for coverage
# and combining the 1k windows into 10k windows.

import gzip
from python_modules import Tabfile
from python_modules import Options

options = Options.parse( { "inputpath"   : str,   "inputsuffix" : str,
                           "output"      : str,   "coverage"    : int,
                           "mappability" : float, "NoOfWindows" : int } )
if options:

	outfile = gzip.open( options["output"], "wb" )

	def  process_accumulated_lines( lines ):

		if len( lines ) >= options["NoOfWindows"]:

			chromo = lines[0]["chr"]
			
			if chromo.startswith( "chr" ):
				chromo = chromo.replace("chr", "")
				
			chromo = chromo.replace( "X", "23" )
			chromo = chromo.replace( "Y", "24" )
			
			position = ( int( lines[0]["pos"] ) // 10000 ) * 10000 + 1
			
			normal = sum( [ float(line["normal"]) for line in lines ] ) * 10 / len( lines )
			tumor  = sum( [ float(line["tumor" ]) for line in lines ] ) * 10 / len( lines )
			
			# This function emulates the exact formating that perl uses for
			# float values. It is of no functional importance for this script.
			def perl_like_float_format( value ):
				int_part = str(int(value))
				format_string="{0:."+str(15-len(int_part))+"f}"
				return format_string.format(value).rstrip('0').rstrip('.')

			outfile.write( chromo+'\t'+str(position)+'\t'+
			               perl_like_float_format(normal)+'\t'+
			               perl_like_float_format(tumor )+'\n' ) 
		

	for chromo in [ str( n ) for n in range( 1, 22+1 ) ] + [ "X", "Y" ]:
		print chromo
		input_filename = options["inputpath"] + chromo + options["inputsuffix"]
		
		infile = Tabfile.Input( gzip.open( input_filename , 'rb') )

		accum_lines = [ infile.readline() ]

		for line in infile:

			if( int(line["normal"]) < options["coverage"   ] or
			    float(line["map" ]) < options["mappability"] ):
				continue
				
			if( line["chr"] == accum_lines[0]["chr"] and
			    int(line["pos"])//10000 == int(accum_lines[0]["pos"])//10000 ):
				accum_lines.append( line )
			else:
				process_accumulated_lines( accum_lines )
				accum_lines = [ line ]

		process_accumulated_lines( accum_lines )

