#!/usr/bin/env python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

# This script merges all segmentation approaches into a final segmentation.

from python_modules import Tabfile
from python_modules import Options

options = Options.parse( { "crest_deldupinv" : str, "crest_tx"   : str,
                           "known_segments"  : str, "output"     : str,
                           "crest_out"       : str, "DDI_length" : int } )
if options:

	crest_ddi_file = Tabfile.Input( open( options["crest_deldupinv"], "r" ) )
	crest_tx_file  = Tabfile.Input( open( options["crest_tx"       ], "r" ) )
	crest_out      = open( options["crest_out"], "w" )
	file_out       = open( options["output"   ], "w" )
	
	
	breakpoints = []

	for line in crest_ddi_file:
	
		line["LENGTH"] = str(int(line["END"])-int(line["POS"])+1)
		
		if( line["SOMATIC_GERMLINE_CLASSIFICATION"] == "somatic" and
		    int(line["LENGTH"]) >= options["DDI_length"] and "CHROM" in line ):

			line["CHROM"] = line["CHROM"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )

			crest_out.write( line["CHROM"]+"\t"+line["POS"]+"\t"+
			                 line["END"]+ "\t"+line["LENGTH"]+"\t"+
			                 line["SV_TYPE"]+"\tNA\tNA\tNA\n" )

			breakpoints += [ ( int(line["CHROM"]), int(line["POS"]) ) ]
			breakpoints += [ ( int(line["CHROM"]), int(line["END"]) ) ]

	for line in crest_tx_file:

		line["CHROM_FROM"] = line["CHROM_FROM"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
		line["CHROM_TO"  ] = line["CHROM_TO"  ].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
			
		if( line["SOMATIC_GERMLINE_CLASSIFICATION"] == "somatic" and
		    line["CHROM_FROM"] in [ str(i) for i in range( 1, 24+1 ) ] and
		    line["CHROM_TO"  ] in [ str(i) for i in range( 1, 24+1 ) ] ):

			if line["SV_TYPE"] == "CTX":

				crest_out.write( line["CHROM_FROM"]+"\t"+line["POS_FROM"]+"\t"+
				                 "NA\tNA"          +"\t"+line["SV_TYPE" ]+"\t"+
				                 line["CHROM_TO"  ]+"\t"+line["POS_TO"  ]+"\tNA\n" )

				crest_out.write( line["CHROM_TO"  ]+"\t"+line["POS_TO"  ]+"\t"+
				                 "NA\tNA"          +"\t"+line["SV_TYPE" ]+"\t"+
				                 line["CHROM_FROM"]+"\t"+line["POS_FROM"]+"\tNA\n" )

				breakpoints += [ ( int(line["CHROM_FROM"]), int(line["POS_FROM"]) ) ]
				breakpoints += [ ( int(line["CHROM_TO"  ]), int(line["POS_TO"  ]) ) ]

			elif line["SV_TYPE"] == "ITX":

				crest_out.write( line["CHROM_FROM"]+"\t"+line["POS_FROM"]+"\t"+
				                 line["POS_TO"    ]+"\t"+"NA"            +"\t"+
				                 line["SV_TYPE"   ]+"\t"+"NA\tNA\tNA\n" )
				                 
				breakpoints += [ ( int(line["CHROM_FROM"]), int(line["POS_FROM"]) ) ]
				breakpoints += [ ( int(line["CHROM_TO"  ]), int(line["POS_TO"  ]) ) ]

	knownseg_file = Tabfile.Input( open( options["known_segments"], "r" ) );

	for line in knownseg_file:
		if line["start"] != "-Inf":
			breakpoints.append( ( int(line["chromosome"]), int(line["start"]) ) )

	breakpoints.sort();
	
	if breakpoints:
		file_out.write( str(breakpoints[0][0]  )+"\t-Inf\t"+
		                str(breakpoints[0][1]-1)+"\tInf\n" )

	for i in range( 0, len(breakpoints)-1 ):
		if breakpoints[i][0] != breakpoints[i+1][0]:
		
			file_out.write( str(breakpoints[i  ][0]  )+"\t"+
			                str(breakpoints[i  ][1]  )+"\tInf\tInf\n" )
			file_out.write( str(breakpoints[i+1][0]  )+"\t-Inf\t"+
			                str(breakpoints[i+1][1]-1)+"\tInf\n" )
		
		elif breakpoints[i][1] != breakpoints[i+1][1]:
		    
			file_out.write( str(breakpoints[i  ][0]  )+"\t"+
			                str(breakpoints[i  ][1]  )+"\t"+
			                str(breakpoints[i+1][1]-1)+"\t"+
			                str(breakpoints[i+1][1]-breakpoints[i][1]-1)+"\n" )

	if breakpoints:
		file_out.write( str(breakpoints[-1][0])+"\t"+
		                str(breakpoints[-1][1])+"\tInf\tInf\n" )
