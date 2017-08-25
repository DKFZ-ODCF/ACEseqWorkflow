#!/usr/bin/python

# This script merges all segmentation approaches into a final segmentation.

from python_modules import Tabfile
import argparse
import sys

parser=argparse.ArgumentParser(description = "Add SVs detected by other algorithm")
parser.add_argument('--variants',  '-v', type = str, help = "structural variations, combination file of del, dup, inv, tra.")
parser.add_argument('--known_segments','-k', type=str, help = "File containing breakpoints")
parser.add_argument('--sv_out',	   '-s', type = str, help = 'Outfile for SVs')
parser.add_argument('--output',	   '-o', type = str, help = 'Outfile for new breakpoints')
parser.add_argument('--selectCol', '-c', default = "eventscore", type = str, help = 'column name of column to annotate in sv_points.txt')
parser.add_argument('--DDI_length','-l', type = int, help= 'minimum length of del,dup,inv to be considered')

args = parser.parse_args()

selCol = args.selectCol

if not args.known_segments:
        sys.stderr.write("Please specify all known_segments file. For more information use -h\n")
        sys.exit(2)

if not args.sv_out or not args.output:
	sys.stderr.write("Please specify all output files. For more information use -h\n")
        sys.exit(2)

if not args.DDI_length:
	sys.stderr.write("Please specify all minimum duplication deletion and inversion (DDI) lengths. For more information use -h\n")
        sys.exit(2)

try:
	sv_file = Tabfile.Input( open( args.variants, "r" ) )
	sv_out     = open( args.sv_out, "w" )
	knownseg_file = Tabfile.Input( open( args.known_segments, "r" ) )
	file_out   = open( args.output, "w" )
	files = [sv_file]
except IOError as (errno, strerr):
	sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
	sys.exit(errno)

breakpoints = []
chromosomes = [ str(a) for a in range(1,24+1) ]

for f in files:
	for line in f:

		if line['svtype'] == 'INV' or line['svtype'] == 'DUP' or line['svtype'] == 'DEL':

			line["LENGTH"] = str(abs(int(line["start2"])-int(line["start1"]))+1)

			if( int(line["LENGTH"]) >= args.DDI_length and "chrom1" in line ):

				line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )

				#ignore decoy sequences
				if line["chrom1"] in chromosomes :
					sv_out.write( "%s\t%i\t%i\t%s\t%s\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["LENGTH"], line['svtype'], line.get(selCol,"NA")) )


					breakpoints += [ ( str(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( str(line["chrom1"]), int(line["start2"])+1 ) ]

		elif line['svtype'] == 'TRA' or line['svtype']=='INS':

			line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
			line["chrom2"  ] = line["chrom2"  ].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
				
#			if( line["chrom1"] in [ str(i) for i in range( 1, 24+1 ) ] ): # and
			if( line["chrom1"] in [ str(i) for i in range( 1, 24+1 ) ] or line["chrom2"] in [ str(i) for i in range( 1, 24+1 ) ] ): # decoy can be first or second chromosome
#			    line["chrom2"  ] in [ str(i) for i in range( 1, 24+1 ) ] ):

				if line["chrom1"] != line["chrom2"]:

					line["SV_TYPE"] = "CTX"


					sv_out.write( "%s\t%i\tNA\tNA\t%s\t%s\t%i\t%s\n"% (line["chrom1"], int(line["start1"])+1, line["SV_TYPE"], line["chrom2"], int(line["start2"])+1, line.get(selCol, "NA")) )
					sv_out.write( "%s\t%i\tNA\tNA\t%s\t%s\t%i\t%s\n"% (line["chrom2"], int(line["start2"])+1, line["SV_TYPE"], line["chrom1"], int(line["start1"])+1, line.get(selCol, "NA")) )


					if ( line["chrom1"] in [ str(i) for i in range( 1, 24+1 ) ] ):
						breakpoints += [ ( str(line["chrom1"]), int(line["start1"])+1 ) ]

					#don't add second breakpoint to breakpoint file if it maps to decoy
					if ( line["chrom2"] in [ str(i) for i in range( 1, 24+1 ) ] ):
						breakpoints += [ ( str(line["chrom2"]), int(line["start2"])+1 ) ]

				elif line["chrom1"] == line["chrom2"]:

					line["SV_TYPE"] = "ITX"

					sv_out.write( "%s\t%i\t%i\tNA\t%s\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["SV_TYPE"], line.get(selCol,"NA") ) )

					breakpoints += [ ( str(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( str(line["chrom2"]), int(line["start2"])+1 ) ]


for line in knownseg_file:
	if line["start"] != "-Inf":
		breakpoints.append( ( str(line["chromosome"]), int(line["start"]) ) )

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
