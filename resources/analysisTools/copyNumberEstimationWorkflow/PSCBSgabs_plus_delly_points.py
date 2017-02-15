#!/usr/bin/python

# This script replaces PSCBSgabs_plus_CRESTpoints.pl.
#
# This script merges all segmentation approaches into a final segmentation.

from python_modules import Tabfile
import argparse
import sys

parser=argparse.ArgumentParser(description = "Add SVs detected by DELLY")
parser.add_argument('--variants',  '-v', type = str, help = "DELLY structural variations, combination file of del, dup, inv, tra... used if single files are not provided")
parser.add_argument('--dele',	   '-d', type = str, help = "DELLY deletions")
parser.add_argument('--dup',	   '-p', type = str, help = "DELLY duplications")
parser.add_argument('--inv',	   '-i', type = str, help = "DELLY inversions")
parser.add_argument('--tx',	   '-t', type = str, help = 'DELLY translocations')
parser.add_argument('--known_segments','-k', type=str, help = "File containing breakpoints")
parser.add_argument('--sv_out',	   '-s', type = str, help = 'Outfile for delly SVs')
parser.add_argument('--output',	   '-o', type = str, help = 'Outfile for new breakpoints')
parser.add_argument('--selectCol', '-c', default = "eventscore", type = str, help = 'column name of column to annotate in sv_points.txt')
parser.add_argument('--DDI_length','-l', type = int, help= 'minimum length of del,dup,inv to be considered')

args = parser.parse_args()
filenumber = 4		#svs split up into 4 bedpe files
if not args.dele  or not args.dup or not args.inv or not args.tx:
	if not args.variants:
	        sys.stderr.write("Please specify all DELLY file(s). For more information use -h\n")
        	sys.exit(2)
	else:
		filenumber = 1	#svs combined in one bedpe file

if not args.known_segments:
        sys.stderr.write("Please specify all known_segments file. For more information use -h\n")
        sys.exit(2)

if not args.sv_out or not args.output:
	sys.stderr.write("Please specify all output files. For more information use -h\n")
        sys.exit(2)

if not args.DDI_length:
	sys.stderr.write("Please specify all minimum duplication deletion and inversion (DDI) lengths. For more information use -h\n")
        sys.exit(2)

if filenumber==4:
	try:
		delly_del_file = Tabfile.Input( open( args.dele, "r" ) )
		delly_dup_file  = Tabfile.Input( open( args.dup, "r" ) )
		delly_inv_file = Tabfile.Input( open( args.inv, "r" ) )
		delly_tx_file = Tabfile.Input( open( args.tx, "r" ) )
		knownseg_file = Tabfile.Input( open( args.known_segments, "r" ) )
		sv_out      = open( args.sv_out, "w" )
		file_out       = open( args.output, "w" )
		files = [ delly_del_file, delly_dup_file, delly_inv_file, delly_tx_file ]
	except IOError as (errno, strerr):
		sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
		sys.exit(errno)
else:
	try:
		delly_file = Tabfile.Input( open( args.variants, "r" ) )
		sv_out     = open( args.sv_out, "w" )
		knownseg_file = Tabfile.Input( open( args.known_segments, "r" ) )
		file_out   = open( args.output, "w" )
		files = [delly_file]
	except IOError as (errno, strerr):
		sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
		sys.exit(errno)

breakpoints = []
chromosomes = [ str(a) for a in range(1,24+1) ]

for f in files:
	for line in f:
		if line['svtype'] == 'DEL':
			line["LENGTH"] = str(int(line["start2"])-int(line["start1"])+1)

			if( int(line["LENGTH"]) >= args.DDI_length and "chrom1" in line ):

				line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
				#ignore decoy sequences
				if line["chrom1"] in chromosomes :

					sv_out.write( "%s\t%i\t%i\t%s\tDEL\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["LENGTH"], line['id']) )

					breakpoints += [ ( int(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( int(line["chrom1"]), int(line["start2"])+1 ) ]

		elif line['svtype'] == 'DUP':

			line["LENGTH"] = str(int(line["start2"])-int(line["start1"])+1)

			if( int(line["LENGTH"]) >= args.DDI_length and "chrom1" in line ):

				line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )

				#ignore decoy sequences
				if line["chrom1"] in chromosomes :


					sv_out.write( "%s\t%i\t%i\t%s\tDUP\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["LENGTH"], line['id']) )

					breakpoints += [ ( int(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( int(line["chrom1"]), int(line["start2"])+1 ) ]

		elif line['svtype'] == 'INV':

			line["LENGTH"] = str(int(line["start2"])-int(line["start1"])+1)

			if( int(line["LENGTH"]) >= args.DDI_length and "chrom1" in line ):

				line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )

				#ignore decoy sequences
				if line["chrom1"] in chromosomes :


					sv_out.write( "%s\t%i\t%i\t%s\tINV\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["LENGTH"], line['id']) )

					breakpoints += [ ( int(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( int(line["chrom1"]), int(line["start2"])+1 ) ]

		elif line['svtype'] == 'TRA':

			line["chrom1"] = line["chrom1"].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
			line["chrom2"  ] = line["chrom2"  ].replace( "chr", "" ).replace( "X", "23" ).replace("Y", "24" )
				
			if( line["chrom1"] in [ str(i) for i in range( 1, 24+1 ) ] and
			    line["chrom2"  ] in [ str(i) for i in range( 1, 24+1 ) ] ):

				if line["chrom1"] != line["chrom2"]:

					line["SV_TYPE"] = "CTX"


					sv_out.write( "%s\t%i\tNA\tNA\t%s\t%s\t%i\t%s\n"% (line["chrom1"], int(line["start1"])+1, line["SV_TYPE"], line["chrom2"], int(line["start2"])+1, line['id']) )

					sv_out.write( "%s\t%i\tNA\tNA\t%s\t%s\t%i\t%s\n"% (line["chrom2"], int(line["start2"])+1, line["SV_TYPE"], line["chrom1"], int(line["start1"])+1, line['id']) )

					breakpoints += [ ( int(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( int(line["chrom2"]), int(line["start2"])+1 ) ]

				elif line["chrom1"] == line["chrom2"]:

					line["SV_TYPE"] = "ITX"

					sv_out.write( "%s\t%i\t%i\tNA\t%s\tNA\tNA\t%s\n"% (line["chrom1"], int(line["start1"])+1, int(line["start2"])+1, line["SV_TYPE"], line['id']) )

					breakpoints += [ ( int(line["chrom1"]), int(line["start1"])+1 ) ]
					breakpoints += [ ( int(line["chrom2"]), int(line["start2"])+1 ) ]


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
