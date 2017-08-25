#!/usr/bin/python

# This script replaces segments_to_data.pl and segments_to_data_2.pl.
#
# usage: segments_to_data.py --pscbs [FILE] --input [FILE] --output [FILE]

from python_modules import Tabfile
from python_modules import Options
import subprocess
import gzip
import sys

options = Options.parse( { "pscbs" : str, "input" : str, "output" : str } )


if options:

	pscbsfile = Tabfile.Input( gzip.open( options["pscbs" ] ) )

	#SNPs could be gzipped or not
	try:
		outfile   = subprocess.Popen("bgzip >%s"% options["output"], shell=True, stdin=subprocess.PIPE)
	except IOError as (errno, strerror):
		syst.stderr.write( "I/O error (%i): %s\n"% (errno,strerror) )

	pscbs_line = pscbsfile.readline()
	while pscbs_line:

		current_chromo = pscbs_line["chromosome"]
		print current_chromo

		#segments might be zipped
		if options["input"][-2:]=='gz':	
			try:	
				infile  = Tabfile.Input( gzip.open( options["input"], "r" ) )
			except IOError as (errno, strerror):
				sys.stderr.write( "I/O error (%i): %s\n"% (errno, strerror) )
				sys.exit(errno)
		else:
			try:	
				infile  = Tabfile.Input(      open( options["input"], "r" ) )
			except IOError as (errno, strerror):
				sys.stderr.write( "I/O error (%i): %s\n"% (errno, strerror) )
				sys.exit(errno)
	
		in_line = infile.readline()
		
		while in_line["chromosome"] != current_chromo:
			in_line = infile.readline()

		while pscbs_line and pscbs_line["chromosome"] == current_chromo:

			while( in_line and
			       in_line["chromosome"] == current_chromo and
			       float(pscbs_line["x"]) > float(in_line["end"]) ):

				in_line = infile.readline()

			if( in_line and in_line["chromosome"] == current_chromo and
			    float(in_line   ["start"]) <= float(pscbs_line["x"  ]) and
			    float(pscbs_line["x"    ]) <= float(in_line   ["end"]) ):

				outfile.stdin.write(
				    pscbs_line["chromosome"]+'\t'+pscbs_line["x"      ]+'\t'+
				    in_line   ["start"     ]+'\t'+in_line   ["end"    ]+'\t'+
				    in_line   ["SV.Type"     ]+'\t'+pscbs_line["CT"     ]+'\t'+
				    pscbs_line["covT"      ]+'\t'+in_line   ["tcnMean"]+'\t'+
				    pscbs_line["betaT"     ]+'\t'+pscbs_line["betaN"  ]+'\t'+
				    pscbs_line["Atumor"    ]+'\t'+pscbs_line["Btumor" ]+'\t'+
				    pscbs_line["Anormal"   ]+'\t'+pscbs_line["Bnormal"]+'\t'+
				    pscbs_line['haplotype' ]+'\t'+in_line   ["map"    ]+'\n' )

			pscbs_line = pscbsfile.readline()

outfile.communicate()


