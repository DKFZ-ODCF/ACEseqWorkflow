#!/usr/bin/python

# This script replaces segments_to_data.pl and segments_to_data_2.pl.
#
# usage: segments_to_data.py --pscbs [FILE] --input [FILE] --output [FILE]

from python_modules import Tabfile
from python_modules import Options
import sys
import gzip

options = Options.parse( { "pscbs" : str, "input" : str, "output" : str } )

if options:

	pscbsfile = Tabfile.Input( gzip.open( options["pscbs" ] ) )
	outfile   =                gzip.open( options["output"], "w" )

	pscbs_line = pscbsfile.readline()
	while pscbs_line:

		current_chromo = pscbs_line["chromosome"]
		print current_chromo

		infile  = Tabfile.Input( open( options["input"], "r" ) )
		in_line = infile.readline()
		
		while in_line["chromosome"] != current_chromo:
			in_line = infile.readline()
			if in_line==None:
				sys.stderr.write("WARNING: Chr %s not listed in %s"% (str(current_chromo), options["input"]))
				break

		while pscbs_line and pscbs_line["chromosome"] == current_chromo:

			while( in_line and
			       in_line["chromosome"] == current_chromo and
			       float(pscbs_line["x"]) > float(in_line["end"]) ):

				in_line = infile.readline()

			if( in_line and in_line["chromosome"] == current_chromo and
			    float(in_line   ["start"]) <= float(pscbs_line["x"  ]) and
			    float(pscbs_line["x"    ]) <= float(in_line   ["end"]) ):

				outfile.write(
				    pscbs_line["chromosome"]+'\t'+pscbs_line["x"      ]+'\t'+
				    in_line   ["start"     ]+'\t'+in_line   ["end"    ]+'\t'+
				    in_line   ["SV.Type"     ]+'\t'+pscbs_line["CT"     ]+'\t'+
				    pscbs_line["covT"      ]+'\t'+in_line   ["tcnMean"]+'\t'+
				    pscbs_line["betaT"     ]+'\t'+pscbs_line["betaN"  ]+'\t'+
				    pscbs_line["Atumor"    ]+'\t'+pscbs_line["Btumor" ]+'\t'+
				    pscbs_line["Anormal"   ]+'\t'+pscbs_line["Bnormal"]+'\t'+
				    pscbs_line['haplotype' ]+'\t'+in_line   ["map"    ]+'\n' )

#				    in_line   ["map"       ]+'\n' )

			pscbs_line = pscbsfile.readline()

