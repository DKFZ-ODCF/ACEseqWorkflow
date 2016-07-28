#!/usr/bin/python

# This script replaces all.snp.pl.
#
# usage: merge_and_filter_snp.py --inputpath [PATH] --inputsuffix [SUFFIX] --output [FILE] --coverage [INT]
#
# This script generates a file name ( inputpath + chromosome + inputsuffix ) for
# each chromosome. 
# For example "~/Data/patient.chr" and ".snp" result in
# "~/Data/patient.chr1.snp", ... "~/Data/patient.chr22.snp",
# "~/Data/patient.chrX.snp", "~/Data/patient.chrY.snp".
# It than merges these files into one output file with some filtering for
# coverage and a randomization of the A/B alleles.
#
# The functionality is described in the Bachelor thesis of Isabell Bludau at
# page 17, lines 11-17.

import gzip
import subprocess

from python_modules import Tabfile
from python_modules import Options

options = Options.parse( { "inputpath" : str, "inputsuffix" : str,
                           "output"    : str, "coverage"    : int } )
if options:

	outfile = subprocess.Popen("bgzip >%s"% options["output"], shell=True, stdin=subprocess.PIPE)
	
	outfile.stdin.write("chr\tstartPos\tAnormal\tBnormal\tAtumor\tBtumor\thaplotype\n")		#header

	for chromo in [ str( n ) for n in range( 1, 22+1 ) ] + [ "X", "Y" ]:

		infile = gzip.open( options["inputpath"] + chromo + options["inputsuffix"], 'rb' )

		for line in Tabfile.Input( infile ):

			if( int(line["An"]) + int(line["Bn"]) >= options["coverage"] ) :

				line['haplotype'] = 0 

				if line["chr"].startswith( "chr" ):
					line["chr"] = line["chr"].replace("chr", "")


				line["chr"] = line["chr"].replace( "X", "23" );
				line["chr"] = line["chr"].replace( "Y", "24" );

				outfile.stdin.write(
				    line["chr"]+'\t'+line["pos"]+'\t'+line["An"]+'\t'+
				    line["Bn" ]+'\t'+line["At" ]+'\t'+line["Bt"]+'\t'+
				    str(line['haplotype'])      +'\n' )

