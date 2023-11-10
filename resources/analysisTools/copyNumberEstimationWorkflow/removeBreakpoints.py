#!/usr/bin/env python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

import numpy
import argparse
import sys
from python_modules import Tabfile

parser = argparse.ArgumentParser(description ="list all genes with their CNV status")

parser.add_argument( '--file',	'-f', type=file, help="segment file with copy number information" )
parser.add_argument( '--out',	'-o', default=sys.stdout, type=str,  help='outputfile' )

args = parser.parse_args()
out = args.out

try:
	infile= Tabfile.Input( args.file )
except IOError as (errno, strerr):
        sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
        sys.exit(2)

if out != sys.stdout:
        try:
                out = open(args.out, 'w')
        except IOError as ( errno, strerr ):
                sys.exit("IOError %i: %s\n"% (errno, strerr) )


line = infile.readline()
out.write(("\t").join(infile.header)+ "\n")

for newline in infile:
	if (int(newline["start"])-int(line["end"]) <= 2
	    and line["chromosome"] == newline["chromosome"] 
	    and line[ "cluster" ] == newline[ "cluster" ] 
	    and line["cluster"] != "NA" ): 

		line["end"] = newline["end"]
		line[ "length" ] = int( line ["end"] ) - int( line["start"] )
	else:
		out.write( "\t".join( str(line[key]) for key in infile.header ) + "\n")
		line = newline

out.write( "\t".join( str(line[key]) for key in infile.header ) + "\n" ) 
