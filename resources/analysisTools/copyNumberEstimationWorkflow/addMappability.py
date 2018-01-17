#!/usr/bin/python

### Kortine Kleinheinz 10.02.2014
### k.kleinheinz@dkfz.de
### This script allows to obtain mappability values from the vcf_annotate.pl script using the report --reportBFeatCoord option

import argparse
import gzip
import sys
from python_modules import Tabfile

parser=argparse.ArgumentParser(description="Sum mappability values as output from vcf_annotate.pl")
parser.add_argument( "--infile",   "-i", type=argparse.FileType('r'), default=sys.stdin, help="Tab-seperatd input file for which mappability values should be added up, (default: STDIN)" )
parser.add_argument( "--out",      "-o", default=sys.stdout, help="Name of Output (gzip format)" )
parser.add_argument( "--startCol", "-s", type=str, default='pos', help="Name of column with start coordinates" ) 
parser.add_argument( "--endCol",   "-e", type=str, default='end', help="Name of column with end coordinates" )
parser.add_argument( "--mapCol",   "-m", type=str, default='map', help="Name of column with end coordinates" )

args=parser.parse_args()

out      = args.out
mapCol   = args.mapCol
endCol   = args.endCol
startCol = args.startCol

try:
	infile= Tabfile.Input( args.infile )
except IOError as (errno, strerr):
	sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
	sys.exit(2)

if not mapCol in infile.header:
	sys.stderr.write( "Error: mapCol %s not found in columnNames!\n"% mapCol )
	sys.exit(2)
if not endCol in infile.header:
	sys.stderr.write( "Error: endCol %s not found in columnNames!\n"% endCol )
	sys.exit(2)
if not startCol in infile.header:
	sys.stderr.write( "Error: startCol %s not found in columnNames!\n"% startCol )
	sys.exit(2)

if out != sys.stdout:
	try:
		out = gzip.open(args.out, 'wb')
	except IOError as ( errno, strerr ):
	#	sys.stderr.write( "IOError %i: %s\n"% (errno, strerr) )
		sys.exit("IOError %i: %s\n"% (errno, strerr) )



def addMap(start, end, maps):
	'''get relevant segment length and adjust mappability value'''
	if maps[0] ==".":
		return(0)
	maps   = maps.split("&")
	values = []
	for m in maps:
		[ m_start, m_end, m_map ] = m.split(";")
		m_start= m_start.replace("POS=", "")
		m_end  = m_end.replace("END=", "")
		if float(m_start) < float(start):
			m_start = start
		if float(m_end) > float(end):
			m_end  = end

		length = float(m_end)-float(m_start)+1
		values.append( length*float( m_map ) )
	return(sum(values))

#print header
out.write("#%s\n"% ("\t".join(infile.header ) ) )

for line in infile:
	start 		= line[startCol]
	end   		= line[endCol]
	mappability 	= line[mapCol]
	line[mapCol] 	= addMap(start, end, mappability)
	out.write( "%s\n"% ( "\t".join( [str(line[key]) for key in infile.header ] ) ) )
