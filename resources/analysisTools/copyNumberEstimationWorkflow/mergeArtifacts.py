#!/usr/bin/python



### idea if segment too small check for crest of prior and following segment if both neighbouring merge with that that has identical crest type
### if copy number closer to one of the segments merge with that one

import numpy
import argparse
import sys
from python_modules import Tabfile

parser = argparse.ArgumentParser(description ="list all genes with their CNV status")

parser.add_argument( '--file',	'-f', type=file, help="segment file with copy number information" )
parser.add_argument( '--out',	'-o', default=sys.stdout, type=str,  help='outputfile' )
parser.add_argument( '--length','-l', default=900, type=str,  help='outputfile' )

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




def merge_lines_CN(prior_line, newline, next_line):
	start = "" 
	end   = "" 
	#check whether segment could be merged
	nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
	#merge to closest segment
	while (start != newline["start"] or end != newline["end"]):
		start = newline["start"]
		end   = newline["end"]
		#if next segment is closest redefined newline and check whether next line could be merged
		if nextTrue:
			if  float(next_line["length"]) > float(newline["length"]):
				next_line["start"] = newline[ "start" ]
			else:
				newline["end"] = next_line["end"]
				next_line = newline
			next_line[ "length" ] = float( next_line ["end"] ) - float( next_line["start"] )
			newline = next_line
			next_line = infile.readline()
			if float(newline["length"]) >= 900:
				break
			else:
				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
		#else if piro segment is close
		elif priorTrue:
			if  float(prior_line["length"]) > float(newline["length"]):
				prior_line[ "end" ] = newline[ "end" ]
			else:
				newline[ "start" ] = prior_line[ "start" ]
				prior_line = newline
			prior_line[ "length" ] = float( prior_line ["end"] ) - float( prior_line["start"] )
			newline = next_line
			next_line = infile.readline()
			if float(newline["length"]) >= 900:
				break
			else:
				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
	out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n")
	prior_line = newline
	newline = next_line
	return ( [prior_line, newline] )

def find_closer_match_CN(prior_line, newline, next_line):
	"Find closest match according to tcn and crest definition"
	neighbours =  [ next_line, prior_line ]
	#considered for merging if segments are directly neighboured and central segment is either homozygous Deletion or less than 0.3 different from neighbour TCN
	#additional criteria: segment is exceptionally short (<51)
	nextTrue =  ( int(next_line["start"])-int(newline["end"]) < 2
		      and next_line["chromosome"] == newline["chromosome"]
		      and ( abs( float(next_line[ "tcnMean" ] ) - float(newline[ "tcnMean" ])) <= 0.5
		      or float(newline["tcnMean"]) == 0
		      or  float(newline["length"]) < 101)  ) 
	priorTrue = ( int(newline["start"])-int(prior_line["end"]) < 2 
			and prior_line["chromosome"] == newline["chromosome"] 
			and ( abs( float(prior_line[ "tcnMean" ]) - float(newline[ "tcnMean" ])) <= 0.5
			or float(newline["tcnMean"]) == 0 ) )

	#additional criterium for small segment with homozygous deletion between two segments of identical TCN
	#this shouldn't apply anymore as homozygous deletions are detected
#	if not nextTrue and not priorTrue and prior_line["tcnMean"] == next_line["tcnMean"]:
#		nextTrue = int(next_line["start"])-int(newline["end"]) < 2 
#				and next_line["chromosome"] == newline["chromosome"] 
#				and float(newline["length"]) < float(next_line["length"])  
#		priorTrue = int(newline["start"])-int(prior_line["end"]) < 2
#				 and prior_line["chromosome"] == newline["chromosome"] 
#				 and  float(newline["length"]) < float(prior_line["length"])
	if nextTrue and priorTrue:
		tcn_dists = [ abs(float(j[ "tcnMean" ]) - float(newline["tcnMean"] ) ) for i, j in enumerate( neighbours ) ]
		same_crest = [ i for i, j in enumerate( neighbours ) if j[ "crest" ] == newline["crest"] and newline["crest"] != "NA" ]
		crest_defined = [ i for i, j in enumerate( neighbours ) if j[ "crest" ]   != "NA" ]
			
		if tcn_dists[0] == tcn_dists[1] or newline["tcnMean"] == 0:
			if len(same_crest) == 1 :
				if same_crest[0]==0:
					priorTrue, nextTrue = [ False, True ]
				elif same_crest[0] == 1:
					priorTrue, nextTrue = [ True, False ]
				elif len(crest_defined) > 0: 
					if crest_defined[0]==0:
						priorTrue, nextTrue = [ False, True ]
					elif crest_defined[0] == 1:
						priorTrue, nextTrue = [ True, False ]
				else:				
					priorTrue, nextTrue = [ True, False ]
		elif tcn_dists[0] > tcn_dists[1]:
			priorTrue, nextTrue = [False, True]
		elif tcn_dists[0] < tcn_dists[1]:
			priorTrue, nextTrue = [True, False]
	return( [nextTrue, priorTrue ])

if __name__ == "__main__":
	#read first two lines of file before looping over all lines
	prior_line = infile.readline()
	newline = infile.readline()
	out.write(("\t").join(infile.header)+ "\n")

	#loop over all lines and check whether newline is shorter than 900kb and might be merged
	for next_line in infile:
		#check for merging
		if float(newline["length"]) <= 900:
			 prior_line, newline = merge_lines_CN(prior_line, newline, next_line)
		#proceed to next line
		else:
			out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n")
			prior_line, newline = [ newline, next_line ]

	#print last line(s)
	out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
	if newline:
		out.write( "\t".join( str(newline[key]) for key in infile.header ) + "\n" ) 
