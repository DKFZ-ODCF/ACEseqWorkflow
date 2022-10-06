#!/usr/bin/env python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).



### idea if segment too small check for sv of prior and following segment if both neighbouring merge with that that has identical sv type
### if copy number closer to one of the segments merge with that one

import argparse
import sys
from python_modules import Tabfile

parser = argparse.ArgumentParser(description ="list all genes with their CNV status")

parser.add_argument( '--file',	'-f', type=file, help="segment file with copy number information" )
parser.add_argument( '--maxlen','-m', type=int,  help='maximum length of segments to be merged', default=50000 )
parser.add_argument( '--maxdist','-d', type=int,  help='maximum distance between segments to be merged', default=50000 )
parser.add_argument( '--cutoff','-c', type=float, default=0.5, help='maximum difference between two segments to be merged (in case they round to the same copy number)')
parser.add_argument( '--out',	'-o', default=sys.stdout, type=str,  help='outputfile' )

args = parser.parse_args()
out = args.out
cutoff = args.cutoff
maxLen = args.maxlen
maxDist = args.maxdist

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

def compare_TCN(TCN1,TCN2,cutoff):
	cond1= abs(TCN1-TCN2) <= cutoff
	cond2= round(TCN1) == round(TCN2)
	if(cond1 and cond2):
		return(True)
	else:
		return(False)

def compare_ACN(a1,a2,cutoff,length1, length2, maxLen=maxLen):
	isna1= a1=="NA"
	isna2= a2=="NA"
	#one or both undefined
	if isna1 or isna2:
		cond1=True
		cond2=True
		if isna1:
			cond1=isna1 and length1 <= maxLen
		if isna2:
			cond2=isna2 and length2 <= maxLen
	#both defined
	else:
		cond1= abs(float(a1)-float(a2) ) <= cutoff
		cond2= round(float(a1)) == round(float(a2))

	if(cond1 and cond2):
		return(True)
	else:
		return(False)


def merge_lines_CN(prior_line, newline, next_line):
	start = "" 
	end   = "" 
	nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)

	while (start != newline["start"] or end != newline["end"]):
		start = newline["start"]
		end   = newline["end"]
		if nextTrue:
			if  float(next_line["length"]) > float(newline["length"]):
				next_line["start"] = newline[ "start" ]
				next_line["tcnMean"] = ( (float(newline["tcnMean"]) * float(newline["length"])) + (float(next_line["tcnMean"]) * float(next_line["length"])) ) /( float(next_line["length"]) + float( newline["length"] ) )
				if( newline["c1Mean"]!= "NA" and next_line["c1Mean"] != "NA" ):
					next_line["c2Mean"] = ( (float(newline["c2Mean"]) * float(newline["length"])) + (float(next_line["c2Mean"]) * float(next_line["length"])) ) /( float(next_line["length"]) + float( newline["length"] ) )
					next_line["c1Mean"] = ( (float(newline["c1Mean"]) * float(newline["length"])) + (float(next_line["c1Mean"]) * float(next_line["length"])) ) /( float(next_line["length"]) + float( newline["length"] ) )
				elif( newline["c1Mean"] == "NA" ):
					newline["c1Mean"] = next_line["c1Mean"]
					newline["c2Mean"] = next_line["c2Mean"]
			else:
				newline["end"] = next_line["end"]
				newline["tcnMean"] = ( (float(newline["tcnMean"]) * float(newline["length"])) + (float(next_line["tcnMean"]) * float(next_line["length"])) ) /(   float(next_line["length"])   + float( newline["length"] ) )
				next_line = newline
			next_line[ "length" ] = int(float( next_line ["end"] ) - float( next_line["start"] ) + 1 )
			newline = next_line
			next_line = infile.readline()
#			if float(newline["length"]) >= 900:
#				break
#			else:
#				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
			if next_line == None:
				nextTrue, priorTrue = [ False, False ]
				break	
			else:
				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)

		elif priorTrue:
			if  float(prior_line["length"]) > float(newline["length"]):
				prior_line[ "end" ] = newline[ "end" ]
				prior_line["tcnMean"] = ( (float(newline["tcnMean"]) * float(newline["length"])) + (float(prior_line["tcnMean"]) * float(prior_line["length"]) )) /( float(prior_line["length"]) + float( newline["length"] ) )
				if( newline["c1Mean"]!= "NA" and prior_line["c1Mean"] != "NA" ):
					prior_line["c2Mean"] = ( (float(newline["c2Mean"]) * float(newline["length"])) + (float(prior_line["c2Mean"]) * float(prior_line["length"])) ) /( float(prior_line["length"]) + float( newline["length"] ) )
					prior_line["c1Mean"] = ( (float(newline["c1Mean"]) * float(newline["length"])) + (float(prior_line["c1Mean"]) * float(prior_line["length"])) ) /( float(prior_line["length"]) + float( newline["length"] ) )
				elif( newline["c1Mean"] == "NA" ):
					newline["c1Mean"] = prior_line["c1Mean"]
					newline["c2Mean"] = prior_line["c2Mean"]
			else:
				newline[ "start" ] = prior_line[ "start" ]
				newline["tcnMean"] = ( (float(newline["tcnMean"]) * float(newline["length"])) + (float(prior_line["tcnMean"]) * float(prior_line["length"]) )) /( float(prior_line["length"]) + float( newline["length"] ) )
				if( newline["c1Mean"]!= "NA" and prior_line["c1Mean"] != "NA" ):
					newline["c2Mean"] = ( (float(newline["c2Mean"]) * float(newline["length"])) + (float(prior_line["c2Mean"]) * float(prior_line["length"])) ) /( float(prior_line["length"]) + float( newline["length"] ) )
					newline["c1Mean"] = ( (float(newline["c1Mean"]) * float(newline["length"])) + (float(prior_line["c1Mean"]) * float(prior_line["length"])) ) /( float(prior_line["length"]) + float( newline["length"] ) )
				elif( newline["c1Mean"] == "NA" ):
					newline["c1Mean"] = prior_line["c1Mean"]
					newline["c2Mean"] = prior_line["c2Mean"]
				prior_line = newline
			prior_line[ "length" ] = int(float( prior_line ["end"] ) - float( prior_line["start"] ))
			newline = next_line
			next_line = infile.readline()
#			if float(newline["length"]) >= 900:
#				break
#			else:
#				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
			if next_line == None:
				nextTrue, priorTrue = [ False, False ]
				break
			else:
				nextTrue, priorTrue = find_closer_match_CN(prior_line, newline, next_line)
	out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n")
	prior_line = newline
	newline = next_line
	return ( [prior_line, newline] )

def find_closer_match_CN(prior_line, newline, next_line):
	"Find closest match according to tcn and sv definition"
	neighbours =  [ next_line, prior_line ]
	#considered for merging if segments are directly neighboured and central segment is either homozygous Deletion or less than 0.3 different from neighbour TCN
	#additional criteria: segment is exceptionally short (<101)
	nextTrue =  ( int(next_line["start"])-int(newline["end"]) < maxDist 
		      and next_line["chromosome"] == newline["chromosome"]
		      and ( ( compare_TCN( float(next_line["tcnMean"]), float(newline[ "tcnMean"]), cutoff) 
		      and   compare_ACN( next_line["c1Mean"], newline["c1Mean"] , cutoff, float(next_line["length"]), float(newline["length"]) ) 
		      and   compare_ACN( next_line["c2Mean"], newline["c2Mean"], cutoff, float(next_line["length"]), float(newline["length"]) ) )
		      or ( float(newline["tcnMean"]) == 0 and float(newline['length']) <= 1000 )
		      or float(newline["length"]) < 101)  )  
	priorTrue = ( int(newline["start"])-int(prior_line["end"]) < maxDist
		      and prior_line["chromosome"] == newline["chromosome"] 
		      and ( ( compare_TCN( float(prior_line["tcnMean"]), float(newline[ "tcnMean"]), cutoff) 
		      and   compare_ACN( prior_line["c1Mean"], newline["c1Mean"] ,cutoff, float(prior_line["length"]), float(newline["length"]) ) 
		      and   compare_ACN( prior_line["c2Mean"], newline["c2Mean"], cutoff, float(prior_line["length"]), float(newline["length"]) ) )
		      or float(newline["tcnMean"]) == 0 ) ) 

#	#additional criterium for small segment with homozygous deletion between two segments of identical TCN
#	if not nextTrue and not priorTrue and prior_line["tcnMean"] == next_line["tcnMean"]:
#		nextTrue = int(next_line["start"])-int(newline["end"]) < 10
#			   and next_line["chromosome"] == newline["chromosome"] 
#			   and float(newline["length"]) < float(next_line["length"])  
#		priorTrue = int(newline["start"])-int(prior_line["end"]) < 10 
#			    and prior_line["chromosome"] == newline["chromosome"] 
#			    and  float(newline["length"]) < float(prior_line["length"])   
#
	if nextTrue and priorTrue:
		
		same_tcn = [ i for i, j in enumerate( neighbours ) if j[ "tcnMean" ] == newline["tcnMean"] ]
		same_sv = [ i for i, j in enumerate( neighbours ) if j[ "SV.Type" ]   == newline["SV.Type"] ]
		sv_defined = [ i for i, j in enumerate( neighbours ) if j[ "SV.Type" ]   != "NA" ]
			
		if len( same_tcn ) == 2:
			if not float( newline["tcnMean"] ) == 0:
				if abs( float(prior_line[ "tcnMean" ]) - float(newline[ "tcnMean" ]) ) > abs( float(next_line[ "tcnMean" ] )  - float(newline[ "tcnMean" ]) ): 
					priorTrue, nextTrue = [ False, True ] 
				else:
					 priorTrue, nextTrue = [ True, False ]
			elif len(same_sv) == 2 or len(same_sv) < 1 or same_sv[0] == 0:
				priorTrue, nextTrue = [ False, True ]
			elif same_sv[0] == 1:
				priorTrue, nextTrue = [ True, False ]
		elif len(same_tcn) > 0:
			if same_tcn[0] == 1:
					priorTrue, nextTrue = [ True, False ]
			elif same_tcn[0] == 0:
					priorTrue, nextTrue = [ False, True ]
		elif len(same_sv) > 0:
			if same_sv[0] == 0:
				priorTrue, nextTrue = [ False, True ]
			elif same_sv[0] == 1:
				priorTrue, nextTrue = [ True, False ]
		elif len(sv_defined) > 0:
			if sv_defined[0] == 0:
				priorTrue, nextTrue = [ False, True ]
			elif sv_defined[0] == 1:
				priorTrue, nextTrue = [ True, False ]

	return( [nextTrue, priorTrue ])

if __name__ == "__main__":
	#read first two lines of file before looping over all lines
	prior_line = infile.readline()
	newline = infile.readline()
	#print header
	out.write(("\t").join(infile.header)+ "\n")
	#loop over all lines and check whether newline is shorter than 900kb and might be merged
	for next_line in infile:
		#check for merging
#		if float(newline["length"]) <= 900:
		if next_line != None:
			prior_line, newline = merge_lines_CN(prior_line, newline, next_line)
		#proceed to next line
#		else:
#			out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n")
#			prior_line, newline = [ newline, next_line ]

	#print last line(s)
	out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
	if newline:
		out.write( "\t".join( str(newline[key]) for key in infile.header ) + "\n" ) 
