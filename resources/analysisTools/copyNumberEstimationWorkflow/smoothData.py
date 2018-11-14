#!/usr/bin/python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).



### idea if segment too small check for sv of prior and following segment if both neighbouring merge with that that has identical sv type
### if copy number closer to one of the segments merge with that one

import argparse
import sys
from python_modules import Tabfile

parser = argparse.ArgumentParser(description ="list all genes with their CNV status")

parser.add_argument( '--file',	'-f', type=file, help="segment file with copy number information" )
parser.add_argument( '--out',	'-o', default=sys.stdout, type=str,  help='outputfile' )
parser.add_argument( '--maxDistToNext',	'-m', default=1000, type=int,  help='maximum allowed distance between segments to be merged' )
# maxDistToNext is used only used for first segments [in chromosome | after large gap]


# comment by warsow on 2018-09-27: TCN levels are never compared when merging segments (compare_TCN/compare_ACN is never used...)
# Schlesner and Warsow agreed on not taking TCN/ACN levels into account as otherwise focal events may introduce gaps which lead
# to separating of otherwise merged segments which in turn increases the HRD score although most likely only 1 LOH event took place

args = parser.parse_args()
maxDistToNext=args.maxDistToNext
out = args.out
cutoff = 0.5
maxLen = 3000000

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

def compare_ACN(a1,a2,cutoff,length1, length2, maxLen=3000000):
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
	if nextTrue:
		if  float(next_line["length"]) > float(newline["length"]):
			next_line["start"] = newline[ "start" ]
		else:
			newline["end"] = next_line["end"]
			next_line = newline
		next_line[ "length" ] = int(float( next_line ["end"] ) - float( next_line["start"] ) + 1 )
		newline = next_line

	elif priorTrue:
		if  float(prior_line["length"]) > float(newline["length"]):
			prior_line[ "end" ] = newline[ "end" ]
		else:
			newline[ "start" ] = prior_line[ "start" ]
			prior_line = newline
		prior_line[ "length" ] = int(float( prior_line ["end"] ) - float( prior_line["start"] ) + 1 )
		newline = next_line

	return ( [prior_line, newline] )

def find_closer_match_CN(prior_line, newline, next_line):
	# function had comment that TCN counts are considered here, but no evidence of TCN consideration was found...
	# in agreement with comments at the beginning of this script, TCN counts will not be taken into account for merging decision
	neighbours =  [ next_line, prior_line ]
	#considered for merging if segments are directly neighboured
	diffPriorTcn = abs( float(prior_line["tcnMean"]) - float(newline["tcnMean"]) )
	diffNextTcn = abs( float(next_line["tcnMean"]) - float(newline["tcnMean"]) )
	if( prior_line["c1Mean"] == "NA" or newline["c1Mean"] == "NA" ):
		diffPriorC1 = 10000	
		diffPriorC2 = 10000
	else:
		diffPriorC1 = abs( float(prior_line["c1Mean"]) - float(newline["c1Mean"]) )
		diffPriorC2 = abs( float(prior_line["c2Mean"]) - float(newline["c2Mean"]) )

	if ( newline["c1Mean"] == "NA" or next_line["c1Mean"] == "NA" ) :
		diffNextC1 = 10000	
		diffNextC2 = 10000
	else:
		diffNextC1 = abs( float(next_line["c1Mean"]) - float(newline["c1Mean"]) )
		diffNextC2 = abs( float(next_line["c2Mean"]) - float(newline["c2Mean"]) )

	lengthPrior = float(prior_line["length"]) > float(next_line["length"])

	nextTrue  =  ( next_line["chromosome"] == newline["chromosome"]  and int(next_line["start"]) - int(newline["end"]) < 2 )
	priorTrue =  ( prior_line["chromosome"] == newline["chromosome"] and int(newline["start"]) - int(prior_line["end"]) < 2 )

	if (priorTrue and nextTrue):
		nextTrue  =  ( diffPriorTcn >= diffNextTcn  )
		priorTrue =  ( diffPriorTcn <= diffNextTcn )

	if (priorTrue and nextTrue):
		nextTrue  =  ( diffPriorC2 >= diffNextC2  )
		priorTrue =  ( diffPriorC2 <= diffNextC2 )
	if (priorTrue and nextTrue):
		nextTrue  =  (  diffPriorC1 >= diffNextC1 )
		priorTrue =  (  diffPriorC1 <= diffNextC1 )
	if (priorTrue and nextTrue):
		nextTrue  =  ( not lengthPrior  )
	if (priorTrue and nextTrue):
		nextTrue = False

	return( [nextTrue, priorTrue ])

def merge_gap( priorGap_line, afterGap_line ):
		distance = (float(afterGap_line["start"]) - float(priorGap_line["end"]))
		if ( distance > maxLen ):
			return( [priorGap_line, afterGap_line] )

		priorGap_line["end"]  = int(int(priorGap_line["end"]) + round( (distance/2)-0.5 ) )
		priorGap_line["length"]  = int(int(priorGap_line["end"]) - int(priorGap_line["start"]) + 1 )
		afterGap_line["start"]= int(int(afterGap_line["start"]) - round( distance/2 ) )
		afterGap_line["length"]  = int(int(afterGap_line["end"]) - int(afterGap_line["start"]) + 1 )
		return( [priorGap_line, afterGap_line] )

def merge_first_segment(prior_line, newline):
		if( float(prior_line["length"]) < float( newline["length"] ) ):
			newline["start"] = prior_line["start"]
			prior_line = newline
		else:
			prior_line["end"] = newline["end"]
		prior_line["length"] = int(int(prior_line["end"]) - int(prior_line["start"]) +1 )
		return(prior_line)

	

if __name__ == "__main__":
	#read first two lines of file before looping over all lines
	prior_line = infile.readline()
	newline = infile.readline()
	priorChrom = 0

	#merge first segment in case it is too short
	while( float(prior_line["length"]) <= maxLen  and int(newline["start"])-int(prior_line["end"]) <= maxDistToNext):
		prior_line=merge_first_segment(prior_line, newline)
		newline = infile.readline()

	#print header
	out.write(("\t").join(infile.header)+ "\n")

	#loop over all lines and check whether newline is shorter than 900kb and might be merged
	for next_line in infile:
		#merge first segments per chromosome in case it is too short
		# also merge first segment after long gap (e.g. centromere)
		# last segments before long gap will be handled via merge_lines_CN
		if ( newline != None and (prior_line["chromosome"] != newline["chromosome"] or (int(newline["start"]) - int(prior_line["end"]) > maxLen)) ):
			out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
			prior_line = newline
			newline = next_line
			while( newline is not None and float(prior_line["length"]) <= maxLen and int(newline["start"])-int(prior_line["end"]) <= maxDistToNext):
				prior_line=merge_first_segment(prior_line, newline)
				newline = infile.readline()
			continue


		#last segment can only be merged to single prior
		if ( next_line != None and next_line["chromosome"] != newline["chromosome"] ) :
				if( float(newline["length"]) <= maxLen and int(newline["start"])-int(prior_line["end"]) <= maxDistToNext):
					prior_line=merge_first_segment(prior_line, newline)
					newline=next_line
					continue
				else:
					out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
					prior_line=newline
					newline = next_line 
					continue

		#check for merging
		if ( next_line != None) :
			if  ( int(newline["start"]) - int(prior_line["end"]) > 1 and prior_line["chromosome"] == newline["chromosome"] ) :
				distanceToPriorSegment = int(newline["start"]) - int(prior_line["end"])
				if(distanceToPriorSegment  <= maxLen):
					prior_line, newline = merge_gap( prior_line, newline )
				else:

					out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
					prior_line=newline
					newline=next_line
					continue
			elif  ( int(next_line["start"])-int(newline["end"]) > 1 and next_line["chromosome"] == newline["chromosome"] ) :
				# was:      if(int(newline["start"]) - int(next_line["end"]) <= maxLen):
				# but must be as follows.
				# otherwise, the difference will be always negative and the if-condition would be ALWAYS fulfilled. thus,this cannot be the intention...
				if(int(next_line["start"]) - int(newline["end"]) <= maxLen):
					newline, next_line = merge_gap( newline, next_line )



			if ( float(newline["length"]) <= maxLen):
				prior_line, newline = merge_lines_CN(prior_line, newline, next_line)
			else:
				out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
				prior_line=newline
				newline=next_line

	if ( newline != None and prior_line["chromosome"] != newline["chromosome"] and float(newline["length"]) <= maxLen ):
			prior_line = merge_gap(prior_line, newline)
			newline=None

	#print last line(s)
	out.write( "\t".join( str(prior_line[key]) for key in infile.header ) + "\n" ) 
	if newline:
		out.write( "\t".join( str(newline[key]) for key in infile.header ) + "\n" ) 
