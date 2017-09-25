#!/usr/bin/python

import argparse
import sys
import json
from python_modules import Tabfile

parser=argparse.ArgumentParser(description='Convert tab seperated CNV calls to vcf file' )
parser.add_argument('--solutionFile',  	'-s', type=str, help='File with ACEseq solutions')
parser.add_argument('--pid', 	  	'-i', type=str, help='Patient identifier printed in json')
parser.add_argument('--path', 	  	'-p', type=str, help='Path to aceseq directory')
parser.add_argument('--out', 	  	'-o', type=str, help='Output file')

args = parser.parse_args()

if __name__=='__main__':

	if not args.pid or not args.solutionFile or not args.path:
		sys.stderr.write("ERROR: Please specify patient id, ploidy_purity file and input path. For more information, use -h.\n\n\n")
		sys.exit(2)

	if not args.out:
		out=sys.stdout

	else:
		try:
			out=open(args.out, 'w')
		except IOError as (errno, strerr ):
			sys.stderr.write("WARNING: Specified outputfile cannot be written. Please check given path.\n")
			sys.exit("IOError %i:%s\n" % (errno, strerr))

	#get best solution
	try:	
			pid = args.pid
			path = args.path		

			ppFile = open( args.solutionFile ) 

			distances = []
			ploidies = []
			entries =[]
			for line in ppFile:
				#colnames ppFile: ploidy ploidy_factor purity distance
				if line.startswith("ploidy"):
					continue	
				fields = line.rstrip("\n").split("\t")
				distances.append(float(fields[3]) )
				ploidies.append(float(fields[1]) )
				entries.append(fields)
			contin=1
	except:
			sys.stderr.write( "FILE for %s does not exist\n"% pid)
			sys.exit(2)

	m = min( [ abs( j-2.0 ) for j in ploidies ] )
	index = [i for i,j in enumerate(ploidies) if abs(j-2.0)==m ]

	solutions={1: "%s_%s"% (entries[index[0]][1], entries[index[0]][2] ) }
	count=2
	for i,j  in enumerate(entries):
		if i==index: continue
		solutions[count]="%s_%s"% ( entries[index[0]][1], entries[index[0]][2] )
		count+=1
	
	jsonMain={}
	for key in solutions.keys():	
		try:
			infile="%s/%s_cnv_parameter_%s.txt"% (path, pid, solutions[key])
			tabfile = Tabfile.Input( open(infile) )
		except IOError as (errno, strerr ):
			sys.exit("IOError %i:%s\n" % (errno, strerr))

		#convert simple tab seperated file wit header 
		jsonMain[key]=tabfile.readline()

	out.write( json.dumps(jsonMain, indent=2, separators=(",",":")) )
