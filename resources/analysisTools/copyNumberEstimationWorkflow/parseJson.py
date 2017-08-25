

import json
import argparse
import sys

parser=argparse.ArgumentParser(description='Convert tab seperated CNV calls to vcf file' )
parser.add_argument('--file', '-f', type=file, help='Tabfile with headerline')

parser.add_argument
args = parser.parse_args()

json_data=args.file.read()

solutions =  json.loads(json_data)

for index in solutions.keys():
	outstring = ""
	for detail in solutions[index].keys():
		outstring = outstring + detail + "=" + solutions[index][detail] + " "
#		sys.stdout.write( " " + detail + "=" + solutions[index][detail] +" ")
#	sys.stdout.write("\n")
	print(outstring)
