#!/usr/bin/python

import sys
import numpy
import matplotlib as mpl
mpl.use( 'Agg' )
import matplotlib.pyplot as plt

points_x = []
points_y = []
points_c = []

for line in sys.stdin:
	
	if line[0] == '#' or len(line.split())<9:
		continue
	
	info_field = [ entry.split("=") for entry in line.split()[7].split(";") ]
	
	reads_entry = None
	
	for entry in info_field:
		if entry[0] == "DP4":
			reads_entry =  [ float(val) for val in entry[1].split(",") ]
	
	position = float( line.split("\t")[1] )
	baf_value = sum( reads_entry[2:4] ) / sum( reads_entry[0:4] )

	points_x.append( position )
	points_y.append( baf_value )
		
	if line.split()[9].split(":")[0] == "0|1":
		points_c.append("r")
	elif line.split()[9].split(":")[0] == "1|0":
		points_c.append("b")
	else:
		points_c.append("k")

plt.figure(figsize=(20,10))

plt.scatter( points_x, points_y, c=points_c, s=(4), linewidths=(0) )

plt.savefig( sys.stdout, format="png", dpi=300 )
