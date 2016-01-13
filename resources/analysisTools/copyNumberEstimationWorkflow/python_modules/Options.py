#!/usr/bin/python

import sys

def parse( mandatory, optional={} ):

	options = {}
	for pairs in zip( sys.argv[1:], sys.argv[2:] ):                             # read command line arguments
		if pairs[0][:2] == '--':
			if pairs[0][2:] in mandatory:
				options[ pairs[0][2:] ] = mandatory[ pairs[0][2:] ]( pairs[1] ) # mandatory[ pairs[0][2:] ] contains the expected data type
			elif pairs[0][2:] in optional:
				options[ pairs[0][2:] ] = optional[ pairs[0][2:] ]( pairs[1] )  # optional[ pairs[0][2:] ] contains the expected data type
			else:
				raise Exception( "ERROR: invalid argument " + pairs[0][2:] )

	missing_args = set(mandatory.keys()) - set(options.keys())

	if missing_args:
		raise Exception( "ERROR: missing arguments "+", ".join( missing_args ) )
	else:
		return options
