#!/usr/bin/python

from python_modules import Options

options = Options.parse( { "vcf_file"   : str,   "outfile" : str  } )
if options:

	vcf_infile = open( options["vcf_file"], "r" )
	outfile    = open( options["outfile" ], "w" )
	
	for line in vcf_infile:
		if line[0] != "#":
			line = line.rstrip("\n").split("\t")
			if len(line) > 9:
				gt_index = line[8].split(":").index("GT")
				genotype = line[9].split(":")[gt_index].split("/")
				
				line[0] = line[0].replace('chr','')
				
				outfile.write( line[0]+"-"+line[1]+" "+line[1]+" "+
				               line[1]+" "+line[3]+" "+line[4]+" " )

				if genotype[0] == "0" and genotype[1] == "0":
					outfile.write( "1 0 0\n" )
				elif genotype[0] == "0" and genotype[1] == "1":
					outfile.write( "0 1 0\n" )
				elif genotype[0] == "1" and genotype[1] == "1":
					outfile.write( "0 0 1\n" )
				else:
					raise(Exception("Invalid genotype"+repr(genotype)))
