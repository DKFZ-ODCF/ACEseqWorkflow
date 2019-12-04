 #!/usr/bin/python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--in_file', help='Input .vcf file')
parser.add_argument('--out_file', help='Out .vcf file')
args = parser.parse_args()

vcf_infile = open( args.in_file, "r" )
outfile    = open( args.out_file, "w" )

for vcf_line in vcf_infile:
    
    if vcf_line[0] == "#":
        if "#CHROM" in vcf_line:
            vcf_line = vcf_line.rstrip().split("\t")
            vcf_line.append('sample0')
            vcf_line = "\t".join( vcf_line )+"\n"

    else:
        vcf_line = vcf_line.rstrip().split("\t")
        vcf_line.append(vcf_line[9])
        vcf_line = "\t".join( vcf_line )+"\n"

    outfile.write( vcf_line )
