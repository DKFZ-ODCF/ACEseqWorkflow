#!/usr/bin/python

import sys
import argparse
import time


parser=argparse.ArgumentParser(description='Convert tab seperated CNV calls to vcf file' )
parser.add_argument('--file',	  	'-f', type=file, help='ACEseq most important info output file')
parser.add_argument('--makehead', 	'-m', type=int, default=1, help='Set to 1 if you want to produce pancancer conform head')
parser.add_argument('--id', 	  	'-i', type=str, default="NA", help='Patient identifier')
parser.add_argument('--out', 	  	'-o', type=str, help='Output file')
parser.add_argument('--center',   	'-c', type=str, default="DKFZ", help='Center of calling pipeline')
parser.add_argument('--url',   	  	'-u', type=str, default="NA", help='The URL for the variant calling workflow')
parser.add_argument('--projectCode',	'-p', type=str, default="NA", help='project code for this donor')
parser.add_argument('--additionalTag',	'-a', type=str, help='additional header strings.')


parser.add_argument
args = parser.parse_args()


makehead = args.makehead
date	 = time.strftime("%Y%m%d")
print_info = 0
refgenome = ["hs37d5","ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"]


#HEADER

def vcf_header(file_hd, sample):
    #GENERAL INFORMATION
    file_hd.write("##fileformat=VCFv4.1\n##CreatedBy=TabToVCF-v1.1\n")
    file_hd.write("##fileDate=%s\n"% date)
    file_hd.write("##pancancerversion=1.0\n")
    file_hd.write("##reference=<ID=%s,Source=%s\n"% ( refgenome[0], refgenome[1] ) )
    file_hd.write("##center=" + args.center + "\n")
    file_hd.write("##filedata=" + sample['ID'] + "\n")
    file_hd.write("##workflowName=DKFZ_CNV_workflow\n")
    file_hd.write("##workflowVersion=1.0.0\n")
    file_hd.write("##workflowURL=" + args.url + "\n")
    file_hd.write("##donor=submitter_donor_id:" + args.id + "\n")
    file_hd.write("##projectCode=dcc_project_code:" + args.projectCode + "\n")
    if (args.additionalTag):
	if (not args.additionalTag.endswith("\n")):
		args.additionalTag = args.additionalTag + "\n"
	file_hd.write(args.additionalTag)
    file_hd.write("##Purity=" + sample['purity'] + "\n")
    file_hd.write("##Ploidy=" + sample['ploidy'] + "\n")
    file_hd.write("##Full ploidy (CN of most segments)=" + sample['fullPloidy'] + "\n")
    file_hd.write("##Quality =%.2f"% float(sample['quality']) + "\n")
    file_hd.write("##Assumed sex=" + sample['assumed sex'] + "\n")
    #INFO FIELDS
    file_hd.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'+"\n")
    file_hd.write('##INFO=<ID=STARTINT,Number=2,Type=String,Description="Intervall for start of segment">' + "\n")
    file_hd.write('##INFO=<ID=ENDINT,Number=2,Type=String,Description="Intervall for end of segment">' + "\n")
#   file_hd.write('##INFO=<ID=CNVTYPE,Number=1,Type=Flag,Description="CNV type">'+"\n")
    file_hd.write('##INFO=<ID=CNVLEN,Number=1,Type=Integer,Description="Length of copy number altered region">'+"\n")
    file_hd.write('##INFO=<ID=TCNexact,Number=1,Type=Float,Description="Precise total copy number of region">' + "\n")
    file_hd.write('##INFO=<ID=CN1,Number=1,Type=Float,Description="Precise copy number allele 1">' + "\n")
    file_hd.write('##INFO=<ID=CN2,Number=1,Type=Float,Description="Precise copy number allele 2">' + "\n")
    file_hd.write('##INFO=<ID=sub,Number=0,Type=Flag,Description="Allele specific copy number indicates subpopulation with aberration in this region">'+"\n")
    #ALT FIELD
    file_hd.write('##ALT=<ID=NEUTRAL,Description="No CN change with regards to diploid genome">' + "\n" )
    file_hd.write('##ALT=<ID=DEL,Description="Deletion with regards to diploid genome">' + "\n" )
    file_hd.write('##ALT=<ID=AMP,Description="Amplification with regards to diploid genome">'+"\n")
    file_hd.write('##ALT=<ID=cnLOH,Description="copy number neutral Loss of Heterozygosity">' + "\n")
    file_hd.write('##ALT=<ID=LOHgain,Description="Gain with loss of heterozygosity with regards to diploid genome">' + "\n")
    file_hd.write('##ALT=<ID=LOH,Description="Deletion with loss of heterozygosity with regards to diploid genome">' + "\n")
    #FORMAT FIELD
    file_hd.write('##FORMAT=<ID=TCN,Number=1,Type=Integer,Description="Total copy number">' + "\n")
    file_hd.write('##FORMAT=<ID=ACN,Number=1,Type=String,Description="Allele specific copy number">' + "\n")
    #SAMPLE FIELD
    file_hd.write('##SAMPLE=<ID=TUMOR,SampleName=tumor_%s,Individual=%s,Description=\"Tumor\">'% (args.id,args.id) + "\n")
    #COLUMN NAMES
    file_hd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "TUMOR" + "\n")

def parse_header(infile):
	"Parse information out of header and get column names"
	sample = {}
	for line in infile:
		if line.startswith('#'):
			line = line.lstrip('#').rstrip('\n')
			line = line.split(':')
			sample[line[0]]=line[1]
		else:
			colnames = line.rstrip('\n').split('\t')
			break
	return(sample, colnames)
	


def vcf_line(file_vcf, line, alias, fullPloidy, sex, refField='N', qualField='.', filterField='PASS', formatField='', idField='.'):
	"Rewrite information to vcf format"
	start = str(int( round(float(line['start']) ) ) ) #should be the higher number
	end = int( round(float(line['end']) - 0.5 ) ) 	  #should be the lower number
	TCN = str( int( round (float(line['TCN']) ) ) )
	intervalls = [ getCoord(start, line['minStart']), getCoord(start, line['maxStart']), getCoord(end, line['minEnd']), getCoord(end, line['maxEnd']) ]
	info = "END=%i;"% (end)
	info = info + "STARTINT=[%s,%s];ENDINT=[%s,%s];"% (intervalls[0], intervalls[1], intervalls[2], intervalls[3])
	info = info + "CNVLEN=%s;TCNexact=%.2f;"% ( int(round(float(line['length']))), float(line['TCN']) )
	GNL = line['GNL']
	if not line['c1Mean']=='NA' and not line['c1Mean']=='NA':
		c1 = str( int( round (float(line['c1Mean']) ) ) )
		c2 = str( int( round (float(line['c2Mean']) ) ) )
	else:
		c1 = "."
		c2 = "."

	if GNL == 'sub' :
		GNL = checkCopyNumberState(2, TCN)
		info = info + "CN1=%.2f;CN2=%.2f;sub;"% ( float(line['c1Mean']), float(line['c2Mean']) )
		c1 = '.'
		c2 = '.'
	#copy number state defined 
	elif GNL in alias:
		GNL = alias[GNL]
		if not c1==".":
			info = info + "CN1=%.2f;CN2=%.2f;"% ( float(line['c1Mean']), float(line['c2Mean']) )
			if  round(round(float(line['c1Mean'])) + round(float(line['c2Mean'])) ) != int(TCN):
				info = info + "sub;"
				c1 = '.'
				c2 = '.'
		#Ploiy not diploid but still defined as neutral == > rewrite to loss or gain
		if ( int(TCN) != 2 ) and GNL == "NEUTRAL":
			GNL = checkCopyNumberState( 2, TCN)
			if ( int(TCN) == 1 and sex=='male' and (line['chromosome']=="23" or line['chromosome'] == "24") ):
				GNL = checkCopyNumberState( 1, TCN)
				
	#imbalance in diploid genome
	elif GNL == 'qmixDH' and int(fullPloidy) == 2:
		if not c1 == ".":
			info = info + "CN1=%.2f;CN2=%.2f;sub;"% ( float(line['c1Mean']), float(line['c2Mean']) )
		GNL = 'NEUTRAL'
	#other CN state such as NA due to to few SNPs in segment or small segments
	elif not GNL in alias:
		GNL = checkCopyNumberState(2, TCN)
		c1 = "."
		c2 = "."

 	form = "TCN:ACN"
	sample = "%s:%s/%s"% ( TCN, c1, c2 )
	file_vcf.write( line['chromosome']+'\t' + start + '\t' + '.' + '\t' + refField + '\t<' + GNL + '>\t' + qualField + '\t' +filterField + '\t' + info + '\t' + form + "\t" + sample + "\n")
	
def getCoord(coord, intervallCoord):
	"If no max or min is set set it to actual start or end coordinate"
	if ("Inf" in intervallCoord):
		return(coord)
	else:
		return(intervallCoord)

def checkCopyNumberState(fullPloidy, TCN):
	'''Check whether copy number indicates gain, loss or neutral copy number state'''
	fullPloidy = int(fullPloidy)
	TCN = int(TCN)
	if (TCN > fullPloidy):
		return "AMP"
	elif TCN < fullPloidy :
		return "DEL"
	else:
		return "NEUTRAL"


if __name__=='__main__':

	if not args.id or not args.file:
		sys.stderr.write("ERROR: Please specify patient id and input file. For more information, use -h.\n\n\n")
		sys.exit(2)

	if not args.out:
		out=sys.stdout
	else:
		try:
			out=open(args.out, 'w')
		except IOError as (errno, strerr ):
			sys.stderr.write("WARNING: Specified outputfile cannot be written. Please check given path.\n")
			sys.exit("IOError %i:%s\n" % (errno, strerr))
	

	try:
		tabfile = args.file
	except IOError as (errno, strerr ):
		sys.exit("IOError %i:%s\n" % (errno, strerr))

	alias = {"gain":"AMP", "loss":"DEL", "cnLOH":"cnLOH", "LOHgain":"LOHgain", "LOH":"LOH", "neutral":"NEUTRAL", "sub":"sub"}

	colnames=[]
	sample, colnames  = parse_header(tabfile)
	sample['ID'] = args.id
	vcf_header(out,sample)

	for line in tabfile:
		line =line.rstrip("\n")
		line = line.split('\t')
		fields = {}
		for i in range(0, len(colnames)):
			fields[colnames[i]] = line[i]
		vcf_line(out, fields, alias, sample['fullPloidy'], sample['assumed sex'])
