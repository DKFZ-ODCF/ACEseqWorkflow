#!/usr/bin/python

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

#######################################
# snp_cnv.py
# v0.0.1
# 2 nov 2011
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# retrieves SNP and CNV information from samtools mpileup output of control AND tumor whole-genome sequencing data.
# to pipe from mpileup, use "-". to only consider specific genomic regions, set "-r" in mpileup.
#
# returns base counts of A and B allele of control and tumor, if there is coverage in either one (A allele is reference base).
# returns coverage in 1kb bins of control and tumor, if there is coverage in either one.
#
# arguments:
#   -Q INT     minimum phred-scaled base quality [13]


def mod_snp_cnv(sysargv):
	import sys
	import pysam
	import gzip	
	import argparse	
	#######################################
	# arguments, filehandles

	#print sysargv
	#print len(sysargv) 
	parser=argparse.ArgumentParser(description = "Parse mpileup output for ACEseq raw data")
	parser.add_argument('--withoutcontrol', '-w', default=False, type = bool, help = "")
	parser.add_argument('--dbsnp',   '-d', required=1, type = str, help = "dbSNP position vcf")
	parser.add_argument('--outsnps', '-s', required=1, help = "Output file for snp coverage")
	parser.add_argument('--outcov',  '-c', required=1, help = 'Output file for coverage window output')
	parser.add_argument('--infile',  '-i', type=argparse.FileType('r'), default='-' ,help="mpileup output, usually given as stdin")
	parser.add_argument('--quality', '-q', default=13 , type=int, help="minimum base quality for base in read to be considered")

	args = parser.parse_args()

	
	if not args.dbsnp or not args.outsnps or not args.outcov:	# at least five arguments
		sys.exit("not all parameters set, please use option -h for mor information")
	
	withoutControl = args.withoutcontrol
	tabixPositions = pysam.Tabixfile( args.dbsnp, "r")
	fin   = args.infile
	fout  = gzip.open( args.outsnps, "wb" )
	fout2 = open( args.outcov, "w" )
	bq = args.quality


	#######################################
	# functions
	
	def mp2bc(mps, mpq, mr, mq):												# count base occurences in mpileup sequence string
		mps = mps.upper()

		mr = mr.upper()
		if mr == "A": mrn = 0
		elif mr == "C": mrn = 1
		elif mr == "G": mrn = 2
		else: mrn = 3
			
		bc = [0, 0, 0, 0]														# [A, C, G, T]
																																									
		i1 = 0																	# index for seq
		i2 = 0																	# index for qual
		
		try:
			while i1 < len(mps):												# parse mpileup ..
				if mps[i1] in [".", ","]:
					if ord(mpq[i2]) >= mq+33:									# only consider bases with baseQ >= mq
						bc[mrn] += 1
					i1 += 1
					i2 += 1
				elif mps[i1] == "A":
					if ord(mpq[i2]) >= mq+33:
						bc[0] += 1
					i1 += 1
					i2 += 1
				elif mps[i1] == "T":
					if ord(mpq[i2]) >= mq+33:
						bc[3] += 1
					i1 += 1
					i2 += 1
				elif mps[i1] == "C":
					if ord(mpq[i2]) >= mq+33:
						bc[1] += 1
					i1 += 1
					i2 += 1
				elif mps[i1] == "G":
					if ord(mpq[i2]) >= mq+33:
						bc[2] += 1
					i1 += 1
					i2 += 1
				elif mps[i1] == "^":											# start of read
					i1 += 2														# .. jump over alignment quality
				elif mps[i1] in ["$", "*"]:										# end of read
					i1 += 1
				elif mps[i1] in ["+", "-"]:										# insertion/deletion
					ni = 1
					nis = ""
					while mps[i1+ni] in [`n` for n in range(10)]:				# .. number of inserted/deleted bases (might be >9)
						nis += mps[i1+ni]
						ni += 1
					i1 += int(nis)+1+len(nis)									# .. baseQ of inserted bases are not in mpileup 
				elif mps[i1] == "N":
					i1 += 1
					i2 += 1
				else:															# should not happen
					sys.exit(mps)
		except IndexError:
			print mps
			print mpq
		return [bc,mrn]																# return base count and index of reference base

	#######################################
	# main

	r = "foo"
	[bl, bsn, bst] = [0, 0, 0]													# [last bin id, bin coverage control, bin coverage tumor]
	fout.write("#chr\tpos\tid\tAn\tBn\tAt\tBt\n")								# snp header
	fout2.write("#chr\tpos\tend\tnormal\ttumor\n")									# cnv header
	
	l = fin.readline()
	while l:
		l = l.rstrip().split("\t")
		
		l[0] = l[0].replace('chr','')

		if l[0] != r:															# new chromosome
			if ( (  ( bsn > 0 or bst > 0 ) and not withoutControl ) or withoutControl ) and r != "foo" :
				 fout2.write("%s\t%i\t%i\t%i\t%i\n" % (r, bl*1000+1, (bl+1)*1000, bsn, bst))
			[bl, bsn, bst] = [0, 0, 0]											# write last cnv bin and reset
				
			r = l[0]
			sys.stderr.write("%s .. " % r)
			dp = dict([(int(t[1]), t[2]) for t in tabixPositions.fetch( region=r, parser=pysam.asTuple())])
			sys.stderr.write("%i SNPs in file\n" % len(dp))						# load positions in dbSNP
		
		l[1] = int(l[1])
		
		if l[1] in dp:
			if not withoutControl:														# if in dbSNP
				[b1, ref] = mp2bc(l[4], l[5], l[2], bq)									# calculate base count
				b2 = mp2bc(l[7], l[8], l[2], bq)[0]
			else:
				b1	 = [ 0 ,0 ,0 ,0 ]
				[b2,ref] = mp2bc(l[4], l[5], l[2], bq)
			
			b_ref = zip( [b1[ref]], [b2[ref]] )
			b1 = [ b1[i] for i in range(len(b1)) if i != ref ]
			b2 = [ b2[i] for i in range(len(b2)) if i != ref ]
			b  = sorted(zip(b1, b2), key=lambda a: a[0])							# determine allele A (higher base count) and B in control, applies to tumor
#			if b[2][0] != 0 and b_ref[0][0] != 0: 		#control has alternative allele count
			fout.write("%s\t%i\t%s\t%i\t%i\t%i\t%i\n" % (r, l[1], dp[l[1]], b_ref[0][0], b[2][0], b_ref[0][1], b[2][1]))
		
		bi = int((l[1]-1)/1000)												# cnv bin id
		if bi == bl:
			if not withoutControl:													# if not new
				bsn += int(l[3])													# .. add to coverage
				bst += int(l[6])
			else:
				bsn = 0
				bst += int(l[3])
		else:																	# else write last bin and make new bin
			if ( ( bsn > 0 or bst > 0 ) and not withoutControl ) or withoutControl :
				 fout2.write("%s\t%i\t%i\t%i\t%i\n" % (r, bl*1000+1, (bl+1)*1000, bsn, bst))
			bl = bi
			if not withoutControl:
				bsn = int(l[3])
				bst = int(l[6])
			else:
				bsn = 0
				bst = int(l[3])
						
		l = fin.readline()
	
	if ( bsn > 0 or bst > 0 and not withoutControl ) or withoutControl :
		 fout2.write("%s\t%i\t%i\t%i\t%i\n" % (r, bl*1000+1, (bl+1)*1000, bsn, bst))
		
	
	fin.close
	fout.close
		

if __name__ == "__main__":
    import sys
    mod_snp_cnv(sys.argv)

