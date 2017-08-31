#script to download dbSNP Variants and split them according to variant type

# CITATION
 
#As a NCBI Resource: "Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11."
#As a whole for a specific build (use this!) : "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. (dbSNP Build ID: 141 ). Available from: http://www.ncbi.nlm.nih.gov/SNP/"
#A single or a range of Submitted SNP (ss) or Reference SNP (rs) entries: "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. dbSNP accession:{ss1 or ss1 – ss100}, (dbSNP Build ID: 141). Available from: http://www.ncbi.nlm.nih.gov/SNP/"

# DOWNLOAD

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00--README.txt 
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi 

# POST PROCESSING

###UNFILTERED
# grep header
# sort and compress using bgzip
# created tabix index

(zcat 00-All.vcf.gz | head -n 1000 | grep ^"#" > 00-All.SNV.vcf; zcat 00-All.vcf.gz | grep "VC=SNV") | #sort -k1,1 -k2,2n -T ./) |
	bgzip > 00-All.SNV.vcf.gz && tabix -p vcf 00-All.SNV.vcf.gz
