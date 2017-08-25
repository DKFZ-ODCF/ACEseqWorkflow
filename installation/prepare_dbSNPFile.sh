#script to download dbSNP Variants and split them according to variant type

# CITATION
 
#As a NCBI Resource: "Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11."
#As a whole for a specific build (use this!) : "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. (dbSNP Build ID: 141 ). Available from: http://www.ncbi.nlm.nih.gov/SNP/"
#A single or a range of Submitted SNP (ss) or Reference SNP (rs) entries: "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. dbSNP accession:{ss1 or ss1 – ss100}, (dbSNP Build ID: 141). Available from: http://www.ncbi.nlm.nih.gov/SNP/"

# DOWNLOAD

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00--README.txt 
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz.tbi 

# POST PROCESSING

###UNFILTERED
# grep header
# grep SNV/INDEL/MNV with weight 1 (= unique) instead of 2 (non unique) or 3 (many) - as done previously by Barbara Hutter for dbSNP 135 (?) and Naveed Isaque dbSNP 141
# add CLN tag if OM tag is present
# sort and compress using bgzip
# created tabix index

zcat 00-All.vcf.gz | head -n 1000 | grep ^"#" > 00-All.SNV.vcf
zcat 00-All.vcf.gz | grep "VC=SNV" | sed 's/;OM/;OM;CLN/' | sort -k1,1 -k2,2n -T /ibios/co01/kleinhei/ >> 00-All.SNV.vcf
#add new line in header to explain CLN tag
sed -i '54i\##INFO=<ID=CLN,Number=0,Type=Flag,Description="CLN=OM">' 00-All.SNV.vcf

zcat 00-All.vcf.gz | head -n 1000 | grep ^"#" > 00-All.MNV.vcf
zcat 00-All.vcf.gz | grep "VC=MNV" | sed 's/;OM/;OM;CLN/' | sort -k1,1 -k2,2n -T /ibios/co01/kleinhei/ >> 00-All.MNV.vcf
#add new line in header to explain CLN tag
sed -i '54i\##INFO=<ID=CLN,Number=0,Type=Flag,Description="CLN=OM">' 00-All.MNV.vcf

zcat 00-All.vcf.gz | head -n 1000 | grep ^"#" > 00-All.INDEL.vcf
zcat 00-All.vcf.gz | grep "VC=MNV" | sed 's/;OM/;OM;CLN/' | sort -k1,1 -k2,2n -T /ibios/co01/kleinhei/ >> 00-All.INDEL.vcf
#add new line in header to explain CLN tag
sed -i '54i\##INFO=<ID=CLN,Number=0,Type=Flag,Description="CLN=OM">' 00-All.INDEL.vcf


### FILTERED
# grep header
# grep SNV/INDEL/MNV with weight 1 (= unique) instead of 2 (non unique) or 3 (many) - as done previously by Barbara Hutter for dbSNP 135 (?) and Naveed Isaque dbSNP 141
# sort and compress using bgzip
# created tabix index

(zcat 00-All.vcf.gz | head -n 1000 | grep ^"#"; zcat 00-All.filtered.vcf.gz | grep "WGT=1" | grep "VC=SNV" | sort -k1,1 -k2,2n) | bgzip > 00-All.filtered..SNV.vcf.gz  && tabix -p vcf 00-All.filtered.SNV.vcf.gz
(zcat 00-All.vcf.gz | head -n 1000 | grep ^"#"; zcat 00-All.filtered.vcf.gz | grep "WGT=1" | grep "VC=DIV" | sort -k1,1 -k2,2n) | bgzip > 00-All.filtered.INDEL.vcf.gz && tabix -p vcf 00-All.filtered.INDEL.vcf.gz
(zcat 00-All.vcf.gz | head -n 1000 | grep ^"#"; zcat 00-All.filtered.vcf.gz | grep "WGT=1" | grep "VC=MNV" | sort -k1,1 -k2,2n) | bgzip > 00-All.filtered.MNV.vcf.gz   && tabix -p vcf 00-All.filtered.MNV.vcf.gz
