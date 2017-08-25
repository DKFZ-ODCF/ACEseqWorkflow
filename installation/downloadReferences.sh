
REFERENCE_GENOME=false
dbSNP_FILE=false
MAPPABILITY_FILE=false
CHROMOSOME_LENGTH_FILE=false
IMPUTE_FILES=false

if [[ $REFERENCE_GENOME != "true" ]] 
then
        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
fi

if [[  $dbSNP_FILE != "true" ]]
then
       bash prepare_dbSNPFile.sh
fi

if [[  $MAPPABILITY_FILE != "true" ]]
then
       wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bw.gz
fi

if [[  $CHROMOSOME_LENGTH_FILE != "true" ]] 
then
       wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz  | zcat | grep -Pv "(_)|(chrM)" | sed -e '1i\#chrom\tsize\tfileName' >chrlengths.txt
fi

##provide online or load into docker?
#       <cvalue name="REPLICATION_TIME_FILE" value="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/ENCODE/ReplicationTime_10cellines_mean_10KB.Rda" type="path"/>
#       <cvalue name="GC_CONTENT_FILE" value="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_GRch37_100genomes_gc_content_10kb.txt" type="path"/>
if [[  $IMPUTE_FILES != "true" ]]
then
	wget https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
	tar -xzvf ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
fi

