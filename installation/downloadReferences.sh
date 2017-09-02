
REFERENCE_GENOME=true
dbSNP_FILE=true
MAPPABILITY_FILE=true
CHROMOSOME_LENGTH_FILE=true
statFiles=true
IMPUTE_FILES=false

SCRIPT=$(readlink -f $0)
scriptdir=`dirname $SCRIPT`

echo $scriptdir

mkdir -p sequence/1KGRef
mkdir databases
mkdir stats
mkdir databases/ENCODE
mkdir -p databases/1000genomes/IMPUTE
mkdir databases/UCSC
mkdir -p databases/dbSNP/dbSNP_135

if [[ $REFERENCE_GENOME != "true" ]] 
then
	echo downloading reference genome///
        wget -P sequence/1KGRef ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gunzip sequence/1KGRef/hs37d5.fa.gz
fi

if [[  $dbSNP_FILE != "true" ]]
then
	echo dbSNP file....
       installationpath=`pwd`
       cd databases/dbSNP/dbSNP_135
       bash $scriptdir/prepare_dbSNPFile.sh
	cd $installationpath
fi

if [[  $MAPPABILITY_FILE != "true" ]]
then
       echo downloading mappability file....
       wget -P databases/UCSC  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
       echo "please convert the bigwig file to bedgraph"
fi

if [[  $CHROMOSOME_LENGTH_FILE != "true" ]] 
then
	echo downloading dbSNP file....
       wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz  | zcat | grep -Pv "(_)|(chrM)" | sed -e '1i\#chrom\tsize\tfileName' >stats/chrlengths.txt
fi

if [[ $statFiles == "true" ]]
then
	echo downloading statsfile
       cp $scriptdir/hg19_GRch37_100genomes_gc_content_10kb.txt stats/
       cp $scriptdir/ReplicationTime_10cellines_mean_10KB.Rda ENCODE/
fi

if [[  $IMPUTE_FILES != "true" ]]
then
	echo downloading impute files....
	wget -P databases/1000genomes/IMPUTE https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
	tar -xzvf databases/1000genomes/IMPUTE/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz -C databases/1000genomes/IMPUTE
	wget -P databases/1000genomes/IMPUTE https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
	tar -xzvf databases/1000genomes/IMPUTE/ALL_1000G_phase1integrated_v3_impute.tgz -C databases/1000genomes/IMPUTE
fi

