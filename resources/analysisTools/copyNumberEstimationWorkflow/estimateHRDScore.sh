#!/usr/bin/sh

path=/icgc/pcawg/analysis/train2full/projects/
ACEseq=
target=/icgc/pcawg/analysis/BRCAness_biomarker/cnv_new/newgenelists/
genelistFile=/icgc/pcawg/analysis/BRCAness_biomarker/gene_lists/all_BRCAness.bed
blacklists=/icgc/dkfzlsdf/analysis/mmml/genome/cohort_analysis/CNV_analysis_SOPHIA_v32/comb_pro_extra_SOPHIA/blacklistRegions
TOOLSDIR=/home/kleinhei/Project/panCan/pancan.git/CNV_analysis/

date=`date +%y%m%d%H`
		 [[  ! -f $path/$project/$pid/$ACEseq/checkpointVcf.txt ]] &&  continue;

##remove artifact regions
echo $pid

#smooth Data
bedtools-2.24.0 intersect -header -v -f 0.7 -a - -b $blacklists/artifact.homoDel1500.txt $file | \
bedtools-2.24.0 intersect -header -v -f 0.7 -a - -b $blacklists/potentialArtifacts.txt >$combProFile.tmp

python ${TOOLSDIR}/removeBreakpoints.py -f $combProFile.tmp -i 4108101 -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
python ${TOOLSDIR}/mergeArtifacts.py -f $combProFile.tmp -i 4108101 -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
python ${TOOLSDIR}/removeBreakpoints.py -f $combProFile.tmp -i 4108101 -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
#this file could be written out and sorted according to chromosomes
python ${TOOLSDIR}/smoothData.py -f $combProFile.tmp  -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp && \
python ${TOOLSDIR}/removeBreakpoints.py -f $combProFile.tmp -i 4108101 -o $combProFile.tmp.tmp && mv $combProFile.tmp.tmp $combProFile.tmp

patientsex=`cat $SEX_FILE`
fullploidy=`grep "fullPloidy" $file | sed 's/.*://'`
purity=`grep "purity" $file | sed 's/.*://'`

Rscript-3.3.1 /home/kleinhei/Project/panCan/pancan.git/CNV_analysis/analyseCNVs.new.R \
        --patientsex $patientsex \
	--ploidy $fullploidy \
	--purity $purity \
	--id $pid \
	--segmentfile $file.tmp \
	--mergedfile $combProFile.tmp 
