#!/bin/bash
set -ue -o pipefail

nr_cores=4
kmer=100
reference_genome="${reference_path}/GRCh38_decoy_ebv_phiX_alt_hla_chr.fa"

mkdir -p ref out
# copy the reference but only chr1-22, X and Y
echo "copy reference genome"
awk '/^>chrM/ {exit} /^>/ {print $1}  !/^>/ {print}' ${reference_genome} > ref/reference.fa
echo "gem-indexer"
./gemtools-1.7.1-i3/bin/gem-indexer -m force-slow-algorithm -T ${nr_cores} -c dna -i ref/reference.fa -o out/index
echo "gem-mappability"
./gemtools-1.7.1-i3/bin/gem-mappability -m2 -e2 -T ${nr_cores} -I out/index.gem -l ${kmer} -o out/outfiles
echo "gem-2-wig"
./gemtools-1.7.1-i3/bin/gem-2-wig -I out/index.gem -i out/outfiles.mappability -o out/outfiles
echo "wigToBigWig"

./tools/wigToBigWig out/outfiles.wig out/outfiles.sizes out/outfiles.bigWig
echo "bigWigToBedGraph"
./tools/bigWigToBedGraph out/outfiles.bigWig out/outfiles.bedGraph

echo "create FAKE scores for HLA and ALT contigs"                                                                                                                                                                                   
echo "Remove primary chrs and phix contigs"                                                                                                                                                                                         
grep ">" GRCh38_decoy_ebv_phiX_alt_hla_chr.fa > contigs                                                                                                                                                                             
cat contigs | grep -v HLA | cut -f1,7 -d " " | sed -e 's/ /\t1\t/g' -e 's/LN://' -e 's/>//' -e 's/$/\t1.0/' > contigs_fake_mappability.bedGraph                                                                                     
cat contigs | grep HLA | cut -f1,2 -d " " | sed -e 's/\t/ /' | cut -f1,3 -d " " | sed -e 's/ /\t1\t/g' -e 's/>//' -e 's/$/\t1.0/' >> contigs_fake_mappability.bedGraph

echo "Filter lines and compress"
(awk '$4 > 0.0 {print $0}' out/outfiles.bedGraph ; cat contigs_fake_mappability.bedGraph ) | ./tools/bgzip > GRCh38_Mappability_Align_100mer_m2e2_ALT_HLA.bedGraph.gz

echo "Create Index with Tabix"
./tools/tabix -p bed GRCh38_Mappability_Align_100mer_m2e2_ALT_HLA.bedGraph.gz

echo "cleaning"
rm -r ref out
