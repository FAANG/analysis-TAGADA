#!/bin/bash
# /work2/project/fragencode/tools/peaks_to_distrib.sh

# usage
# peaks_to_distrib.sh peaks.bed annot.tsv

# example
# cd /work2/project/fragencode/workspace/sdjebali/atacseq/pig/liver/ATAC14_ACCACT_L002.q10.rmdup.rmmt.macs_peaks.narrowPeak
# peaks=~kmunyard/fragencode/results/atacseq/sus_scrofa/liver/peaks/ATAC14_ACCACT_L002.q10.rmdup.rmmt.macs_peaks.narrowPeak
# annot=/work2/project/fragencode/workspace/sdjebali/atacseq/pig/annotelt_file.tsv
# /work2/project/fragencode/tools/peaks_to_distrib.sh $peaks $annot

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: peaks_to_distrib.sh peaks.bed annot.tsv > dist_stats.tsv" >&2
    echo "where:" >&2
    echo "- peaks.bed is the list of peaks in bed format" >&2
    echo "- annot.tsv is a 2 column tsv file with no header with the name of the annotated element and the file containing these elements" >&2
    echo "  this should include at least the following categories: prom1Kb, prom5Kb, five_prime_utr, three_prime_utr, gene and intron" >&2
    echo "  for example for pig: /work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/gx.domains.tsv" >&2
    echo "- dist_stats.tsv is a tsv file containing the distribution of peaks into the annotation" >&2
    echo "BE CAREFUL: do not run several times in the same directory since it uses fixed names for outputs" >&2
    echo "" >&2
    exit 1
fi

# To improve
# do not use fixed names for output

peaks=$1
meta=$2

# Intersect the peaks with each kind of element from the annotation
###################################################################
# assuming the following categories
###################################
# cds	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.cds.positions.bed
# exon	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.exon.positions.bed
# five_prime_utr	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.five_prime_utr.positions.bed
# gene	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.gene.positions.bed
# intron	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.intron.positions.bed
# prom1Kb	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.prom1Kb.positions.bed
# prom5Kb	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.prom5Kb.positions.bed
# three_prime_utr	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.three_prime_utr.positions.bed
# transcript	/work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/sus_scrofa.transcript.positions.bed
# like in this file /work/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.84/gx.domains.tsv 
cat $meta | while read elt f
do
intersectBed -a $peaks -b $f -u > peaks_over_$elt.bed
done

# Make the distribution of peaks with 2 different promoter definitions (1Kb and 5KB) and using this order in case of conflict
#############################################################################################################################
# - prom
# - utr 5'
# - utr 3'
# - intron
# - ig
# which means the following rules:
#################################
# - prom = over prom
# - utr 5' = over utr 5' but not in the above categories
# - utr 3' = over utr 3' but not in the above categories
# - intron = over intron but not in the above categories
# - ig = not overlapping genes and not in the above categories
# - rest = within genes but not in the above categories
tot=`wc -l $peaks | awk '{print $1}'`
for d in 1Kb 5Kb
do
prom=`wc -l peaks_over_prom$d.bed | awk '{print $1}'`
utr5=`awk -v fileRef=peaks_over_prom$d.bed 'BEGIN{while (getline < fileRef >0){ko[$0]=1}} ko[$0]!=1' peaks_over_five_prime_utr.bed | wc -l | awk '{print $1}'`
utr3=`awk -v fileRef=peaks_over_prom$d.bed 'BEGIN{while (getline < fileRef >0){ko[$0]=1}} ko[$0]!=1' peaks_over_three_prime_utr.bed | wc -l | awk '{print $1}'`
intron=`awk -v fileRef1=peaks_over_prom$d.bed -v fileRef2=peaks_over_five_prime_utr.bed -v fileRef3=peaks_over_three_prime_utr.bed 'BEGIN{while (getline < fileRef1 >0){ko[$0]=1}while (getline < fileRef2 >0){ko[$0]=1}while (getline < fileRef3 >0){ko[$0]=1}} ko[$0]!=1' peaks_over_intron.bed | wc -l | awk '{print $1}'`
ig=`awk -v fileRef1=peaks_over_prom$d.bed -v fileRef2=peaks_over_five_prime_utr.bed -v fileRef3=peaks_over_three_prime_utr.bed -v fileRef4=peaks_over_gene.bed 'BEGIN{while (getline < fileRef1 >0){ko[$1"_"$2"_"$3]=1}while (getline < fileRef2 >0){ko[$1"_"$2"_"$3]=1}while (getline < fileRef3 >0){ko[$1"_"$2"_"$3]=1} while (getline < fileRef4 >0){ko[$1"_"$2"_"$3]=1}} ko[$1"_"$2"_"$3]!=1' $peaks | wc -l | awk '{print $1}'`
rest=$((tot-prom-utr5-utr3-intron-ig))
echo $prom $utr5 $utr3 $intron $ig $rest | awk -v d=$d -v tot=$tot 'BEGIN{OFS="\t"}{print d, $1, $1/tot*100, $2, $2/tot*100, $3, $3/tot*100, $4, $4/tot*100, $5, $5/tot*100, $6, $6/tot*100}'
done | awk 'BEGIN{OFS="\t"; print "dist", "prom_nb", "prom_pcent", "utr5_nb", "utr5_pcent", "utr3_nb", "utr3_pcent", "intron_nb", "intron_pcent", "ig_nb", "ig_pcent", "rest_nb", "rest_pcent"}{print}' 
# dist	prom_nb	prom_pcent	utr5_nb	utr5_pcent	utr3_nb	utr3_pcent	intron_nb	intron_pcent	ig_nb	ig_pcent	rest_n
# b	rest_pcent
# 1Kb	6228	41.0818	102	0.672823	44	0.290237	1934	12.7573	6804	44.8813	48	0.316623
# 5Kb	6896	45.4881	93	0.613456	37	0.244063	1865	12.3021	6226	41.0686	43	0.283641
