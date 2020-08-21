#!/bin/bash

# analyse_transcript_models.sh
# script to make several analyses on a set of transcripts models predicted (usually from rnaseq)
# and produce several files that can then be used by the script compare_multiple_trsets_wrt_reftrset.sh
# to generate a comparison of several sets of predicted transcript models with respect to a given reference
# (or set of transcripts that we expect to find (usually in the rnaseq sample in question))
# this script is expected to use a lot of resources such as many tens of G of ram and also time
# so it is better to launch it on a cluster

# on Dec 15th 2016 changed the _ delimitor by _ in order to accept ncbi annotation
# on March 27th 2020 changed the file _nbex.tsv into _nbex_intermclass.tsv due to change in refine_comptr script

# takes as input
################
# - a gtf or gff version 2 file containing the transcript models to analyse and including at least exon rows
#   and with gene_id followed by transcript_id as the first two (key,value) pairs
# - a gtf or gff version 2 file containing a set of reference transcripts one expect to find in the previous set
#   (typically annotated transcripts that are thought to be expressed in the rnaseq sample under study) 
# produces as output
####################
# - a set of intermediate files to be used by another script that will gather info from many predictions
#   and then make plots and tables to be used in a presentation (odp, ppt or latex)

# example
#########
# cd ~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie/WithoutAnnot/Test
# annot=~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Annot/gencode.v19.annotation.gtf
# time analyse_transcript_models.sh K562_polya+_wcell_biorep1_star2_stringtie1.0.4.gtf $annot 2> analyse_transcript_models.err
# real    16m0.409s

# Check inputs are provided
###########################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: analyse_transcript_models.sh predtr.gff annot.gtf >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a gtf or gff version 2 file containing the transcript models to analyse and including at least exon rows" >&2
    echo "  and with gene_id followed by transcript_id as the first two (key,value) pairs" >&2
    echo "- a gtf or gff version 2 file containing a set of reference transcripts one expect to find in the previous set" >&2
    echo "  (typically annotated transcripts that are thought to be expressed in the rnaseq sample under study) " >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- a set of intermediate files to be used by another script that will gather info from many predictions" >&2
    echo "  and then make plots and tables to be used in a presentation (odp, ppt or latex)" >&2
    echo "" >&2
    echo "Note: since this script needs many tens of Gigabytes of ram and also some time (not parallelized though), it is better to use a cluster to run it" >&2
    echo "Note: Needs bedtools, comptr and overlap to be installed" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
pred=$1
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
annot=$2
annotbasetmp=`basename ${annot%.gtf}`
annotbase=${annotbasetmp%.gff}
annotdir=`dirname $annot`
anntss=$annotdir/$annotbase\_capped_sites.gff
annrtss=$annotdir/$annotbase\_capped_sites_nr.gff
anntts=$annotdir/$annotbase\_tts_sites.gff
annrtts=$annotdir/$annotbase\_tts_sites_nr.gff

# Programs
##########
MAKETSS=$rootDir/make_TSS_file_from_annotation_simple.sh
MAKETTS=$rootDir/make_TTS_file_from_annotation_simple.sh
MAKESUMMARY=$rootDir/make_summary_stat_from_annot.sh
REFINECOMPTR=$rootDir/refine_comptr_output.sh
EXTRACT5P=$rootDir/extract_most_5p.awk 
EXTRACT3P=$rootDir/extract_most_3p.awk
GFF2GFF=$rootDir/gff2gff.awk


# Start of the script
#####################

# Make the TSS, nrTSS, TTS, nrTTS files for the annotation where the annotation is if they are not already present there
########################################################################################################################
echo "Making TSS from annotated transcripts if not already present where the annotation is" >&2
if [ ! -f $anntss ] || [ ! -f $annrtss ] 
then  
$MAKETSS $annot
mv $annotbase\_capped_sites.gff $annotbase\_capped_sites_nr.gff $annotdir
fi
echo "done" >&2
echo "Making TTS from annotated transcripts if not already present where the annotation is" >&2
if [ ! -f $anntts ] || [ ! -f $annrtts ]
then 
$MAKETTS $annot
mv $annotbase\_tts_sites.gff $annotbase\_tts_sites_nr.gff $annotdir
fi
echo "done" >&2

# Make a complete gff file for the predicted transcripts and compute some basic statistics
##########################################################################################
echo "I am making a complete gff file out of the predicted transcript file and computing some basic statistics" >&2
$MAKESUMMARY $pred > make_summary_stat_from_annot_$base.out 
echo "done" >&2
# nbex nbdistinctex nbtr nbgn nbintrons nbdistinctintrons
# 639561 335376 163349 80844 476212 159060 *** real    7m13.891s

# Refine the output of comptr in order to have more precise transcript classes with respect to the annotation
#############################################################################################################
echo "I am refining the output of comptr in order to have more precise transcript classes with respect to the annotation" >&2
$REFINECOMPTR $base\_complete.gff $annot
echo "done" >&2
# trid  comptrclass     annottrlist     refinedtrclass  nbex
# tSTRG.80746.1   Monoexonic      .       intergenic      n1
# 163350 (5 fields)  *** some minutes


# Make small intermediate 1 column files to be used by plotting script in subsequent action (length of objects from predictions)
#################################################################################################################################
echo "I am making small intermediate 1 column files to be used by plotting script for length of objects from predictions" >&2
# exon length
awk '$3=="exon"{print $5-$4+1}' $base\_complete.gff > exon_length_$base.txt
# distinct exon length
awk '$3=="exon"{print $1":"$4":"$5":"$7}' $base\_complete.gff | sort | uniq | awk '{split($1,a,":"); print a[3]-a[2]+1}' > distinct_exon_length_$base.txt
# spliced and stranded transcript 5' exon length 
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} (($3=="exon")&&(ko[$12]!=1))' $base\_complete.gff | awk -v fldno=12 -f $EXTRACT5P > tr_5p_exon_$base.gff 
awk '{print $5-$4+1}' tr_5p_exon_$base.gff > tr_5p_exon_$base\_length.txt
# spliced and stranded tr gene 5' exon length 
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} (($3=="exon")&&(ko[$12]!=1))' $base\_complete.gff | awk -v fldno=10 -f $EXTRACT5P  > gn_5p_exon_$base.gff 
awk '{print $5-$4+1}' gn_5p_exon_$base.gff > gn_5p_exon_$base\_length.txt
# spliced and stranded transcript 3' exon length
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} (($3=="exon")&&(ko[$12]!=1))' $base\_complete.gff | awk -v fldno=12 -f $EXTRACT3P > tr_3p_exon_$base.gff 
awk '{print $5-$4+1}' tr_3p_exon_$base.gff > tr_3p_exon_$base\_length.txt
# spliced and stranded transcript 3' exon length
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} (($3=="exon")&&(ko[$12]!=1))' $base\_complete.gff | awk -v fldno=10 -f $EXTRACT3P > gn_3p_exon_$base.gff 
awk '{print $5-$4+1}' gn_3p_exon_$base.gff > gn_3p_exon_$base\_length.txt
# spliced and stranded tr internal exon length
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} ($3=="exon")&&(ko[$12]!=1)&&(($7=="+")||($7=="-"))' $base\_complete.gff | awk -v fileRef1=tr_5p_exon_$base.gff -v fileRef2=tr_3p_exon_$base.gff 'BEGIN{while(getline < fileRef1 >0){ko[$1":"$4":"$5":"$7":"$12]=1} while(getline < fileRef2 >0){ko[$1":"$4":"$5":"$7":"$12]=1}} ko[$1":"$4":"$5":"$7":"$12]!=1{print}' | awk -f $GFF2GFF > internal_exons_$base.gff
awk '{print $5-$4+1}' internal_exons_$base.gff > internal_exons_$base\_length.txt
# spliced and stranded transcript 3' exon length
awk '{print $1":"$4":"$5":"$7}' internal_exons_$base.gff | sort | uniq | awk '{split($1,a,":"); print a[3]-a[2]+1}' > distinct_internal_exon_length_$base.txt
# spliced and stranded transcript 3' exon length
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ok["\""substr($1,2)"\";"]=1}}} ($3=="exon")&&(ok[$12]==1){print}' $base\_complete.gff | awk -f $GFF2GFF > monoextr_exons_$base.gff
awk '{print $5-$4+1}' monoextr_exons_$base.gff > monoextr_exons_$base\_length.txt
# spliced and stranded transcript 3' exon length
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} (($3=="exon")&&(ko[$12]!=1))' $base\_complete.gff | sort -k12,12 -k4,4n -k5,5n | awk '{nbex[$12]++; ex[$12,nbex[$12]]=$4"_"$5;} END{for(t in nbex){split(ex[t,1],a,"_"); split(ex[t,nbex[t]],b,"_"); print t, b[2]-a[1];}}' > tr_lgwithexintr_$base.txt
# spliced tr cDNA length (only exons) - here tr can be unstranded
awk -v fileRef=$base\_complete_comp_refinedclass_nbex_intermclass.tsv 'BEGIN{while(getline < fileRef >0){if($2=="Monoexonic"){ko["\""substr($1,2)"\";"]=1}}} ($3=="exon")&&(ko[$12]!=1){lg[$12]+=($5-$4+1)}END{for(t in lg){print t, lg[t]}}' $base\_complete.gff > tr_cumulexlg_$base.txt
echo "done" >&2

# Make TSS and TTS from the predicted transcripts
#################################################
echo "I am making TSS from the predicted transcripts" >&2
$MAKETSS $pred
echo "done" >&2
echo "I am making TTS from the predicted transcripts" >&2
$MAKETTS $pred
echo "done" >&2

# Add to the exact predicted transcripts the following information
##################################################################
# - the coord of the tss of the predicted tr
# - the list of coord of tss of the corresponding annotated tr
# - the coord of the tts of the precicted tr
# - the list of coord of tts of the corresponding annotated tr
echo "To each exact predicted transcript I am adding the coord of the tss/tts of the predicted tr as well as the list of coord of tss/tts of the corresponding annotated tr" >&2
awk -v fileRef1=$base\_capped_sites.gff -v fileRef2=$anntss -v fileRef3=$base\_tts_sites.gff -v fileRef4=$anntts 'BEGIN{while(getline < fileRef1 >0){split($12,a,"\""); tss1[a[2]]=$1"_"$4"_"$7} while(getline < fileRef2 >0){split($12,a,"\""); tss2[a[2]]=$1"_"$4"_"$7} while(getline < fileRef3 >0){split($12,a,"\""); tts1[a[2]]=$1"_"$4"_"$7} while(getline < fileRef4 >0){split($12,a,"\""); tts2[a[2]]=$1"_"$4"_"$7}} {if($2=="Exact"){s1=""; s2=""; t=substr($1,2); gsub(/\"/,"",$3); gsub(/\;/,"",$3); split($3,b,","); k=1; while(b[k]!=""){s1=(s1)(tss2[b[k]])(","); s2=(s2)(tts2[b[k]])(","); k++} print t, $3, tss1[t], s1, tts1[t], s2}}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv > $base\_predtr_spliced_stranded_exact_$annotbase\_correstrlist_TSS_AnnotTSSlist_TTS_AnnotTTSlist.txt
echo "done" >&2

# Add the smaller distance to a matching annotated tr for tss and for tts
##########################################################################
echo "I am adding the smaller distance to a matching annotated tr for tss and for tts" >&2
awk '{dtss=1000000000; dtts=1000000000; split($4,a,","); k=1; while(a[k]!=""){split($3,b,"_"); split(a[k],c,"_"); if(abs(b[2]-c[2])<dtss){dtss=abs(b[2]-c[2]);} k++} split($6,a,","); k=1; while(a[k]!=""){split($5,b,"_"); split(a[k],c,"_"); if(abs(b[2]-c[2])<dtts){dtts=abs(b[2]-c[2]);} k++} print $0, dtss, dtts} function abs(x){return x >= 0 ? x : -x;}' $base\_predtr_spliced_stranded_exact_$annotbase\_correstrlist_TSS_AnnotTSSlist_TTS_AnnotTTSlist.txt > $base\_predtr_spliced_stranded_exact_$annotbase\_correstrlist_TSS_AnnotTSSlist_TTS_AnnotTTSlist_smallerdist_TSS_TTS.txt
echo "done" >&2

# To assess sn and precision when detecting annotated tss and tts, the annotated tss and tts (nr per gene) 
##########################################################################################################
# and the predicted tss (nr per gene) are compared using bedtools windowBed function
####################################################################################
echo "I am comparing the annotated tss and tts (nr per gene) and the predicted tss (nr per gene) using bedtools windowBed function" >&2
for dist in 50 100 500
do
# for tss
windowBed -a $annrtss -b $base\_capped_sites_nr.gff -w $dist -sm -u > $annotbase\_capped_sites_nr_with_pred_less$dist\bp.gff
windowBed -a $base\_capped_sites_nr.gff -b $annrtss -w $dist -sm -u > $base\_capped_sites_nr_with_annottss_less$dist\bp.gff
# for tts
windowBed -a $annrtts -b $base\_tts_sites_nr.gff -w $dist -sm -u > $annotbase\_tts_sites_nr_with_pred_less$dist\bp.gff
windowBed -a $base\_tts_sites_nr.gff -b $annrtts -w $dist -sm -u > $base\_tts_sites_nr_with_annottts_less$dist\bp.gff
done
echo "done" >&2


# Clean
#######
echo "I am cleaning" >&2
rm $base\_complete.gff 
rm tr_5p_exon_$base.gff tr_3p_exon_$base.gff 
rm gn_5p_exon_$base.gff gn_3p_exon_$base.gff 
rm internal_exons_$base.gff monoextr_exons_$base.gff
echo "done" >&2
