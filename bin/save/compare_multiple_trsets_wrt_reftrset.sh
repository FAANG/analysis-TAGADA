#!/bin/bash

# compare_multiple_trsets_wrt_reftrset.sh
# script to gather the analysis results of the script analyse_transcript_models.sh in order to generate
# several tables and plots for a presentation that aims at comparing several transcript sets (typically
# predicted transcripts from a givem rnaseq sample) to a reference set of transcripts (typically the 
# transcripts that we expect to be present in this rnaseq sample)

# On March 27th 2020, all calls to _nbex.tsv file replaced by _nbex_intermclass.tsv file
# because of a change in refine_comptr script

# takes as input
################
# - a gtf or gff version 2 file containing the reference set of transcripts (typically the transcripts
#   we expect to be present in the rnaseq sample on which many transcript modellers have been run)
#   This file has to contain at least exon rows and should have the gene_id and transcript_id first in 
#   the 9th field. This file should be the same as the one given as input to analyse_transcript_models.sh
#   for each prediction set.
# - a 2 column tsv file with header that lists several predicted transcript sets with the first column
#   indicating the source (name of the program for example) and the second column indicating the absolute
#   path to the transcript model gtf or gff version 2 file, and where we expect the outputs of 
#   analyse_transcript_models.sh are with respect to this transcript set. Note that this file could include
#   the reference set itself since we usually want to compare the predicted transcripts with respect to it

# produces as output
####################
# many tables and plots for a presentation

# Note: expects the output files from the script analyse_transcript_models.sh for each prediction set in their directory
# Note: cannot be run twice in the same dir since uses files with fixed names
# Note: Needs R to be installed as well as the R libraries ggplot2 and reshape2


# example
# cd ~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie/WithoutAnnot/Test
# annot=~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Annot/gencode.v19.annotation.gtf
# time compare_multiple_trsets_wrt_reftrset.sh $annot annot_to_eval_small.tsv 2> compare_multiple_trsets_wrt_reftrset.err
# real    0m50.384s



# Check the inputs do exist
###########################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: compare_multiple_trsets_wrt_reftrset.sh reftrset.gff predtrsets.tsv >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a gtf or gff version 2 file containing the reference set of transcripts (typically the transcripts" >&2
    echo "  we expect to be present in the rnaseq sample on which many transcript modellers have been run)." >&2
    echo "  This file has to contain at least exon rows and should have the gene_id and transcript_id first" >&2
    echo "  in the 9th field. This file should be the same as the one given as input to analyse_transcript_models.sh" >&2
    echo "  for each prediction set." >&2
    echo "- a 2 column tsv file with header that lists several predicted transcript sets with the first column" >&2
    echo "  indicating the source (name of the program for example) and the second column indicating the absolute" >&2
    echo "  path to the transcript model gtf or gff version 2 file, and where we expect the outputs of " >&2
    echo "  analyse_transcript_models.sh are with respect to this transcript set. Note that this file could include" >&2
    echo "  the reference set itself since we usually want to compare the predicted transcripts with respect to it" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- some tables and plots about the success of the different predictions that can be used for a ppt" >&2
    echo "" >&2
    echo "Note: expects the output files from the script analyse_transcript_models.sh for each prediction set in their directory" >&2
    echo "Note: cannot be run twice in the same directory since would erase the output of the previous run" >&2
    echo "Note: Needs R to be installed as well as the R libraries ggplot2 and reshape2" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
ref=$1
predtrsets=$2
refbasetmp=`basename ${ref%.gtf}`
refbase=${refbasetmp%.gff}
refdir=`dirname $ref`
reftss=$refdir/$refbase\_capped_sites.gff
refnrtss=$refdir/$refbase\_capped_sites_nr.gff
reftts=$refdir/$refbase\_tts_sites.gff
refnrtts=$refdir/$refbase\_tts_sites_nr.gff

# Programs
##########
GFFOK=$rootDir/../Awk/make_gff_ok.awk
REFINEOUTPUT=$rootDir/../Awk/refine_comptr_to_table_stats.awk
BOXPLOT=$rootDir/boxplots.sh

# Start of the script
#####################
# Organize the directory structure
##################################
echo "I am organizing the directory structure" >&2
mkdir -p Plots Tables
mkdir -p Plots/ExonPerTranscript
mkdir -p Plots/TranscriptPerGene
mkdir -p Plots/ExonLength
mkdir -p Plots/DistinctExonLength
mkdir -p Plots/5pExonLength_Tr
mkdir -p Plots/5pExonLength_Gn
mkdir -p Plots/3pExonLength_Tr
mkdir -p Plots/3pExonLength_Gn
mkdir -p Plots/InternalExonLength
mkdir -p Plots/DistinctInternalExonLength
mkdir -p Plots/MonoExTrExLength
mkdir -p Plots/TrLength
mkdir -p Plots/cDNALength
mkdir -p Plots/Exact_tr_dist_to_Genc_TSS
echo "done" >&2


# List the ids of the spliced transcripts in the reference, making sure the transcript_id is in the 12th field by applying $GFFOK before
########################################################################################################################################
echo "I am listing the ids of the spliced transcripts in the reference" >&2
awk -f $GFFOK $ref | awk '$3=="exon"{split($12,a,"\""); nbex[a[2]]++}END{for(t in nbex){if(nbex[t]>=2){print t}}}' > $refdir/ref_spliced_tr.txt
echo "done" >&2
# ENST00000000233.5
# 171656 (1 fields)

# Compute some basic statistics about the basic categories of spliced transcripts with respect to the reference transcripts
###########################################################################################################################
echo "I am computing some basic statistics about the basic categories of spliced transcripts with respect to the reference transcripts" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
echo $src
awk 'NR>=2{print $2}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv | sort | uniq -c | sort -k1,1nr
done > Tables/basic_sumstats_from_comptr.txt
cat Tables/basic_sumstats_from_comptr.txt >&2
echo "done" >&2
# Cuff_annot_all
#  171130 Exact
#   80763 Monoexonic
#    4073 Unstranded
#    3056 Extension
#    1046 Overlap
#     984 Intergenic_or_antisense
#     771 Inclusion
# ...


# Make a table with the distribution of predicted transcripts in the different classes and with respect to the reference set
############################################################################################################################
echo "I am making a table with the distribution of predicted transcripts in the different classes and with respect to the reference set" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
awk -v lid=$src -f $REFINEOUTPUT $base\_complete_comp_refinedclass_nbex_intermclass.tsv
done | awk 'BEGIN{OFS="\t"} NR==1||NR%2==0' > Tables/refine_comptr_prediction_sets_for_table.txt
cat Tables/refine_comptr_prediction_sets_for_table.txt >&2
echo "done" >&2
# Cuff_annot_all  261824  232400  88.7619 171130  65.3607 771     0.294473        423     0.161559        60076   22.9452 3802    1.45212 3056    1.1672  623     0.237946        123 0.0469781       19686   7.51879 689     0.263154        18997   7.25564 1862    0.711165        295     0.112671        1567    0.598494                1.55563


# Make a table of basic summary statistics for the different prediction sets
#############################################################################
echo "I am making a table of basic summary statistics for the different prediction sets" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
awk 'NR==2' make_summary_stat_from_annot_$base.out | awk -v set=$src '{OFS="\t"; print set, $0, $1/$3, $1/$4, $2/$4, $3/$4, $5/$3, $5/$4, $6/$4}' 
done | awk 'BEGIN{print "prediction\tnbex\tnnbdistinctex\tnbtr\tnbgn\tnbintrons\tnbdistinctintrons\tnbexpertr\tnbexpergn\tnbdistinctexpergn\tnbtrpergn\tnbintronpertr\tnbintronpergn\tnbdistinctintronpergn"}{print}' > Tables/predtrsets_basic_sumstats.tsv
cat Tables/predtrsets_basic_sumstats.tsv >&2
echo "done" >&2
# predotation      nbex    nnbdistinctex   nbtr    nbgn    nbintrons       nbdistinctintrons       nbexpertr       nbexpergn       nbdistinctexpergn       nbtrpergn       nbintronpertr      nbintronpergn   nbdistinctintronpergn
# Gencv19_all     1196293 562680 196520 57820 999773 344651       6.08739 20.69   9.73158 3.39882 5.08739 17.2911 5.96076

# Gather across all prediction sets the lengths of several objects in order to have one plot per object type
############################################################################################################
# making sure transcript_id is in the 12th field by applying $GFFOK before
##########################################################################
echo "I am gathering across all prediction sets the lengths of several objects in order to have one plot per object type" >&2

awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
    WORKDIR=`dirname $pred`
    basepred=`basename $pred`
    cd $WORKDIR
    lid=$no\_$src
    awk -f $GFFOK $pred | awk -v lid=$lid 'BEGIN{OFS="\t"} $3=="exon"{nbex[$12]++}END{for(t in nbex){print lid, nbex[t]}}' 
done | awk 'BEGIN{OFS="\t"; print "prediction", "nb_ex_in_transcript"}{print}' > Plots/ExonPerTranscript/prediction_nbexintr_forggplot.tsv
echo "  done" >&2

echo "  2. number of transcripts per gene" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
lid=$no\_$src
awk -f $GFFOK $pred | awk -v lid=$lid 'BEGIN{OFS="\t"} $3=="exon"{seen[$12]++; if(seen[$12]==1){nbtr[$10]++}}END{for(g in nbtr){print lid, nbtr[g]}}' 
done | awk 'BEGIN{OFS="\t"; print "prediction", "nb_tr_in_gene"}NF==2{print}' > Plots/TranscriptPerGene/prediction_nbtringn_forggplot.tsv
echo "  done" >&2

echo "  3. exon length" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' exon_length_$base.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "exon_length"}{print}' > Plots/ExonLength/prediction_exlg_forggplot.tsv
echo "  done" >&2

echo "  4. distinct exon length" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' distinct_exon_length_$base.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "distinct_exon_length"}{print}' > Plots/DistinctExonLength/prediction_distexlg_forggplot.tsv
echo "  done" >&2

echo "  5. transcript 5' exon length (for spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' tr_5p_exon_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "fivep_exon_length_tr"}{print}' > Plots/5pExonLength_Tr/prediction_5pexlgtr_forggplot.tsv
echo "  done" >&2

echo "  6. gene 5' exon length (for genes with spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' gn_5p_exon_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "fivep_exon_length_gn"}{print}' > Plots/5pExonLength_Gn/prediction_5pexlggn_forggplot.tsv
echo "  done" >&2

echo "  7. transcript 3' exon length (for spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' tr_3p_exon_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "threep_exon_length_tr"}{print}' > Plots/3pExonLength_Tr/prediction_3pexlgtr_forggplot.tsv
echo "  done" >&2

echo "  8. gene 3' exon length (for genes with spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' gn_3p_exon_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "threep_exon_length_gn"}{print}' > Plots/3pExonLength_Gn/prediction_3pexlggn_forggplot.tsv
echo "  done" >&2

echo "  9. internal exon length (for spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
lid=$no\_$src
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' internal_exons_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "internal_exon_length"}{print}' > Plots/InternalExonLength/prediction_internexlg_forggplot.tsv
echo "  done" >&2

echo "  10. distinct internal exon length (for spliced and stranded tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' distinct_internal_exon_length_$base.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "distinct_internal_exon_length"}{print}' > Plots/DistinctInternalExonLength/prediction_distinternexlg_forggplot.tsv
echo "  done" >&2

echo "  11. monoexonic transcript exon length" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $1}' monoextr_exons_$base\_length.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "monoextr_exon_length"}{print}' > Plots/MonoExTrExLength/prediction_monoextrexlg_forggplot.tsv
echo "  done" >&2

echo "  12. transcript length (exons + introns) (for spliced transcripts)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $2}' tr_lgwithexintr_$base.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "transcript_length"}{print}' > Plots/TrLength/prediction_trlg_forggplot.tsv
echo "  done" >&2

echo "  13. cDNA length (exons only) (for spliced tr)" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid 'BEGIN{OFS="\t"}{print lid, $2}' tr_cumulexlg_$base.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "cDNA_length"}{print}' > Plots/cDNALength/prediction_cdnalg_forggplot.tsv
echo "  done" >&2

# Compute Sn and precision using 3 different defintions of TP
#############################################################
# trid  comptrclass     annottrlist     refinedtrclass  nbex
# tSTRG.293.3     Extension       ENST00000457540.1,ENST00000414273.1,    extension       n2
# 163350 (5 fields)
echo "I am computing the sensitivity and precision of each prediction set with respect to the reference transcripts" >&2
echo "using three definitions of true positives (Exact, Exact+Extension, Exact+Extension+Inclusion)" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
tot_ref=`wc -l $refdir/ref_spliced_tr.txt | awk '{print $1}'`
tp1_ref=`awk '$2=="Exact"{split($3,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv | sort | uniq | wc -l | awk '{print $1}'`
tp2_ref=`awk -v fileRef=$refdir/ref_spliced_tr.txt 'BEGIN{while (getline < fileRef >0){ok[$1]=1}} $2=="Exact"||$2=="Extension"{split($3,a,","); k=1; while(a[k]!=""){if(ok[a[k]]==1){print a[k]} k++}}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv | sort | uniq | wc -l | awk '{print $1}'`
tp3_ref=`awk -v fileRef=$refdir/ref_spliced_tr.txt 'BEGIN{while (getline < fileRef >0){ok[$1]=1}} $2=="Exact"||$2=="Extension"||$2=="Inclusion"{split($3,a,","); k=1; while(a[k]!=""){if(ok[a[k]]==1){print a[k]} k++}}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv | sort | uniq | wc -l | awk '{print $1}'`
awk -v lid=$src -v tp1_ref=$tp1_ref -v tp2_ref=$tp2_ref -v tp3_ref=$tp3_ref -v tot_ref=$tot_ref 'BEGIN{OFS="\t"}{split($0,a,"\t"); tot++; if(a[2]=="Monoexonic"){mono++}else{if(a[2]=="Unstranded"){unstr++}else{if(a[2]=="Exact"){exact++}else{if(a[2]=="Extension"){ext++}else{if(a[2]=="Inclusion"){incl++}}}}}}END{splstr=tot-(mono+unstr); tp1_pred=exact; tp2_pred=exact+ext; tp3_pred=exact+ext+incl; sn1=tp1_ref/tot_ref*100; sp1=tp1_pred/splstr*100; sn2=tp2_ref/tot_ref*100; sp2=tp2_pred/splstr*100; sn3=tp3_ref/tot_ref*100; sp3=tp3_pred/splstr*100; print lid, tot, splstr, tp1_ref, sn1, tp1_pred, sp1, (sn1+sp1)/2, tp2_ref, sn2, tp2_pred, sp2, (sn2+sp2)/2, tp3_ref, sn3, tp3_pred, sp3, (sn3+sp3)/2}' $base\_complete_comp_refinedclass_nbex_intermclass.tsv
done | awk 'BEGIN{OFS="\t"; print "lid", "tot_tr", "str_spl_tr", "tp1_exact_ref", "sn1_over_ref", "tp1_exact_pred", "sp1_over_pred", "avg1_sn1_sp1", "tp2_exact_ext_ref", "sn2_over_ref", "tp2_exact_ext_pred", "sp2_over_pred", "avg2_sn2_sp2", "tp3_exact_ext_incl_ref", "sn3_over_ref", "tp3_exact_ext_incl_pred", "sp3_over_pred", "avg3_sn3_sp3"}{print}' > Tables/prediction_sets_eval_wrt_ref_for_table.txt
cat Tables/prediction_sets_eval_wrt_ref_for_table.txt >&2
echo "done" >&2
# lid     tot_tr  str_spl_tr      tp1_exact_ref sn1_over_ref  tp1_exact_pred  sp1_over_pred   avg1_sn1_sp1    tp2_exact_ext_ref     sn2_over_ref  tp2_exact_ext_pred      sp2_over_pred      avg2_sn2_sp2    tp3_exact_ext_incl_ref        sn3_over_ref  tp3_exact_ext_incl_pred sp3_over_pred   avg3_sn3_sp3
# Cuff_annot_all  261824  176988  170973  99.6021 171130  96.6902 98.1461 170973  99.6021 174186  98.4168 99.0095 171078  99.6633 174957  98.8525 99.2579

# Make table saying for different min dist x (50, 100, 500) the number and prop of exact tr of each prediction set
##################################################################################################################
# that has its TSS closer than x to one of the matching reference tr, and the same for tts
###########################################################################################
echo "I am making a table saying for different min dist x (50, 100, 500) the number and prop of exact tr of each prediction set that has its TSS closer than x to one of the matching reference tr, and the same for tts" >&2
for dist in 50 100 500
do
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
echo $dist >&2
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
awk -v lid=$src -v dist=$dist 'BEGIN{OFS="\t"}{tot++; if($7<=dist){tss++;} if($8<=dist){tts++;}}END{print lid, tot, tss, tss/tot*100, tts, tts/tot*100}' $base\_predtr_spliced_stranded_exact_$refbase\_correstrlist_TSS_AnnotTSSlist_TTS_AnnotTTSlist_smallerdist_TSS_TTS.txt
done | awk -v dist=$dist 'BEGIN{OFS="\t"; print "pred_set", "nb_exact_tr", "hitting_ref_tss_"dist"bp_nb", "hitting_ref_tss_"dist"bp_pcent", "hitting_ref_tts_"dist"bp_nb", "hitting_ref_tts_"dist"bp_pcent"}{print}' > Tables/prediction_sets_nbexacttr_hitting_ref_tss_$dist\bp_nb_pcent_hitting_ref_tts_$dist\bp_nb_pcent.tsv 
cat Tables/prediction_sets_nbexacttr_hitting_ref_tss_$dist\bp_nb_pcent_hitting_ref_tts_$dist\bp_nb_pcent.tsv >&2
done
echo "done" >&2
# pred_set        nb_exact_tr     hitting_ref_tss_50bp_nb       hitting_ref_tss_50bp_pcent    hitting_ref_tts_50bp_nb       hitting_ref_tts_50bp_pcent
# Gencv19_all     171656  171656  100     171656  100

# Distribution of the smaller distance from predicted TSS to its matching refated transcript TSS for all prediction sets
########################################################################################################################
echo "I am computing the distribution of the smaller distance from predicted TSS to its matching refated transcript TSS for all prediction sets" >&2
awk 'NR>=2{print $1, $2, ((NR-2)<=9 ? "0"(NR-2) : (NR-2))}' $predtrsets | while read src pred no
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
lid=$no\_$src
awk -v lid=$lid -v dist=$dist 'BEGIN{OFS="\t"}{print lid, $7}' $base\_predtr_spliced_stranded_exact_$refbase\_correstrlist_TSS_AnnotTSSlist_TTS_AnnotTTSlist_smallerdist_TSS_TTS.txt
done | awk 'BEGIN{OFS="\t"; print "prediction", "dist_to_tss"}{print}' > Plots/Exact_tr_dist_to_Genc_TSS/prediction_tssdist_forggplot.tsv
echo "done" >&2

# Make a simple table with number of stranded predicted transcripts, their genes, their individual and nr TSS and TTS
#####################################################################################################################
echo "I am making a simple table with number of stranded predicted transcripts, their genes, their individual and nr TSS and TTS" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
ntng=`awk '{seeng[$10]++; if(seeng[$10]==1){ng++} split($12,a,"\""); split(a[2],b,","); k=1; while(b[k]!=""){seent[b[k]]++; if(seent[b[k]]==1){nt++} k++}}END{print nt, ng}' $base\_capped_sites_nr.gff`
tss=`wc -l $base\_capped_sites.gff | awk '{print $1}'`
tssnr=`wc -l $base\_capped_sites_nr.gff | awk '{print $1}'`
tts=`wc -l $base\_tts_sites.gff | awk '{print $1}'`
ttsnr=`wc -l $base\_tts_sites_nr.gff | awk '{print $1}'`
echo $b $ntng $tss $tssnr $tts $ttsnr
done | awk 'BEGIN{OFS="\t"; print "pred_set", "nbtr", "nbgn", "ntss", "ntssnr", "ntts", "nttsnr"}{print $1, $2, $3, $4, $5, $6, $7}' > Tables/prediction_sets_nbstrtr_gn_nbtss_nbtssnr_nbtts_nbttsnr.tsv
cat Tables/prediction_sets_nbstrtr_gn_nbtss_nbtssnr_nbtts_nbttsnr.tsv >&2
echo "done" >&2
# pred_set        nbtr    nbgn    ntss    ntssnr  ntts    nttsnr
# Gencv19_all     196520  57820   196520  179185  196520  120235

# Make a sensitivity table for nr reference tss and how they are hit by any predicted tss
#########################################################################################
echo "I am making a sensitivity table for nr reference tss and how they are hit by any predicted tss" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
distreftss=`wc -l $refnrtss | awk '{print $1}'`
tsshit1=`wc -l $refbase\_capped_sites_nr_with_pred_less50bp.gff | awk '{print $1}'`
tsshit2=`wc -l $refbase\_capped_sites_nr_with_pred_less100bp.gff | awk '{print $1}'`
tsshit3=`wc -l $refbase\_capped_sites_nr_with_pred_less500bp.gff | awk '{print $1}'`
echo $src $distreftss $tsshit1 $tsshit2 $tsshit3 
done | awk 'BEGIN{OFS="\t"; print "pred_set", "dist_ref_tss", "tss_hit_by_pred_50bp_nb", "tss_hit_by_pred_50bp_pcent", "tss_hit_by_pred_100bp_nb", "tss_hit_by_pred_100bp_pcent", "tss_hit_by_pred_500bp_nb", "tss_hit_by_pred_500bp_pcent"} {print $1, $2, $3, $3/$2*100, $4, $4/$2*100, $5, $5/$2*100}' > Tables/prediction_sets_nrreftss_hitbypredtss_50bp_100bp_500bp_nb_pcent.tsv
cat Tables/prediction_sets_nrreftss_hitbypredtss_50bp_100bp_500bp_nb_pcent.tsv >&2
echo "done" >&2
# pred_set        dist_ref_tss   tss_hit_by_pred_50bp_nb tss_hit_by_pred_50bp_pcent      tss_hit_by_pred_100bp_nb        tss_hit_by_pred_100bp_pcent     tss_hit_by_pred_500bp_nb  tss_hit_by_pred_500bp_pcent
# Cuff_annot_all  179185  179013  99.904  179013  99.904  179014  99.9046

# Make a sensitivity table for nr reference tts and how they are hit by any predicted tts
#########################################################################################
echo "I am making a sensitivity table for nr reference tts and how they are hit by any predicted tts" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
distreftts=`wc -l $refnrtts | awk '{print $1}'`
ttshit1=`wc -l $refbase\_tts_sites_nr_with_pred_less50bp.gff | awk '{print $1}'`
ttshit2=`wc -l $refbase\_tts_sites_nr_with_pred_less100bp.gff | awk '{print $1}'`
ttshit3=`wc -l $refbase\_tts_sites_nr_with_pred_less500bp.gff | awk '{print $1}'`
echo $src $distreftts $ttshit1 $ttshit2 $ttshit3
done | awk 'BEGIN{OFS="\t"; print "pred_set", "dist_ref_tts", "tts_hit_by_pred_50bp_nb", "tts_hit_by_pred_50bp_pcent", "tts_hit_by_pred_100bp_nb", "tts_hit_by_pred_100bp_pcent", "tts_hit_by_pred_500bp_nb", "tts_hit_by_pred_500bp_pcent"} {print $1, $2, $3, $3/$2*100, $4, $4/$2*100, $5, $5/$2*100}' > Tables/prediction_sets_nrreftts_hitbypredtts_50bp_100bp_500bp_nb_pcent.tsv
cat Tables/prediction_sets_nrreftts_hitbypredtts_50bp_100bp_500bp_nb_pcent.tsv >&2
echo "done" >&2
# pred_set        dist_ref_tts   tts_hit_by_pred_50bp_nb tts_hit_by_pred_50bp_pcent      tts_hit_by_pred_100bp_nb        tts_hit_by_pred_100bp_pcent     tts_hit_by_pred_500bp_nb  tts_hit_by_pred_500bp_pcent
# Cuff_annot_all  120235  119763  99.6074 119821  99.6557 120022  99.8228

# Make a precision table for all distinct predicted tss with respect to the reference tss
#########################################################################################
echo "I am making a precision table for all distinct predicted tss with respect to the reference tss" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
distpredtss=`wc -l $base\_capped_sites_nr.gff | awk '{print $1}'`
tsshit1=`wc -l $base\_capped_sites_nr_with_annottss_less50bp.gff | awk '{print $1}'`
tsshit2=`wc -l $base\_capped_sites_nr_with_annottss_less100bp.gff | awk '{print $1}'`
tsshit3=`wc -l $base\_capped_sites_nr_with_annottss_less500bp.gff | awk '{print $1}'`
echo $src $distpredtss $tsshit1 $tsshit2 $tsshit3 
done | awk 'BEGIN{OFS="\t"; print "pred_set", "dist_pred_tss", "tss_hit_by_pred_50bp_nb", "tss_hit_by_pred_50bp_pcent", "tss_hit_by_pred_100bp_nb", "tss_hit_by_pred_100bp_pcent", "tss_hit_by_pred_500bp_nb", "tss_hit_by_pred_500bp_pcent"} {print $1, $2, $3, $3/$2*100, $4, $4/$2*100, $5, $5/$2*100}' > Tables/prediction_sets_nrpredtss_hitbypredtss_50bp_100bp_500bp_nb_pcent.tsv
cat Tables/prediction_sets_nrpredtss_hitbypredtss_50bp_100bp_500bp_nb_pcent.tsv >&2
echo "done" >&2
# pred_set        dist_pred_tss   tss_hit_by_pred_50bp_nb tss_hit_by_pred_50bp_pcent      tss_hit_by_pred_100bp_nb        tss_hit_by_pred_100bp_pcent     tss_hit_by_pred_500bp_nb  tss_hit_by_pred_500bp_pcent
# Cuff_annot_all  201752  179345  88.8938 179549  88.9949 180065  89.2507

# Make a precision table for all distinct predicted tts with respect to the reference tts
#########################################################################################
echo "I am making a precision table for all distinct predicted tts with respect to the reference tts" >&2
awk 'NR>=2{print $1, $2}' $predtrsets | while read src pred
do
WORKDIR=`dirname $pred` 
cd $WORKDIR
basetmp=`basename ${pred%.gtf}`
base=${basetmp%.gff}
distpredtts=`wc -l $base\_tts_sites_nr.gff | awk '{print $1}'`
ttshit1=`wc -l $base\_tts_sites_nr_with_annottts_less50bp.gff | awk '{print $1}'`
ttshit2=`wc -l $base\_tts_sites_nr_with_annottts_less100bp.gff | awk '{print $1}'`
ttshit3=`wc -l $base\_tts_sites_nr_with_annottts_less500bp.gff | awk '{print $1}'`
echo $src $distpredtts $ttshit1 $ttshit2 $ttshit3 
done | awk 'BEGIN{OFS="\t"; print "pred_set", "dist_pred_tts", "tts_hit_by_pred_50bp_nb", "tts_hit_by_pred_50bp_pcent", "tts_hit_by_pred_100bp_nb", "tts_hit_by_pred_100bp_pcent", "tts_hit_by_pred_500bp_nb", "tts_hit_by_pred_500bp_pcent"} {print $1, $2, $3, $3/$2*100, $4, $4/$2*100, $5, $5/$2*100}' > Tables/prediction_sets_nrpredtts_hitbypredtts_50bp_100bp_500bp_nb_pcent.tsv
cat Tables/prediction_sets_nrpredtts_hitbypredtts_50bp_100bp_500bp_nb_pcent.tsv >&2
echo "done" >&2
# pred_set        dist_pred_tts   tts_hit_by_pred_50bp_nb tts_hit_by_pred_50bp_pcent      tts_hit_by_pred_100bp_nb        tts_hit_by_pred_100bp_pcent     tts_hit_by_pred_500bp_nb  tts_hit_by_pred_500bp_pcent
# Cuff_annot_all  176744  119547  67.6385 119757  67.7573 122787  69.4717

# Make the plots for numbers and lengths of several objects across all the predictions
######################################################################################
echo "I am making the plots for numbers and lengths of several objects across all the predictions" >&2
echo "  1. number of exons per transcript" >&2
$BOXPLOT Plots/ExonPerTranscript/prediction_nbexintr_forggplot.tsv prediction nb_ex_in_transcript "Number of exons per transcript" 0 20 Plots/ExonPerTranscript/prediction_nbexintr_forggplot.png
echo "  done" >&2
echo "  2. number of transcripts per gene" >&2
$BOXPLOT Plots/TranscriptPerGene/prediction_nbtringn_forggplot.tsv prediction nb_tr_in_gene "Number of transcripts per gene" 0 20 Plots/TranscriptPerGene/prediction_nbtringn_forggplot.png
echo "  done" >&2
echo "  3. exon length" >&2
$BOXPLOT Plots/ExonLength/prediction_exlg_forggplot.tsv prediction exon_length "Exon length" 0 500 Plots/ExonLength/prediction_exlg_forggplot.png
echo "  done" >&2
echo "  4. distinct exon length" >&2
$BOXPLOT Plots/DistinctExonLength/prediction_distexlg_forggplot.tsv prediction distinct_exon_length "Distinct exon length" 0 500 Plots/DistinctExonLength/prediction_distexlg_forggplot.png
echo "  done" >&2
echo "  5. transcript 5' exon length (for spliced and stranded tr)" >&2
$BOXPLOT Plots/5pExonLength_Tr/prediction_5pexlgtr_forggplot.tsv prediction fivep_exon_length_tr "5p exon length (spliced, stranded tr)" 0 500 Plots/5pExonLength_Tr/prediction_5pexlgtr_forggplot.png
echo "  done" >&2
echo "  6. gene 5' exon length (for genes with spliced and stranded tr)" >&2
$BOXPLOT Plots/5pExonLength_Gn/prediction_5pexlggn_forggplot.tsv prediction fivep_exon_length_gn "5p exon length of genes (with spl str tr)" 0 500 Plots/5pExonLength_Gn/prediction_5pexlggn_forggplot.png
echo "  done" >&2
echo "  7. transcript 3' exon length (for spliced and stranded tr)" >&2
$BOXPLOT Plots/3pExonLength_Tr/prediction_3pexlgtr_forggplot.tsv prediction threep_exon_length_tr "3p exon length (spliced, stranded tr)" 0 500 Plots/3pExonLength_Tr/prediction_3pexlgtr_forggplot.png
echo "  done" >&2
echo "  8. gene 3' exon length (for genes with spliced and stranded tr)" >&2
$BOXPLOT Plots/3pExonLength_Gn/prediction_3pexlggn_forggplot.tsv prediction threep_exon_length_gn "3p exon length of genes (with spl str tr)" 0 500 Plots/3pExonLength_Gn/prediction_3pexlggn_forggplot.png
echo "  done" >&2
echo "  9. internal exon length (for spliced and stranded tr)" >&2
$BOXPLOT Plots/InternalExonLength/prediction_internexlg_forggplot.tsv prediction internal_exon_length "Internal exon length (spliced, stranded tr)" 0 260 Plots/InternalExonLength/prediction_internexlg_forggplot.png
echo "  done" >&2
echo "  10. distinct internal exon length (for spliced and stranded tr)" >&2
$BOXPLOT Plots/DistinctInternalExonLength/prediction_distinternexlg_forggplot.tsv prediction distinct_internal_exon_length "Distinct internal exon length (spliced, stranded tr)" 0 260 Plots/DistinctInternalExonLength/prediction_distinternexlg_forggplot.png
echo "  done" >&2
echo "  11. monoexonic transcript exon length" >&2
$BOXPLOT Plots/MonoExTrExLength/prediction_monoextrexlg_forggplot.tsv prediction monoextr_exon_length "Monoexonic transcript exon length" 0 500 Plots/MonoExTrExLength/prediction_monoextrexlg_forggplot.png
echo "  done" >&2
echo "  12. transcript length (exons + introns) (for spliced transcripts)" >&2
$BOXPLOT Plots/TrLength/prediction_trlg_forggplot.tsv prediction transcript_length "Transcript length (exons + introns) (spliced tr)" 0 50000 Plots/TrLength/prediction_trlg_forggplot.png
echo "  done" >&2
echo "  13. cDNA length (exons only) (for spliced tr)" >&2
$BOXPLOT Plots/cDNALength/prediction_cdnalg_forggplot.tsv prediction cDNA_length "cDNA length (only exons) (for spliced tr)" 0 4000 Plots/cDNALength/prediction_cdnalg_forggplot.png
echo "  done" >&2
echo "done" >&2


# Make the plot for distance between predicted and reference TSS for exact transcripts
######################################################################################
echo "I am making the plot for distance between predicted and reference TSS for exact transcripts" >&2
$BOXPLOT Plots/Exact_tr_dist_to_Genc_TSS/prediction_tssdist_forggplot.tsv prediction dist_to_tss "Distance to Gencode TSS (exact tr)" 0 10 Plots/Exact_tr_dist_to_Genc_TSS/prediction_tssdist_forggplot.png
echo "done" >&2

# Clean
########
echo "I am cleaning" >&2
rm Rplots.pdf
echo "done" >&2
