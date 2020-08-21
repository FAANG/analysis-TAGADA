#!/bin/bash

# ref_pred_all_expr_2_tr_ratio_size_pos.sh
##########################################
# takes as input:
#################
# - a reference gene annotation file in gtf (such as the ones provided by ensembl) or gff2 format 
# - a predicted gene annotation file in gtf or gff2 format
# - a matrix of transcript TPM values in a set of samples (tsv file with header whose columns are transcript_id, gene_id and TPM values in a set of samples)
#   this file should contain transcript ids from the reference and from the predicted gene annotation
# - an optional output directory (otherwise will write in $PWD)
# outputs:
##########
# - the reference and predicted transcript gff file filtered for TPM above 0.1 in at least 2 samples"
# - a Plots and a Tables directory with plots and tables about ratios, sizes and positions of transcripts from 4 different sets"
#   (ref, ref_expr, string, string_expr) with respect to the reference annotation"
# - some .err files about the different processing steps"

# TO DO
#######
# make the script more general for
# - format of transcript expression matrix, being able to specify the fld no from which the first TPM occurs in the file
# - min expression threshold
# - min number of samples where the min expression has to occur

# example
#########
# cd /work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/elements
# module load bioinfo/bedtools2-2.29.0
# module load system/R-3.6.2
# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/ref_pred_all_expr_2_tr_ratio_size_pos.sh
# ref=/work/project/fragencode/workspace/geneswitch/pipelines/rnaseq/tests/sus_scrofa/allbig/data/species/sus_scrofa.gtf
# pred=/work/project/fragencode/workspace/geneswitch/results/counts/rnaseq.sus_scrofa.liver.pig1.gff
# expr=/work/project/fragencode/workspace/geneswitch/results/counts/transcripts_TPM.tsv
# time $pgm $ref $pred $expr 2> ref_pred_all_expr_2_tr_ratio_size_pos.err

# $expr
# transcript_id	gene_id	TPM_rnaseq.sus_scrofa.cd8.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig4.R1	TPM_rnaseq.sus_scrofa.cd4.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig2.R1	TPM_rnaseq.sus_scrofa.liver.pig3.R1	TPM_rnaseq.sus_scrofa.cd8.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig1	TPM_rnaseq.sus_scrofa.liver.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig1.R1
# ENSSSCT00000000003	MSTRG.11569	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# 88267 (13 fields)

# check all the needed inputs are there
########################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo "Usage: ref_pred_annot2stats_number_size_position.sh ref_annot.gtf pred_annot.gff tr_TPM.tsv <outdir>" &2
    echo "" >&2
    echo where >&2
    echo "- ref_annot.gtf and pred_annot.gff are two gtf or gff2 files with at least exons rows and with gene_id and transcript_id in the 9th field" >&2
    echo "- tr_TPM.tsv is a tsv file with header that has transcript id, gene id and then TPM expression in a set of samples" >&2
    echo "- outdir is the optional output directory (default current working directory)" >&2
    echo "and will write to output directory" >&2
    echo "- the reference and predicted transcript gff file filtered for TPM above 0.1 in at least 2 samples" >&2
    echo "- a Plots and a Tables directory with plots and tables about ratios, sizes and positions of transcripts from 4 different sets" >&2
    echo "  (ref, ref_expr, string, string_expr) with respect to the reference annotation" >&2
    echo "- some .err files about the different processing steps" >&2
    echo "" >&2
    echo "Needs bedtools and R to be installed (tested with versions 2.29.0 and 3.6.2 resp) and to have ggplot2 library available in R" >&2
    exit 1
fi

# Set general variables
#######################
path="`dirname \"$0\"`" # relative path
rootdir="`( cd \"$path\" && pwd )`" # absolute path
ref=$1
refbase=`basename $ref`
str=$2
strbase=`basename $str`
tr=$3
if [ -n "$4" ]
then
    outdir=$4
else
    outdir=$PWD
fi

# Set variables of programs
###########################
EXPRFILTER=$rootdir/exprfilter_annot.awk
GFFOK=$rootdir/make_gff_ok.awk
ANALYSE=$rootdir/analyse_transcript_models.sh
COMPARE=$rootdir/compare_multiple_trsets_wrt_reftrset.sh


# 1. Make the two expression filtered gff files and hard copy the 4*2 files given as input to analyse_transcript_models in the 4 directories where it is run
#############################################################################################################################################################
cd $outdir
awk -v mexpr=0.1 -v msamp=2 -v fstexpr=3 -v fileRef=$tr -f $EXPRFILTER $ref | awk -f $GFFOK > ref.annot.tpm0.1.2samples.exons.gff
awk -v mexpr=0.1 -v msamp=2 -v fstexpr=3 -v fileRef=$tr -f $EXPRFILTER $str | awk -f $GFFOK > stringtie.annot.tpm0.1.2samples.exons.gff
mkdir -p $outdir/ref
mkdir -p $outdir/ref_expr
mkdir -p $outdir/string
mkdir -p $outdir/string_expr
cp $ref $outdir/ref
cp ref.annot.tpm0.1.2samples.exons.gff $ref $outdir/ref_expr
cp $str $ref $outdir/string
cp stringtie.annot.tpm0.1.2samples.exons.gff $ref $outdir/string_expr

# 2. Run the analyse transcript script on the 4 sets (ref gene, ref gene tpm filter, string gene, string gene tpm filter)
##########################################################################################################################
# but each one in a different directory and one after the other and also put a symb link to the gff file in each dir
####################################################################################################################
# !!! these scripts cannot be launched at the same time since produce files relative to ref that will be the same !!!
# !!! each of them take 3 minutes !!!
# a. for ref
############
cd $outdir/ref
$ANALYSE $refbase $refbase 2> analyse_transcript_models_ref.err

# b. for ref expr
#################
cd $outdir/ref_expr
$ANALYSE ref.annot.tpm0.1.2samples.exons.gff $refbase 2> analyse_transcript_models_ref_expr.err

# c. for stringtie
##################
cd $outdir/string
$ANALYSE $strbase $refbase 2> analyse_transcript_models_string.err

# d. for stringtie expr
########################
cd $outdir/string_expr
$ANALYSE stringtie.annot.tpm0.1.2samples.exons.gff $refbase 2> analyse_transcript_models_string_expr.err

# 3. produce tsv file of 4 sets to be used by the gathering info script
#######################################################################
cd $outdir
printf "lid\tfile\n" > trsets.tsv
printf "ref\t"$outdir"/ref/"$refbase"\n" >> trsets.tsv
printf "ref_expr\t"$outdir"/ref_expr/ref.annot.tpm0.1.2samples.exons.gff\n" >> trsets.tsv
printf "string\t"$outdir"/string/"$strbase"\n" >> trsets.tsv
printf "string_expr\t"$outdir"/string_expr/stringtie.annot.tpm0.1.2samples.exons.gff\n" >> trsets.tsv
# lid	file
# ref	/work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/elements/ref/sus_scrofa.gtf
# 5 (2 fields)

# 4. run the gather info script
################################
$COMPARE $ref trsets.tsv 2> compare_multiple_trsets_wrt_reftrset.err


# 5. Erase all intermediate directories
#######################################
# rm -r $outdir/ref
# rm -r $outdir/ref_expr
# rm -r $outdir/string
# rm -r $outdir/string_expr
