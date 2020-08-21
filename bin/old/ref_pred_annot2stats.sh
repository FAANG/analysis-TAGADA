#!/bin/bash

# ref_pred_annot2stats.sh
#########################
# takes as input:
#################
# - a reference gene annotation file in gtf (such as the ones provided by ensembl) or gff2 format 
# - a predicted gene annotation file in gtf or gff2 format
# - a matrix of transcript TPM values in a set of samples (tsv file with header whose columns are transcript_id, gene_id and TPM values in a set of samples)
#   this file should contain transcript ids from the reference and from the predicted gene annotation
# outputs:
##########
# - a tsv file with header called make_summary_stat_from_4annots.tsv that has:
###############################################################################
#   * 5 rows, one for the header, one for the reference gene annotation, one for the ref gene annotation where only
#     transcripts with tpm above 0.1 in at least 2 samples were retained, one for the predicted gene annotation and one for the pred gene annotation+
#     where only transcripts with tpm above 0.1 in at least 2 samples were retained
#   * 14 columns with number and ratios of elements in each of the 4 gene prediction sets

# TO DO
#######
# make the script more general for
# - format of transcript expression matrix, being able to specify the fld no from which the first TPM occurs in the file
# - min expression threshold
# - min number of samples where the min expression has to occur

# example
#########
# cd /work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/elements
# module load system/R-3.6.2
# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/ref_pred_annot2stats.sh
# ref=/work/project/fragencode/workspace/geneswitch/pipelines/rnaseq/tests/sus_scrofa/allbig/data/species/sus_scrofa.gtf
# pred=/work/project/fragencode/workspace/geneswitch/results/counts/rnaseq.sus_scrofa.liver.pig1.gff
# expr=/work/project/fragencode/workspace/geneswitch/results/counts/transcripts_TPM.tsv
# time $pgm $ref $pred $expr > ref_pred_annot2stats.tsv 2> ref_pred_annot2stats.err
# real	2m40.725s

# $expr
# transcript_id	gene_id	TPM_rnaseq.sus_scrofa.cd8.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig4.R1	TPM_rnaseq.sus_scrofa.cd4.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig2.R1	TPM_rnaseq.sus_scrofa.liver.pig3.R1	TPM_rnaseq.sus_scrofa.cd8.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig1	TPM_rnaseq.sus_scrofa.liver.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig1.R1
# ENSSSCT00000000003	MSTRG.11569	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# 88267 (13 fields)


# check all the inputs are there
################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage: ref_pred_annot2stats.sh ref_annot.gtf pred_annot.gff tr_TPM.tsv >&2
    echo "" >&2
    echo where >&2
    echo "- ref_annot.gtf and pred_annot.gff are two gtf or gff2 files with at least exons rows and with gene_id and transcript_id in the 9th field" >&2
    echo "- tr_TPM.tsv is a tsv file with header that has transcript id, gene id and then TPM expression in a set of samples" >&2
    echo "and will write to standard output" >&2
    echo "- a tsv file with header with number and ratio of elements in 4 gene annotations (ref annot, ref annot with transcripts whose TPM is at least 0.1" >&2
    echo "  in at least 2 samples, stringtie annot and stringtie annot with transcript whose TPM is at least 0.1 in at least 2 samples)" >&2
    echo "Needs R to be installed (tested with version 3.6.2)" >&2
    exit 1
fi

# Set general variables
#######################
path="`dirname \"$0\"`" # relative path
rootdir="`( cd \"$path\" && pwd )`" # absolute path
ref=$1
str=$2
expr=$3
basereftmp=`basename $ref`
basereftmp2=${basereftmp%.gtf}
baseref=${basereftmp2%.gff}
basestrtmp=`basename $str`
basestrtmp2=${basestrtmp%.gtf}
basestr=${basestrtmp2%.gff}

# Set variables of programs
###########################
gffok=$rootdir/../Awk/make_gff_ok.awk
exprfilter=$rootdir/../Awk/exprfilter_annot.awk
makesum=$rootdir/make_summary_stat_from_annot.sh

# 1. Make the exon gff file of reference transcripts with tpm 0.1 in 2 samples
##############################################################################
awk -v mexpr=0.1 -v msamp=2 -v fstexpr=3 -v fileRef=$expr -f $exprfilter $ref | awk -f $gffok > ref.annot.tpm0.1.2samples.exons.gff

# 2. Make the exon gff file of all stringtie transcripts with tpm 0.1 in 2 samples
##################################################################################
awk -v mexpr=0.1 -v msamp=2 -v fstexpr=3 -v fileRef=$expr -f $exprfilter $str | awk -f $gffok > stringtie.annot.tpm0.1.2samples.exons.gff

# 3. Apply the make_summary_stat_from_annot.sh script to the ref annot, the ref annot with tr with tpm 0.1 in 2 samples
#######################################################################################################################
# the stringtie annot, the stringtie annot with tr with tpm 0.1 in 2 samples
############################################################################
$makesum $ref > make_summary_stat_from_annot_ref.out
$makesum ref.annot.tpm0.1.2samples.exons.gff > make_summary_stat_from_annot_ref_expr.out
$makesum $str > make_summary_stat_from_annot_str.out
$makesum stringtie.annot.tpm0.1.2samples.exons.gff > make_summary_stat_from_annot_str_expr.out

# 4. Gather all the information in the wanted tsv file with header that is the output of the script
###################################################################################################
for lid in ref ref_expr str str_expr
do
    awk -v lid=$lid 'NR==2{OFS="\t"; print lid, $0}' make_summary_stat_from_annot_$lid.out
done | awk 'BEGIN{OFS="\t"; print "trset\tnbex\tnbdistinctex\tnbtr\tnbgn\tnbintrons\tnbdistinctintrons\tnbexpertr\tnbexpergn\tnbdistinctexpergn\tnbtrpergn\tnbintronpertr\tnbintronpergn\tnbdistinctintronpergn"}{print $0, $2/$4, $2/$5, $3/$5, $4/$5, $6/$4, $6/$5, $7/$5}' 
# nbex      nbdistinctex  nbtr    nbgn   nbintrons  nbdistinctintrons  nbexpertr  nbexpergn  nbdistinctexpergn  nbtrpergn  nbintronpertr  nbintronpergn  nbdistinctintronpergn
# ref       536859        262635  49448  25880      487411             217796     10.857     20.7442            10.1482    1.91066        9.85704        18.8335                8.41561
# ref_expr  357234        194945  29715  15646      327519             167254     12.022     22.8323            12.4597    1.89921        11.022         20.9331                10.6899
# str       1009341       315571  88266  26559      921075             256531     11.4352    38.0037            11.8819    3.32339        10.4352        34.6803                9.65891
# str_expr  805853        254904  66446  17277      739407             212563     12.1279    46.6431            14.754     3.84592        11.1279        42.7972                12.3032

# 5. Remove unuseful files
##########################
rm ref.annot.tpm0.1.2samples.exons.gff
rm stringtie.annot.tpm0.1.2samples.exons.gff
rm $baseref\_trid_nbex.txt
rm $baseref\_gnid_nbtr.txt
rm $baseref\_complete.gff
rm ref.annot.tpm0.1.2samples.exons_trid_nbex.txt
rm ref.annot.tpm0.1.2samples.exons_gnid_nbtr.txt
rm ref.annot.tpm0.1.2samples.exons_complete.gff
rm $basestr\_trid_nbex.txt
rm $basestr\_gnid_nbtr.txt
rm $basestr\_complete.gff
rm stringtie.annot.tpm0.1.2samples.exons_trid_nbex.txt
rm stringtie.annot.tpm0.1.2samples.exons_gnid_nbtr.txt
rm stringtie.annot.tpm0.1.2samples.exons_complete.gff
