#!/bin/bash

# compute_detected_tr_gn.sh
############################
# takes as input:
#################
# - a reference gene annotation file in gtf format (such as the ones provided by ensembl)
# - a matrix of transcript TPM values in a set of samples (tsv file with header whose columns are transcript_id, gene_id and TPM values in a set of samples)
# - a matrix of gene TPM values in a set of samples (tsv file with header whose columns are gene_id and TPM values in a set of samples)
# outputs in the current working directory:
###########################################
# - a tsv file with header that has:
####################################
#   * 3 rows: one for the header, one for transcripts and one for genes
#   * 9 columns: one for the element type, one for the number of elements in the reference, the number of elements output by stringtie and its % of the total,
#     then the number and % of total of the ref elements that have tpm 0.1 in at least 2 samples, the number of novel (not in ref) elements output by stringtie,
#     and the number and % of total novel elements that have tpm 0.1 in at least 2 samples

# TO DO
#######
# make the script more general for
# - format of expression matrices, being able to specify the fld no from which TPM occur in the two files
# - min expression threshold
# - min number of samples where the min expression has to occur


# example
#########
# cd /work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/elements
# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/compute_detected_tr_gn.sh
# ref=/work/project/fragencode/workspace/geneswitch/pipelines/rnaseq/tests/sus_scrofa/allbig/data/species/sus_scrofa.gtf
# tr=/work/project/fragencode/workspace/geneswitch/results/counts/transcripts_TPM.tsv
# gn=/work/project/fragencode/workspace/geneswitch/results/counts/genes_TPM.tsv
# time $pgm $ref $tr $gn > compute_detected_tr_gn.tsv 2> compute_detected_tr_gn.err
# Element     ref_nb  ref_in_file_nb  ref_in_file_pcent  ref_det_nb  ref_det_pcent  nov_nb  nov_det_nb  nov_det_pcent
# transcript  49448   49442           99.9879            29715       60.0934        38824   36731       94.609
# gene        25880   8492            32.813             1905        7.3609         18067   15414       85.3158
# real	0m30.100s


# check all the inputs are there
################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage: compute_detected_tr_gn.sh ref_annot.gtf tr_TPM.tsv gn_TPM.tsv >&2
    echo "" >&2
    echo "where" >&2
    echo "- ref_annot.gtf is a reference gene annotation in gtf format (ensembl-like)" >&2
    echo "- tr_TPM.tsv is a transcript expression matrix with transcript id, gene id and then TPM values in a set of samples" >&2
    echo "- gn_TPM.tsv is a gene expression matrix with gene id and then TPM values in a set of samples" >&2
    echo "produces in the standard output a tsv file with header that has the number of initial reference transcripts and genes" >&2
    echo "as well as those in the expression matrices and those with expression above 0.1 in at least 2 samples" >&2
    exit 1
fi


# Set general variables
#######################
path="`dirname \"$0\"`" # relative path
rootdir="`( cd \"$path\" && pwd )`" # absolute path
ref=$1
tr=$2
gn=$3

# Set variables of programs and of palettes
###########################################
extract=$rootdir/extract.gtf.tags.sh

# Extract the transcript ids from the reference gtf file
########################################################
awk '$3=="exon"' $ref | $extract - transcript_id | sort | uniq > ref_tr_id.txt
# ENSSSCT00000000003
# 49448 (1 fields)

# Extract the gene ids from the reference gtf file
##################################################
awk '$3=="exon"' $ref | $extract - gene_id | sort | uniq > ref_gn_id.txt
# ENSSSCG00000000002
# 25880 (1 fields) 

# Make the header of the output file
####################################
printf "Element\tref_nb\tref_in_file_nb\tref_in_file_pcent\tref_det_nb\tref_det_pcent\tnov_nb\tnov_det_nb\tnov_det_pcent\n"
# Element	ref_nb	ref_in_file_nb	ref_in_file_pcent	ref_det_nb	ref_det_pcent	nov_nb	nov_det_nb	nov_det_pcent

# Compute the stats for the transcript elements
###############################################
awk -v fileRef=ref_tr_id.txt 'BEGIN{OFS="\t"; while (getline < fileRef >0){refelt++; inref[$1]=1}} NR>=2{ok=0; k=3; while(ok<=1&&k<=NF){if($k>=0.1){ok++} k++} if(inref[$1]==1){refeltinfile++; if(ok==2){refeltdet++}} else{noveltinfile++; if(ok==2){noveltdet++}}} END{print "transcript", refelt, refeltinfile, refeltinfile/refelt*100, refeltdet, refeltdet/refelt*100, noveltinfile, noveltdet, noveltdet/noveltinfile*100}' $tr 
# transcript 49448	49442	99.9879	29715	60.0934	38824	36731	94.609

# Compute the stats for the gene elements
#########################################
awk -v fileRef=ref_gn_id.txt 'BEGIN{OFS="\t"; while (getline < fileRef >0){refelt++; inref[$1]=1}} NR>=2{ok=0; k=2; while(ok<=1&&k<=NF){if($k>=0.1){ok++} k++} if(inref[$1]==1){refeltinfile++; if(ok==2){refeltdet++}} else{noveltinfile++; if(ok==2){noveltdet++}}} END{print "gene", refelt, refeltinfile, refeltinfile/refelt*100, refeltdet, refeltdet/refelt*100, noveltinfile, noveltdet, noveltdet/noveltinfile*100}' $gn
# gene 25880	8492	32.813	1905	7.3609	18067	15414	85.3158 
