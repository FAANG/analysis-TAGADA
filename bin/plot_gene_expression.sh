#!/bin/bash

# plot_gene_expression.sh
#########################
# takes as input:
#################
# - a gene TPM expression tsv file output by the gene-switch rnaseq pipeline (as many columns in header and in body)
# - a gene read count expression tsv file output by the gene-switch rnaseq pipeline (as many columns in header and in body)
# - a metadata tsv file with labExpId as a minimum
# - a factor for colouring the histograms (the palettes are available in palettes)
# - an output directory where the plots will be put
# outputs:
##########
# - 5 plots in pdf format in the output directory
#   * histogram of gene TPM with zeros for each sample in facet wrap
#   * histogram of gene TPM without zeros for each sample in facet wrap
#   * histogram of gene read count with zeros for each sample in facet wrap
#   * histogram of gene read count without zeros for each sample in facet wrap
#   * fraction of total TPM captured by N most expressed gene for all samples

# example
#########
# module load system/R-3.6.2
# tpm=/work/project/fragencode/workspace/geneswitch/results/counts/genes_TPM.tsv
# count=/work/project/fragencode/workspace/geneswitch/results/counts/genes_count.tsv
# meta=/work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/rnaseq.fastq.R1.metadata.tsv
# outdir=/work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/gnexpr/ref
# time plot_gene_expression.sh $tpm $count $meta labExpId $outdir 
# real	0m32.161s


# inputs
########
# gene TPM file
###############
# gene_id	TPM_rnaseq.sus_scrofa.cd8.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig4.R1	TPM_rnaseq.sus_scrofa.cd4.pig3.R1	TPM_rnaseq.sus_scrofa.cd4.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig2.R1	TPM_rnaseq.sus_scrofa.liver.pig3.R1	TPM_rnaseq.sus_scrofa.cd8.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig1	TPM_rnaseq.sus_scrofa.liver.pig4.R1	TPM_rnaseq.sus_scrofa.liver.pig2.R1	TPM_rnaseq.sus_scrofa.cd8.pig1.R1
# ENSSSCG00000000026	0.0	0.0	0.0	0.0	0.0	0.0	0.019718	0.0	0.0	0.0	0.0
# 26560 (12 fields)
# gene read count file
######################
# gene_id	cov_rnaseq.sus_scrofa.cd8.pig3.R1	cov_rnaseq.sus_scrofa.cd4.pig4.R1	cov_rnaseq.sus_scrofa.cd4.pig3.R1	cov_rnaseq.sus_scrofa.cd4.pig2.R1	cov_rnaseq.sus_scrofa.cd8.pig2.R1	cov_rnaseq.sus_scrofa.liver.pig3.R1	cov_rnaseq.sus_scrofa.cd8.pig4.R1	cov_rnaseq.sus_scrofa.liver.pig1	cov_rnaseq.sus_scrofa.liver.pig4.R1	cov_rnaseq.sus_scrofa.liver.pig2.R1	cov_rnaseq.sus_scrofa.cd8.pig1.R1
# ENSSSCG00000000026	0.0	0.0	0.0	0.0	0.0	0.0	0.228127	0.0	0.0	0.0	0.0
# 26560 (12 fields)
# metadata file
###############
# gene_id	cov_rnaseq.sus_scrofa.cd8.pig3.R1	cov_rnaseq.sus_scrofa.cd4.pig4.R1	cov_rnaseq.sus_scrofa.cd4.pig3.R1	cov_rnaseq.sus_scrofa.cd4.pig2.R1	cov_rnaseq.sus_scrofa.cd8.pig2.R1	cov_rnaseq.sus_scrofa.liver.pig3.R1	cov_rnaseq.sus_scrofa.cd8.pig4.R1	cov_rnaseq.sus_scrofa.liver.pig1	cov_rnaseq.sus_scrofa.liver.pig4.R1	cov_rnaseq.sus_scrofa.liver.pig2.R1	cov_rnaseq.sus_scrofa.cd8.pig1.R1
# ENSSSCG00000000026	0.0	0.0	0.0	0.0	0.0	0.0	0.228127	0.0	0.0	0.0	0.0
# 26560 (12 fields)

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] 
then
    echo "" >&2
    echo Usage: plot_gene_expression.sh genes_TPM.tsv genes_count.tsv metadata.tsv colorfactor outdir >&2
    echo "" >&2
    echo "where colorfactor is a factor from metadata.tsv which is used for colouring the plot. You can use labExpId" >&2
    echo "Needs R to be installed (tested with version 3.6.2)" >&2
    exit 1
fi

# Set general variables
#######################
path="`dirname \"$0\"`" # relative path
rootdir="`( cd \"$path\" && pwd )`" # absolute path
tpmfile=$1
covfile=$2
meta=$3
factor=$4
outdir=$5

# Set variables of programs and of palettes
###########################################
PLOT_EXPR=$rootdir/rpkm_distribution.R
PLOT_CUMUL=$rootdir/rpkm_fraction.R
palette=$rootdir/palettes/cbbPalette.8.txt

# 0. Go to output directory to produce the intermediate text files and the plots
#################################################################################
cd $outdir

# 1. make an ok TPM matrix file (with header with 1 column less than body)
##########################################################################
awk 'NR==1{OFS="\t"; gsub(/TPM_/,"",$0); gsub(/.R1/,"",$0); for(i=2; i<=NF-1; i++){s=(s)($i)("\t")} print (s)($i)} NR>=2{print}' $tpmfile > genes_TPM.tsv
# rnaseq.sus_scrofa.cd8.pig3	rnaseq.sus_scrofa.cd4.pig4	rnaseq.sus_scrofa.cd4.pig3	rnaseq.sus_scrofa.cd4.pig2	rnaseq.sus_scrofa.cd8.pig2	rnaseq.sus_scrofa.liver.pig3	rnaseq.sus_scrofa.cd8.pig4	rnaseq.sus_scrofa.liver.pig1	rnaseq.sus_scrofa.liver.pig4	rnaseq.sus_scrofa.liver.pig2	rnaseq.sus_scrofa.cd8.pig1
# ENSSSCG00000000026	0.0	0.0	0.0	0.0	0.0	0.0	0.019718	0.0	0.0	0.0	0.0
# 1 (11 fields)
# 26559 (12 fields)

# 2. make an ok read count matrix file
######################################
awk 'NR==1{OFS="\t"; gsub(/cov_/,"",$0); gsub(/.R1/,"",$0); for(i=2; i<=NF-1; i++){s=(s)($i)("\t")} print (s)($i)} NR>=2{print}' $covfile > genes_readcount.tsv
# rnaseq.sus_scrofa.cd8.pig3	rnaseq.sus_scrofa.cd4.pig4	rnaseq.sus_scrofa.cd4.pig3	rnaseq.sus_scrofa.cd4.pig2	rnaseq.sus_scrofa.cd8.pig2	rnaseq.sus_scrofa.liver.pig3	rnaseq.sus_scrofa.cd8.pig4	rnaseq.sus_scrofa.liver.pig1	rnaseq.sus_scrofa.liver.pig4	rnaseq.sus_scrofa.liver.pig2	rnaseq.sus_scrofa.cd8.pig1
# ENSSSCG00000000026	0.0	0.0	0.0	0.0	0.0	0.0	0.228127	0.0	0.0	0.0	0.0
# 1 (11 fields)
# 26559 (12 fields)

# 3. run the first expression script for TPM with zeros
########################################################
$PLOT_EXPR -i genes_TPM.tsv -r histogram -w -l -v "log10(TPM)" -m $meta -P $palette -f $factor --verbose -o genes_TPM

# 4. run the first expression script for TPM without zeros
###########################################################
$PLOT_EXPR -i genes_TPM.tsv -r histogram -w -l -p 0 -v "log10(TPM+0)" -m $meta -P $palette -f $factor --verbose -o genes_TPM

# 5. run the first expression script for read counts with zeros
################################################################
$PLOT_EXPR -i genes_readcount.tsv -r histogram -w -l -v "log10(readcount)" -m $meta -P $palette -f $factor --verbose -o genes_readcount

# 6. run the first expression script for read counts without zeros
####################################################3#############
$PLOT_EXPR -i genes_readcount.tsv -r histogram -w -l -p 0 -v "log10(readcount+0)" -m $meta -P $palette -f $factor --verbose -o genes_readcount

# 7. run the second expression script on the TPM file  
######################################################
$PLOT_CUMUL -i genes_TPM.tsv -m $meta -c $factor -o genes_totalTPM_captured_by_top_genes
