#!/bin/bash


# compute_gene_rpkm_from_simple_counts.sh
# for pe stranded rnaseq data and gencode like annotation file
# Takes as input the 4 compulsory arguments:
############################################
# - an exon gff file associated to gene id
# - a bam file of mapped reads from pe stranded rnaseq data
# - a parameter for the orientation of the paired end 
# - an optional library id for the rnaseq experiment (by default taken from the bam file basename with no extension)
# Provides as output in the working directory:
##############################################
# - a three column file with the gene id, the number of reads and the RPKM of the gene based on simple count of mappings totally 
#   and strandedly included in the projected exons of the gene. This file is called $b12\_genes_readcount_RPKM_simple_counts_$lid.txt

# modified on nov 7th 2015 to divide by number of reads in exons rather than total number of reads in the bam file (samtools view -c $bam)

# usage:
########
# /users/rg/sdjebali/bin/compute_gene_rpkm_from_simple_counts.sh annot.gff mappings.bam mate_strand [exp_lib_id] 2> compute_gene_rpkm_from_simple_counts.err

# !!! be careful: assumes that gene id is in field no 10 and transcript id is in field no 12 !!!
# !!! be careful: should not be run in parallel in the same directory !!!


# Test on input files
######################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage: compute_gene_rpkm_from_simple_counts.sh annot.gff mappings.bam mate_strand [exp_lib_id] >&2
    echo "       will provide a 3 column file \in the working directory that has gene id, gene mapping count, gene RPKM" >&2
    echo BE CAREFUL: for exon row \in the gff file I am expecting gene id \in field no 10 and tr id \in field no 12 >&2
    echo BE CAREFUL: should not be run \in parallel \in the same directory >&2
    echo "" >&2
    exit 0
fi

# Set variables
###############
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
b1=`basename $annot`
b12tmp=${b1%.gff}
b12=${b12tmp%.gtf}
mappings=$2
b2=`basename $mappings`
b22=${b2%.bam}
if [ "$3" != "MATE1_SENSE" ] &&  [ "$3" != "MATE2_SENSE" ] && [ "$3" != "MATE_STRAND_CSHL" ]
then 
    echo "" >&2
    echo Usage: compute_gene_rpkm_from_simple_counts.sh annot.gff mappings.bam mate_strand [exp_lib_id] >&2
    echo "       will provide a 3 column file \in the working directory that has gene id, gene mapping count, gene RPKM" >&2
    echo Wrong value for mate_strand: it can only be MATE1_SENSE, MATE2_SENSE, or MATE_STRAND_CSHL not anything else >&2
    echo "" >&2
    exit 0 
else
    mate_strand=$3
fi
if [ -n "$4" ]
then
    lid=$4
else
    lid=$b22
fi


# Programs
##########
CUTGFF=$rootDir/cutgff.awk 
GFF2GFF=$rootDir/gff2gff.awk
MAKESP=$rootDir/../bin/makeSP
COUNT=$rootDir/add_sense_antisense_read_counts_to_segments_frombam.sh

# Call the projected exons of each gene 
#######################################
echo I am making the projected exons of each gene >&2
echo "    first make a proper exon file for makeSP..." >&2
awk '$3=="exon"&&($7=="+"||$7=="-"){split($0,a,"\t"); split(a[9],b,"; "); k=1; while(b[k]!=""){split(b[k],c," "); if(c[1]=="gene_id"){gn=c[2]} if(c[1]=="transcript_id"){tr=c[2]} k++;} $10=tr";"; $12=gn";"; $9="transcript_id"; $11="gene_id"; print}' $annot | awk -v to=12 -f $CUTGFF | awk -f $GFF2GFF > $b12.formakeSP.gff
echo "    second apply makeSP..." >&2
$MAKESP $b12.formakeSP.gff -f gff
echo "    third make the projected exon file..." >&2
awk '{split($10,a,":"); print a[1], "makeSP", "projex", a[2], a[3], ".", (a[4]=="p" ? "+" : "-"), ".", "gene_id", "\""$12"\"\;";}' $b12.formakeSP.gff_segproj.gff | sort | uniq | awk -f $GFF2GFF > $b12.exonproj.gff
echo done >&2

# Compute the cumulative projected exon length for each gene
############################################################
echo I am computing the cumulative projected exon length for each gene >&2
awk '{split($10,a,"\""); cumullg[a[2]]+=($5-$4+1)}END{for(g in cumullg){print g, cumullg[g]}}' $b12.exonproj.gff > $b12.gn.with.projexlg.txt
echo done >&2

# Compute the RPKM of each gene based on simple count of mappings totally and strandedly included in the projected exons of the gene
####################################################################################################################################
echo I am computing the actual RPKM of each gene based on counts of mappings totally and strandedly included in the projected exons of the gene >&2
echo "    first count the reads" >&2
$COUNT $b12.exonproj.gff $mappings $mate_strand 0
echo "    second compute RPKM" >&2
awk -v fileRef=$b12.gn.with.projexlg.txt 'BEGIN{while(getline < fileRef >0){lg[$1]=$2}} {split($10,a,"\""); split($12,b,"\""); nbreads[a[2]]+=(b[2]); tot+=(b[2])}END{for(g in lg){print g, nbreads[g], ((nbreads[g])*10^9)/(lg[g]*tot)}}' $b12.exonproj.ext0.withreadcount.of.$b22.gff > $b12\_genes_readcount_RPKM_simple_counts_$lid.txt
echo done >&2

# Remove intermediate files
###########################
echo I am cleaning >&2
rm $b12.formakeSP.gff $b12.formakeSP.gff_segproj.gff $b12.exonproj.gff $b12.gn.with.projexlg.txt
rm $b12.exonproj.ext0.withreadcount.of.$b22.gff
echo done >&2
