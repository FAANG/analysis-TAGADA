#!/bin/bash

#make_TSS_file_from_annotation_with_confidence_better.sh
# this script takes as input:
# - a Gencode annotation file in GTF format (mandatory)
# - a list of transcript biotypes to consider in this file (optional, default: no filtering)
# and provides as output the following files in the directory where it is executed:
# - a file of most 5' exons from each transcript in GFF2 format
# - a file of TSS from each transcript in GFF2 format
# - a file of distinct TSS from each gene in GFF2 format, associated to a low confidence level in case the cds_start_NF 
#   tag was present in any transcript with this TSS

# Be careful:
#############
# makes the assumption that gene id and transcript id come as first and second (tag,value) pairs in the annotation GTF file

# Usage:
########
# make_TSS_file_from_annotation_with_confidence_better.sh annot.gtf [tr_biotypes.txt] 2> make_TSS_file_with_confidence_better_from_annotation.err &


# Check that the annotation file exists, otherwise exit
#######################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_TSS_file_from_annotation_with_confidence_better.sh annot.gtf [tr_biotypes.txt] >&2
    echo Be careful: this script makes the assumption that gene id and transcript id come as first >&2
    echo and second \(tag\,value\) pairs \in the GTF file >&2
    echo In case no list of transcript biotypes is specified there will be no filtering done >&2
    echo "" >&2
    exit 1
fi

# Initialize variables
######################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annotation=$1
annotbase=`basename ${annotation%.gtf}`
if [ -n "$2" ] 
then
    trbiotypes=$2
fi
    
# Programs
##########
EXTRACT5p=$rootDir/../Awk/extract_most_5p.awk
CUTGFF=$rootDir/../Awk/cutgff.awk
GFF2GFF=$rootDir/../Awk/gff2gff.awk


##########################################################
# Make the TSS file for the asked transcript biotypes    #
##########################################################

# a. Extract most 5' exons of transcripts from specified transcript biotypes
#############################################################################
echo I am extracting the most 5\' exons of transcripts from the biotypes specified by the user \(default is no selection\) >&2
if [ ! -n "$trbiotypes" ] || [ ! -e "$trbiotypes" ]
then
awk '$3=="exon"' $annotation | awk -v fldno=12 -f $EXTRACT5p > $annotbase\_exons_most5p.gff
else
awk -v fileRef=$trbiotypes 'BEGIN{while(getline < fileRef >0){biotype["\""$1"\"\;"]=1;}} $3=="exon"{i=9; while($i!="transcript_type"){i+=2}if(($i=="transcript_type")&&(biotype[$(i+1)]==1)){print}}' $annotation | awk -v fldno=12 -f $EXTRACT5p > $annotbase\_exons_most5p.gff
fi

# b. Then the most 5' bp of each transcript for each gene (not that all tss are said to come from gencode)
##########################################################################################################
#    with associated confidence level= low confidence level whenever the tss belongs to a tr where the CDS start
#################################################################################################################
#    not found tag was set, not_low otherwise (syntax was CDS start not found for v3c)
######################################################################################
echo I am extracting the most 5\' bp of each transcript for each gene, >&2
echo associating a low confidence level when the tss comes from a tr where the cds_start_NF tag was set, >&2
echo and adding the list of tr and of tr biotypes the tss comes from. >&2
awk '{i=9; while($i!="transcript_type"){i++}if($i=="transcript_type"){split($(i+1),c,"\""); trbiot=c[2];} ($0~/cds_start_NF/) ? confidence="low" : confidence="not_low"; ($7=="+") ? tsspos=$4 : tsspos=$5; split($10,a,"\""); split($12,b,"\""); print $1, "Gencode", "CapSite", tsspos, tsspos, ".", $7, ".", "gene_id", a[2], "tr", b[2], "trbiot", trbiot, "confidence", confidence;}' $annotbase\_exons_most5p.gff | awk -f $GFF2GFF > $annotbase\_capped_sites.gff

# c. Finally collapse TSS per gene and put a low confidence level when the collapsed tss is composed of at least one low tss
#############################################################################################################################
#    also add the gene biotype at the end
##########################################
echo I am collapsing all TSSs per gene, put a low confidence level whenever one of the indiv tss has a low confidence level, >&2
echo and add the gene biotype. >&2
cat $annotbase\_capped_sites.gff | awk -v to=10 -f $CUTGFF | sort -n | uniq -c | awk '{$1=""; print $0}' | awk -f $GFF2GFF | awk -v fileRef=$annotbase\_capped_sites.gff 'BEGIN{while (getline < fileRef >0){trlist[$1"_"$4"_"$5"_"$7,$10]=(trlist[$1"_"$4"_"$5"_"$7,$10])($12)(","); trbiotlist[$1"_"$4"_"$5"_"$7,$10]=(trbiotlist[$1"_"$4"_"$5"_"$7,$10])($14)(","); if($16=="low"){low[$1"_"$4"_"$5"_"$7,$10]=1;}}} {$11="trlist"; $12=trlist[$1"_"$4"_"$5"_"$7,$10]; $13="trbiotlist"; $14=trbiotlist[$1"_"$4"_"$5"_"$7,$10]; $15="confidence"; (low[$1"_"$4"_"$5"_"$7,$10]==1) ? $16="low" : $16="not_low"; print $0}' | awk -v fileRef=$annotation 'BEGIN{while(getline < fileRef >0){if($3=="gene"){split($10,a,"\""); i=9; while($i!="gene_type"){i+=2} if($i=="gene_type"){split($(i+1),b,"\"");} biotype[a[2]]=b[2]}}}{print $0, "gene_biotype", biotype[$10]}' | awk -f $GFF2GFF > $annotbase\_capped_sites_nr_with_confidence.gff


