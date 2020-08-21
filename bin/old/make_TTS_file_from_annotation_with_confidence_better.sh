#!/bin/bash

# make_TTS_file_from_annotation_with_confidence_better.sh
# improved on 06/05/2012 replacing i++ by i+=2 when looking for tr biotype in step1
# and more importantly replacing the assumption of gene_type being in 14th field
# in gene rows by looking for the gene_type key
# same as make_TSS_file_from_annotation_with_confidence.sh except but for tts (termination site of tr)

# improved on dec 20th 2013 in order to be able to specify the gtf file on the command line

# Be careful: this script cannot be launched several times in the same directory
# because it is writing files not indexed by the input file like tmp and exp_fld1_fld2.txt for example
# This script supposes that the annotation file contains gene name of exon features in $10 and tr name in $12

# - uses awk scripts

# Usage:
########
# make_TTS_file_from_annotation_with_confidence_better.sh annot.gff [tr_bt_list] 2> make_TTS_file_with_confidence_better_from_annotation.err &

# Check it has all it needs
###########################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_TTS_file_from_annotation_with_confidence_better.sh annot.gtf [tr_biotypes.txt] >&2
    echo Be careful: it is not possible to run two instances at the same time \in the same directory since produces intermediate files not indexed >&2
    echo Be careful: the gene id and the transcript id of exon features have to be \in column no 10 and 12 respectively >&2
    echo The transcript biotypes must be surrounded by double quotes and end with semi colons as \in the gtf file >&2
    echo In case no list of transcript biotypes is specified it will take all of them >&2
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
###########
EXTRACT3p=$rootDir/../Awk/extract_most_3p.awk
CUTGFF==$rootDir/../Awk/cutgff.awk
GFF2GFF==$rootDir/../Awk/gff2gff.awk


##########################################################
# Make the TTS file for the asked transcript biotypes    #
##########################################################

# a. Extract most 3' exons of transcripts from specified transcript biotypes
#############################################################################
echo I am extracting the most 3\' exons of transcripts from the biotypes specified by the user \(default is no selection\) >&2
if [ ! -n "$trbiotypes" ] || [ ! -e "$trbiotypes" ]
then
awk '$3=="exon"' $annotation | awk -v fldno=12 -f $EXTRACT3p > $annotbase\_exons_most3p.gff
else
awk -v fileRef=$trbiotypes 'BEGIN{while(getline < fileRef >0){biotype[$1]=1;}} $3=="exon"{i=9; while($i!="transcript_type"){i+=2}if(($i=="transcript_type")&&(biotype[$(i+1)]==1)){print}}' $annotation | awk -v fldno=12 -f $EXTRACT3p > $annotbase\_exons_most3p.gff
fi


# b. Then the most 3' bp of each transcript for each gene (not that all tts are said to come from gencode)
##########################################################################################################
#    with associated confidence level= low confidence level whenever the tts belongs to a tr where the CDS end
#################################################################################################################
#    not found tag was set, not_low otherwise (syntax was "CDS end not found" for v3c)
###################################################################################
echo I am extracting the most 3\' bp of each transcript for each gene, >&2
echo associating a low confidence level when the tts comes from a tr where the CDS end not found tag was set, >&2
echo and adding the list of tr and of tr biotypes the tts comes from. >&2
awk '{i=9; while($i!="transcript_type"){i++}if($i=="transcript_type"){split($(i+1),c,"\""); trbiot=c[2];} ($0~/cds_end_NF/) ? confidence="low" : confidence="not_low"; ($7=="+") ? ttspos=$5 : ttspos=$4; split($10,a,"\""); split($12,b,"\""); print $1, "Gencode", "TTS", ttspos, ttspos, ".", $7, ".", "gene_id", a[2], "tr", b[2], "trbiot", trbiot, "confidence", confidence;}' $annotbase\_exons_most3p.gff | awk -f $GFF2GFF > $annotbase\_tts_sites.gff

# c. Finally collapse per gene and put a low confidence level when the collapsed tts has at least low inside it
###############################################################################################################
#    also add the gene biotype at the end
##########################################
echo I am collapsing all TTSs per gene, put a low confidence level whenever one of the indiv tts has a low confidence level, >&2
echo and add the gene biotype. >&2
cat $annotbase\_tts_sites.gff | awk -v to=10 -f $CUTGFF | sort -n | uniq -c | awk '{$1=""; print $0}' | awk -f $GFF2GFF | awk -v fileRef=$annotbase\_tts_sites.gff 'BEGIN{while (getline < fileRef >0){trlist[$1"_"$4"_"$5"_"$7,$10]=(trlist[$1"_"$4"_"$5"_"$7,$10])($12)(","); trbiotlist[$1"_"$4"_"$5"_"$7,$10]=(trbiotlist[$1"_"$4"_"$5"_"$7,$10])($14)(","); if($16=="low"){low[$1"_"$4"_"$5"_"$7,$10]=1;}}} {$11="trlist"; $12=trlist[$1"_"$4"_"$5"_"$7,$10]; $13="trbiotlist"; $14=trbiotlist[$1"_"$4"_"$5"_"$7,$10]; $15="confidence"; (low[$1"_"$4"_"$5"_"$7,$10]==1) ? $16="low" : $16="not_low"; print $0}' | awk -v fileRef=$annotation 'BEGIN{while(getline < fileRef >0){if($3=="gene"){split($10,a,"\""); i=9; while($i!="gene_type"){i+=2} if($i=="gene_type"){split($(i+1),b,"\"");} biotype[a[2]]=b[2]}}}{print $0, "gene_biotype", biotype[$10]}' | awk -f $GFF2GFF > $annotbase\_tts_sites_nr_with_confidence.gff


