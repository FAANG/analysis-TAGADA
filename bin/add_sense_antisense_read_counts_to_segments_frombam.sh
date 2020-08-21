#!/bin/bash

# add_sense_antisense_read_counts_to_segments_frombam.sh
# Works for stranded pe rnaseq data.
# Counts the number of mappings falling in some segments extended on each side,
# in the sense and the antisense orientation. Designed for paired end RNAseq mappings 
# with the two mates pointing towards each other, or in the same orientation as the tr, 
# and was mainly used for gene elements. 

# on 12/03/2013 I added the -f 1 option to only have features that are totally included
# in our features and removed it on April28th 2020, and on 03/17/2014 I added as an argument the output dir and indexed
# the annotation file by the bam file name in order to avoid conflicts.

# This script takes as input:
#############################
# - a gff file of element with the element id in field no 10, 
# - a bam file of paired end mappings, 
# - the convention for the mate strand with respect to the transcript, which can be either MATE1_SENSE (mate 1
#   on the same strand as the transcript, mate 2 pointing towards mate 1), MATE2_SENSE (mate 2 on the same strand
#   as the transcript, mate 1 pointing towards mate 2), or MATE_SENSE_CSHL (two mates on same strand as the transcript),
# - a compulsory parameter saying how many nucleotides to extend each object on each side for read overlap (default 0).
# - an optional parameter for the output directory (default directory where it is run).
# This script provides as output:
#################################
# - the same gff file but extended by a number given in parameter no 4 and with two additional (key,value) pairs
#   representing the number of sense and antisense reads falling in this region extended by argument 4 on each side
# Improvements would be:
########################
# 1. put a quantity that is number of read weighted by the prop of the total length of the read that 
#    falls there (but then the quantity will not be an integer).
# 2. weight the multimaps to reflect the uncertainty of mapping (but again the the quantity will not 
#    be an integer).
# 3. have an awk script that can work for all possible configurations of the reads, in order to have a single command here
# 4. make as an option to have -f 1 or not in the intersectBed

# usage
########
# add_sense_antisense_read_counts_to_segments_frombam.sh annot.gff mappings.bam mate_strand nb_bp_to_extend_annot [outputdir]

# Test on input files
######################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ]
then
    echo "" >&2
    echo Usage: add_sense_antisense_read_counts_to_segments_frombam.sh annot.gff mappings.bam mate_strand nb_bp_to_extend_annot [outputdir] >&2
    echo "      - the id of the element has to be in column no 10 of annot.gff" >&2
    echo "      - mate_strand should be one of the following values: MATE1_SENSE, MATE2_SENSE or MATE_SENSE_CSHL (for old star v1 cshl files)" >&2
    echo "!!! Note: requires bedtools !!!"  >&2
    echo "" >&2
    exit -1
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
    echo Usage: add_sense_antisense_read_counts_to_segments_frombam.sh annot.gff mappings.bam mate_strand nb_bp_to_extend_annot [outputdir] >&2
    echo Wrong value for mate_strand: it can only be MATE1_SENSE, MATE2_SENSE, or MATE_STRAND_CSHL, and not anything else >&2
    echo the id of the element has to be in column no 10 of annot.gff >&2
    echo "!!! Note: requires bedtools !!!"  >&2
    echo "" >&2
    exit 1
else
    mate_strand=$3
fi
toextend=$4
if [ -n "$5" ]
then
    outdir=$5
else
    outdir=.
fi

# Will exit if there is an error or in a pipe
#############################################
set -e -o pipefail

# Programs
##########
GFF2GFF=$rootDir/gff2gff.awk
INTER2GFF=$rootDir/intersectBed_pestrandedbam_with_eltgff_to_gff_with_s_as_mappings.awk

# Extend the annotation by the number of bp required
####################################################
echo I am extending the annotation by the number of bp required as 4th argument >&2
awk -v toadd=$toextend '{$4=$4-toadd; $5=$5+toadd; print $0}' $annot | awk -f $GFF2GFF > $outdir/$b12.ext$toextend.$b22.gff 
echo done >&2

# Compute the overlap between this annotation and the reads and divide the overlapping reads into sense and antisense
#####################################################################################################################
echo I am computing the number of reads overlapping the extended annotation and dividing the reads into sense and antisense >&2
intersectBed -abam $mappings -b $outdir/$b12.ext$toextend.$b22.gff -bed -wo | awk -v fileRef=$outdir/$b12.ext$toextend.$b22.gff -v mate_strand=$mate_strand -f $INTER2GFF | awk -f $GFF2GFF > $outdir/$b12.ext$toextend.withreadcount.of.$b22.gff
echo done >&2
# totally included \in
# -f 1

# Remove intermediate files
###########################
echo I am cleaning >&2
rm $outdir/$b12.ext$toextend.$b22.gff
echo done >&2
