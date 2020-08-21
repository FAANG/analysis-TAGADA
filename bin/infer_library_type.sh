#!/bin/bash

# Description
##############
# Takes as input a bam file with a set of mapped reads, a reference gene annotation and infers the sequencing library protocol (Unstranded, Mate2_sense & Mate1_sense) used to generate the rna-seq data. It does it by comparing the mapping strand in 1% of the aligments with the strand of the gene the read maps. Finally it produces three numbers: 

# 1) Fraction of reads explained by "1++,1--,2+-,2-+"
# 2) Fraction of reads explained by "1+-,1-+,2++,2--"
# 3) Fraction of reads explained by other combinations

# They give information regarding the library. They contain several strings of three characters, i.e. 1+-, where:
#   Character 1. 1 and 2 are mate1 and mate2 respectively.
#   Character 2. + and - is the strand where the read maps.
#   Character 3. + and - is the strand where the gene in which the read overlaps is annotated.

# You can apply the following rules to infer the used library from this information:

#    NONE. Not strand-specific protocol (unstranded data). Fraction of reads explained by “1++,1–,2+-,2-+” and “1+-,1-+,2++,2–” close to 0.5000 in both cases.

# Strand-specific protocols (stranded data):
#    MATE1_SENSE. Fraction of reads explained by “1++,1–,2+-,2-+” close to 1.0000.
#    MATE2_SENSE. Fraction of reads explained by “1+-,1-+,2++,2–” close to 1.0000.


# usage
#######
# infer_library_type.sh alignments.bam annotation.gff

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts
# - uses bedtools
# - comes from chimpipe but changed a bit on July 3rd 2020 to call awk script at the same level as the bash script
# - could be made faster by parallelizing the samtools part

# will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
echo "" >&2
echo "infer_library_type.sh"
echo "" >&2
echo "Takes as input a bam file with a set of mapped reads, a reference gene annotation and infers the sequencing library protocol (Unstranded, Mate2_sense & Mate1_sense) used to generate the rna-seq data. It does it by comparing the mapping strand in 1% of the aligments with the strand of the gene the read maps."
echo "" >&2
echo Usage:  infer_library_type.sh alignments.bam annotation.gff >&2
echo "" >&2
echo "" >&2
exit 1
fi

# GETTING INPUT ARGUMENTS
#########################
bamfile=$1
annot=$2

# Directories 
#############
## Set root directory
path="`dirname \"$0\"`"              # relative path
rootdir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootdir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi

# PROGRAMS
##########
cutgff=$rootdir/cutgff.awk

# START
########
## Comment: samtools view -F 260 ** filter out unmapped reads (4) + secondary alignments (254) = 260
## Note about random sampling with samtools view -s:
# Current versions of samtools (ab)use the integer part of the -s value to set the seed:
# -s FLOAT    Integer part is used to seed the random number generator. Part after the decimal point sets the fraction of templates/pairs to subsample
# So two runs with e.g. view -s 123.1 / view -s 456.1 should select two distinct randomly-selected 10% subsets.
# I will take a random number between 0 and 10 as a seed
randomSeed=`shuf -i1-10 -n1`

awk -v elt='exon' '$3==elt' $annot | awk -v to=8 -f $cutgff | sort -k1,1 -k4,4 -k5,5 | uniq | bedtools intersect -abam <(samtools view -b -s ${randomSeed}.01 -F 260 $bamfile) -b stdin -split -bed -wo | awk '{print $4, $6, $19;}' | uniq | awk '{split($1,a,"/"); readCount["total"]++; readCount[a[2]":"$2":"$3]++;}END{fraction1=(readCount["1:+:+"]+readCount["1:-:-"]+readCount["2:+:-"]+readCount["2:-:+"]); fraction2=(readCount["1:+:-"]+readCount["1:-:+"]+readCount["2:+:+"]+readCount["2:-:-"]); other=(readCount["total"]-(fraction1+fraction2)); print (fraction1/readCount["total"]*100), (fraction2/readCount["total"]*100), other;}' 



