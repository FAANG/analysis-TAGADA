#!/bin/bash

# Shell script
###############
# ~sdjebali/bin/extract_coord_from_maf_totally_included.sh
# - takes as input a gff file of segments in a given species
# - provides as output a set of maf files in the MafPrecise directory 
#   that it is creating where the script is launched, and that correspond
#   to the multiple alignments of those segments with other species.
#   Note for the ma of a segment to be reported there must exist a MA
#   that totally includes the segment in the maf file.
# NOTES: 
# - two instances of this script should not be run at the same time
#   on the same directory
# - even if a wanted segment is on the - strand it will be provided the MA regions
#   of the + strand including it, since in a MA the ref species is always on the +.

# Usage
#######
# extract_coord_from_maf_totally_included.sh coord.gff 

# Example of usage
###################
# extract_coord_from_maf_totally_included.sh data/exons.v4.hg19.10fst.gtf

# Makes several assumptions:
############################
# On the system
###############
# - Made for using on a 64 bit linux architecture
# - the script needs to be able to see the home directories, in particular ~sdjebali
# - the script needs to be able to see the /users/rg/projects/ directory
# - the user needs to have write access in the directory where the script is launched
# - the script uses awk scripts
# On the data
#############
# - MAF files are located in $MAFDIR (below) and named chr$i.maf
# - the human regions to which these maf file represent are located in $dirmareg (below)
#   and named chr$i_human_regions_hg19_with_MA_46species.gff
#   (note that they could be deduced from the MAF files but this is to gain time)
# - the assembly should be provided below in the assembly variable



# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage:    extract_coord_from_maf_totally_included.sh coord.gff" >&2
    echo "          coord.gff is a gff or a gtf file with segments for which" >&2
    echo "          the user wants to get the multiple alignments they correspond to" >&2
    echo "Example:  extract_coord_from_maf_totally_included.sh data/exons.v4.hg19.10fst.gtf" >&2
    echo "" >&2
    echo "Takes a gff or gtf file of segments and retrieves the multiple alignments they correspond to, provided that for a segment there is one multiple alignment (MA) that totally include it." >&2
    echo "" >&2
    exit 0
fi

# Set general variables
#######################
path="`dirname \"$0\"`" # relative path
rootdir="`( cd \"$path\" && pwd )`" # absolute path

# Input files
#############
# we suppose that there is no strand in the input file and that we do not want to look at the strand
input=$1
basetmp=${input%.gff}
base=${basetmp%.gtf}

# Programs
##########
GFF2GFF=$rootdir/gff2gff.awk
OVERLAP=$rootdir/overlap
MAKEMULTIPLEMAF=$rootdir/make_multiple_maf_out_of_one_maf.awk
EXTRACTCOORD=$rootdir/xtract_coord_from_maf.awk 

# Data
#######
MAFDIR=/seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf
assembly=hg19
# Directory of per chromosome gff files of coordinates of alignments of human present in maf files
dirmareg=/users/rg/projects/encode/scaling_up/whole_genome/HumanRegWhereMA46species


# Creates the MafPrecise and MaReg directories where the script is launched
############################################################################
# The first one will contain the maf files of interest and the other one will 
# contain intermediate files and will thus be removed afterwards
mkdir MaReg
if [[ ! -e MafPrecise ]]; then
mkdir MafPrecise 
fi

# 0) per chromosome gff files of coordinates of alignments of human present in maf files
########################################################################################
# $dirmareg/chr$i_human_regions_hg19_with_MA_46species.gff
# chr1  .       .       10918   11396   .       +       .       . .
# 2723537 (10 fields)


# Actions
##########
# 1) Make gff file out of coordinates of ref species corresponding to the wanted submafs  
########################################################################################
echo "I am making the list of chromosomes corresponding to the wanted segments ..." >&2
awk '$1!~/#/{print $1}' $1 | sort -n | uniq -c | awk '{print $2}' > chr_names.txt 
echo "I am making the list of non redundant wanted segments ..." >&2
cat chr_names.txt | while read c; do awk -v searchedchr=$c '$1!~/#/&&($1==searchedchr)' $1 | sort -n | uniq -c | awk '{print $2, $3, $4, $5, $6, $7, $8, $9, ".", "."}' | awk -f $GFF2GFF > $base\_$c\_nr.gff; done


# 2) For each segment in 1) find the segment in 0) that totally includes it (there should be only one) 
######################################################################################################
echo "For each wanted segment, I am computing the MA region that totally includes it ..." >&2
cat chr_names.txt | while read c; do $OVERLAP $base\_$c\_nr.gff $dirmareg/$c\_human_regions_hg19_with_MA_46species.gff -m -1 -i 1 -f ma -o $base\_$c\_nr_withMA_includingthem.gff; done
# time consuming AND linear in the number of initial segments

# 3) For each segment in 1) with a segment in 0) that totally includes it,
##########################################################################
# extract the submaf region corresponding to it
###############################################
# a) Make the file of unique submafs (we have exact beg of big maf file) in bed format
######################################################################################
# Note that here we have only one submaf in the list since total inclusion (not the case when overlap)
#######################################################################################################
echo "I am making the file of distinct MA regions that totally include the wanted segments ..." >&2
cat chr_names.txt | while read c; do awk '$12>=1{split($14,a,","); k=1; while(a[k]!=""){split(a[k],b,"_"); print b[1], b[2]-1, b[3], ".", ".", b[4]; k++;}}' $base\_$c\_nr_withMA_includingthem.gff | sort -n | uniq -c | awk 'BEGIN{OFS="\t";}{print $2, $3, $4, $5, $6, $7}' > $base\_$c\_nr_withMA_includingthem_onlymareg_includingthem.bed; done
# b) then actually extract the submafs corresponding to them reading the maf file once
#######################################################################################
echo "I am making the maf files corresponding to those distinct MA regions that totally include the wanted segments ..." >&2
cat chr_names.txt | while read c; do awk -v refspecies=$assembly.$c -v searched_str=+ -v fileRef=$base\_$c\_nr_withMA_includingthem_onlymareg_includingthem.bed -f $MAKEMULTIPLEMAF $MAFDIR/$c.maf; done
# very time consuming BUT not linear in the number of initial segments
# c) finally run the script extract_coord_from_maf.awk on those intermediate files 
###################################################################################
# to get the _precise.maf file corresponding to each wanted segment
###################################################################
echo "I am actually extracting the wanted segments from those MA regions ..." >&2
cat chr_names.txt | while read c; do awk '$12>=1{split($14,a,","); split(a[1],b,"_"); print $1"_"($4-1)"_"$5"_"$7, b[1]"_"(b[2]-1)"_"b[3]"_"b[4]}' $base\_$c\_nr_withMA_includingthem.gff | sort -n | uniq -c | awk '{print $2, $3}' | while read regseg regma; do awk -v assembly=$assembly -v coord=$regseg -f $EXTRACTCOORD MaReg/$regma\_precise.maf > MafPrecise/$regseg\_precise.maf; done; done
# linear in the number of initial segments
# d) and remove at the end the _precise.maf file in MaReg/ representing the MA regions
#######################################################################################
echo "I am removing the maf files corresponding to the MA regions including the wanted segments as well as some additional intermediate files ..." >&2
cat chr_names.txt | while read c; do awk '$12>=1{split($14,a,","); split(a[1],b,"_"); print b[1]"_"(b[2]-1)"_"b[3]"_"b[4]}' $base\_$c\_nr_withMA_includingthem.gff | sort -n | uniq -c | while read i regma; do rm MaReg/$regma\_precise.maf; done; done
rmdir MaReg
rm chr_names.txt
