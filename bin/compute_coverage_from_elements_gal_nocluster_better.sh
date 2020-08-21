#!/bin/bash

# compute_coverage_from_elements_gal_nocluster_better.sh
# same as ~/bin/compute_coverage_from_elements_gal_nocluster.sh except that the projex and projgn bed files
# this script takes as input any projection bed file (only 3 elements separated by tabs) and compare them to projexons and projgenes
# in order to provide as output 3 numbers (unstranded):
#######################################################
# - total nt covered
# - exonic nt
# - genic nt
# These three numbers are all we need to compute:
#################################################
# - proportion of total genome covered
# - proportion of detected nt that are exonic / intronic / intergenic
# - proportion of exonic / intronic / intergenic nt that are detected 
# (for boxplots suggested by Roderic)

# Note:
#######
# - needs to be able to see the home directories, in particular ~sdjebali
# - needs to be able to see the /users/rg/projects/ directory
# - uses awk scripts

# Usage:
########
# compute_coverage_from_elements_gal_nocluster_better.sh absolute_path_to_element_file.bed projected_exon.bed projected_gene.bed [absolute_path_to_output_file]

# Example
#########
# compute_coverage_from_elements_gal_nocluster_better.sh ~/ENCODE_AWG/Analyses/SummaryStatistics/Projections/rnafrac_longPolyA_cellcompart_cell_proj.bed /users/rg/projects/encode/scaling_up/whole_genome/Gencode/OfficialPartition/v10_hg19/gen10.partition.unstr.exonic.proj.bed /users/rg/projects/encode/scaling_up/whole_genome/Gencode/OfficialPartition/v10_hg19/gen10.partition.unstr.genic.proj.bed
# produces the file rnafrac_longPolyA_cellcompart_cell_proj_coverage.txt
# named after the input proj file
# 1115344347 81718300 835183850



# In case the user does not provide one of the 3 mandatory files
################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage:    compute_coverage_from_elements_gal_nocluster_better.sh absolute_path_to_element_file.bed projected_exon.bed projected_gene.bed [absolute_path_to_output_file] >&2
    echo "         absolute_path_to_element_file.bed is a bed file of elements from which we want the common intersection with the exons and the genes" >&2
    echo "         projected_exon.bed and projected_gene.bed are the unstranded projections of all exons and genes respectively" >&2
    echo "         needs bedtools" >&2
    echo "         " >&2
    echo Example:  compute_coverage_from_elements_gal_nocluster_better.sh ~/ENCODE_AWG/Analyses/SummaryStatistics/Projections/rnafrac_longPolyA_cellcompart_cell_proj.bed /users/rg/projects/encode/scaling_up/whole_genome/Gencode/OfficialPartition/v10_hg19/gen10.partition.unstr.exonic.proj.bed /users/rg/projects/encode/scaling_up/whole_genome/Gencode/OfficialPartition/v10_hg19/gen10.partition.unstr.genic.proj.bed >&2
    echo "" >&2
    exit -1 
fi

# Assigns the 3 mandatory files to variables
############################################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
infile=$1
ex=$2
gn=$3 

# Determine the basename according to whether the file is gzipped or not
########################################################################
test=`file $infile | awk '{print $2}'`
if [ $test == gzip ]
then
inbase=`basename ${infile%.bed.gz}`
else
inbase=`basename ${infile%.bed}`
fi

# in case the user specifies the location of the output file (without coverage.txt) 
###################################################################################
# we use it otherwise we put the file in current dir and name it after the input file
#####################################################################################
if [ ! -n "$4" ]
then
output=$inbase\_coverage.txt
else
output=$4\_coverage.txt
fi

# According to gzip status will make the 3 nt essential numbers
################################################################
if [ $test == gzip ]
then
for i in 1
do
echo Computing total nt covered >&2
total=`zcat $infile | awk '{s+=$3-$2}END{print s}'`
echo Computing exonic nt covered >&2
zcat $infile | coverageBed -a stdin -b $ex > projex_with_coverage_from_$inbase.txt
exonic=`awk '{s+=$(NF-2)}END{print s}' projex_with_coverage_from_$inbase.txt`
echo Computing genic nt covered >&2
zcat $infile | coverageBed -a stdin -b $gn > projgn_with_coverage_from_$inbase.txt
genic=`awk '{s+=$(NF-2)}END{print s}' projgn_with_coverage_from_$inbase.txt`
echo $total $exonic $genic 
echo Removing intermediate files >&2
rm projex_with_coverage_from_$inbase.txt projgn_with_coverage_from_$inbase.txt
done > $output

else
for i in 1
do
echo Computing total nt covered >&2
total=`awk '{s+=$3-$2}END{print s}' $infile`
echo Computing exonic nt covered >&2
coverageBed -a $infile -b $ex > projex_with_coverage_from_$inbase.txt
exonic=`awk '{s+=$(NF-2)}END{print s}' projex_with_coverage_from_$inbase.txt`
echo Computing genic nt covered >&2
coverageBed -a $infile -b $gn > projgn_with_coverage_from_$inbase.txt
genic=`awk '{s+=$(NF-2)}END{print s}' projgn_with_coverage_from_$inbase.txt`
echo $total $exonic $genic 
echo Removing intermediate files >&2
rm projex_with_coverage_from_$inbase.txt projgn_with_coverage_from_$inbase.txt
done > $output
fi
