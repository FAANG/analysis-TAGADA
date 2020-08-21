#!/bin/bash

# make_sample_equal_distrib.sh
# takes a set of elements associated to a class and a value and returns a sample of those elements where the distribution of values 
# are the same for all classes

# This script takes as input:
#############################
# - a tsv file with header and no space in any column that has at least 3 columns: 
#  * the element id
#  * the element class
#  * the element value (mandatory)
# - the class key in header (mandatory)
# - the value key in header (mandatory)
# - the number of bins used for the sampling (optional) (default = 25)
# This script provides as output:
#################################
# - a tsv file that is a subset of the input tsv file where the distribution of values is equal for all classes
# - a log file with some summary stat information

# Example of usage
##################
# cd ~sdjebali/ENCODE_AWG/Analyses/Mouse_Human/GeneExpression/ConstrainedGenes/SameExpressionUnconstrained/AvgNoZerohs
# time make_sample_equal_distrib.sh merged_orthologs_HumanMouse_1_1_genes_nozero_notissuespec_eachspec_min_max_diff_bothspecies_constrclass_matchexpr_gnname_top1000_log10avgnozeroexpr_hs_mm_hsmm.tsv constr_class log10_avgnozero_hs 2> make_sample_equal_distrib_log10avgnozerohs.err > make_sample_equal_distrib_log10avgnozerohs.out
# real	0m10.411s
# input is
# hs_gn	mm_gn	log10min_rpkm	log10max_rpkm	dynrange_both	constr_class	matched_set	gene_name	top1000	log10_avgnozero_hs	log10_avgnozero_mm	log10_avgnozero_hsmm
# 14364 (12 fields)
# output is
# hs_gn	mm_gn	log10min_rpkm	log10max_rpkm	dynrange_both	constr_class	matched_set	gene_name	top1000	log10_avgnozero_hs	log10_avgnozero_mm	log10_avgnozero_hsmmbin_no
# 10537 (13 fields)


# Check that we have the necesary arguments
###########################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage: make_sample_equal_distrib.sh elements.tsv class_key value_key [nb_bins_for_sampling] >&2
    echo where: >&2
    echo - elements.tsv is a tsv file with header and no space in any column that has at least 3 columns: the element id, class and value >&2
    echo - class_key \in header >&2
    echo - value key \in header >&2
    echo Will output a tsv file which is a sample of the input file where the distribution of the values is equal for each class >&2
    echo as well as a log file with the distribution before and after the equalization >&2
    echo "" >&2
    exit 1
fi

# And if we have the 4th optional argument we use it otherwise we use default value
####################################################################################
if [ -n "$4" ]
then
nobin=$4
else
nobin=25
fi

# Variable from input
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

input=$1
b=`basename $input`
b2=${b%.tsv}
classkey=$2
valkey=$3

# Programs
##########
STAT=$rootDir/stats.sh
SAMPLE=$rootDir/sample_lines_better.sh

# Start the computation
#######################
# Column no for class and value in input file
echo I am computing the column number for the class and the value \in the input file >&2
ck=`head -1 $input | awk -v classkey=$classkey '{n=split($0,a,"\t"); found=0; i=1; while((found==0)&&(i<=n)){if(a[i]==classkey){found=1; ck=i} i++}}END{print ck}'`
vk=`head -1 $input | awk -v valkey=$valkey '{n=split($0,a,"\t"); found=0; i=1; while((found==0)&&(i<=n)){if(a[i]==valkey){found=1; vk=i} i++}}END{print vk}'`
echo done >&2

# List of classes from input file
echo I am making a list of classes from the input file >&2
awk -v ck=$ck 'NR>=2{print $ck}' $input | sort | uniq > $b2\_classes.txt
echo done >&2

# Distribution of values for each class
echo I am computing the distribution of values for each class >&2
cat $b2\_classes.txt | while read c
do
echo class $c >&2
awk -v ck=$ck -v vk=$vk -v c=$c '((NR>=2)&&($ck==c)){print $vk}' $input > $b2\_class$c\_values.txt
$STAT $b2\_class$c\_values.txt 1 > $b2\_class$c\_values_stats.txt
cat $b2\_class$c\_values_stats.txt >&2
done
echo done >&2

# Extract the min and the max of all distributions
echo I am extracting the min and max values of all distributions >&2
m=`cat $b2\_classes.txt | while read c; do awk 'NR==3{print $2}' $b2\_class$c\_values_stats.txt; done | awk 'NR==1{m=$1}NR>=2{if($1<m){m=$1}}END{print m}'`
M=`cat $b2\_classes.txt | while read c; do awk 'NR==3{print $NF}' $b2\_class$c\_values_stats.txt; done | awk 'NR==1{M=$1}NR>=2{if($1>M){M=$1}}END{print M}'`
echo done >&2

# Assign a bin from 0 to no_bin-1 to each element by considering the min, the max and the number of bins
echo I am assigning a bin from 0 to no_bin-1 to each element of the input file \(using min, max and nb bins\) >&2
awk -v vk=$vk -v m=$m -v M=$M -v nobin=$nobin 'BEGIN{OFS="\t"; binlg=(M-m)/nobin}NR==1{print $0, "bin_no"}NR>=2{print $0, int(($vk-m)/binlg)}' $input > $b2\_binno.tsv
echo done >&2

# Make Aux dir and divide the input file in as many classes * bin files
echo I am dividing the input file in as many files as nb classes times nb bins >&2
mkdir -p Aux
cat $b2\_classes.txt | while read c
do
awk -v b2=$b2 -v ck=$ck -v c=$c '((NR>=2)&&($ck==c)){print $0 > "Aux/"b2"_class_"c"_binno_"$NF".tsv"}' $b2\_binno.tsv
done
echo done >&2

# Make a new tsv file with header that has as many rows as bins and providing for each bin the number of elements
# of each class from this bin as well as the min number of elements and the class where it was found
echo I am making a new tsv file with header that has for each bin the number of elements of each class for this bin as well as the min and the class of the min >&2
awk -v ck=$ck -v fileRef=$b2\_classes.txt 'BEGIN{shead="binno\t"; while (getline < fileRef >0){n++; class[n]=$1; shead=(shead)($1)("\t")} print (shead)("min_elt\tclass_min_elt")} NR>=2{nb[$ck,$NF]++; if($NF>M){M=$NF}} END{OFS="\t"; for(i=0; i<=M; i++){s=i"\t"; for(j=1; j<=n; j++){if(nb[class[j],i]==""){nb[class[j],i]=0} if(j==1){m=nb[class[j],i]; cm=class[j]}else{if(nb[class[j],i]<m){m=nb[class[j],i]; cm=class[j]}} s=(s)(nb[class[j],i])("\t")} print (s)(m)("\t")(cm);}}' $b2\_binno.tsv > binno_nbelt_eachclass_min_classofmin.tsv
echo done >&2

# Make Res dir and for each class go along the binned tsv file and for each bin if the min value is not zero 
# and the last column is not the class then sample a number of elements corresponding to this min in the class elt file
echo I am sampling for each class and each bin number when necesary >&2
mkdir -p Res
awk '{print $1, NR+1}' $b2\_classes.txt | while read c nc
do
awk -v c=$c -v nc=$nc '((NR>=2)&&($(NF-1)>0)&&($NF!=c)){print $1, $(NF-1), c}' binno_nbelt_eachclass_min_classofmin.tsv | while read bin m cl
do
$SAMPLE Aux/$b2\_class_$cl\_binno_$bin.tsv $m
mv $m.sample Res/genes_$cl\_binno_$bin\_$m.tsv
done
done
echo done >&2   

# For each class and for each bin where the min value was not zero, concatenate one of the two following files
# - the elt file if the class was the one with the min for this bin
# - the sampled elt file otherwise
echo I am concatenating the different bin files for each class >&2
awk '{print $1, NR+1}' $b2\_classes.txt | while read c nc
do
awk -v c=$c -v nc=$nc '((NR>=2)&&($(NF-1)>0)){print $1, $(NF-1), c, $NF}' binno_nbelt_eachclass_min_classofmin.tsv | while read bin m cl clm
do
if [ "$clm" = "$cl" ]
then
cat Aux/$b2\_class_$cl\_binno_$bin.tsv
else
cat Res/genes_$cl\_binno_$bin\_$m.tsv
fi 
done > $c\_elt_with_matched_values.tsv
done
echo done >&2

# Compute the distribution of values for each class after the equalization
echo I am computing the distribution of values for each class after the equalization process >&2
cat $b2\_classes.txt | while read c
do
echo class $c >&2
awk -v vk=$vk 'NR>=2{print $vk}' $c\_elt_with_matched_values.tsv > $c\_elt_with_matched_values_temp.tsv
$STAT $c\_elt_with_matched_values_temp.tsv 1 > $b2\_class$c\_values_after_equalization_stats.txt
cat $b2\_class$c\_values_after_equalization_stats.txt >&2
done
echo done >&2

# Make a final tsv file with equalized distributions and a header
echo I am making a final tsv file with equalized value distributions for each class >&2
awk 'NR==1' $b2\_binno.tsv > $b2\_samples_with_equal_distrib_each_class.tsv
cat $b2\_classes.txt | while read c
do
cat $c\_elt_with_matched_values.tsv
done >> $b2\_samples_with_equal_distrib_each_class.tsv
awk '{s=""; for(i=1; i<=NF; i++){s=(s)($i)("\t")} print s}' $b2\_samples_with_equal_distrib_each_class.tsv > $b2\_samples_with_equal_distrib_each_class_temp.tsv
mv $b2\_samples_with_equal_distrib_each_class_temp.tsv $b2\_samples_with_equal_distrib_each_class.tsv
echo done >&2

# Clean
echo I am cleaning >&2
cat $b2\_classes.txt | while read c
do
rm $b2\_class$c\_values.txt 
rm $b2\_class$c\_values_stats.txt 
rm $c\_elt_with_matched_values_temp.tsv
rm $c\_elt_with_matched_values.tsv
done
echo done >&2






