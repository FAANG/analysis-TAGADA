#!/bin/bash

# make_hist_rna_in_ref_better_until_gff.sh 
###########################################
# - same as make_hist_rna_in_ref_better.sh but stops after making the gff file (which
# is the action that takes more time). Then we can choose to make dist bins for a barplot
# that is either a reverse plot (number of features hit by elements in y axis) or a normal 
# (but unbugged) plot (number of elements in the feature in y axis)
#   * it uses the strand option of overlap instead of dividing the files
#   * it supposes the two input files are gff and have gff extension
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

# Be careful: you should not launch two instances of this script in the same directory
# since it is using files generally named tmp...

# Takes as input:
#################
# - a set of reference segments in gtf format (eg a kind of exons of a particular chr)
# - a set of nr/distinct RNA elements in gff format with information about redundancy in the 10th field
# - a flag about the strandedness of the nr/distinct RNA elements (= whether there is a strand 
#   in the gff file AND this strand actually MEANS something)
# - optionally a short name for the RNA experiment
# Provides as output:
#####################
# - the gff file of the initial features witgh 4 additional columns for the list of elements overlapping it
#   and the associated redundancy list
# - in the error file it provides the proportion of features hit by elements

# Usage:
########
# make_hist_rna_in_ref_better_until_gff.sh refseg.gtf nrrnaelt.gff strandedness cpkc
# where:
# - refseg.gtf is a set of reference segments in gtf format (eg a kind of exons of a particular chr)
# - nrrnaelt.gff is a set of nr/distinct RNA elements in gff format with information about redundancy in the 10th field
# - strandedness is a boolean for the strandedness of the nr/distinct RNA elements 
#  (= whether there is a strand in the gff file AND this strand actually MEANS something)

# Example:
##########
# make_hist_rna_in_ref_better_until_gff.sh /projects/encode/scaling_up/whole_genome/Gencode/Feb2009/gencode_data.rel2_exons_pcgtr_pcggn_of_longestcodingtr_withmorethan3ex_most5pex.gtf /projects/encode/scaling_up/whole_genome/CAGE/CAGE_PolyA+_K562_Cy_1_1_hg18.bed.gz.nr.gtf 1 


# Programs
###########
OVERLAP=../bin/overlap
GFF2GFF=gff2gff.awk

# If the user does not provide any reference segment file or any nr/distinct RNA element file or any strandedness boolean
#################################################################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo  >&2
    echo Usage:    make_hist_rna_in_ref_better_until_gff.sh refseg.gtf nrrnaelt.gtf strandedness expname >&2
    echo "         where:" >&2
    echo "         - refseg.gtf is a set of reference segments in gtf format" >&2
    echo "         - nrrnaelt.gff is a set of nr/distinct RNA elements in gff format with information about redundancy in the 10th field" >&2
    echo "         - strandedness is a boolean for the strandedness of the nr/distinct RNA elements" >&2
    echo "        - expname is a short name for the experiment (optional) " >&2
    echo "!!! Do not launch several times in the same dir since uses common file names !!! " >&2
    echo Example: make_hist_rna_in_ref_better_until_gff.sh /projects/encode/scaling_up/whole_genome/Gencode/Feb2009/gencode_data.rel2_exons_pcgtr_pcggn_of_longestcodingtr_withmorethan3ex_most5pex.gtf /projects/encode/scaling_up/whole_genome/CAGE/CAGE_PolyA+_K562_Cy_1_1_hg18.bed.gz.nr.gff 1 cpkc >&2
    echo >&2
    echo From a set of reference elements and a set of nr/distinct RNA elements, produces the gff of the features with list of overlapping elements and list of redundancy values. >&2
    echo >&2
    exit 1
fi

if [ ! -n "$4" ]
then
    exp="e"
else
    exp=$4
fi


# Input files and base names
#############################
REFSEG=$1
NRRNAELT=$2
refseg=${REFSEG%.gff}
refsegbase=`basename $refseg`
nrrnaelt=${NRRNAELT%.gff}


# Compute the overlap between the reference segments and the nr/distinct RNA elements
############################################################################# 
# (twice to keep the redundancy info (even though not used here)), in a way given 
#################################################################################
# by the strandedness parameter flag
#####################################
echo Computing the actual overlaps >&2
#######################################
echo First the list of coordinates of the overlapping nr/distinct RNA elements... >&2
######################################################################################
$OVERLAP $refseg.gff $nrrnaelt.gff -m -1 -st $3 -f $exp -v | awk -v e=$exp '{$(NF-3)=""; $(NF-2)=""; $(NF-1)=e; if($NF!="."){split($NF,a,","); k=1; s=""; while(a[k]!=""){split(a[k],b,"_"); s=(s)(b[2]"_"b[3])","; k++;} $NF=s;} print}' | awk -f $GFF2GFF > $refsegbase\_$exp.gff
echo Then the list of corresponding redundancy values... >&2
#############################################################
$OVERLAP $refsegbase\_$exp.gff $nrrnaelt.gff -m 10 -st $3 -f $exp -v | awk -v e=$exp '{$(NF-3)=""; $(NF-2)=""; $(NF-1)=e; print}' | awk -f $GFF2GFF > tmp
mv tmp $refsegbase\_$exp.gff
echo done >&2


# Compute the number and proportion of reference segments covered by nr/distinct RNA elements
###############################################################################################
echo I am computing the number and proportion of reference segments covered by nr/distinct RNA elements... >&2
awk -v e=$exp '{n++; if($NF!="."){n1++;}}END{print "#", e, n, (n1!="" ? n1 : 0), (n1!="" ? n1/n*100 : 0)}' $refsegbase\_$exp.gff >&2







