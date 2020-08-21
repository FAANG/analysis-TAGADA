#!/bin/bash

# add_overlap_info_to_junctions.sh
# Takes as input
# - a tsv file with header with at least 1 column and where the 1st column has to be junctions in chimpipe's format (donchr_donpos_donstr:accchr_accpos_accstr)
# - an annotation file in gtf or gff2 format that has to contain at least the features specified in the next argument 
# - the features from the annotation one wants to intersect with each part of each junction (e.g. exon)
# - the key from the annotation corresponding to the information to be reported for each part of the junctions (e.g. gene_name)
# Provides as output a tsv file of junctions that has exactly the same information as in the input tsv file with two additional columns
# corresponding to the information wanted by the user for the overlap of the features with each part of the junction

# BE careful:
# should not be run in parallel in the same working dir since produces files with fixed names

# example
# cd ~/Chimeras/Benchmark/Data
# junctions=/users/rg/sdjebali/Chimeras/Benchmark/Data/Berger/16berger_junctions_ok.txt
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/gencode.v19.annotation.gtf
# time add_overlap_info_to_junctions.sh $junctions $annot exon gene_name
# real    0m21.347s


# Check the 4 arguments are present
###################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] 
then
    echo "" >&2
    echo Usage: add_overlap_info_to_junctions.sh junctions.tsv annot.gtf feature gtfkey >&2
    echo "where:" >&2
    echo "- junctions.tsv is a junction file in tsv format with header and at least 1 column corresponding to junctions chimpipe's format (donchr_donpos_donstr:accchr_accpos_accstr)" >&2
    echo "- annot.gtf is an annotation file in gtf or gff2 format with at least rows corresponding to the feature name provided in next argument" >&2
    echo "- feature is the name of the features one wants to intersect with each part of each junction and report information about" >&2
    echo "- gtf key is the key in the gtf file corresponding to the information the user wants to obtain for each part of the junction" >&2
    echo "It will provide as output a tsv file of junctions that has exactly the same information as in the input tsv file with two additional columns" >&2
    echo "corresponding to the information wanted by the user for the overlap of the features with each part of the junction" >&2
    echo "!!! Be careful: should not be run in parallel in the same working directory since it produces files with fixed names!!!" >&2
    echo "" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
junctions=$1
annot=$2
feat=$3
key=$4

basetmp=`basename $junctions`
base=${basetmp%.tsv}

i=$key\_from_$feat

# Programs
###########
GFF2GFF=$rootDir/gff2gff.awk
OVERLAP=$rootDir/../bin/overlap

# Make a gff file for the junctions
###################################
awk 'NR>=2{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); print a1[1], ".", "d", a1[2], a1[2], ".", a1[3], ".", "juncid", "\""$1"\"\;"; print a2[1], ".", "a", a2[2], a2[2], ".", a2[3], ".", "juncid", "\""$1"\"\;"}' $junctions | awk -f $GFF2GFF > junctions.gff

# Make an annotation file that only has the features indicated by the user and only the key wanted by the user
##############################################################################################################
awk -v feat=$feat -v key=$key '$3==feat{s=""; for(i=1; i<=8; i++){s=(s)($i)("\t")} for(i=9; i<=(NF-1); i+=2){if($i==key){s=(s)(key)(" ")($(i+1))}} print s}' $annot | awk -f $GFF2GFF > annot.gff

# Overlap the two last gff file together, the key should now be in field no 10 in the annotation file
#####################################################################################################
$OVERLAP junctions.gff annot.gff -f $key -m 10 -nr -o junctions_over_annot.gff

# Reports the information back to the initial tsv file
######################################################
awk -v fileRef=junctions_over_annot.gff -v i=$i 'BEGIN{OFS="\t"; while (getline < fileRef >0){split($10,a,"\""); info[a[2],$3]=$14}} NR==1{print $0, i"_A", i"_B"} NR>=2{gsub(/\"/,"",info[$1,"d"]);gsub(/\;/,"",info[$1,"d"]); gsub(/\"/,"",info[$1,"a"]);gsub(/\;/,"",info[$1,"a"]); print $0, info[$1,"d"], info[$1,"a"]}' $junctions > $base\_with_$i.tsv

# Clean
#######
rm junctions.gff annot.gff junctions_over_annot.gff