#!/bin/bash

# sample_lines_better.sh

# usage: 
# sample_lines_better.sh [DATAPATH]inputfile numberlines 

# example:
# ~/People/Chays/sample_lines_better.sh /home/ug/cmanicha/METAGENOME/Sargasso_sea/sargasso.trimmed.fasta_random.tbl 40

# In case the user does not specify any input file
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
echo "*** you must provide an input file and a number of lines to randomly extract from it ***" >&2
echo "usage: sample_lines.sh [DATAPATH]inputfile numberlines" >&2
echo "example: ~/People/Chays/sample_lines.sh /home/ug/cmanicha/METAGENOME/Sargasso_sea/sargasso.trimmed.fasta_random.tbl 40" >&2
exit 1
fi


# Retrieve the program path and the data path from the executable name and the gff file name
# this in order to be able to launch this script from anywhere and have the output at the same
# place as the input data
PGMPATH=`dirname $0`
DATAPATH=`dirname $1`

echo 'data=read.table("'$1'"); write.table(data[sample(dim(data)[1],'$2'),],"'$2'.sampletmp")' | R --vanilla

awk 'BEGIN{OFS="\t";}NR>=2{s=""; for(i=2; i<=NF; i++){if($i~/"/){split($i,a,"\""); s=(s)(a[2])" ";}else{s=(s)($i)" ";}} print s}' $2\.sampletmp > $2\.sample
# $DATAPATH\/$2\.sample
#$2\.sample

rm $2\.sampletmp