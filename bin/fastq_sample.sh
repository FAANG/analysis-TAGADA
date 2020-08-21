#!/bin/bash

# usage
# fastq_sample.sh R1.fastq.gz R2.fastq.gz pcent

# Improved on June 15th 2020 so that intermediate files are named after the first fastq.gz file and uses wc -l instead of an awk script to count rows in count file

# example
#########
# cd ~/fragencode/workspace/sdjebali/irsd/courses/rnaseq/SIB_august2020/analysis/sample_1pcent
# datadir=~/fragencode/workspace/sdjebali/irsd/courses/rnaseq/SIB_august2020/data
# pgm=~/fragencode/tools/multi/Scripts/Bash/fastq_sample.sh
# time $pgm $datadir/ENCFF000EWJ.fastq.gz $datadir/ENCFF000EWX.fastq.gz 1 2> ENCFF000EWJ.fastq_sample.err
# real	13m39.501s  *** makes 2*85M files of 1 million reads each from two 8.1G files of 120 million reads each 



if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo "Usage: fastq_sample.sh R1.fastq.gz R2.fastq.gz pcent" >&2
    echo "" >&2
    exit -1
else
    R1=$1
    R2=$2
    pcent=$3
fi

basetmp1=`basename $R1`
base1=${basetmp1%.fastq.gz}

basetmp2=`basename $R2`
base2=${basetmp2%.fastq.gz}

zcat $R1 | awk 'NR%4==1{n++; print n}' > $base1.count.txt
nb=`wc -l $base1.count.txt | awk '{print $1}'`
shuf $base1.count.txt | head -n $((nb*$pcent/100)) > $base1.sample_count.txt
zcat $R1 | awk -v fileRef=$base1.sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base1.fastq.gz
zcat $R2 | awk -v fileRef=$base1.sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base2.fastq.gz

rm $base1.count.txt
