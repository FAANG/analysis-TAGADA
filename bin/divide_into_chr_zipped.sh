#!/bin/bash

# divide_into_chr_zipped.sh 

# takes as input a gzipped file and its format (gff (default), gtf, bed) and makes as many files as chr, ie from 1 to 22 and then X, Y and M

# usage:
# divide_intro_chr_zipped.sh file format

if [ ! -n "$1" ]
then
echo "usage: divide_intro_chr_zipped.sh file format" >&2
echo "       where format is gff (default), gtf or bed" >&2
exit 0
fi

if [ ! -n "$2" ]
then
format="gff"
else
format=$2
fi

zcat $1 | awk -v f=$format -v base=${1%"."$format".gz"} '$1!~/#/{print > base"_"$1"."f}' 