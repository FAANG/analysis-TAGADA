#!/bin/bash

# divide_into_chr.sh 

# takes as input a file and its format (gff (default), gtf, bed) and makes as many files as chr, ie from 1 to 22 and then X, Y and M
# note: removes the comment lines

# usage:
#  divide_into_chr.sh file format

if [ ! -n "$1" ]
then
echo "usage: divide_into_chr.sh file format" >&2
echo "       where format is gff, gtf or bed (default is gff)" >&2
exit 0
fi

if [ ! -n "$2" ]
then
format="gff"
else
format=$2
fi

input=$1
awk -v format=$format -v base=`basename ${input%"."$format}` '$1!~/#/{print > base"_"$1"."format}' $input