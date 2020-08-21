#!/bin/bash

# divide_intro_strand.sh 

# takes as input a file and its format (gff, gtf, bed, default is gff) and makes as many files as strands, usually two (+ and -)

# usage:
# divide_intro_strand.sh file format

if [ ! -n "$1" ]
then
echo "usage: divide_intro_strand.sh file format" >&2
echo "       where format is gff (default), gtf or bed" >&2
exit 0
fi

if [ ! -n "$2" ]
then
format="gff"
else
format=$2
fi

base=`basename ${1%%"."$format}`
awk -v f=$format -v base=$base '$1!~/#/{if(f=="bed"){str=$6}else{str=$7} print > base"_"((str=="+") ? "plus" : ((str=="-") ? "minus" : str))"."f}' $1