#!/bin/bash

awk -v no=$2 '{if(no=="") print $1; else print $no}' $1 > /tmp/coucoufromsarah.txt
echo 'd<-read.table("/tmp/coucoufromsarah.txt")[,1]; title <- "'$1'"; postscript(file="'$1'.eps",onefile=FALSE,paper="a4",pointsize=12,horizontal=TRUE); hist(d,main=title,breaks=50,col=("lightblue")); dev.off()' | R --vanilla
for i in $1.eps ; do { f=${i%.eps}.pdf ; convert -rotate 90 $i $f ; } ; done
rm $1.eps
