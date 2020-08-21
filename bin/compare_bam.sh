#!/bin/bash


# make sam file of the reads that map
samtools view -F 260 $1 | sort -T ~/Tmp > file1.sam
samtools view -F 260 $2 | sort -T ~/Tmp > file2.sam

# make the list of common rows
##############################
cat file1.sam file2.sam | sort -T ~/Tmp | uniq -c | awk '$1==2{k=2; s=""; while($k!=""){s=(s)($k)("\t"); k++} print s}' > same.sam

# make the list of different rows from her file and my file
###########################################################
awk -v fileRef=same.sam 'BEGIN{while (getline < fileRef >0){ok[$0]=1}} ok[$0"\t"]!=1' file1.sam > file1_diff_file2.sam
awk -v fileRef=same.sam 'BEGIN{while (getline < fileRef >0){ok[$0]=1}} ok[$0"\t"]!=1' file2.sam > file2_diff_file1.sam
wc -l same.sam file1_diff_file2.sam file2_diff_file1.sam
