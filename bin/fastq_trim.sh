#!/bin/bash
myname=`basename $0`
usage() { echo "###########################################################################################################
                   $myname: trim fastq paired reads

Usage: $myname <reads_R1.fastq(.gz)> <site>

 reads_R1.fastq(.gz) : Input reads in fastq format (suffix '_R1.fastq' is required ; 'R2' will be processed too)
                       The file name has to end by _R1.fastq (can be compressed though). Full path expected.
                       A similar file with R2 in the name has to be present in the same directory.
                       The fastq file has to contain one entry every 4 lines: NO MULTI-LINE SEQUENCES ALLOWED!

 site: motif to be trimmed (adapter or restriction site etc). For instance: AAGCTAGCTT for HindIII religation site (Hi-C)

Output (in the same directory as fastq input): 
 reads_R1_trim.fastq.gz
 reads_R2_trim.fastq.gz

WARNING: 
=> no more than 4 lines per read (no multi-line sequence is allowed)
###########################################################################################################" 1>&2; exit 1;}

fastq1=$1 # reads (first in pair)
if [ ! -r "$fastq1" ]
 then
  echo "ERROR: could not read input sequences $fastq1"
  usage
fi

if [ ! "$2" ]
then  
    echo "Please provide a subsequence to look for (adapter or alike)"
  usage
fi

echo "#### Starting $0 - `date`"
fastq2=${fastq1/R1.fastq/R2.fastq} # reads (second in pair)
trimmed1=${fastq1%_R1*}_trimmed_R1.fastq
trimmed2=${fastq2%_R2*}_trimmed_R2.fastq
trim1=${trimmed1%_trimmed_R1*}_trim_R1.fastq
trim2=${trimmed2%_trimmed_R2*}_trim_R2.fastq
bin=/work2/project/fragencode/tools
seq=$2 

echo "## Looking for site in reads..."

cutadapt -a $seq -m 20 -e 0 $fastq1 2> $trim1.log > $trimmed1
echo "## ...done for R1 in $trim1.log"

cutadapt -a $seq -m 20 -e 0 $fastq2 2> $trim2.log > $trimmed2
echo "## ...done for R2 in $trim2.log"
   
echo "## Matching read pair IDs after trimming..."
   
id1=$trimmed1.fastqIDlist1
id2=$trimmed2.fastqIDlist2
ids=$trimmed1.fastqIDs
   
##  ONLY FOR DIXON DATA
#    awk 'NR%4==1' $trimmed1 | sort -u | sed 's/^@//' > $id1
#    awk 'NR%4==1' $trimmed2 | sort -u | sed 's/^@//' > $id2
#    comm -12 $id1 $id2 > $ids
#    cat $trimmed1 | $bin/fastq_extract $ids > $trim1
#    cat $trimmed2 | $bin/fastq_extract $ids > $trim2
##

## FRAGENCODE DATA
awk 'NR%4==1{print $1}' $trimmed1 | sort -u | sed 's/^@//' > $id1
awk 'NR%4==1{print $1}' $trimmed2 | sort -u | sed 's/^@//' > $id2
comm -12 $id1 $id2 > $ids
suffix1=`head -1 $trimmed1 | awk '{print $2}'`
suffix2=`head -1 $trimmed2 | awk '{print $2}'`
cat $ids | awk -v sfx=$suffix1 '{print $1,sfx}' > $ids.full1
cat $ids | awk -v sfx=$suffix2 '{print $1,sfx}' > $ids.full2
cat $trimmed1 | $bin/fastq_extract $ids.full1 > $trim1
cat $trimmed2 | $bin/fastq_extract $ids.full2 > $trim2
##
# clean up
rm $id1
rm $id2
rm $ids
rm $ids.full1
rm $ids.full2
rm $trimmed1
rm $trimmed2
# gzip  reads
echo "## gzipping the new fastq..."
gzip $trim1
gzip $trim2
echo "#### Completed $0 - `date`"
