#!/bin/bash

# gtf2fasta.sh
# Takes as input an annotation file in gtf format with gene id followed by transcript id in the 9th field
# as well as a genome index file in gem format, and outputs a fasta file with the nucleotide sequences
# of all the transcripts present in the annotation file
# Note: bedtools getfasta could be able to do the same, need to test it and see what takes less time (since we use an index it might be the present one)

# Usage
# gtf2fasta.sh annot.gtf genome_index.gem

# Example
# cd /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/seq
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/gen10.gtf
# genome=/users/rg/projects/references/Genome/H.sapiens/hg19/gem2/Homo_sapiens.GRCh37.chromosomes.chr.M.fa.gem
# time gtf2fasta.sh $annot $genome 2> gtf2fasta.err > gtf2fasta.out
# real	3m18.403s  *** and checked using fastlength that the length is the same as when obtained with julien's perl script extract_spliced_transc_seqs.pl

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: gtf2fasta.sh annot.gtf genome_index.gem >&2
    echo This script taks as input an annotation file \in gtf format with gene id followed by transcript id \in the 9th field >&2
    echo as well as a genome index file \in gem format, and outputs a fasta file with the nucleotide sequences of all the transcripts >&2
    echo present \in the annotation file >&2
    echo "!!! Needs gem-retriever software installed !!!" >&2
    echo "" >&2
    exit 1
else
annot=$1
genome=$2
fi

# Variable from input
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

b=`basename $annot`
b2tmp=${b%.gtf}
b2=${b2tmp%.gff}


# Programs
##########
RETRIEVER=$rootDir/../bin/gem-retriever 
# for example in /users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/ChimPipe/bin/gemtools-1.7.1-i3/


# Make the list of distinct exon coordinates (fast)
###################################################
echo I am making the list of distinct exon coordinates >&2
awk 'BEGIN{OFS="\t"} $3=="exon"{print $1, $7, $4, $5}' $annot | sort | uniq > $b2\_distinct_exon_coord.tsv 
echo done >&2
# chr10	+	100003848	100004106
# 518897 (4 fields)


# Retrieve the exon sequences (less than 2 minutes)
###################################################
echo I am retrieving the exon sequences >&2
cat $b2\_distinct_exon_coord.tsv | $RETRIEVER $genome > $b2\_distinct_exon_coord.seq
echo done >&2
# AGAGAAAGCGGTTGGAAGCCAAGCAACGGGAAGACATCTGGGAAGGCAGAGACCAGTCTACAGTTTGAACATCACTCAATGAAAGGGATAATTCCATGAATCAGAAAATGTTTCCATAGCCTTCAGATAAGATGATCCTTCCAGAGCTCTATGTACATGCAGATGTGCATGTTAAAGAGATAAAGTGATCGAGACAAGGACTGACTGGGTATAGAAGGAAGACAGACTCCTGTCTTCACTCCTAAATGCAGTTCTTTGG
# 518897 (1 fields)

# Make a file that both has the exon coordinates and sequence (quite fast)
##########################################################################
echo I am making a file that both has the exon coordinates and sequence >&2
paste $b2\_distinct_exon_coord.tsv $b2\_distinct_exon_coord.seq | awk '{print $1"_"$3"_"$4"_"$2, $5}' > $b2\_distinct_exon_coord_seq.txt
echo done >&2
# chr10_100003848_100004106_+ AGAGAAAGCGGTTGGAAGCCAAGCAACGGGAAGACATCTGGGAAGGCAGAGACCAGTCTACAGTTTGAACATCACTCAATGAAAGGGATAATTCCATGAATCAGAAAATGTTTCCATAGCCTTCAGATAAGATGATCCTTCCAGAGCTCTATGTACATGCAGATGTGCATGTTAAAGAGATAAAGTGATCGAGACAAGGACTGACTGGGTATAGAAGGAAGACAGACTCCTGTCTTCACTCCTAAATGCAGTTCTTTGG
# 518897 (2 fields)


# For each transcript make a list of exon coordinates from 5' to 3' (a bit slow)
################################################################################
echo For each transcript I am making a list of exon coordinates from 5\' to 3\' >&2
awk '$3=="exon"' $annot | sort -k12,12 -k4,4n -k5,5n | awk '{nbex[$12]++; strand[$12]=$7; ex[$12,nbex[$12]]=$1"_"$4"_"$5"_"$7}END{for(t in nbex){s=""; split(t,a,"\""); if(strand[t]=="+"){for(i=1; i<=nbex[t]; i++){s=(s)(ex[t,i])(",")}} else{if(strand[t]=="-"){for(i=nbex[t]; i>=1; i--){s=(s)(ex[t,i])(",")}}} print a[2], s}}' > $b2\_trid_exonlist_5pto3p.txt
echo done >&2
# ENST00000556118.1 chr15_91573118_91573226_+,
# 173599 (2 fields)


# For each transcript make its sequence by concatenating the sequences of its exons from 5' to 3' (15 seconds)
##############################################################################################################
echo For each transcript I am making its sequence by concatenating the sequences of its exons from 5\' to 3\'  >&2
awk -v fileRef=$b2\_distinct_exon_coord_seq.txt 'BEGIN{while (getline < fileRef >0){seqex[$1]=$2}} {split($2,a,","); k=1; while(a[k]!=""){seqtr[$1]=(seqtr[$1])(seqex[a[k]]); k++}} END{for(t in seqtr){s=seqtr[t]; print ">"t; n=length(s); n2=int(n/60); for(i=0; i<=(n2-1); i++){print substr(s,i*60+1,60)} if(n>n2*60){print substr(s,n2*60+1,n-n2*60)}}}' $b2\_trid_exonlist_5pto3p.txt > $b2\_tr.fasta
echo done >&2


# Clean
#######
echo I am cleaning >&2
rm $b2\_distinct_exon_coord.tsv $b2\_distinct_exon_coord.seq 
rm $b2\_distinct_exon_coord_seq.txt $b2\_trid_exonlist_5pto3p.txt
echo done >&2
