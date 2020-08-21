#!/bin/bash

# junc_txt_to_canonical.sh
# Take as input:
################
# - a tsv file of general junctions (could be chimeric) with header where the first column is the junction identifier in chimpipe's format 
#   (donchr_donpos_donstrand:accchr_accpos_accstrand)
# - a genome index in gem format made with gemtools-1.7.1-i3
# Provides as output in the working directory:
##############################################
# - a tsv file of 6 columns with no header that has that only has the junctions with non NA coordinates:
#   * in column no 1 the same information as the input file
#   * in column no 2 the 24 bp sequence surrounding the donor site (the donor splice site being in pos 13 and 14 of this sequence)
#   * in column no 3 the 24 bp sequence surrounding the acceptor site (the acceptor splice site is always on pos 11 and 12 of this sequence)
#   * in column no 4 the 2 bp sequence corresponding to the donor splice site
#   * in column no 5 the 2 bp sequence corresponding to the acceptor splice site
#   * in column no 6 a boolean (0 or 1) saying whether the junction is canonical

# example
# cd ~/Chimeras/Benchmark/Data/Edgren
# genome=/users/rg/projects/references/Genome/H.sapiens/hg19/gemtools1.7.1-i3/Homo_sapiens.GRCh37.chromosomes.chr.M.gem
# time ~sdjebali/bin/junc_txt_to_canonical.sh 46edgren_junctions.txt $genome
# input: 46edgren_junctions.txt
# juncid	gnpair
# chr17_35479453_-:chr17_37374426_-	ACACA-STAC2
# 47 (2 fields)
# output: 46edgren_junctions_ok_donacc_seq_can.tsv
# chr17_35479453_-:chr17_37374426_-	TCTGAAGCCAAGGTAGGCCACCTC	TTTTTCTCCCAGCTCCAGCGATTC	1
# 45 (4 fields)

if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: junc_txt_to_canonical.sh junctions.txt genome.gem >&2
    echo "Take as input:" >&2
    echo "   - a tsv file of general junctions (could be chimeric) with header where the first column is the junction identifier" >&2
    echo "     in chimpipe's format (donchr_donpos_donstrand:accchr_accpos_accstrand)" >&2 
    echo "   - a genome index in gem format made with gemtools-1.7.1-i3" >&2
    echo "Provides as output in the working directory:" >&2
    echo "   - a tsv file of 6 columns with no header that only has the junctions with non NA coordinates:" >&2
    echo "     * in column no 1 the same information as the input file" >&2
    echo "     * in column no 2 the 24 bp sequence surrounding the donor site (12 bp in exon, 12 bp in intron)" >&2
    echo "     * in column no 3 the 24 bp sequence surrounding the acceptor site (12 bp in intron, 12 bp in exon)" >&2
    echo "     * in column no 4 the 2 bp sequence corresponding to the donor splice site" >&2
    echo "     * in column no 5 the 2 bp sequence corresponding to the acceptor splice site" >&2
    echo "     * in column no 6 a boolean (0 or 1) saying whether the junction is canonical" >&2
    echo "!!! needs gem-retriever to be installed !!!" >&2
    echo "" >&2
    exit 1
fi

# Variables
###########
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
junc=$1
genome=$2
btmp=${junc%.txt}
btmp2=${btmp%.tsv}
b=`basename $btmp2`

# Programs
##########
CHIMTOBED=$rootDir/chim_txt_to_bedpe.awk 
BEDTODONACC=$rootDir/chim_bedpe_to_don_acc_tsv.awk
RETRIEVER=$rootDir/../bin/gem-retriever
CANONICAL=$rootDir/intron_is_canonical_gal.awk 

awk 'NR>=2&&$1!~/NA/{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); print $1, $2, ".", a1[2], a2[2]}' $junc | awk -f $CHIMTOBED > $b\_ok.bedpe
# chr17	35479452	35479453	chr17	37374425	37374426	chr17_35479453_-:chr17_37374426_-	.	-	-
# 45 (10 fields)

awk -f $BEDTODONACC $b\_ok.bedpe > $b\_ok_donacc.bedpe
# chr17	-	35479441	35479464
# 90 (4 fields)

cat $b\_ok_donacc.bedpe | $RETRIEVER $genome > $b\_ok_donacc_seq.txt
# TCTGAAGCCAAGGTAGGCCACCTC
# 90 (1 fields)

paste $b\_ok_donacc.bedpe $b\_ok_donacc_seq.txt | awk -v fileRef=$b\_ok.bedpe 'BEGIN{while (getline < fileRef >0){n++; id[n]=$7}} NR%2==1{mem=$0}NR%2==0{print mem, $0, id[NR/2]}' | awk -v flddon=5 -v fldacc=10 -f $CANONICAL | awk 'BEGIN{OFS="\t"} {split($5,a,""); split($10,b,""); print $11, $5, $10, a[13]a[14], b[11]b[12], $12}' > $b\_ok_donacc_seq_can.tsv
